import pysam
import pytest

from sbam.base.base import BaseAccuracy
from sbam.cli import main as sbam_main
from sbam.struct.eval import JunctionEvaluator
from sbam.struct.ori import BioVal, _find_dnaa_box_positions
from sbam.struct.struct import GenomeStructure


def _write_bam(bam_path, header_refs, reads_spec):
    """Build a coordinate-sorted, indexed BAM using pysam."""
    unsorted_path = str(bam_path) + ".unsorted.bam"
    header = {"HD": {"VN": "1.6"}, "SQ": [{"SN": sn, "LN": ln} for sn, ln in header_refs]}
    with pysam.AlignmentFile(unsorted_path, "wb", header=header) as outf:
        for spec in reads_spec:
            seg = pysam.AlignedSegment()
            seg.query_name = spec["name"]
            seg.reference_id = spec.get("ref_id", 0)
            seg.reference_start = spec["start"]
            seg.flag = spec.get("flag", 0)
            seg.mapping_quality = spec.get("mapq", 60)
            seg.cigartuples = spec["cigar"]
            seg.query_sequence = spec["seq"]
            seg.query_qualities = pysam.qualitystring_to_array("I" * len(spec["seq"]))
            outf.write(seg)
    pysam.sort("-o", str(bam_path), unsorted_path)
    pysam.index(str(bam_path))


def test_cigar_eqx_and_m_and_supplementary_exclusion_and_tracks(tmp_path):
    """
    Integration test with real pysam BAMs:
    - Tests `--eqx` (=/X/I/D) and `M` CIGAR sweeps
    - Tests supplementary alignment (0x800) exclusion (Bug 5)
    - Tests SNP/Indel stratification and read-strand bias Fisher test
    - Tests --bam CLI execution and track exports (.bed, .tsv, .vcf)
    """
    seq = "ACGT" * 250  # 1000 bp
    fasta_path = tmp_path / "asm.fa"
    fasta_path.write_text(f">chr1\n{seq}\n>empty_contig\n\n", encoding="utf-8")

    # Verify GenomeStructure drops 0-length contig (Bug 7)
    gs = GenomeStructure(str(fasta_path), str(tmp_path))
    _cyclic_fa, orig_lengths = gs.create_cyclic_buffer(buffer_size=100)
    assert "empty_contig" not in orig_lengths
    assert orig_lengths == {"chr1": 1000}

    bam_path = tmp_path / "test.bam"
    reads = []
    # 20 forward reads + 20 reverse reads with =/X/I/D
    # At pos 100 ('A'): forward reads match (20/20), reverse reads mismatch 'C' (0/20) -> 50% conc, strand-biased!
    # Wait: make 12 reverse reads mismatch and 8 forward match -> 8/20? Better:
    # 15 forward reads: 100= + 1X (pos 100) + 899=  -> wait: let's have 10 fwd match, 10 rev mismatch? Then 10/20 = 0.50 (not < 0.50)
    # Use 5 fwd match + 15 rev mismatch at pos 100 -> 5/20 = 0.25 (<0.50 consensus error, strand p < 0.01)!
    # At pos 200: 15 reads have 1I insertion after pos 200 -> 5/20 = 0.25 (<0.50 consensus error, dominant_error_type=INS)!
    # At pos 300: 15 reads have 1D deletion -> 5/20 = 0.25 (<0.50 consensus error, dominant_error_type=DEL)!
    for i in range(5):
        # Clean forward read covering 0..1000 using '=' (op 7)
        reads.append(
            {
                "name": f"fwd_clean_{i}",
                "start": 0,
                "flag": 0,
                "cigar": [(7, 1000)],
                "seq": seq,
            }
        )
    for i in range(15):
        # Reverse read with mismatch at pos 100, insertion after pos 200, deletion at pos 300
        # Reference consumed: 100 + 1 + 99 + 99 + 1 + 700 = 1000
        # Query consumed:     100 + 1 + 99 + 1(I) + 99 + 700 = 1000
        q_seq = list(seq)
        q_seq[100] = "C"  # mismatch at 100
        q_str = "".join(q_seq[:200]) + "T" + "".join(q_seq[200:299]) + "".join(q_seq[300:])
        reads.append(
            {
                "name": f"rev_err_{i}",
                "start": 0,
                "flag": 16,  # reverse
                "cigar": [(7, 100), (8, 1), (7, 99), (1, 1), (7, 99), (2, 1), (7, 700)],
                "seq": q_str,
            }
        )

    # Inject a supplementary alignment (0x800) with mismatch at pos 500; must be ignored (Bug 5)!
    reads.append(
        {
            "name": "supp_ignore",
            "start": 490,
            "flag": 0x800,
            "cigar": [(8, 20)],
            "seq": "T" * 20,
        }
    )

    # Also add a spanning read across cyclic seam 1000 (spanning 500..1100, i.e., buffer 100 bp)
    #Wait: for 1000 bp contig, min_overlap is min(500, 1000//2) = 500; buffer=100 -> linear mode or 2000 bp contig
    _write_bam(bam_path, [("chr1_cyclic", 1100)], reads)

    acc = BaseAccuracy(str(bam_path), str(fasta_path), threads=1, k_list=[4])
    acc.build_fidelity_profile()

    assert acc.consensus_error_bases == 3  # pos 100 (SNP), pos 199 (INS anchor), pos 299 (DEL)
    assert acc.consensus_error_snp_sites == 1
    assert acc.consensus_error_indel_sites == 2
    assert acc.read_strand_biased_sites >= 1

    # Also test M (op 0) fallback on linear contig
    bam_m_path = tmp_path / "test_m.bam"
    m_reads = [
        {"name": f"m_{i}", "start": 0, "flag": 0, "cigar": [(0, 1000)], "seq": seq}
        for i in range(10)
    ]
    _write_bam(bam_m_path, [("chr1", 1000)], m_reads)
    acc_m = BaseAccuracy(str(bam_m_path), str(fasta_path), threads=1, k_list=[4], linear_mode=True)
    acc_m.build_fidelity_profile()
    assert acc_m.consensus_error_bases == 0
    assert acc_m.evaluated_bases == 1000

    # Test CLI --bam end-to-end with --linear and track emission
    outdir = tmp_path / "cli_out"
    rc = sbam_main(
        [
            "-a", str(fasta_path),
            "--bam", str(bam_m_path),
            "-o", str(outdir),
            "--prefix", "lin_sample",
            "--linear",
            "--offline",
        ]
    )
    assert rc == 0
    assert (outdir / "lin_sample_sbam_report.html").exists()
    assert (outdir / "lin_sample_sbam_metrics.json").exists()
    assert (outdir / "lin_sample_masked_regions.bed").exists()
    assert (outdir / "lin_sample_consensus_errors.tsv").exists()
    assert (outdir / "lin_sample_consensus_errors.vcf").exists()


def test_even_palindrome_gatc_position1_and_noncanonical_5mer_merging(tmp_path):
    """
    Bug 6: GATC position 1 ('A') error is detected!
    Bug 9: canonical=False does not crash _cluster_into_families and 5-mers with 0% overlap don't collapse.
    """
    unit = "TTTTTGATCTTTTT"  # GATC at indices 5..9; 'A' is at index 6 (offset 1 = k//2 - 1)
    seq = unit * 200
    seq_len = len(seq)
    fasta_path = tmp_path / "gatc.fa"
    fasta_path.write_text(f">chr1\n{seq}\n", encoding="utf-8")

    # Build BAM where base at index 6 ('A') of every repeat has 4/10 mismatches (60% concordance)
    bam_path = tmp_path / "gatc.bam"
    reads = []
    for i in range(10):
        q_list = list(seq)
        if i < 4:
            for r_i in range(200):
                q_list[r_i * len(unit) + 6] = "G"
        reads.append(
            {
                "name": f"r_{i}",
                "start": 0,
                "flag": 0,
                "cigar": [(0, seq_len)],
                "seq": "".join(q_list),
            }
        )
    _write_bam(bam_path, [("chr1", seq_len)], reads)

    acc = BaseAccuracy(
        str(bam_path),
        str(fasta_path),
        threads=1,
        k_list=[4],
        canonical=True,
        linear_mode=True,
    )
    families = acc.analyze_motifs(min_occurrences=100, min_discordant_count=10)
    gatc = [f for f in families if f["motif"] == "GATC/GATC"]
    assert len(gatc) == 1
    assert gatc[0]["occurrences"] == 200
    assert gatc[0]["discordant_count"] == 200

    # Non-canonical mode must not crash on split("/") and detects TGAT (offset 2 -> pos 6)
    acc_nc = BaseAccuracy(
        str(bam_path),
        str(fasta_path),
        threads=1,
        k_list=[4],
        canonical=False,
        linear_mode=True,
    )
    families_nc = acc_nc.analyze_motifs(min_occurrences=100, min_discordant_count=10)
    assert any(f["motif"] == "TGAT" for f in families_nc)


def test_20_to_49bp_blind_spot_fixed_and_more_than_50_masked_regions(tmp_path):
    """
    Bug 4:
    - 30 bp low-concordance block (20-49 bp) is captured in masked_regions & low_concordance_blocks.
    - >50 high-depth blocks do NOT drop a shorter <90% low-concordance block from low_concordance_blocks.
    - Two adjacent low-concordance blocks separated by <500 bp are NOT diluted by intervening 100% bases.
    """
    # Build a 120 kb linear chromosome
    seq_len = 120_000
    seq = ("ACGT" * (seq_len // 4))
    fa = tmp_path / "mask.fa"
    fa.write_text(f">chr1\n{seq}\n", encoding="utf-8")

    # 10 baseline reads covering whole contig
    reads = [
        {"name": f"bg_{i}", "start": 0, "cigar": [(7, seq_len)], "seq": seq}
        for i in range(10)
    ]
    # Add 55 separate 80 bp high-depth repeat spikes (30x extra -> 40x > 2.5*10=25x), spaced 1000 bp apart
    for idx in range(55):
        pos = 2000 + idx * 1500
        for r in range(30):
            reads.append(
                {
                    "name": f"hd_{idx}_{r}",
                    "start": pos,
                    "cigar": [(7, 80)],
                    "seq": seq[pos : pos + 80],
                }
            )
    # Add a 30 bp unpolished block at 110,000..110,030 (8/10 reads have mismatches -> 20% concordance!)
    # AND another 30 bp block 200 bp later at 110,230..110,260 (within merge_gap_bp=500!)
    # Wait: remove 8 background reads at those 30 bp? Easier: add 30 mismatching reads? No, that changes depth!
    # Instead, let background reads be 2 reads + 8 reads? Let's just use 2 clean bg reads + 8 mismatching bg reads? No,
    # bg reads above were 10 clean reads; let's replace 8 of them with reads that mismatch at 110000..110030 and 110230..110260!
    for i in range(8):
        # CIGAR over 120,000: 110000= + 30X + 200= + 30X + 9740=
        reads[i]["cigar"] = [(7, 110000), (8, 30), (7, 200), (8, 30), (7, 9740)]

    bam_path = tmp_path / "mask.bam"
    _write_bam(bam_path, [("chr1", seq_len)], reads)

    acc = BaseAccuracy(str(bam_path), str(fa), threads=1, k_list=[4], linear_mode=True)
    acc.build_fidelity_profile()

    assert len(acc.low_concordance_blocks) >= 1
    assert any(b["mean_concordance"] == pytest.approx(20.0, abs=0.5) for b in acc.low_concordance_blocks)
    assert any(b["mean_concordance"] < 90.0 for b in acc.masked_regions)


def test_5mer_non_overmerging_and_exclusive_consensus_errors(tmp_path):
    """
    Bug 9: Two distinct 5-mers sharing a 4-bp suffix (e.g. ACCAG vs TCCAG) with 0% overlap do not merge.
    Bug 1: Overlapping motifs share consensus errors exclusively so sum(consensus_error_count) <= consensus_error_bases.
    """
    # Part A: Two distinct 5-mers sharing 4-bp suffix 'CCAG' (A[CCAG] vs T[CCAG]) with disjoint coords
    # Note: CCAGG/CCTGG vs ACCAG/CTGGT: ACCAG (pos 2 = C) vs TCCAG (pos 2 = C)
    unit_a = "AAAAAACCAGAAAAA"  # ACCAG at 5..10, mid=7 ('C')
    unit_b = "AAAAATCCAGAAAAA"  # TCCAG at 5..10, mid=7 ('C')
    seq = (unit_a * 150) + (unit_b * 150)
    seq_len = len(seq)
    fa = tmp_path / "five.fa"
    fa.write_text(f">chr1\n{seq}\n", encoding="utf-8")

    # 10 reads: 4 reads mismatch at index 7 of every 15-bp repeat (60% concordance)
    reads = []
    for i in range(10):
        ql = list(seq)
        if i < 4:
            for r in range(300):
                ql[r * 15 + 7] = "T"
        reads.append({"name": f"r_{i}", "start": 0, "cigar": [(0, seq_len)], "seq": "".join(ql)})
    bam = tmp_path / "five.bam"
    _write_bam(bam, [("chr1", seq_len)], reads)

    acc = BaseAccuracy(str(bam), str(fa), threads=1, k_list=[5], canonical=True, linear_mode=True)
    fams = acc.analyze_motifs(min_occurrences=100, min_discordant_count=10)
    # Must remain 2 distinct families (ACCAG/CTGGT and CTGGA/TCCAG), NOT collapsed via 4-bp core!
    assert len(fams) == 2

    # Part B: Overlapping 5-mer & 6-mer sharing consensus error sites -> exclusive partitioning (Bug 1)
    unit_c = "TTTTACCAGGTTTT"  # ACCAG (5..10, mid 7) and CCAGG (6..11, mid 8) and ACCAGG (5..11, mid 8)
    seq_c = unit_c * 200
    fa_c = tmp_path / "cons.fa"
    fa_c.write_text(f">chr1\n{seq_c}\n", encoding="utf-8")
    reads_c = []
    for i in range(10):
        ql = list(seq_c)
        if i < 6:
            # 6/10 mismatch at pos 7 and pos 8 -> 40% concordance (<50% consensus errors!)
            for r in range(200):
                ql[r * len(unit_c) + 7] = "T"
                ql[r * len(unit_c) + 8] = "T"
        reads_c.append({"name": f"r_{i}", "start": 0, "cigar": [(0, len(seq_c))], "seq": "".join(ql)})
    bam_c = tmp_path / "cons.bam"
    _write_bam(bam_c, [("chr1", len(seq_c))], reads_c)

    acc_c = BaseAccuracy(str(bam_c), str(fa_c), threads=1, k_list=[5, 6], canonical=True, linear_mode=True)
    fams_c = acc_c.analyze_motifs(min_occurrences=100, min_discordant_count=10)
    assert sum(f["consensus_error_count"] for f in fams_c) <= acc_c.consensus_error_bases


def test_junction_clip_window_and_small_contig_and_linear_termini(tmp_path):
    """
    Bug 2: soft-clipped reads ending 200 bp before junction are caught; adapter clips 700 bp away are neutral.
    Bug 7: <500 bp circular contigs can pass with scaled min_overlap.
    1.1: linear mode evaluates termini cleanly.
    """
    fasta_path = tmp_path / "contigs.fa"
    fasta_path.write_text(
        ">tiny\n" + ("A" * 300) + "\n>chr\n" + ("C" * 30000) + "\n", encoding="utf-8"
    )
    bam_path = tmp_path / "junc.bam"
    reads = []
    # Tiny contig (300 bp, cyclic buffer 300 bp -> 600 bp): reads spanning 50..550 cover 250 bp each side (>150)
    for i in range(10):
        reads.append(
            {
                "name": f"tiny_span_{i}",
                "ref_id": 0,
                "start": 50,
                "cigar": [(7, 500)],
                "seq": "A" * 500,
            }
        )
    # Chr (30,000 bp, junction=30000, min_overlap=1000):
    # - 10 spanning reads (28500..31500)
    # - 5 reads ending/clipped 200 bp from seam (29000..29800 + 100S) -> MUST be broken_reads!
    # - 5 reads with 80 bp leading adapter clip starting at 29200..30200 (800 bp from seam, right_extension=200 < 1000? wait: left_extension=800 < 1000, abs(29200-30000)=800 > clip_window=500 -> neutral!)
    for i in range(10):
        reads.append(
            {
                "name": f"chr_span_{i}",
                "ref_id": 1,
                "start": 28500,
                "cigar": [(7, 3000)],
                "seq": "C" * 3000,
            }
        )
    for i in range(5):
        reads.append(
            {
                "name": f"chr_broken200_{i}",
                "ref_id": 1,
                "start": 29000,
                "cigar": [(7, 800), (4, 100)],
                "seq": "C" * 900,
            }
        )
    for i in range(5):
        reads.append(
            {
                "name": f"chr_adapter800_{i}",
                "ref_id": 1,
                "start": 29200,
                "cigar": [(4, 80), (7, 1000)],  # 29200..30200, overlaps junction, adapter clip 800 bp away
                "seq": "C" * 1080,
            }
        )
    _write_bam(bam_path, [("tiny_cyclic", 600), ("chr_cyclic", 35000)], reads)

    ev = JunctionEvaluator(str(bam_path), {"tiny": 300, "chr": 30000})
    res = ev.evaluate_junctions()
    assert res["tiny"]["status"] == "PASS"
    assert res["tiny"]["spanning_reads"] == 10
    assert res["chr"]["spanning_reads"] == 10
    assert res["chr"]["broken_reads"] == 5
    assert res["chr"]["neutral_reads"] == 5


def test_dnaa_seam_wrapping_and_linear_bioval(tmp_path):
    """Bug 11 & 3.3: DnaA box spanning coordinate 0/L seam is detected in circular mode, not in linear."""
    # Place TTATCCACA across the seam: last 4 bp 'TTAT' + first 5 bp 'CCACA'
    seq = "CCACA" + ("A" * 599991) + "TTAT"
    assert len(seq) == 600000
    circ_hits = _find_dnaa_box_positions(seq, max_mismatches=1, circular=True)
    lin_hits = _find_dnaa_box_positions(seq, max_mismatches=1, circular=False)
    assert 599996 in circ_hits
    assert 599996 not in lin_hits

    # Linear BioVal symmetry when oriC is central (minimum cumulative skew at 300 kb)
    lin_seq = ("C" * 300000) + ("G" * 300000)
    fa = tmp_path / "lin.fa"
    fa.write_text(f">lin_chr\n{lin_seq}\n", encoding="utf-8")
    bv = BioVal(str(fa), linear_mode=True)
    res = bv.calculate_symmetry()
    assert res["lin_chr"]["rotation_status"] == "LINEAR"
    assert res["lin_chr"]["viability"] == "ACCEPTABLE"
    assert res["lin_chr"]["symmetry"] >= 175.0
