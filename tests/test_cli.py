import json
import pytest

from sbam import __version__
from sbam.cli import parse_run_args, run_summarize


def test_cli_version_flag(capsys):
    with pytest.raises(SystemExit) as exc:
        parse_run_args(["--version"])
    assert exc.value.code == 0
    captured = capsys.readouterr()
    assert __version__ in captured.out


def test_cli_input_validation_errors():
    base = ["-a", "asm.fa", "-r", "reads.fq", "-o", "out"]

    with pytest.raises(SystemExit):
        parse_run_args(base + ["--depth-variance", "0"])
    with pytest.raises(SystemExit):
        parse_run_args(base + ["--depth-variance", "1.0"])
    with pytest.raises(SystemExit):
        parse_run_args(base + ["--min-concordance", "0.50"])
    with pytest.raises(SystemExit):
        parse_run_args(base + ["--min-concordance", "1.2"])
    with pytest.raises(SystemExit):
        parse_run_args(base + ["--min-ori-length", "0"])

    # Mutually exclusive -r and --bam
    with pytest.raises(SystemExit):
        parse_run_args(["-a", "asm.fa", "-r", "reads.fq", "--bam", "aln.bam", "-o", "out"])


def test_cli_bam_and_prefix_and_no_emit_tracks():
    args = parse_run_args(
        [
            "-a", "asm.fa",
            "--bam", "aln.bam",
            "-o", "out",
            "--prefix", "sample_A",
            "--no-emit-tracks",
            "--linear",
        ]
    )
    assert args.bam == "aln.bam"
    assert args.reads is None
    assert args.prefix == "sample_A"
    assert args.emit_tracks is False
    assert args.linear is True


def test_summarize_respects_sample_name_and_exact_counts(tmp_path):
    j1 = tmp_path / "s1_sbam_metrics.json"
    j2 = tmp_path / "s2_sbam_metrics.json"

    payload_base = {
        "schema_version": 2,
        "sample_name": "iso_01",
        "assembly": "assembly.fasta",
        "verdict": "ACCEPT",
        "completeness_pct": 100.0,
        "junction_metrics": {
            "chr1": {"length": 1_000_000, "avg_depth": 50.0, "spanning_score": 0.99, "status": "PASS"}
        },
        "replication_metrics": {
            "chr1": {"symmetry": 178.0, "oric_offset_kb": 5.0, "rotation_status": "ROTATED", "dnaa_fold_enrichment": 3.0, "dnaa_poisson_pval": 0.001, "viability": "ACCEPTABLE"}
        },
        "base_fidelity": {
            "masked_bases": 100,
            "masked_pct": 0.01,
            "evaluated_bases": 999900,
            "discordant_bases": 25,
            "consensus_error_bases": 2,
            "subconsensus_bias_bases": 23,
            "low_concordance_blocks": [{"start": 0, "end": 60, "length": 60, "mean_concordance": 80.0}],
            "masked_regions": [],
            "systematic_motifs": [
                {
                    "motif": "CCAGG/CCTGG",
                    "assigned_unique_sites": 15,
                    "consensus_error_count": 2,
                    "fold_increase": 10.0,
                    "artifact_class": "Consensus Error (2 site(s) <50%)",
                },
                {
                    "motif": "ACCAGG/CCTGGT",
                    "assigned_unique_sites": 5,
                    "consensus_error_count": 2,  # Legacy/overlapping sum must NOT inflate consensus_error_bases
                    "fold_increase": 8.0,
                    "artifact_class": "Consensus Error (2 site(s) <50%)",
                },
            ],
        },
    }
    j1.write_text(json.dumps(payload_base), encoding="utf-8")

    payload2 = json.loads(json.dumps(payload_base))
    payload2["sample_name"] = "iso_02"
    j2.write_text(json.dumps(payload2), encoding="utf-8")

    out_tsv = tmp_path / "cohort.tsv"
    rc = run_summarize([str(j1), str(j2), "-o", str(out_tsv)])
    assert rc == 0

    lines = out_tsv.read_text(encoding="utf-8").strip().splitlines()
    assert len(lines) == 3
    header = lines[0].split("\t")
    row1 = lines[1].split("\t")
    row_dict = dict(zip(header, row1))
    assert row_dict["sample"] == "iso_01"
    assert row_dict["consensus_error_bases_bp"] == "2"
    assert row_dict["low_concordance_blocks_count"] == "1"

    motifs_tsv = tmp_path / "cohort_recurrent_motifs.tsv"
    assert motifs_tsv.exists()
    assert "Recurrent Batch/Basecaller Signature" in motifs_tsv.read_text(encoding="utf-8")
