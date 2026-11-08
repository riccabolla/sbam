from __future__ import annotations

import argparse
import json
import os
import sys

from sbam.base.base import BaseAccuracy
from sbam.report.report import DashboardBuilder
from sbam.struct.eval import JunctionEvaluator
from sbam.struct.mapping import ReadAligner
from sbam.struct.ori import BioVal
from sbam.struct.struct import GenomeStructure


def parse_size(size_str: str | None) -> int | None:
    if not size_str:
        return None
    s = size_str.lower().strip().replace(",", "").replace("_", "")
    # Strip trailing 'bp' or 'b' if preceded by another unit or digit
    if s.endswith("bp"):
        s = s[:-2].strip()

    multipliers = [
        ("gb", 1_000_000_000),
        ("g", 1_000_000_000),
        ("mb", 1_000_000),
        ("m", 1_000_000),
        ("kb", 1_000),
        ("k", 1_000),
        ("b", 1),
    ]

    try:
        for suffix, mult in multipliers:
            if s.endswith(suffix):
                num_part = s[: -len(suffix)].strip()
                value = int(float(num_part) * mult)
                if value <= 0:
                    raise ValueError
                return value
        value = int(float(s))
        if value <= 0:
            raise ValueError
        return value
    except ValueError:
        raise argparse.ArgumentTypeError(
            f"Invalid expected genome size format: '{size_str}' (examples: '4.8m', '4800kb', '4800000')"
        )


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        prog="sbam",
        description="SBAM: Structural and Base-Level analysis of long-read microbial assemblies",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    # Required inputs
    parser.add_argument(
        "-a",
        "--assembly",
        required=True,
        help="Path to the consensus assembly (FASTA)",
    )
    parser.add_argument(
        "-r",
        "--reads",
        required=True,
        help="Path to the raw long reads (FASTQ/FASTA)",
    )
    parser.add_argument(
        "-o",
        "--outdir",
        required=True,
        help="Directory to save the BAM files and HTML report",
    )

    # Optional parameters
    parser.add_argument(
        "-t",
        "--threads",
        type=int,
        default=4,
        help="Number of CPU threads for minimap2 and base pileup",
    )
    parser.add_argument(
        "-g",
        "--genome_size",
        "--genome-size",
        dest="genome_size",
        type=parse_size,
        default=None,
        help="Expected genome size (e.g., 4.8m, 5000kb, 4800000)",
    )
    parser.add_argument(
        "--buffer-size",
        type=int,
        default=50000,
        help="Size of the cyclic buffer edge (bp) for mapping",
    )
    parser.add_argument(
        "--read-type",
        choices=["map-ont", "lr:hq", "map-pb", "map-hifi"],
        default="map-ont",
        help="minimap2 preset for sequencing technology",
    )
    parser.add_argument(
        "--min-concordance",
        type=float,
        default=0.75,
        help="Minimum read-to-assembly base concordance threshold (0.0 - 1.0)",
    )
    parser.add_argument(
        "--depth-variance",
        type=float,
        default=2.5,
        help="Fold-change limit from median depth before masking a region",
    )

    args = parser.parse_args(argv)

    if args.threads < 1:
        parser.error("--threads must be >= 1")
    if args.buffer_size < 1000:
        parser.error("--buffer-size must be >= 1000 bp")
    if not (0.0 < args.min_concordance <= 1.0):
        parser.error("--min-concordance must be between 0.0 and 1.0")

    return args


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)

    # Fail fast if input files are missing
    if not os.path.isfile(args.assembly):
        print(f"[ERROR] Assembly file not found: {args.assembly}", file=sys.stderr)
        return 1
    if not os.path.isfile(args.reads):
        print(f"[ERROR] Reads file not found: {args.reads}", file=sys.stderr)
        return 1

    os.makedirs(args.outdir, exist_ok=True)
    print("[INIT] Starting SBAM analysis...")
    print(f"[INIT] Assembly: {args.assembly}")
    print(f"[INIT] Output directory: {args.outdir}\n")

    # 1. Structural validation step
    print("Structural Validation")
    genome = GenomeStructure(args.assembly, args.outdir)
    if not genome.records:
        print(
            f"[ERROR] No valid FASTA records found in {args.assembly}.",
            file=sys.stderr,
        )
        return 1

    print(f" > Found {len(genome.records)} contig(s).")

    cyclic_fasta, orig_lengths = genome.create_cyclic_buffer(
        buffer_size=args.buffer_size
    )
    print(f" > Cyclic buffer applied. Saved to: {cyclic_fasta}")

    aligner = ReadAligner(args.reads, cyclic_fasta, args.outdir, args.threads)
    bam_path = aligner.run_mapping(read_type=args.read_type)

    evaluator = JunctionEvaluator(bam_path, orig_lengths)
    junction_metrics = evaluator.evaluate_junctions()

    print("\n Junction Structural Scores")
    for contig, metrics in junction_metrics.items():
        length_kb = metrics["length"] / 1000
        depth = metrics["avg_depth"]
        print(
            f"   - {contig} ({length_kb:.1f} kb | Depth: {depth}x): "
            f"Score {metrics['spanning_score']:.2f} [{metrics['status']}] "
            f"({metrics['spanning_reads']} spanning / {metrics['broken_reads']} broken)"
        )

    # 2. Replication architecture validation
    physics = BioVal(args.assembly)
    physics_metrics = physics.calculate_symmetry()

    print("\n Biological Scores")
    if not physics_metrics:
        print("   - No contigs >= 100 kb available for GC skew symmetry analysis.")
    else:
        for contig, metrics in physics_metrics.items():
            print(
                f"   - {contig}: Symmetry {metrics['symmetry']}° [{metrics['viability']}] "
                f"(oriC ≈ {metrics['oric_pos']} bp, ter ≈ {metrics['ter_pos']} bp)"
            )

    # 3. Base-level validation
    print("\n Base-Level Validation")
    fidelity = BaseAccuracy(
        bam_path=bam_path,
        fasta_path=args.assembly,
        threads=args.threads,
        min_concordance=args.min_concordance,
        depth_variance_limit=args.depth_variance,
    )

    motif_results = fidelity.analyze_motifs()

    print("\n > Top Systematic Error Motifs (by k-mer size)")
    if not motif_results:
        print(
            "   No significant systematic error motifs detected. Basecalling is highly concordant."
        )
    else:
        for k in fidelity.k_list:
            print(f"\n   [{k}-mers]")
            k_results = [res for res in motif_results if res["kmer_size"] == k]
            if not k_results:
                print("     No significant motifs detected.")
            else:
                for res in k_results[:3]:
                    print(
                        f"     {res['motif']} | Discordance Rate: {res['error_rate']:.1%} "
                        f"({res['discordant_count']}/{res['occurrences']} hits) | "
                        f"Fold Increase: {res['fold_increase']:.1f}x"
                    )

    # 4. Reporting
    print("\n Generating Dashboard")
    builder = DashboardBuilder(args.outdir, args.assembly)
    masked_pct = getattr(fidelity, "masked_pct", 0.0)
    masked_bases = getattr(fidelity, "masked_bases", 0)
    masked_regions = getattr(fidelity, "masked_regions", [])

    report_file = builder.generate_report(
        junction_metrics=junction_metrics,
        physics_metrics=physics_metrics,
        motif_results=motif_results,
        masked_bases=masked_bases,
        masked_pct=masked_pct,
        masked_regions=masked_regions,
        genome_size=args.genome_size,
    )

    # Save machine-readable JSON summary (excluding large plot arrays)
    summary_json_path = os.path.join(
        args.outdir, f"{os.path.basename(args.assembly)}_sbam_metrics.json"
    )
    compact_physics = {
        cid: {k: v for k, v in m.items() if k not in ("plot_theta", "plot_r")}
        for cid, m in physics_metrics.items()
    }
    with open(summary_json_path, "w") as jf:
        json.dump(
            {
                "assembly": os.path.basename(args.assembly),
                "expected_genome_size": args.genome_size,
                "junction_metrics": junction_metrics,
                "replication_metrics": compact_physics,
                "base_fidelity": {
                    "masked_bases": masked_bases,
                    "masked_pct": round(masked_pct, 4),
                    "masked_regions": masked_regions,
                    "systematic_motifs": motif_results,
                },
            },
            jf,
            indent=2,
        )

    print("\n[DONE] SBAM execution finished successfully.")
    print(f"View your HTML report at: {report_file}")
    print(f"View your JSON metrics at: {summary_json_path}")
    return 0