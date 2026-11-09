from __future__ import annotations

import argparse
import csv
import glob
import json
import os
import sys
from contextlib import contextmanager
from datetime import datetime

from sbam.base.base import BaseAccuracy
from sbam.report.report import DashboardBuilder
from sbam.struct.eval import JunctionEvaluator
from sbam.struct.mapping import ReadAligner
from sbam.struct.ori import BioVal
from sbam.struct.struct import GenomeStructure


class _TeeStream:
    """Mirrors stdout/stderr writes to both the terminal and a unified log file."""

    def __init__(self, ConsoleStream, log_file):
        self.console = ConsoleStream
        self.log_file = log_file

    def write(self, message):
        self.console.write(message)
        self.log_file.write(message)

    def flush(self):
        self.console.flush()
        self.log_file.flush()

    def isatty(self):
        return getattr(self.console, "isatty", lambda: False)()


@contextmanager
def unified_logger(log_path: str | None):
    if not log_path:
        yield None
        return

    log_dir = os.path.dirname(os.path.abspath(log_path))
    if log_dir:
        os.makedirs(log_dir, exist_ok=True)

    with open(log_path, "w", encoding="utf-8") as lf:
        lf.write(f"# SBAM Unified Execution Log — Started {datetime.now().isoformat()}\n")
        lf.flush()
        orig_out, orig_err = sys.stdout, sys.stderr
        sys.stdout = _TeeStream(orig_out, lf)
        sys.stderr = _TeeStream(orig_err, lf)
        try:
            yield log_path
        finally:
            lf.write(f"\n# SBAM Execution Ended {datetime.now().isoformat()}\n")
            lf.flush()
            sys.stdout = orig_out
            sys.stderr = orig_err


def parse_size(size_str: str | None) -> int | None:
    if not size_str:
        return None
    s = size_str.lower().strip().replace(",", "").replace("_", "")
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


def parse_run_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        prog="sbam",
        description="SBAM: Structural and Base-Level analysis of long-read microbial assemblies "
        "(use 'sbam summarize *.json -o cohort_summary.tsv' for multi-sample cohort aggregation)",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

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
        help="Directory to save outputs and the HTML report",
    )

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
    parser.add_argument(
        "--min-ori-length",
        type=int,
        default=500000,
        help="Minimum contig length (bp) for chromosomal GC-skew replication analysis",
    )
    parser.add_argument(
        "--clean",
        "--remove-intermediates",
        dest="clean",
        action="store_true",
        help="Remove intermediate files (cyclic FASTA, BAM, and BAI) upon completion",
    )
    parser.add_argument(
        "--log",
        "--log-file",
        dest="log_file",
        nargs="?",
        const="sbam.log",
        default=None,
        help="Save a complete unified SBAM log file (defaults to <outdir>/sbam.log if flag is passed without a path)",
    )

    args = parser.parse_args(argv)

    if args.threads < 1:
        parser.error("--threads must be >= 1")
    if args.buffer_size < 1000:
        parser.error("--buffer-size must be >= 1000 bp")
    if not (0.0 < args.min_concordance <= 1.0):
        parser.error("--min-concordance must be between 0.0 and 1.0")

    return args


def parse_summarize_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        prog="sbam summarize",
        description="Aggregate multiple SBAM JSON metric files into a cohort TSV/CSV summary table",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "inputs",
        nargs="+",
        help="Input *_sbam_metrics.json files or glob patterns (e.g., results/*/*_sbam_metrics.json)",
    )
    parser.add_argument(
        "-o",
        "--output",
        default="cohort_summary.tsv",
        help="Output summary table path (.tsv or .csv)",
    )
    return parser.parse_args(argv)


def run_summarize(argv: list[str]) -> int:
    args = parse_summarize_args(argv)

    json_files: list[str] = []
    for pattern in args.inputs:
        expanded = glob.glob(pattern)
        if expanded:
            json_files.extend(expanded)
        elif os.path.isfile(pattern):
            json_files.append(pattern)

    # Deduplicate while preserving order
    seen = set()
    unique_files = []
    for f in json_files:
        abs_f = os.path.abspath(f)
        if abs_f not in seen:
            seen.add(abs_f)
            unique_files.append(f)

    if not unique_files:
        print("[ERROR] No valid SBAM JSON files found to summarize.", file=sys.stderr)
        return 1

    fieldnames = [
        "sample",
        "verdict",
        "num_contigs",
        "total_length_bp",
        "completeness_pct",
        "primary_contig_id",
        "primary_length_bp",
        "primary_depth_x",
        "primary_junction_score",
        "primary_circularity_status",
        "primary_symmetry_deg",
        "primary_oric_offset_kb",
        "primary_replication_status",
        "circular_plasmids_count",
        "failed_contigs_count",
        "masked_bases_bp",
        "masked_pct",
        "systematic_error_motifs_count",
        "top_error_motifs",
    ]

    rows = []
    for jpath in unique_files:
        try:
            with open(jpath, "r", encoding="utf-8") as jf:
                data = json.load(jf)
        except Exception as e:
            print(f"[WARNING] Skipping unreadable JSON {jpath}: {e}", file=sys.stderr)
            continue

        junction = data.get("junction_metrics", {})
        replication = data.get("replication_metrics", {})
        fidelity = data.get("base_fidelity", {})
        motifs = fidelity.get("systematic_motifs", [])

        total_length = sum(m.get("length", 0) for m in junction.values())
        primary_id = (
            max(junction.keys(), key=lambda k: junction[k].get("length", 0))
            if junction
            else "N/A"
        )
        p_junc = junction.get(primary_id, {})
        p_rep = replication.get(primary_id, {})

        circular_plasmids = sum(
            1
            for cid, m in junction.items()
            if cid != primary_id and m.get("status") == "PASS"
        )
        failed_contigs = sum(
            1 for m in junction.values() if m.get("status") == "FAIL"
        )

        top_motifs_str = (
            ",".join(
                f"{m['motif']}({m['fold_increase']:.1f}x)" for m in motifs[:5]
            )
            if motifs
            else "None"
        )

        rows.append(
            {
                "sample": data.get("assembly", os.path.basename(jpath)),
                "verdict": data.get("verdict", "UNKNOWN"),
                "num_contigs": len(junction),
                "total_length_bp": total_length,
                "completeness_pct": (
                    data.get("completeness_pct")
                    if data.get("completeness_pct") is not None
                    else "N/A"
                ),
                "primary_contig_id": primary_id,
                "primary_length_bp": p_junc.get("length", 0),
                "primary_depth_x": p_junc.get("avg_depth", 0.0),
                "primary_junction_score": p_junc.get("spanning_score", 0.0),
                "primary_circularity_status": p_junc.get("status", "UNKNOWN"),
                "primary_symmetry_deg": p_rep.get("symmetry", "N/A"),
                "primary_oric_offset_kb": p_rep.get("oric_offset_kb", "N/A"),
                "primary_replication_status": p_rep.get("viability", "N/A"),
                "circular_plasmids_count": circular_plasmids,
                "failed_contigs_count": failed_contigs,
                "masked_bases_bp": fidelity.get("masked_bases", 0),
                "masked_pct": fidelity.get("masked_pct", 0.0),
                "systematic_error_motifs_count": len(motifs),
                "top_error_motifs": top_motifs_str,
            }
        )

    delimiter = "," if args.output.lower().endswith(".csv") else "\t"
    out_dir = os.path.dirname(os.path.abspath(args.output))
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    with open(args.output, "w", newline="", encoding="utf-8") as out_f:
        writer = csv.DictWriter(out_f, fieldnames=fieldnames, delimiter=delimiter)
        writer.writeheader()
        writer.writerows(rows)

    print(
        f"[DONE] Summarized {len(rows)} assembly report(s) to: {args.output}"
    )
    return 0


def main(argv: list[str] | None = None) -> int:
    effective_argv = sys.argv[1:] if argv is None else argv

    # Route 'sbam summarize ...' to the cohort summarizer
    if effective_argv and effective_argv[0] == "summarize":
        return run_summarize(effective_argv[1:])
    if effective_argv and effective_argv[0] == "run":
        effective_argv = effective_argv[1:]

    args = parse_run_args(effective_argv)

    if not os.path.isfile(args.assembly):
        print(f"[ERROR] Assembly file not found: {args.assembly}", file=sys.stderr)
        return 1
    if not os.path.isfile(args.reads):
        print(f"[ERROR] Reads file not found: {args.reads}", file=sys.stderr)
        return 1

    os.makedirs(args.outdir, exist_ok=True)

    resolved_log_path = None
    if args.log_file:
        resolved_log_path = (
            os.path.join(args.outdir, args.log_file)
            if not os.path.isabs(args.log_file)
            and os.path.dirname(args.log_file) == ""
            else args.log_file
        )

    with unified_logger(resolved_log_path):
        print("[INIT] Starting SBAM analysis...")
        print(f"[INIT] Assembly: {args.assembly}")
        print(f"[INIT] Output directory: {args.outdir}")
        if resolved_log_path:
            print(f"[INIT] Unified log file: {resolved_log_path}")
        print()

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

        aligner = ReadAligner(
            args.reads,
            cyclic_fasta,
            args.outdir,
            threads=args.threads,
            log_path=resolved_log_path,
        )
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
                f"({metrics['spanning_reads']} spanning / {metrics['broken_reads']} clipped, "
                f"{metrics.get('neutral_reads', 0)} natural ends ignored)"
            )

        # 2. Replication architecture validation
        physics = BioVal(args.assembly, min_contig_len=args.min_ori_length)
        physics_metrics = physics.calculate_symmetry()

        print("\n Biological Scores")
        if not physics_metrics:
            print(
                f"   - No chromosomal contigs >= {args.min_ori_length // 1000} kb available for GC skew symmetry analysis."
            )
        else:
            for contig, metrics in physics_metrics.items():
                print(
                    f"   - {contig}: Symmetry {metrics['symmetry']}° [{metrics['viability']}] "
                    f"(oriC ≈ {metrics['oric_pos']} bp, ter ≈ {metrics['ter_pos']} bp, "
                    f"oriC offset to 0 bp ≈ {metrics['oric_offset_kb']} kb / {metrics['oric_offset_deg']}°)"
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

        summary_json_path = os.path.join(
            args.outdir, f"{os.path.basename(args.assembly)}_sbam_metrics.json"
        )
        compact_physics = {
            cid: {k: v for k, v in m.items() if k not in ("plot_theta", "plot_r")}
            for cid, m in physics_metrics.items()
        }
        with open(summary_json_path, "w", encoding="utf-8") as jf:
            json.dump(
                {
                    "assembly": os.path.basename(args.assembly),
                    "verdict": builder.last_verdict,
                    "expected_genome_size": args.genome_size,
                    "completeness_pct": builder.last_completeness_pct,
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

        # 5. Optional cleanup of intermediate files
        if args.clean:
            print("\n Cleaning Intermediate Files")
            intermediates = [
                cyclic_fasta,
                bam_path,
                f"{bam_path}.bai",
            ]
            for fpath in intermediates:
                if fpath and os.path.exists(fpath):
                    try:
                        os.remove(fpath)
                        print(f" > Removed: {fpath}")
                    except OSError as e:
                        print(
                            f"[WARNING] Failed to remove {fpath}: {e}",
                            file=sys.stderr,
                        )

        print("\n[DONE] SBAM execution finished successfully.")
        print(f"View your HTML report at: {report_file}")
        print(f"View your JSON metrics at: {summary_json_path}")
        if resolved_log_path:
            print(f"View your unified log at: {resolved_log_path}")

    return 0