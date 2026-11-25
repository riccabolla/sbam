"""Command-line interface and pipeline orchestration for SBAM."""

from __future__ import annotations

import argparse
import csv
import glob
import json
import logging
import os
import subprocess
import sys
from collections import defaultdict
from contextlib import contextmanager
from datetime import datetime
from typing import Any

import pysam

from sbam import __version__
from sbam.base.base import BaseAccuracy
from sbam.models import SCHEMA_VERSION
from sbam.report.report import DashboardBuilder
from sbam.report.tracks import (
    export_consensus_errors_tsv,
    export_consensus_errors_vcf,
    export_masked_regions_bed,
)
from sbam.struct.eval import JunctionEvaluator
from sbam.struct.mapping import AlignmentError, ReadAligner
from sbam.struct.ori import BioVal
from sbam.struct.struct import GenomeStructure
from sbam.utils import json_serializable, load_assembly


@contextmanager
def configure_logging(log_path: str | None):
    """Configure the `sbam` logger for console output and optional unified log file."""
    sbam_logger = logging.getLogger("sbam")
    prev_level = sbam_logger.level
    prev_propagate = sbam_logger.propagate
    prev_handlers = list(sbam_logger.handlers)

    for handler in prev_handlers:
        sbam_logger.removeHandler(handler)

    sbam_logger.setLevel(logging.DEBUG)
    sbam_logger.propagate = False

    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setLevel(logging.INFO)
    console_handler.setFormatter(logging.Formatter("%(message)s"))
    sbam_logger.addHandler(console_handler)

    file_handler = None
    if log_path:
        log_dir = os.path.dirname(os.path.abspath(log_path))
        if log_dir:
            os.makedirs(log_dir, exist_ok=True)
        file_handler = logging.FileHandler(log_path, mode="w", encoding="utf-8")
        file_handler.setLevel(logging.DEBUG)
        file_handler.setFormatter(logging.Formatter("%(message)s"))
        sbam_logger.addHandler(file_handler)
        sbam_logger.info("# SBAM Unified Execution Log — Started %s", datetime.now().isoformat())

    try:
        yield log_path
    finally:
        if file_handler is not None:
            sbam_logger.info("\n# SBAM Execution Ended %s", datetime.now().isoformat())
        for handler in list(sbam_logger.handlers):
            sbam_logger.removeHandler(handler)
            handler.close()
        for handler in prev_handlers:
            sbam_logger.addHandler(handler)
        sbam_logger.setLevel(prev_level)
        sbam_logger.propagate = prev_propagate


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
            f"Invalid expected genome size format: '{size_str}' "
            "(examples: '4.8m', '4800kb', '4800000')"
        )


def parse_k_list(k_str: str) -> list[int]:
    try:
        vals = [int(x.strip()) for x in k_str.split(",") if x.strip()]
        if not vals or any(k < 2 or k > 10 for k in vals):
            raise ValueError
        return sorted(set(vals))
    except ValueError:
        raise argparse.ArgumentTypeError(
            f"Invalid --k-list '{k_str}'. Provide comma-separated integers between 2 and 10 (e.g., '4,5,6')."
        )


def parse_run_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        prog="sbam",
        description="SBAM: Structural and Base-Level analysis of long-read microbial assemblies",
        epilog=(
            "Subcommands:\n"
            "  sbam run [options]            Run full assembly QC (default when no subcommand is given)\n"
            "  sbam summarize <json...>      Aggregate *_sbam_metrics.json files into a cohort TSV/CSV\n"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "-V", "--version", action="version", version=f"%(prog)s {__version__}"
    )

    parser.add_argument(
        "-a", "--assembly", required=True, help="Path to the consensus assembly (FASTA)"
    )
    input_group = parser.add_mutually_exclusive_group(required=True)
    input_group.add_argument(
        "-r", "--reads", default=None, help="Path to the raw long reads (FASTQ/FASTA)"
    )
    input_group.add_argument(
        "--bam",
        default=None,
        help="Path to a pre-aligned coordinate-sorted BAM (auto-indexed if .bai is missing)",
    )

    parser.add_argument(
        "-o", "--outdir", required=True, help="Directory to save outputs and the HTML report"
    )
    parser.add_argument(
        "--prefix",
        "--sample-name",
        dest="prefix",
        default=None,
        help="Sample prefix for output filenames and report title (defaults to assembly filename)",
    )

    parser.add_argument("-t", "--threads", type=int, default=4, help="Number of CPU threads (default: 4)")
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
        help="Size of the cyclic buffer edge in bp (default: 50000)",
    )
    parser.add_argument(
        "--read-type",
        choices=["map-ont", "lr:hq", "map-pb", "map-hifi"],
        default="map-ont",
        help="minimap2 preset for sequencing technology (default: map-ont)",
    )
    parser.add_argument(
        "--min-concordance",
        type=float,
        default=0.85,
        help="Minimum read-to-assembly concordance before flagging a site (default: 0.85; >0.50 and <=1.0)",
    )
    parser.add_argument(
        "--depth-variance",
        type=float,
        default=2.5,
        help="Fold-change limit from median depth before structurally masking (default: 2.5; >1.0)",
    )
    parser.add_argument(
        "--k-list",
        type=parse_k_list,
        default=[4, 5, 6],
        help="Comma-separated k-mer sizes to test for systematic motif discordance (default: 4,5,6)",
    )
    parser.add_argument(
        "--min-ori-length",
        type=int,
        default=500000,
        help="Minimum contig length (bp) for chromosomal GC-skew replication analysis (default: 500000)",
    )
    parser.add_argument(
        "--linear",
        action="store_true",
        help="Treat replicons as naturally linear (e.g., Streptomyces, Borrelia) rather than circular",
    )
    parser.add_argument(
        "--offline",
        action="store_true",
        help="Render HTML report with self-contained inline SVG plots (no external Plotly CDN)",
    )
    parser.add_argument(
        "--emit-tracks",
        dest="emit_tracks",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Export BED/TSV/VCF diagnostic tracks for IGV and polishing (default: True)",
    )
    parser.add_argument(
        "--clean",
        "--remove-intermediates",
        dest="clean",
        action="store_true",
        help="Remove intermediate alignment files created by SBAM upon completion",
    )
    parser.add_argument(
        "--log",
        "--log-file",
        dest="log_file",
        nargs="?",
        const="sbam.log",
        default=None,
        help="Save a complete unified SBAM log file (defaults to <outdir>/sbam.log)",
    )

    args = parser.parse_args(argv)
    if args.threads < 1:
        parser.error("--threads must be >= 1")
    if args.buffer_size < 1000:
        parser.error("--buffer-size must be >= 1000 bp")
    if args.depth_variance <= 1.0:
        parser.error("--depth-variance must be > 1.0")
    if not (0.50 < args.min_concordance <= 1.0):
        parser.error("--min-concordance must be > 0.50 and <= 1.0")
    if args.min_ori_length <= 0:
        parser.error("--min-ori-length must be > 0")
    return args


def parse_summarize_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        prog="sbam summarize",
        description="Aggregate multiple SBAM JSON metric files into a cohort TSV/CSV summary table and recurrent motif report",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("inputs", nargs="+", help="Input *_sbam_metrics.json files or glob patterns")
    parser.add_argument("-o", "--output", default="cohort_summary.tsv", help="Output summary table path (.tsv or .csv)")
    return parser.parse_args(argv)


def _ensure_bam_index_and_detect_cyclic(
    bam_path: str, contig_ids: list[str], threads: int
) -> bool:
    """Ensure a BAM index exists and auto-detect whether reference names carry `_cyclic`."""
    stem, _ = os.path.splitext(bam_path)
    has_index = any(
        os.path.isfile(cand)
        for cand in (f"{bam_path}.bai", f"{bam_path}.csi", f"{stem}.bai", f"{stem}.csi")
    )
    if not has_index:
        logging.getLogger("sbam").info(" > [samtools] Indexing pre-aligned BAM: %s", bam_path)
        proc = subprocess.run(
            ["samtools", "index", "-@", str(threads), bam_path],
            capture_output=True,
            text=True,
            check=False,
        )
        if proc.returncode != 0:
            raise RuntimeError(
                f"Failed to index BAM '{bam_path}': {proc.stderr.strip()}"
            )

    with pysam.AlignmentFile(bam_path, "rb") as bam:
        refs = set(bam.references)

    if any(f"{cid}_cyclic" in refs for cid in contig_ids):
        return True
    if any(cid in refs for cid in contig_ids):
        return False
    raise ValueError(
        f"None of the assembly contigs ({', '.join(contig_ids[:5])}) match references in BAM '{bam_path}'."
    )


def run_summarize(argv: list[str]) -> int:
    args = parse_summarize_args(argv)

    json_files: list[str] = []
    for pattern in args.inputs:
        expanded = glob.glob(pattern)
        if expanded:
            json_files.extend(expanded)
        elif os.path.isfile(pattern):
            json_files.append(pattern)

    seen: set[str] = set()
    unique_files: list[str] = []
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
        "primary_rotation_status",
        "primary_dnaa_fold_enrichment",
        "primary_dnaa_poisson_pval",
        "primary_replication_status",
        "circular_plasmids_count",
        "failed_contigs_count",
        "masked_bases_bp",
        "masked_pct",
        "low_concordance_blocks_count",
        "evaluated_bases_bp",
        "consensus_error_bases_bp",
        "subconsensus_bias_bases_bp",
        "motif_families_count",
        "explained_discordant_pct",
        "top_motif_families",
    ]

    rows: list[dict[str, Any]] = []
    motif_cohort_tracker: dict[str, dict[str, Any]] = defaultdict(
        lambda: {
            "samples": [],
            "total_unique_sites": 0,
            "total_consensus_errors": 0,
            "folds": [],
            "classes": set(),
        }
    )

    for jpath in unique_files:
        try:
            with open(jpath, "r", encoding="utf-8") as jf:
                data = json.load(jf)
        except (OSError, ValueError) as e:
            print(f"[WARNING] Skipping unreadable JSON {jpath}: {e}", file=sys.stderr)
            continue

        sample_name = data.get("sample_name") or data.get("assembly") or os.path.basename(jpath)
        junction = data.get("junction_metrics", {})
        replication = data.get("replication_metrics", {})
        fidelity = data.get("base_fidelity", {})
        motifs = fidelity.get("systematic_motifs", [])
        masked_regions = fidelity.get("masked_regions", [])

        if "low_concordance_blocks" in fidelity and fidelity["low_concordance_blocks"] is not None:
            low_concord_blocks = fidelity["low_concordance_blocks"]
        else:
            low_concord_blocks = [
                r for r in masked_regions if r.get("mean_concordance", 100.0) < 90.0
            ]

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
            if cid != primary_id and m.get("status") in ("PASS", "LINEAR")
        )
        failed_contigs = sum(1 for m in junction.values() if m.get("status") == "FAIL")

        for m in motifs:
            m_name = m["motif"]
            entry = motif_cohort_tracker[m_name]
            entry["samples"].append(sample_name)
            entry["total_unique_sites"] += m.get(
                "assigned_unique_sites",
                m.get("family_unique_sites", m.get("discordant_count", 0)),
            )
            entry["total_consensus_errors"] += m.get("consensus_error_count", 0)
            entry["folds"].append(m.get("fold_increase", 0.0))
            if m.get("artifact_class"):
                entry["classes"].add(m["artifact_class"])

        top_motifs_str = (
            ",".join(
                f"{m['motif']}({m.get('assigned_unique_sites', m.get('family_unique_sites', m.get('discordant_count', 0)))}bp)"
                for m in motifs[:5]
            )
            if motifs
            else "-"
        )

        disc_bp = fidelity.get("discordant_bases", 0)
        if fidelity.get("consensus_error_bases") is not None:
            cons_err_bp = fidelity["consensus_error_bases"]
        else:
            cons_err_bp = sum(m.get("consensus_error_count", 0) for m in motifs)

        sub_bias_bp = fidelity.get(
            "subconsensus_bias_bases", max(0, disc_bp - cons_err_bp)
        )
        expl_bp = fidelity.get(
            "explained_discordant_bases",
            min(
                disc_bp,
                sum(
                    m.get(
                        "assigned_unique_sites",
                        m.get("family_unique_sites", m.get("discordant_count", 0)),
                    )
                    for m in motifs
                ),
            ),
        )
        expl_pct = fidelity.get(
            "explained_discordant_pct",
            round((expl_bp / disc_bp * 100.0), 2) if disc_bp > 0 else 0.0,
        )

        rows.append(
            {
                "sample": sample_name,
                "verdict": data.get("verdict", "UNKNOWN"),
                "num_contigs": len(junction),
                "total_length_bp": total_length,
                "completeness_pct": (
                    data.get("completeness_pct")
                    if data.get("completeness_pct") is not None
                    else "-"
                ),
                "primary_contig_id": primary_id,
                "primary_length_bp": p_junc.get("length", 0),
                "primary_depth_x": p_junc.get("avg_depth", 0.0),
                "primary_junction_score": p_junc.get("spanning_score", 0.0),
                "primary_circularity_status": p_junc.get("status", "UNKNOWN"),
                "primary_symmetry_deg": p_rep.get("symmetry", "-"),
                "primary_oric_offset_kb": p_rep.get("oric_offset_kb", "-"),
                "primary_rotation_status": p_rep.get("rotation_status", "-"),
                "primary_dnaa_fold_enrichment": p_rep.get("dnaa_fold_enrichment", "-"),
                "primary_dnaa_poisson_pval": p_rep.get("dnaa_poisson_pval", "-"),
                "primary_replication_status": p_rep.get("viability", "-"),
                "circular_plasmids_count": circular_plasmids,
                "failed_contigs_count": failed_contigs,
                "masked_bases_bp": fidelity.get("masked_bases", 0),
                "masked_pct": fidelity.get("masked_pct", 0.0),
                "low_concordance_blocks_count": len(low_concord_blocks),
                "evaluated_bases_bp": fidelity.get("evaluated_bases", 0),
                "consensus_error_bases_bp": cons_err_bp,
                "subconsensus_bias_bases_bp": sub_bias_bp,
                "motif_families_count": len(motifs),
                "explained_discordant_pct": expl_pct,
                "top_motif_families": top_motifs_str,
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

    print(f"[DONE] Summarized {len(rows)} assembly report(s) to: {args.output}")

    if motif_cohort_tracker:
        base_out, ext = os.path.splitext(args.output)
        motif_out_path = f"{base_out}_recurrent_motifs{ext or '.tsv'}"
        motif_fields = [
            "lead_motif",
            "samples_count",
            "cohort_prevalence_pct",
            "batch_classification",
            "total_unique_discordant_sites",
            "total_consensus_errors",
            "mean_fold_enrichment",
            "signal_classes",
            "samples_list",
        ]
        num_samples = max(1, len(rows))
        motif_rows = []
        for m_name, info in motif_cohort_tracker.items():
            s_count = len(info["samples"])
            prev_pct = round((s_count / num_samples) * 100.0, 1)
            batch_class = (
                "Recurrent Batch/Basecaller Signature"
                if s_count >= 2
                else "Sample-Specific Motif"
            )
            mean_fold = round(sum(info["folds"]) / len(info["folds"]), 1)
            motif_rows.append(
                {
                    "lead_motif": m_name,
                    "samples_count": s_count,
                    "cohort_prevalence_pct": prev_pct,
                    "batch_classification": batch_class,
                    "total_unique_discordant_sites": info["total_unique_sites"],
                    "total_consensus_errors": info["total_consensus_errors"],
                    "mean_fold_enrichment": mean_fold,
                    "signal_classes": "; ".join(sorted(info["classes"])) or "-",
                    "samples_list": ",".join(info["samples"]),
                }
            )
        motif_rows.sort(
            key=lambda x: (x["samples_count"], x["total_unique_discordant_sites"]),
            reverse=True,
        )
        with open(motif_out_path, "w", newline="", encoding="utf-8") as mf:
            m_writer = csv.DictWriter(mf, fieldnames=motif_fields, delimiter=delimiter)
            m_writer.writeheader()
            m_writer.writerows(motif_rows)
        print(f"[DONE] Saved cross-sample motif aggregation to: {motif_out_path}")

    return 0


def main(argv: list[str] | None = None) -> int:
    effective_argv = sys.argv[1:] if argv is None else argv

    if effective_argv and effective_argv[0] == "summarize":
        return run_summarize(effective_argv[1:])
    if effective_argv and effective_argv[0] == "run":
        effective_argv = effective_argv[1:]

    args = parse_run_args(effective_argv)

    if not os.path.isfile(args.assembly):
        print(f"[ERROR] Assembly file not found: {args.assembly}", file=sys.stderr)
        return 1
    if args.reads and not os.path.isfile(args.reads):
        print(f"[ERROR] Reads file not found: {args.reads}", file=sys.stderr)
        return 1
    if args.bam and not os.path.isfile(args.bam):
        print(f"[ERROR] BAM file not found: {args.bam}", file=sys.stderr)
        return 1

    os.makedirs(args.outdir, exist_ok=True)
    sample_name = args.prefix or os.path.basename(args.assembly)

    resolved_log_path = None
    if args.log_file:
        resolved_log_path = (
            os.path.join(args.outdir, args.log_file)
            if not os.path.isabs(args.log_file) and os.path.dirname(args.log_file) == ""
            else args.log_file
        )

    sbam_logger = logging.getLogger("sbam")

    with configure_logging(resolved_log_path):
        sbam_logger.info("[INIT] Starting SBAM v%s analysis...", __version__)
        sbam_logger.info("[INIT] Sample: %s | Assembly: %s", sample_name, args.assembly)
        sbam_logger.info("[INIT] Output directory: %s\n", args.outdir)

        assembly = load_assembly(args.assembly)
        if not assembly.records:
            sbam_logger.error("[ERROR] No valid FASTA records found in %s.", args.assembly)
            return 1

        sbam_logger.info("Structural Validation")
        created_intermediates: list[str] = []

        if args.bam:
            try:
                is_cyclic = _ensure_bam_index_and_detect_cyclic(
                    args.bam, list(assembly.lengths.keys()), args.threads
                )
            except (RuntimeError, ValueError) as exc:
                sbam_logger.error("[ERROR] %s", exc)
                return 1
            bam_path = args.bam
            orig_lengths = assembly.lengths
            sbam_logger.info(
                " > [BAM] Using pre-aligned BAM (%s, topology=%s)",
                bam_path,
                "cyclic" if is_cyclic else "linear",
            )
        else:
            genome = GenomeStructure(args.assembly, args.outdir, records=assembly.records)
            ref_fasta, orig_lengths, is_cyclic = genome.prepare_reference(
                buffer_size=args.buffer_size, linear=args.linear
            )
            if is_cyclic:
                created_intermediates.append(ref_fasta)

            aligner = ReadAligner(
                args.reads,
                ref_fasta,
                args.outdir,
                threads=args.threads,
                log_path=resolved_log_path,
                is_cyclic=is_cyclic,
            )
            try:
                bam_path = aligner.run_mapping(read_type=args.read_type)
            except AlignmentError as exc:
                sbam_logger.error("[ERROR] %s", exc)
                return 1
            created_intermediates.extend([bam_path, f"{bam_path}.bai"])

        evaluator = JunctionEvaluator(
            bam_path,
            orig_lengths,
            linear_mode=args.linear,
            is_cyclic=is_cyclic,
        )
        junction_metrics = evaluator.evaluate_junctions()

        physics = BioVal(
            args.assembly,
            min_contig_len=args.min_ori_length,
            records=assembly.records,
            linear_mode=args.linear,
        )
        physics_metrics = physics.calculate_symmetry()

        sbam_logger.info("\n Base-Level Validation")
        fidelity = BaseAccuracy(
            bam_path=bam_path,
            fasta_path=args.assembly,
            threads=args.threads,
            k_list=args.k_list,
            min_concordance=args.min_concordance,
            depth_variance_limit=args.depth_variance,
            records=assembly.records,
            linear_mode=args.linear,
            is_cyclic=is_cyclic,
        )
        motif_families = fidelity.analyze_motifs()

        track_paths: dict[str, str] = {}
        if args.emit_tracks:
            bed_path = os.path.join(args.outdir, f"{sample_name}_masked_regions.bed")
            tsv_path = os.path.join(args.outdir, f"{sample_name}_consensus_errors.tsv")
            vcf_path = os.path.join(args.outdir, f"{sample_name}_consensus_errors.vcf")
            export_masked_regions_bed(bed_path, fidelity.all_masked_regions)
            export_consensus_errors_tsv(tsv_path, fidelity.consensus_error_sites)
            export_consensus_errors_vcf(vcf_path, fidelity.consensus_error_sites, orig_lengths)
            track_paths = {
                "masked_regions_bed": bed_path,
                "consensus_errors_tsv": tsv_path,
                "consensus_errors_vcf": vcf_path,
            }

        sbam_logger.info("\n Generating Dashboard")
        builder = DashboardBuilder(
            args.outdir,
            args.assembly,
            offline=args.offline,
            sample_name=sample_name,
        )
        report_file = builder.generate_report(
            junction_metrics=junction_metrics,
            physics_metrics=physics_metrics,
            motif_results=motif_families,
            masked_bases=fidelity.masked_bases,
            masked_pct=fidelity.masked_pct,
            masked_regions=fidelity.masked_regions,
            low_concordance_blocks=fidelity.low_concordance_blocks,
            genome_size=args.genome_size,
            min_concordance=args.min_concordance,
            min_ori_length=args.min_ori_length,
            linear_mode=args.linear,
            fidelity=fidelity,
        )

        summary_json_path = os.path.join(
            args.outdir, f"{sample_name}_sbam_metrics.json"
        )
        compact_physics = {
            cid: {k: v for k, v in m.items() if k not in ("plot_theta", "plot_r", "plot_x_kb")}
            for cid, m in physics_metrics.items()
        }
        with open(summary_json_path, "w", encoding="utf-8") as jf:
            json.dump(
                {
                    "schema_version": SCHEMA_VERSION,
                    "sbam_version": __version__,
                    "sample_name": sample_name,
                    "assembly": os.path.basename(args.assembly),
                    "linear_mode": args.linear,
                    "min_ori_length": args.min_ori_length,
                    "verdict": builder.last_verdict,
                    "expected_genome_size": args.genome_size,
                    "completeness_pct": builder.last_completeness_pct,
                    "junction_metrics": junction_metrics,
                    "replication_metrics": compact_physics,
                    "tracks": track_paths,
                    "base_fidelity": {
                        "masked_bases": fidelity.masked_bases,
                        "masked_pct": round(fidelity.masked_pct, 4),
                        "evaluated_bases": fidelity.evaluated_bases,
                        "discordant_bases": fidelity.discordant_bases,
                        "consensus_error_bases": fidelity.consensus_error_bases,
                        "consensus_errors_per_mb": fidelity.consensus_errors_per_mb,
                        "empirical_qv": fidelity.empirical_qv,
                        "empirical_qv_capped": fidelity.empirical_qv_capped,
                        "subconsensus_bias_bases": fidelity.subconsensus_bias_bases,
                        "subconsensus_bias_per_mb": fidelity.subconsensus_bias_per_mb,
                        "baseline_discordance_rate": round(fidelity.baseline_discordance_rate, 6),
                        "consensus_error_rate": round(fidelity.consensus_error_rate, 6),
                        "snp_rate": fidelity.snp_rate,
                        "indel_rate": fidelity.indel_rate,
                        "total_mismatches": fidelity.total_mismatches,
                        "total_insertions": fidelity.total_insertions,
                        "total_deletions": fidelity.total_deletions,
                        "consensus_error_snp_sites": fidelity.consensus_error_snp_sites,
                        "consensus_error_indel_sites": fidelity.consensus_error_indel_sites,
                        "read_strand_biased_sites": fidelity.read_strand_biased_sites,
                        "explained_discordant_bases": fidelity.explained_discordant_bases,
                        "explained_discordant_pct": round(fidelity.explained_discordant_pct, 2),
                        "low_concordance_blocks": fidelity.low_concordance_blocks,
                        "masked_regions": fidelity.masked_regions,
                        "contig_base_metrics": fidelity.contig_base_metrics,
                        "systematic_motifs": motif_families,
                    },
                },
                jf,
                indent=2,
                default=json_serializable,
            )

        if args.clean:
            for fpath in created_intermediates:
                if fpath and os.path.exists(fpath):
                    try:
                        os.remove(fpath)
                    except OSError:
                        pass

        sbam_logger.info("\n[DONE] HTML report: %s", report_file)
        sbam_logger.info("[DONE] JSON metrics: %s", summary_json_path)
        if track_paths:
            sbam_logger.info("[DONE] Diagnostic tracks: %s", ", ".join(track_paths.values()))

    return 0
