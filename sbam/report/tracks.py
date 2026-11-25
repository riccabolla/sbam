"""Diagnostic genomic track exporters (BED6, TSV, VCFv4.2) for IGV and polishing."""

from __future__ import annotations

import csv
import os
from typing import Any


def export_masked_regions_bed(path: str, masked_regions: list[dict[str, Any]]) -> str:
    """Write all masked structural regions in BED6 format."""
    sorted_regions = sorted(
        masked_regions, key=lambda r: (r["contig"], int(r["start"]), int(r["end"]))
    )
    out_dir = os.path.dirname(os.path.abspath(path))
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    with open(path, "w", encoding="utf-8") as fh:
        fh.write(
            "#chrom\tchromStart\tchromEnd\tname\tscore_mean_concordance\tstrand_or_depth\n"
        )
        for reg in sorted_regions:
            name = (
                f"{reg.get('category', 'MASKED')}|{reg.get('reason', '')}"
                .replace("\t", " ")
            )
            score = round(float(reg.get("mean_concordance", 0.0)) * 10)
            depth_label = f"depth={reg.get('mean_depth', 0.0):.1f}x;concord={reg.get('mean_concordance', 0.0):.1f}%"
            fh.write(
                f"{reg['contig']}\t{int(reg['start'])}\t{int(reg['end'])}\t"
                f"{name}\t{score}\t{depth_label}\n"
            )
    return path


def export_consensus_errors_tsv(
    path: str, consensus_error_sites: list[dict[str, Any]]
) -> str:
    """Write all <50% consensus error sites with SNP/indel breakdown and read-strand stats."""
    fieldnames = [
        "contig",
        "pos_1based",
        "ref_base",
        "depth",
        "concordance",
        "matches",
        "mismatches",
        "insertions",
        "deletions",
        "dominant_error_type",
        "fwd_depth",
        "fwd_concordance",
        "rev_depth",
        "rev_concordance",
        "strand_pvalue",
        "context",
    ]
    sorted_sites = sorted(
        consensus_error_sites, key=lambda s: (s["contig"], int(s["pos_1based"]))
    )
    out_dir = os.path.dirname(os.path.abspath(path))
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    with open(path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for site in sorted_sites:
            writer.writerow({k: site.get(k, "") for k in fieldnames})
    return path


def export_consensus_errors_vcf(
    path: str,
    consensus_error_sites: list[dict[str, Any]],
    contig_lengths: dict[str, int] | None = None,
) -> str:
    """Write a minimal VCFv4.2 file of <50% majority consensus disagreement positions."""
    sorted_sites = sorted(
        consensus_error_sites, key=lambda s: (s["contig"], int(s["pos_1based"]))
    )
    out_dir = os.path.dirname(os.path.abspath(path))
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    with open(path, "w", encoding="utf-8") as fh:
        fh.write("##fileformat=VCFv4.2\n")
        fh.write("##source=sbam\n")
        if contig_lengths:
            fh.writelines(f"##contig=<ID={cid},length={int(clen)}>\n" for cid, clen in contig_lengths.items())
        fh.write('##INFO=<ID=DP,Number=1,Type=Integer,Description="Total read depth">\n')
        fh.write('##INFO=<ID=CONC,Number=1,Type=Float,Description="Read-to-assembly concordance">\n')
        fh.write('##INFO=<ID=ETYPE,Number=1,Type=String,Description="Dominant disagreement type">\n')
        fh.write('##INFO=<ID=SPVAL,Number=1,Type=Float,Description="Read-strand Fisher exact p-value">\n')
        fh.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for site in sorted_sites:
            info = (
                f"DP={site['depth']};CONC={site['concordance']:.4f};"
                f"ETYPE={site['dominant_error_type']};SPVAL={site['strand_pvalue']:.6g}"
            )
            fh.write(
                f"{site['contig']}\t{int(site['pos_1based'])}\t.\t"
                f"{site['ref_base']}\tN\t.\tPASS\t{info}\n"
            )
    return path
