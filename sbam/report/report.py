"""HTML dashboard builder rendering SBAM QC results via Jinja2."""

from __future__ import annotations

import logging
import os
from datetime import datetime
from typing import Any

from jinja2 import Environment, PackageLoader

from sbam.verdict import evaluate_verdict

logger = logging.getLogger("sbam")


class DashboardBuilder:
    """Generates the self-contained HTML assembly quality report."""

    def __init__(
        self,
        outdir: str,
        assembly_path: str,
        offline: bool = False,
        sample_name: str | None = None,
    ):
        self.outdir = outdir
        self.assembly_name = os.path.basename(assembly_path)
        self.sample_name = sample_name or self.assembly_name
        self.offline = offline
        self.last_verdict = "UNKNOWN"
        self.last_completeness_pct: float | None = None

    def generate_report(
        self,
        junction_metrics: dict[str, dict[str, Any]],
        physics_metrics: dict[str, dict[str, Any]],
        motif_results: list[dict[str, Any]],
        masked_bases: int,
        masked_pct: float,
        masked_regions: list[dict[str, Any]] | None = None,
        low_concordance_blocks: list[dict[str, Any]] | None = None,
        genome_size: int | None = None,
        evaluated_bases: int | None = None,
        discordant_bases: int | None = None,
        consensus_error_bases: int | None = None,
        subconsensus_bias_bases: int | None = None,
        baseline_discordance_rate: float | None = None,
        consensus_error_rate: float | None = None,
        explained_discordant_bases: int | None = None,
        explained_discordant_pct: float | None = None,
        min_concordance: float = 0.85,
        min_ori_length: int = 500000,
        linear_mode: bool = False,
        fidelity: Any = None,
    ) -> str:
        logger.info(" > [Report] Compiling HTML dashboard...")

        masked_regions = masked_regions or []
        motif_results = motif_results or []

        verdict_res = evaluate_verdict(
            junction_metrics=junction_metrics,
            physics_metrics=physics_metrics,
            motif_results=motif_results,
            masked_bases=masked_bases,
            masked_pct=masked_pct,
            masked_regions=masked_regions,
            low_concordance_blocks=low_concordance_blocks,
            genome_size=genome_size,
            evaluated_bases=evaluated_bases,
            discordant_bases=discordant_bases,
            consensus_error_bases=consensus_error_bases,
            subconsensus_bias_bases=subconsensus_bias_bases,
            baseline_discordance_rate=baseline_discordance_rate,
            consensus_error_rate=consensus_error_rate,
            explained_discordant_bases=explained_discordant_bases,
            explained_discordant_pct=explained_discordant_pct,
            min_concordance=min_concordance,
            min_ori_length=min_ori_length,
            fidelity=fidelity,
        )

        self.last_verdict = verdict_res.verdict
        self.last_completeness_pct = verdict_res.completeness_pct

        snp_rate = getattr(fidelity, "snp_rate", 0.0) if fidelity is not None else 0.0
        indel_rate = getattr(fidelity, "indel_rate", 0.0) if fidelity is not None else 0.0
        cons_snp_sites = (
            getattr(fidelity, "consensus_error_snp_sites", 0)
            if fidelity is not None
            else 0
        )
        cons_indel_sites = (
            getattr(fidelity, "consensus_error_indel_sites", 0)
            if fidelity is not None
            else 0
        )
        strand_biased_sites = (
            getattr(fidelity, "read_strand_biased_sites", 0)
            if fidelity is not None
            else 0
        )

        has_genome_size = genome_size is not None and genome_size > 0
        total_length_mb = verdict_res.total_length / 1_000_000

        env = Environment(
            loader=PackageLoader("sbam.report", "templates"),
            autoescape=True,
        )
        template = env.get_template("report.html.j2")
        html_content = template.render(
            sample_name=self.sample_name,
            assembly_name=self.assembly_name,
            linear_mode=linear_mode,
            min_ori_length_kb=max(1, min_ori_length // 1000),
            date=datetime.now().strftime("%Y-%m-%d %H:%M"),
            offline=self.offline,
            result=verdict_res.verdict,
            result_color=verdict_res.color,
            result_msg=verdict_res.message,
            info_notes=verdict_res.info_notes,
            num_contigs=len(junction_metrics) if junction_metrics else 0,
            contig_summaries=verdict_res.contig_summaries,
            primary_contig_id=verdict_res.primary_contig_id,
            j_status=verdict_res.j_status,
            p_status=verdict_res.p_status,
            failed_contigs_count=verdict_res.failed_contigs_count,
            primary_junc_score=verdict_res.primary_junc_score,
            primary_symmetry=verdict_res.primary_symmetry,
            motif_count=len(motif_results),
            depth_sufficient_for_base_qc=verdict_res.depth_sufficient_for_base_qc,
            evaluated_bases=verdict_res.evaluated_bases,
            discordant_bases=verdict_res.discordant_bases,
            discordant_per_mb=verdict_res.discordant_per_mb,
            consensus_error_bases=verdict_res.consensus_error_bases,
            consensus_errors_per_mb=verdict_res.consensus_errors_per_mb,
            qv_display=verdict_res.qv_display,
            subconsensus_bias_bases=verdict_res.subconsensus_bias_bases,
            subconsensus_bias_per_mb=verdict_res.subconsensus_bias_per_mb,
            low_concord_block_count=verdict_res.low_concord_block_count,
            low_concord_total_bp=verdict_res.low_concord_total_bp,
            baseline_discordance_rate=verdict_res.baseline_discordance_rate,
            consensus_error_rate=verdict_res.consensus_error_rate,
            subconsensus_bias_rate=verdict_res.subconsensus_bias_rate,
            explained_discordant_bases=verdict_res.explained_discordant_bases,
            explained_discordant_pct=verdict_res.explained_discordant_pct,
            min_concordance=min_concordance,
            snp_rate=snp_rate,
            indel_rate=indel_rate,
            consensus_error_snp_sites=cons_snp_sites,
            consensus_error_indel_sites=cons_indel_sites,
            read_strand_biased_sites=strand_biased_sites,
            junction=junction_metrics,
            physics=physics_metrics,
            motifs=motif_results,
            masked_bases=masked_bases,
            masked_pct=masked_pct,
            evaluated_pct=100.0 - masked_pct,
            masked_regions=masked_regions,
            has_genome_size=has_genome_size,
            total_length_mb=total_length_mb,
            completeness_pct=verdict_res.completeness_pct or 0.0,
            completeness_status=verdict_res.completeness_status,
            completeness_msg=verdict_res.completeness_msg,
        )

        report_path = os.path.join(self.outdir, f"{self.sample_name}_sbam_report.html")
        with open(report_path, "w", encoding="utf-8") as fh:
            fh.write(html_content)

        logger.info(" > [Report] Dashboard saved to: %s", report_path)
        return report_path
