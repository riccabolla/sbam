"""Pure verdict engine for SBAM assembly quality classification."""

from __future__ import annotations

from typing import Any

from sbam.models import Verdict, VerdictResult
from sbam.utils import compute_empirical_qv

_ESCALATABLE_ACCEPT = {
    Verdict.ACCEPT.value,
    Verdict.ACCEPT_WITH_NOTES.value,
    Verdict.ACCEPT_WITH_CAUTION.value,
    Verdict.ACCEPT_WITH_CAUTION_LOW_DEPTH.value,
}


def evaluate_verdict(
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
    contig_base_metrics: dict[str, dict[str, Any]] | None = None,
    fidelity: Any = None,
) -> VerdictResult:
    """Evaluate all structural, replication, base-fidelity, and completeness rules."""
    masked_regions = masked_regions or []
    motif_results = motif_results or []
    contig_base_metrics = contig_base_metrics or {}

    if fidelity is not None:
        evaluated_bases = getattr(fidelity, "evaluated_bases", evaluated_bases)
        discordant_bases = getattr(fidelity, "discordant_bases", discordant_bases)
        consensus_error_bases = getattr(
            fidelity, "consensus_error_bases", consensus_error_bases
        )
        subconsensus_bias_bases = getattr(
            fidelity, "subconsensus_bias_bases", subconsensus_bias_bases
        )
        baseline_discordance_rate = getattr(
            fidelity, "baseline_discordance_rate", baseline_discordance_rate
        )
        consensus_error_rate = getattr(
            fidelity, "consensus_error_rate", consensus_error_rate
        )
        explained_discordant_bases = getattr(
            fidelity, "explained_discordant_bases", explained_discordant_bases
        )
        explained_discordant_pct = getattr(
            fidelity, "explained_discordant_pct", explained_discordant_pct
        )
        contig_base_metrics = getattr(fidelity, "contig_base_metrics", contig_base_metrics) or {}
        if low_concordance_blocks is None:
            low_concordance_blocks = getattr(fidelity, "low_concordance_blocks", None)

    primary_contig_id = None
    if junction_metrics:
        primary_contig_id = max(
            junction_metrics.keys(), key=lambda k: junction_metrics[k]["length"]
        )

    primary_length = (
        junction_metrics[primary_contig_id]["length"]
        if primary_contig_id and primary_contig_id in junction_metrics
        else 0
    )
    primary_depth = (
        junction_metrics[primary_contig_id]["avg_depth"]
        if primary_contig_id and primary_contig_id in junction_metrics
        else 0.0
    )

    if evaluated_bases is None:
        evaluated_bases = max(0, primary_length - int(masked_bases or 0))
    if discordant_bases is None:
        discordant_bases = 0

    # Never overwrite true pileup/CIGAR consensus_error_bases with motif sums
    if consensus_error_bases is None:
        consensus_error_bases = sum(
            m.get("consensus_error_count", 0) for m in motif_results
        )

    if subconsensus_bias_bases is None:
        subconsensus_bias_bases = max(0, int(discordant_bases) - int(consensus_error_bases))
    else:
        subconsensus_bias_bases = int(subconsensus_bias_bases)

    if not baseline_discordance_rate and evaluated_bases > 0:
        baseline_discordance_rate = discordant_bases / evaluated_bases
    baseline_discordance_rate = float(baseline_discordance_rate or 0.0)

    if not consensus_error_rate and evaluated_bases > 0:
        consensus_error_rate = consensus_error_bases / evaluated_bases
    consensus_error_rate = float(consensus_error_rate or 0.0)

    subconsensus_bias_rate = (
        subconsensus_bias_bases / evaluated_bases if evaluated_bases > 0 else 0.0
    )

    eval_mb = evaluated_bases / 1_000_000 if evaluated_bases > 0 else 0.0
    consensus_errors_per_mb = consensus_error_bases / eval_mb if eval_mb > 0 else 0.0
    subconsensus_bias_per_mb = subconsensus_bias_bases / eval_mb if eval_mb > 0 else 0.0
    discordant_per_mb = discordant_bases / eval_mb if eval_mb > 0 else 0.0

    if evaluated_bases > 0:
        qv_val, qv_capped = compute_empirical_qv(consensus_error_bases, evaluated_bases)
        qv_display = f">{qv_val:.1f}" if qv_capped else f"{qv_val:.1f}"
    else:
        qv_display = "N/A"

    if explained_discordant_bases is None and motif_results:
        explained_discordant_bases = min(
            int(discordant_bases),
            sum(
                m.get("assigned_unique_sites", m.get("family_unique_sites", 0))
                for m in motif_results
            ),
        )
    explained_discordant_bases = int(explained_discordant_bases or 0)

    if explained_discordant_pct is None:
        explained_discordant_pct = (
            (explained_discordant_bases / discordant_bases) * 100.0
            if discordant_bases > 0
            else 0.0
        )
    explained_discordant_pct = float(explained_discordant_pct or 0.0)

    if low_concordance_blocks is None:
        low_concord_blocks = [
            r for r in masked_regions if r.get("mean_concordance", 100.0) < 90.0
        ]
    else:
        low_concord_blocks = list(low_concordance_blocks)
    low_concord_block_count = len(low_concord_blocks)
    low_concord_total_bp = sum(r.get("length", 0) for r in low_concord_blocks)

    primary_junc_score = "N/A"
    j_status = "UNKNOWN"
    if primary_contig_id and primary_contig_id in junction_metrics:
        j_status = junction_metrics[primary_contig_id]["status"]
        primary_junc_score = f"{junction_metrics[primary_contig_id]['spanning_score']:.2f}"

    primary_symmetry = "N/A"
    p_status = "N/A"
    primary_rotation = "ROTATED"
    if primary_contig_id and primary_contig_id in physics_metrics:
        p_status = physics_metrics[primary_contig_id]["viability"]
        primary_symmetry = f"{physics_metrics[primary_contig_id]['symmetry']}°"
        primary_rotation = physics_metrics[primary_contig_id].get("rotation_status", "ROTATED")

    motif_count = len(motif_results)
    depth_sufficient_for_base_qc = primary_depth >= 30.0
    failed_contigs_count = (
        sum(1 for m in junction_metrics.values() if m["status"] == "FAIL")
        if junction_metrics
        else 0
    )

    has_genome_size = genome_size is not None and genome_size > 0
    total_length = (
        sum(m["length"] for m in junction_metrics.values()) if junction_metrics else 0
    )
    total_length_mb = total_length / 1_000_000
    genome_size_mb = (genome_size / 1_000_000) if has_genome_size else 0.0
    completeness_pct = (total_length / genome_size * 100.0) if has_genome_size else 0.0
    primary_pct = (primary_length / genome_size * 100.0) if has_genome_size else 0.0

    completeness_status = "UNKNOWN"
    completeness_msg = ""
    if has_genome_size:
        if completeness_pct < 90.0:
            completeness_status = "FAIL"
            completeness_msg = (
                f"Assembly is missing >10% of expected sequence. "
                f"Total size is {total_length_mb:.1f} Mb vs expected {genome_size_mb:.1f} Mb."
            )
        elif completeness_pct > 110.0:
            completeness_status = "WARNING"
            completeness_msg = (
                "Total assembly is >10% larger than expected genome size "
                "(may reflect accessory plasmids, uncollapsed repeats, or contamination)."
            )
        elif primary_pct < 85.0:
            completeness_status = "WARNING"
            completeness_msg = (
                f"Total size is within range, but primary chromosome is only {primary_pct:.1f}% "
                "of expected size. Chromosome is likely fragmented."
            )
        else:
            completeness_status = "PASS"
            completeness_msg = (
                "Total assembly size (90%–110%) and primary chromosome length (>=85%) match expectations."
            )

    contig_summaries: list[dict[str, Any]] = []
    secondary_high_err_contigs: list[str] = []
    secondary_high_mask_contigs: list[str] = []
    secondary_low_concord_contigs: list[str] = []

    if junction_metrics:
        for cid, j_metrics in junction_metrics.items():
            length_kb = j_metrics["length"] / 1000.0
            c_j_stat = j_metrics["status"]
            raw_class = j_metrics.get("classification", "Plasmid/Fragment")
            c_p_stat = physics_metrics.get(cid, {}).get("viability", "N/A")
            c_rot_stat = physics_metrics.get(cid, {}).get("rotation_status", "ROTATED")

            c_base = contig_base_metrics.get(cid)
            if c_base:
                c_err = c_base["consensus_error_bases"]
                c_err_mb = c_base["consensus_errors_per_mb"]
                c_mask_pct = c_base["masked_pct"]
                base_qc_str = f"{c_err} err ({c_err_mb:.1f}/Mb) | {c_mask_pct:.2f}% masked"
                if cid != primary_contig_id:
                    if c_err_mb > 5.0 and c_err >= 2:
                        secondary_high_err_contigs.append(
                            f"contig {cid} ({c_err} sites <50%, {c_err_mb:.1f}/Mb)"
                        )
                    if c_mask_pct > 5.0:
                        secondary_high_mask_contigs.append(
                            f"contig {cid} ({c_mask_pct:.2f}% masked)"
                        )
                    c_low_blocks = c_base.get("low_concordance_blocks", [])
                    c_low_bp = sum(b.get("length", 0) for b in c_low_blocks)
                    if c_low_bp > 1000 or len(c_low_blocks) >= 2:
                        secondary_low_concord_contigs.append(
                            f"contig {cid} ({len(c_low_blocks)} block(s), {c_low_bp} bp <90%)"
                        )
            elif cid == primary_contig_id:
                base_qc_str = (
                    f"{consensus_error_bases} err ({consensus_errors_per_mb:.2f}/Mb) | "
                    f"{masked_pct:.2f}% masked"
                )
            else:
                base_qc_str = "Not profiled"

            if c_j_stat in ("PASS", "LINEAR"):
                if c_p_stat == "ACCEPTABLE":
                    evaluation = "Reliable Replicon Graph"
                    if c_rot_stat == "UNROTATED":
                        evaluation += " (Unrotated oriC >200 kb from 0 bp)"
                    elif c_rot_stat == "NEAR_ORIC":
                        evaluation += " (oriC within 200 kb of 0 bp)"
                elif c_p_stat == "ATYPICAL":
                    evaluation = "Structurally Intact, Atypical replication symmetry"
                elif c_p_stat == "WARNING":
                    evaluation = "Structurally Intact, Distorted replication symmetry"
                elif "Sub-stoichiometric" in raw_class:
                    evaluation = "Structurally Intact (<0.7c: partial plasmid loss or extraction bias)"
                else:
                    evaluation = "Structurally Intact"
            elif c_j_stat == "NO_DATA":
                evaluation = "Insufficient Read Coverage (<5x)"
            else:
                if "Misassembled" in raw_class:
                    evaluation = "Contradicted Junction"
                elif "Debris" in raw_class:
                    evaluation = "Likely Low-Coverage Artifact"
                else:
                    evaluation = "Unbridged / Linear Fragment"

            contig_summaries.append(
                {
                    "id": cid,
                    "type": raw_class,
                    "length_kb": length_kb,
                    "j_stat": c_j_stat,
                    "base_qc_str": base_qc_str,
                    "evaluation": evaluation,
                }
            )

    info_notes: list[str] = []
    review_triggers: list[str] = []
    min_ori_kb = max(1, min_ori_length // 1000)

    if j_status == "NO_DATA":
        result = Verdict.INSUFFICIENT_COVERAGE.value
        result_color = "#7f8c8d"
        result_msg = (
            "Average read depth on the primary contig is below 5x. "
            "Cannot reliably evaluate junction spanning or base-level concordance."
        )
    elif j_status == "FAIL":
        result = Verdict.FAIL.value
        result_color = "#c0392b"
        result_msg = (
            "Primary contig fails junction spanning thresholds. High probability of misassembly, "
            "incomplete circularization, or unresolved terminal repeats."
        )
    elif depth_sufficient_for_base_qc and consensus_errors_per_mb > 50.0:
        result = Verdict.FAIL_GROSS_BASE_DISCORDANCE.value
        result_color = "#c0392b"
        result_msg = (
            f"Primary chromosome exhibits a gross majority consensus error rate of "
            f"{consensus_errors_per_mb:.1f}/Mb ({consensus_error_bases} sites <50% support; QV ≈ {qv_display})."
        )
    else:
        if p_status == "ACCEPTABLE":
            result = Verdict.ACCEPT.value
            result_color = "#27ae60"
            result_msg = (
                f"Strong evidence of structural continuity and bidirectional replication architecture "
                f"(empirical read-concordance QV ≈ {qv_display}, {consensus_errors_per_mb:.2f} consensus err/Mb)."
            )
        elif p_status in ("ATYPICAL", "WARNING"):
            result = Verdict.REVIEW.value
            result_color = "#f39c12"
            result_msg = (
                "Primary junction is intact, but replication symmetry deviates from standard "
                "bidirectional geometry."
            )
        else:
            result = Verdict.ACCEPT_WITH_CAUTION.value
            result_color = "#f39c12"
            result_msg = (
                f"Primary contig is structurally intact, but smaller than the {min_ori_kb} kb threshold "
                "for chromosomal replication analysis."
            )

        if not depth_sufficient_for_base_qc:
            if result == Verdict.ACCEPT.value:
                result = Verdict.ACCEPT_WITH_CAUTION_LOW_DEPTH.value
                result_color = "#f39c12"
            info_notes.append(
                f"Primary depth ({primary_depth:.1f}x) is below 30x; insufficient depth for high-confidence base-level grading."
            )

        if failed_contigs_count > 0:
            review_triggers.append(
                f"{failed_contigs_count} secondary contig(s) failed circularity"
            )
        if masked_pct > 5.0:
            review_triggers.append(
                f"{masked_pct:.2f}% of the primary contig is structurally masked (>5% limit)"
            )

        if depth_sufficient_for_base_qc:
            if consensus_errors_per_mb > 5.0:
                review_triggers.append(
                    f"consensus error rate is {consensus_errors_per_mb:.2f}/Mb "
                    f"({consensus_error_bases} sites <50% support, QV ≈ {qv_display}; threshold <=5.0/Mb)"
                )
            elif consensus_error_bases > 0:
                info_notes.append(
                    f"{consensus_error_bases} site(s) with <50% read support "
                    f"({consensus_errors_per_mb:.2f}/Mb, empirical QV ≈ {qv_display}) — within <=5.0/Mb acceptance threshold."
                )

            if low_concord_total_bp > 1000 or low_concord_block_count >= 2:
                worst_blk = min(
                    low_concord_blocks,
                    key=lambda b: b.get("mean_concordance", 100.0),
                )
                review_triggers.append(
                    f"{low_concord_block_count} excluded block(s) with <90% concordance totaling {low_concord_total_bp} bp "
                    f"(worst: {worst_blk['start']}–{worst_blk['end']} at {worst_blk['mean_concordance']}%)"
                )
            elif low_concord_block_count == 1:
                blk = low_concord_blocks[0]
                info_notes.append(
                    f"1 contiguous excluded block ({blk['length']} bp at {blk['start']}–{blk['end']}) has "
                    f"{blk['mean_concordance']}% mean concordance ({blk['mean_depth']}x depth) — inspect in IGV for local repeat/variant."
                )

            if secondary_high_err_contigs:
                review_triggers.append(
                    "elevated consensus errors (>5/Mb) on plasmid/secondary "
                    + ", ".join(secondary_high_err_contigs)
                )
            if secondary_high_mask_contigs:
                review_triggers.append(
                    "high structural masking (>5%) on secondary "
                    + ", ".join(secondary_high_mask_contigs)
                )
            if secondary_low_concord_contigs:
                review_triggers.append(
                    "low-concordance excluded block(s) on secondary "
                    + ", ".join(secondary_low_concord_contigs)
                )

        if motif_count > 0:
            has_motif_cons_err = any(
                m.get("consensus_error_count", 0) > 0 for m in motif_results
            )
            motif_desc = (
                "enriched read discordance (including <50% consensus errors and sub-consensus bias)"
                if has_motif_cons_err
                else "sub-consensus read disagreement"
            )
            info_notes.append(
                f"{motif_count} motif family(ies) show {motif_desc} "
                f"(explaining {explained_discordant_bases}/{discordant_bases} [{explained_discordant_pct:.1f}%] "
                f"of sites <{min_concordance*100:.0f}% support)."
            )

        if primary_rotation == "UNROTATED":
            info_notes.append(
                "Primary chromosome oriC is >200 kb from coordinate 0 bp (consider reorienting with dnaapler if required)."
            )

        if review_triggers:
            if result in _ESCALATABLE_ACCEPT:
                result = Verdict.REVIEW.value
                result_color = "#f39c12"
            result_msg += " REVIEW FLAGGED: " + "; ".join(review_triggers) + "."
        elif result == Verdict.ACCEPT.value and (
            consensus_error_bases > 0 or low_concord_block_count > 0
        ):
            result = Verdict.ACCEPT_WITH_NOTES.value

    if has_genome_size and j_status != "NO_DATA":
        if completeness_status == "FAIL":
            result = Verdict.FAIL_INCOMPLETE.value
            result_color = "#c0392b"
            result_msg += (
                f" Assembly is critically incomplete ({completeness_pct:.1f}% of expected size)."
            )
        elif completeness_status == "WARNING":
            if result in _ESCALATABLE_ACCEPT:
                result = Verdict.REVIEW_SIZE_ANOMALY.value
                result_color = "#f39c12"
            if primary_pct < 85.0 and 90.0 <= completeness_pct <= 110.0:
                result_msg += f" WARNING: Primary contig is only {primary_pct:.1f}% of expected genome size."
            else:
                result_msg += f" WARNING: Total assembly size is {completeness_pct:.1f}% of expected genome size."

    return VerdictResult(
        verdict=result,
        color=result_color,
        message=result_msg,
        info_notes=info_notes,
        review_triggers=review_triggers,
        completeness_status=completeness_status,
        completeness_msg=completeness_msg,
        completeness_pct=round(completeness_pct, 2) if has_genome_size else None,
        contig_summaries=contig_summaries,
        primary_contig_id=primary_contig_id,
        primary_length=primary_length,
        primary_depth=primary_depth,
        total_length=total_length,
        failed_contigs_count=failed_contigs_count,
        low_concord_blocks=low_concord_blocks,
        low_concord_block_count=low_concord_block_count,
        low_concord_total_bp=low_concord_total_bp,
        evaluated_bases=evaluated_bases,
        discordant_bases=discordant_bases,
        consensus_error_bases=consensus_error_bases,
        subconsensus_bias_bases=subconsensus_bias_bases,
        baseline_discordance_rate=baseline_discordance_rate,
        consensus_error_rate=consensus_error_rate,
        subconsensus_bias_rate=subconsensus_bias_rate,
        consensus_errors_per_mb=consensus_errors_per_mb,
        subconsensus_bias_per_mb=subconsensus_bias_per_mb,
        discordant_per_mb=discordant_per_mb,
        qv_display=qv_display,
        explained_discordant_bases=explained_discordant_bases,
        explained_discordant_pct=explained_discordant_pct,
        depth_sufficient_for_base_qc=depth_sufficient_for_base_qc,
        primary_junc_score=primary_junc_score,
        j_status=j_status,
        primary_symmetry=primary_symmetry,
        p_status=p_status,
        primary_rotation=primary_rotation,
    )
