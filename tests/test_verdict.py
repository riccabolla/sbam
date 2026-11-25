from sbam.models import Verdict
from sbam.verdict import evaluate_verdict


def _base_inputs(primary_depth=80.0, primary_len=5_000_000, j_status="PASS", p_viability="ACCEPTABLE"):
    junction_metrics = {
        "chr1": {
            "length": primary_len,
            "avg_depth": primary_depth,
            "copy_number": 1.0,
            "spanning_score": 0.98 if j_status in ("PASS", "LINEAR") else 0.10,
            "status": j_status,
            "classification": "Primary Chromosome",
        }
    }
    physics_metrics = {
        "chr1": {
            "length": primary_len,
            "symmetry": 176.0 if p_viability == "ACCEPTABLE" else 130.0,
            "viability": p_viability,
            "rotation_status": "ROTATED",
        }
    }
    return junction_metrics, physics_metrics


def test_verdict_clean_accept():
    jm, pm = _base_inputs()
    res = evaluate_verdict(
        junction_metrics=jm,
        physics_metrics=pm,
        motif_results=[],
        masked_bases=0,
        masked_pct=0.0,
        evaluated_bases=5_000_000,
        discordant_bases=0,
        consensus_error_bases=0,
    )
    assert res.verdict == Verdict.ACCEPT.value


def test_verdict_low_depth_escalates_on_failed_plasmid_or_masking_or_size():
    """Bug 3 regression: LOW DEPTH <30x must escalate to REVIEW / REVIEW (SIZE ANOMALY)."""
    jm, pm = _base_inputs(primary_depth=20.0)
    jm["plasmid_1"] = {
        "length": 50_000,
        "avg_depth": 15.0,
        "copy_number": 0.75,
        "spanning_score": 0.1,
        "status": "FAIL",
        "classification": "Misassembled Junction",
    }
    res = evaluate_verdict(
        junction_metrics=jm,
        physics_metrics=pm,
        motif_results=[],
        masked_bases=0,
        masked_pct=0.0,
        evaluated_bases=5_000_000,
        consensus_error_bases=0,
    )
    assert res.verdict == Verdict.REVIEW.value
    assert "failed circularity" in res.message

    # Clean plasmid but >5% masked at <30x depth -> REVIEW
    jm2, pm2 = _base_inputs(primary_depth=20.0)
    res2 = evaluate_verdict(
        junction_metrics=jm2,
        physics_metrics=pm2,
        motif_results=[],
        masked_bases=300_000,
        masked_pct=6.0,
        evaluated_bases=4_700_000,
        consensus_error_bases=0,
    )
    assert res2.verdict == Verdict.REVIEW.value

    # Clean at <30x depth but >110% size anomaly -> REVIEW (SIZE ANOMALY)
    jm3, pm3 = _base_inputs(primary_depth=20.0)
    res3 = evaluate_verdict(
        junction_metrics=jm3,
        physics_metrics=pm3,
        motif_results=[],
        masked_bases=0,
        masked_pct=0.0,
        genome_size=4_000_000,
        evaluated_bases=5_000_000,
        consensus_error_bases=0,
    )
    assert res3.verdict == Verdict.REVIEW_SIZE_ANOMALY.value


def test_verdict_does_not_overwrite_consensus_errors_with_motif_sum():
    """Bug 1 regression: consensus_error_bases=2 must not be overwritten by motif sums."""
    jm, pm = _base_inputs()
    motifs = [
        {"motif": "AAAA/TTTT", "consensus_error_count": 3, "assigned_unique_sites": 20},
        {"motif": "CCCC/GGGG", "consensus_error_count": 3, "assigned_unique_sites": 20},
    ]
    res = evaluate_verdict(
        junction_metrics=jm,
        physics_metrics=pm,
        motif_results=motifs,
        masked_bases=0,
        masked_pct=0.0,
        evaluated_bases=5_000_000,
        discordant_bases=40,
        consensus_error_bases=2,
    )
    assert res.consensus_error_bases == 2
    assert res.subconsensus_bias_bases == 38
    assert any("including <50% consensus errors" in n for n in res.info_notes)


def test_verdict_secondary_contig_mask_and_low_concord_triggers():
    """Bug 4 regression: secondary contig masking >5% and low-concordance blocks trigger REVIEW."""
    jm, pm = _base_inputs()
    jm["plasmid_1"] = {
        "length": 100_000,
        "avg_depth": 80.0,
        "copy_number": 1.0,
        "spanning_score": 0.99,
        "status": "PASS",
        "classification": "Circular Plasmid",
    }
    contig_base = {
        "chr1": {"consensus_error_bases": 0, "consensus_errors_per_mb": 0.0, "masked_pct": 0.0},
        "plasmid_1": {
            "consensus_error_bases": 0,
            "consensus_errors_per_mb": 0.0,
            "masked_pct": 6.0,
            "low_concordance_blocks": [
                {"contig": "plasmid_1", "start": 100, "end": 1600, "length": 1500, "mean_concordance": 75.0}
            ],
        },
    }
    res = evaluate_verdict(
        junction_metrics=jm,
        physics_metrics=pm,
        motif_results=[],
        masked_bases=0,
        masked_pct=0.0,
        evaluated_bases=5_000_000,
        consensus_error_bases=0,
        contig_base_metrics=contig_base,
    )
    assert res.verdict == Verdict.REVIEW.value
    assert "plasmid_1" in res.message


def test_verdict_dynamic_min_ori_length_and_incomplete():
    jm, _ = _base_inputs(primary_len=150_000)
    res = evaluate_verdict(
        junction_metrics=jm,
        physics_metrics={},
        motif_results=[],
        masked_bases=0,
        masked_pct=0.0,
        min_ori_length=250_000,
        evaluated_bases=150_000,
        consensus_error_bases=0,
    )
    assert res.verdict == Verdict.ACCEPT_WITH_CAUTION.value
    assert "250 kb threshold" in res.message

    res_inc = evaluate_verdict(
        junction_metrics=jm,
        physics_metrics={},
        motif_results=[],
        masked_bases=0,
        masked_pct=0.0,
        genome_size=500_000,
        evaluated_bases=150_000,
        consensus_error_bases=0,
    )
    assert res_inc.verdict == Verdict.FAIL_INCOMPLETE.value
