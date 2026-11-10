from unittest.mock import MagicMock, patch
import numpy as np
import pytest
from sbam.base.base import BaseAccuracy


@patch("sbam.base.base.mp.Pool")
@patch("sbam.base.base.SeqIO.parse")
def test_p0_end_to_end_motif_recovery_regression(mock_seqio, mock_pool):
    """
    P0 Regression Test: Runs analyze_motifs() end-to-end on a synthetic pileup
    with injected errors at 'CCAGG'. Verifies that point motif discordances are
    not masked out and that majority/minority counters reconcile.
    """
    unit = "TTTTTCCAGG"
    repeats = 200
    seq = unit * repeats
    seq_len = len(seq)

    mock_record = MagicMock()
    mock_record.id = "contig_1"
    mock_record.seq = seq
    mock_seqio.return_value = [mock_record]

    raw_depths = np.full(seq_len, 100, dtype=np.int32)
    raw_matches = np.full(seq_len, 99, dtype=np.int32)

    # Inject 60% concordance (sub-consensus bias) at 195 sites and 40% (<50% consensus error) at 5 sites
    for i in range(repeats):
        raw_matches[i * len(unit) + 7] = 40 if i < 5 else 60

    pool_instance = MagicMock()
    pool_instance.__enter__.return_value = pool_instance
    pool_instance.map.return_value = [(0, seq_len, raw_depths, raw_matches)]
    mock_pool.return_value = mock_instance = pool_instance

    engine = BaseAccuracy(
        "dummy.bam",
        "dummy.fasta",
        threads=1,
        k_list=[5],
        min_concordance=0.85,
        canonical=True,
    )
    with patch.object(engine, "_get_buffer_length", return_value=0):
        families = engine.analyze_motifs(min_occurrences=100, min_discordant_count=10)

    assert engine.evaluated_bases == seq_len
    assert engine.discordant_bases == 200
    assert engine.consensus_error_bases == 5
    assert engine.subconsensus_bias_bases == 195
    assert engine.explained_discordant_bases == 200
    assert engine.explained_discordant_pct == pytest.approx(100.0)

    assert len(families) == 1
    assert families[0]["motif"] == "CCAGG/CCTGG"
    assert families[0]["occurrences"] == 200
    assert families[0]["discordant_count"] == 200
    assert families[0]["consensus_error_count"] == 5


@patch("sbam.base.base.mp.Pool")
@patch("sbam.base.base.SeqIO.parse")
def test_palindromic_kmer_counted_once_per_site(mock_seqio, mock_pool):
    """
    Regression Test: Ensures even-length palindromic k-mers (e.g., 'TCGA') are
    counted exactly once per genomic occurrence (never double-counted across strands)
    and that fold_increase equals error_rate / baseline_discordance_rate.
    """
    unit = "AAAATCGAAAAA"  # 12 bp unit containing one 'TCGA' at indices 4..8
    repeats = 200
    seq = unit * repeats
    seq_len = len(seq)

    mock_record = MagicMock()
    mock_record.id = "contig_1"
    mock_record.seq = seq
    mock_seqio.return_value = [mock_record]

    raw_depths = np.full(seq_len, 100, dtype=np.int32)
    raw_matches = np.full(seq_len, 99, dtype=np.int32)

    # Inject discordance at offset k//2 = 2 of 'TCGA' ('G' at index 6 of each unit)
    for i in range(repeats):
        raw_matches[i * len(unit) + 6] = 60

    pool_instance = MagicMock()
    pool_instance.__enter__.return_value = pool_instance
    pool_instance.map.return_value = [(0, seq_len, raw_depths, raw_matches)]
    mock_pool.return_value = pool_instance

    engine = BaseAccuracy(
        "dummy.bam",
        "dummy.fasta",
        threads=1,
        k_list=[4],
        min_concordance=0.85,
        canonical=True,
    )
    with patch.object(engine, "_get_buffer_length", return_value=0):
        families = engine.analyze_motifs(min_occurrences=100, min_discordant_count=10)

    tcga_hits = [f for f in families if f["motif"] == "TCGA/TCGA"]
    assert len(tcga_hits) == 1
    hit = tcga_hits[0]
    # Must be 200 occurrences (once per site), NOT 400!
    assert hit["occurrences"] == 200
    assert hit["discordant_count"] == 200
    assert hit["is_palindrome"] is True
    assert hit["fold_increase"] == pytest.approx(
        hit["error_rate"] / engine.baseline_discordance_rate
    )