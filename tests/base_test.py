from unittest.mock import MagicMock, patch
import numpy as np
import pytest
from sbam.base.base import BaseAccuracy


@patch("sbam.base.base.SeqIO.parse")
@patch("sbam.base.base.pysam.AlignmentFile")
def test_motif_analysis_logic(mock_bam, mock_seqio):
    """
    Test that analyze_motifs() detects a 5-mer motif ('CCAGG') with systematic
    basecalling discordance while ignoring concordant background sequence.
    """
    mock_record = MagicMock()
    mock_record.id = "contig_1"

    # Build a sequence with 600 occurrences of 'CCAGG' (>=500 threshold)
    # separated by concordant 'TTTTT' spacers
    unit = "TTTTTCCAGG"
    repeats = 600
    seq = unit * repeats
    mock_record.seq = seq
    mock_seqio.return_value = [mock_record]

    engine = BaseAccuracy("dummy.bam", "dummy.fasta", threads=1, k_list=[5])

    concord_arr = np.ones(len(seq), dtype=float)
    # Drop concordance to 0.40 at the center 'A' (index 7 of each 10-bp unit)
    for i in range(repeats):
        center_pos = i * len(unit) + 7
        concord_arr[center_pos] = 0.40

    mask_arr = np.zeros(len(seq), dtype=bool)

    with patch.object(
        engine, "build_fidelity_profile", return_value=(concord_arr, mask_arr)
    ):
        results = engine.analyze_motifs()

    assert len(results) > 0
    flagged_motifs = {r["motif"]: r for r in results}
    assert "CCAGG" in flagged_motifs
    assert flagged_motifs["CCAGG"]["occurrences"] == 600
    assert flagged_motifs["CCAGG"]["discordant_count"] == 600
    assert flagged_motifs["CCAGG"]["fold_increase"] > 2.0
    assert "adj_p_value" in flagged_motifs["CCAGG"]


@patch("sbam.base.base.mp.Pool")
@patch("sbam.base.base.SeqIO.parse")
def test_cyclic_buffer_depth_folding_and_masking(mock_seqio, mock_pool):
    """
    Verify that build_fidelity_profile folds cyclic buffer depths/matches back
    into [0, buffer_len) and populates mean_depth and mean_concordance on masked regions.
    """
    seq_len = 1000
    buffer_len = 200
    total_len = seq_len + buffer_len

    mock_record = MagicMock()
    mock_record.id = "contig_1"
    mock_record.seq = "A" * seq_len
    mock_seqio.return_value = [mock_record]

    # Simulate split coverage in the first 200 bp (50x at start + 50x in cyclic tail = 100x folded)
    # and a 100-bp low-depth dropout at [500:600] (depth = 10x)
    raw_depths = np.full(total_len, 100, dtype=np.int32)
    raw_depths[:buffer_len] = 50
    raw_depths[seq_len : seq_len + buffer_len] = 50
    raw_depths[500:600] = 10

    raw_matches = raw_depths.copy()

    pool_instance = MagicMock()
    pool_instance.__enter__.return_value = pool_instance
    pool_instance.map.return_value = [(0, total_len, raw_depths, raw_matches)]
    mock_pool.return_value = pool_instance

    engine = BaseAccuracy("dummy.bam", "dummy.fasta", threads=1)
    with patch.object(engine, "_get_buffer_length", return_value=buffer_len):
        concord_arr, mask_arr = engine.build_fidelity_profile()

    # First 200 bp should NOT be masked because 50x + 50x = 100x matches median depth
    assert not np.any(mask_arr[:buffer_len])
    # Region [500:600] has 10x depth (< 100 / 2.5 = 40x) and should be masked
    assert np.all(mask_arr[500:600])
    assert engine.masked_bases == 100
    assert len(engine.masked_regions) == 1
    assert engine.masked_regions[0]["start"] == 500
    assert engine.masked_regions[0]["end"] == 600
    assert engine.masked_regions[0]["mean_depth"] == pytest.approx(10.0)
    assert engine.masked_regions[0]["mean_concordance"] == pytest.approx(100.0)