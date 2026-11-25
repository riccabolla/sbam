from unittest.mock import MagicMock, patch
import pytest
from sbam.struct.eval import JunctionEvaluator


def make_read(length):
    """Creates a mock primary read for depth calculations."""
    r = MagicMock()
    r.is_unmapped = False
    r.is_secondary = False
    r.is_supplementary = False
    r.reference_length = length
    r.query_alignment_length = length
    return r


def make_junction(start, end, clipped=False):
    """
    Creates a mock read crossing or stopping near the junction.
    If clipped=True, attaches 100 bp of soft-clipping (CIGAR op 4) at both ends.
    """
    r = MagicMock()
    r.is_unmapped = False
    r.is_secondary = False
    r.is_supplementary = False
    r.reference_start = start
    r.reference_end = end
    aligned_span = max(1, end - start)
    if clipped:
        r.cigartuples = [(4, 100), (0, aligned_span), (4, 100)]
    else:
        r.cigartuples = [(0, aligned_span)]
    return r


@pytest.fixture
def mock_bam():
    with patch("sbam.struct.eval.pysam.AlignmentFile") as mock_pysam:
        bam_instance = MagicMock()
        bam_instance.__enter__.return_value = bam_instance
        mock_pysam.return_value = bam_instance
        yield bam_instance


def test_junction_scoring_chromosome(mock_bam):
    """Test standard chromosome (>50kb) passes at >0.60 when clipped reads are minority."""
    orig_lengths = {"contig_1": 5000000}
    evaluator = JunctionEvaluator("dummy.bam", orig_lengths)

    depth_reads = [make_read(5000000) for _ in range(100)]

    junction_reads = (
        [make_junction(4998000, 5002000, clipped=False) for _ in range(70)]  # 70 Spanning
        + [make_junction(4999500, 5000500, clipped=True) for _ in range(30)]  # 30 Clipped/Broken
    )

    def fetch_side_effect(contig, start=None, end=None):
        if start is None:
            return depth_reads
        return junction_reads

    mock_bam.fetch.side_effect = fetch_side_effect

    results = evaluator.evaluate_junctions()

    assert "contig_1" in results
    assert results["contig_1"]["status"] == "PASS"
    assert results["contig_1"]["spanning_score"] == pytest.approx(0.70)
    mock_bam.__exit__.assert_called_once()


def test_junction_scoring_plasmid_fail(mock_bam):
    """Test small plasmid (<50kb) fails if score is <0.40."""
    orig_lengths = {"plasmid_1": 15000}
    evaluator = JunctionEvaluator("dummy.bam", orig_lengths)

    depth_reads = [make_read(15000) for _ in range(100)]

    junction_reads = (
        [make_junction(14000, 16000, clipped=False) for _ in range(20)]  # 20 Spanning
        + [make_junction(14900, 15100, clipped=True) for _ in range(80)]  # 80 Clipped/Broken
    )

    def fetch_side_effect(contig, start=None, end=None):
        if start is None:
            return depth_reads
        return junction_reads

    mock_bam.fetch.side_effect = fetch_side_effect

    results = evaluator.evaluate_junctions()

    assert "plasmid_1" in results
    assert results["plasmid_1"]["status"] == "FAIL"
    assert results["plasmid_1"]["spanning_score"] == pytest.approx(0.20)
    mock_bam.__exit__.assert_called_once()


def test_natural_unclipped_read_ends_ignored(mock_bam):
    """
    Verify that unclipped reads ending naturally within min_overlap of the junction
    are recorded as neutral_reads and do not lower the spanning_score.
    """
    orig_lengths = {"plasmid_2": 250000}
    evaluator = JunctionEvaluator("dummy.bam", orig_lengths)

    depth_reads = [make_read(250000) for _ in range(100)]

    junction_reads = (
        [make_junction(248000, 252000, clipped=False) for _ in range(97)]  # 97 Spanning
        + [make_junction(249500, 250500, clipped=False) for _ in range(66)]  # 66 Natural unclipped ends
        + [make_junction(249500, 250500, clipped=True) for _ in range(7)]  # 7 True clipped reads
    )

    def fetch_side_effect(contig, start=None, end=None):
        if start is None:
            return depth_reads
        return junction_reads

    mock_bam.fetch.side_effect = fetch_side_effect

    results = evaluator.evaluate_junctions()

    assert results["plasmid_2"]["status"] == "PASS"
    assert results["plasmid_2"]["spanning_reads"] == 97
    assert results["plasmid_2"]["broken_reads"] == 7
    assert results["plasmid_2"]["neutral_reads"] == 66
    assert results["plasmid_2"]["spanning_score"] == pytest.approx(97 / 104, abs=1e-3)