import pytest
from unittest.mock import MagicMock, patch
from sbam.struct.eval import JunctionEvaluator

def make_read(length):
    """Creates a mock read for depth calculations."""
    r = MagicMock()
    r.is_unmapped = False
    r.is_secondary = False
    r.query_alignment_length = length
    return r

def make_junction(start, end):
    """Creates a mock read simulating a junction crossing."""
    r = MagicMock()
    r.is_unmapped = False
    r.is_secondary = False
    r.is_supplementary = False
    r.reference_start = start
    r.reference_end = end
    return r

@pytest.fixture
def mock_bam():
    with patch('sbam.struct.eval.pysam.AlignmentFile') as mock_pysam:
        bam_instance = MagicMock()
        mock_pysam.return_value = bam_instance
        yield bam_instance


# Tests

def test_junction_scoring_chromosome(mock_bam):
    """Test standard chromosome (>50kb) passes at >0.60."""
    orig_lengths = {"contig_1": 5000000}
    evaluator = JunctionEvaluator("dummy.bam", orig_lengths)
    
    # Simulate 100x depth
    depth_reads = [make_read(5000000) for _ in range(100)]
    
    junction_reads = (
        [make_junction(4998000, 5002000) for _ in range(70)] + # 70 Spanning
        [make_junction(4999500, 5000500) for _ in range(30)]   # 30 Broken
    )
    
    # Mock the fetch method to return depth reads for the whole contig and junction reads for the junction region
    def fetch_side_effect(contig, start=None, end=None):
        if start is None:
            return depth_reads
        return junction_reads
        
    mock_bam.fetch.side_effect = fetch_side_effect
    
    results = evaluator.evaluate_junctions()
    
    assert "contig_1" in results
    assert results["contig_1"]["status"] == "PASS"
    assert results["contig_1"]["spanning_score"] == pytest.approx(0.70)
    
    # Verify resource cleanup
    mock_bam.close.assert_called_once()


def test_junction_scoring_plasmid_fail(mock_bam):
    """Test small plasmid (<50kb) fails if score is <0.40."""
    orig_lengths = {"plasmid_1": 15000}
    evaluator = JunctionEvaluator("dummy.bam", orig_lengths)
    
    # Simulate 100x depth (Total bases needed: 1,500,000)
    depth_reads = [make_read(15000) for _ in range(100)]
    
    junction_reads = (
        [make_junction(14000, 16000) for _ in range(20)] + # 20 Spanning
        [make_junction(14900, 15100) for _ in range(80)]   # 80 Broken
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
    
    # Verify resource cleanup
    mock_bam.close.assert_called_once()