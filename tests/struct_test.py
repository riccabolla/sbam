import pytest
from unittest.mock import patch, MagicMock
from sbam.struct.ori import BioVal

@patch('sbam.struct.ori.SeqIO.parse')
def test_calculate_symmetry_viable(mock_seqio):
    """Test that a perfectly symmetrical GC skew yields a VIABLE status."""
    # Create a 400kb sequence. 
    # First 200kb is mostly G, next 200kb is mostly C.
    half_1 = "G" * 200000
    half_2 = "C" * 200000
    dummy_seq = half_1 + half_2
    
    mock_record = MagicMock()
    mock_record.id = "contig_1"
    mock_record.seq = dummy_seq
    mock_seqio.return_value = [mock_record]
    
    engine = BioVal("dummy.fasta")
    results = engine.calculate_symmetry(window_size=1000)
    
    assert "contig_1" in results
    assert results["contig_1"]["viability"] == "VIABLE"
    assert 170 <= results["contig_1"]["symmetry"] <= 180

@patch('sbam.struct.ori.SeqIO.parse')
def test_plasmid_cutoff(mock_seqio):
    """Test that sequences under 100kb are skipped."""
    mock_record = MagicMock()
    mock_record.id = "plasmid_1"
    mock_record.seq = "G" * 50000 # 50kb
    mock_seqio.return_value = [mock_record]
    
    engine = BioVal("dummy.fasta")
    results = engine.calculate_symmetry()
    
    assert "plasmid_1" not in results