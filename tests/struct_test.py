from unittest.mock import MagicMock, patch
import pytest
from sbam.struct.ori import BioVal


@patch("sbam.struct.ori.SeqIO.parse")
def test_calculate_symmetry_viable(mock_seqio):
    """Test that a symmetrical >=500kb chromosome yields ACCEPTABLE status and oriC offset metrics."""
    half_1 = "G" * 300000
    half_2 = "C" * 300000
    dummy_seq = half_1 + half_2  # 600 kb chromosome

    mock_record = MagicMock()
    mock_record.id = "contig_1"
    mock_record.seq = dummy_seq
    mock_seqio.return_value = [mock_record]

    engine = BioVal("dummy.fasta")
    results = engine.calculate_symmetry(window_size=1000)

    assert "contig_1" in results
    assert results["contig_1"]["viability"] == "ACCEPTABLE"
    assert 170 <= results["contig_1"]["symmetry"] <= 180
    assert "oric_offset_kb" in results["contig_1"]
    assert "oric_offset_deg" in results["contig_1"]
    assert results["contig_1"]["oric_offset_kb"] == pytest.approx(0.0, abs=1.0)


@patch("sbam.struct.ori.SeqIO.parse")
def test_plasmid_and_megaplasmid_cutoff(mock_seqio):
    """Test that sequences under the 500kb chromosomal threshold (e.g., 250kb plasmids) are skipped."""
    mock_record = MagicMock()
    mock_record.id = "plasmid_250kb"
    mock_record.seq = "G" * 250000
    mock_seqio.return_value = [mock_record]

    engine = BioVal("dummy.fasta")
    results = engine.calculate_symmetry()

    assert "plasmid_250kb" not in results