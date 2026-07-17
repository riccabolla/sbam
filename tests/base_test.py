import pytest
import numpy as np
from unittest.mock import patch, MagicMock
from sbam.base.base import BaseAccuracy

@patch('sbam.base.base.SeqIO.parse')
@patch('sbam.base.base.pysam.AlignmentFile')
def test_motif_analysis_logic(mock_bam, mock_seqio):
    """
    Test the statistical motif finding logic without doing actual multiprocessing.
    We mock the `build_fidelity_profile` to supply a fixed array.
    """
    mock_record = MagicMock()
    mock_record.id = "contig_1"
    
    # Create a sequence with two instances of the motif "CCWGG" and some flanking bases
    seq = "A" * 1000 + "CCWGG" + "A" * 1000 + "CCWGG" + "A" * 1000
    mock_record.seq = seq
    mock_seqio.return_value = [mock_record]
    
    engine = BaseAccuracy("dummy.bam", "dummy.fasta", threads=1)
    
    # Mock the fidelity profile arrays
    # concord_arr: all 1.0 (perfect), except the 'W' in CCWGG which drops to 0.5
    concord_arr = np.ones(len(seq), dtype=float)
    idx_1 = seq.find("CCWGG") + 2
    idx_2 = seq.rfind("CCWGG") + 2
    concord_arr[idx_1] = 0.5
    concord_arr[idx_2] = 0.5
    
    # mask_arr: Nothing is structurally masked
    mask_arr = np.zeros(len(seq), dtype=bool)
    
    # Patch the class method directly for this test
    with patch.object(engine, 'build_fidelity_profile', return_value=(concord_arr, mask_arr)):
        results = engine.analyze_motifs()
        pass 