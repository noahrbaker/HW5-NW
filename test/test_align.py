# Importing Dependencies
import pytest
from align import NeedlemanWunsch, read_fasta
import numpy as np

def test_nw_alignment():
    """
    Unit test for NW alignment using test_seq1.fa and test_seq2.fa by
    asserting that the 3 alignment matrices are correctly filled out.
    Uses the BLOSUM62 matrix and a gap open penalty of -10 and a gap extension penalty of -1.
    """
    seq1, _ = read_fasta("./data/test_seq1.fa")
    seq2, _ = read_fasta("./data/test_seq2.fa")
    
    nw = NeedlemanWunsch("./substitution_matrices/BLOSUM62.mat", -10, -1)
    nw.align(seq1, seq2)
    
    # Expected final values for each matrix
    expected_align = np.array([
        [  0., -np.inf, -np.inf, -np.inf],
        [-np.inf,   5., -11., -13.],
        [-np.inf, -12.,   4.,  -8.],
        [-np.inf, -12.,  -1.,   5.],
        [-np.inf, -14.,  -6.,   4.]
    ])
    
    expected_gapA = np.array([
        [-10., -np.inf, -np.inf, -np.inf],
        [-11., -12.,  -6.,  -7.],
        [-12., -13., -14.,  -7.],
        [-13., -14., -15., -12.],
        [-14., -15., -16., -17.]
    ])
    
    expected_gapB = np.array([
        [-10., -11., -12., -13.],
        [-np.inf, -12., -13., -14.],
        [-np.inf,  -6., -14., -15.],
        [-np.inf,  -7.,  -7., -16.],
        [-np.inf,  -8.,  -8.,  -6.]
    ])
    
    # Assert that the matrices are filled correctly
    np.testing.assert_array_equal(nw._align_matrix, expected_align)
    np.testing.assert_array_equal(nw._gapA_matrix, expected_gapA)
    np.testing.assert_array_equal(nw._gapB_matrix, expected_gapB)

def test_nw_backtrace():
    """
    Unit test for NW backtracing using test_seq3.fa and test_seq4.fa by
    asserting that the backtrace is correct.
    Uses the BLOSUM62 matrix with a gap open penalty of -10 and a gap extension penalty of -1.
    """
    seq3, _ = read_fasta("./data/test_seq3.fa")
    seq4, _ = read_fasta("./data/test_seq4.fa")
    
    nw = NeedlemanWunsch("./substitution_matrices/BLOSUM62.mat", -10, -1)
    score, align1, align2 = nw.align(seq3, seq4)
    
    # Expected results
    expected_score = 17
    expected_align1 = "MAVHQLIRRP"
    expected_align2 = "M---QLIRHP"
    
    # Assert that the alignment score and sequences are correct
    assert score == expected_score, f"Expected score {expected_score}, but got {score}"
    assert align1 == expected_align1, f"Expected seq1 alignment {expected_align1}, but got {align1}"
    assert align2 == expected_align2, f"Expected seq2 alignment {expected_align2}, but got {align2}"


