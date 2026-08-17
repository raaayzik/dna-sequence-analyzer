"""
test_orf.py - Unit tests for Open Reading Frame (ORF) detection.
"""

import unittest
from orf import find_orfs


class TestORFModule(unittest.TestCase):
    """Test suite for ORF detection functions."""

    def test_find_orfs_obvious(self):
        """Test detection of an obvious forward-strand ORF."""
        # Long sequence (ATG + 100 bases of repeating codons + TAA) to satisfy min_length=30
        seq = "ATG" + ("AAACCC" * 20) + "TAA"
        orfs = find_orfs(seq, min_length=30)
        self.assertEqual(len(orfs), 1)
        self.assertEqual(orfs[0].strand, "+")
        self.assertTrue(len(orfs[0].protein_sequence) > 0)

    def test_find_orfs_none(self):
        """Test sequence with no valid ORFs."""
        seq = "CGATCGATCGATCGATCGATCGATCGATCGATCG"
        orfs = find_orfs(seq, min_length=30)
        self.assertEqual(len(orfs), 0)

    def test_find_orfs_min_length_filter(self):
        """Test that min_length filters out short ORFs."""
        # Short ORF (9 nucleotides)
        seq = "ATGAAATAG"
        orfs_strict = find_orfs(seq, min_length=30)
        self.assertEqual(len(orfs_strict), 0)

        orfs_loose = find_orfs(seq, min_length=3)
        self.assertEqual(len(orfs_loose), 1)

    def test_find_orfs_reverse_strand(self):
        """Test detection of reverse-complement strand ORFs."""
        # Reverse complement sequence containing a valid long ORF
        seq = "CGATCGATCGATCGATCGATCGATCGATCTTTA" + ("AAACCC" * 20) + "CAT"
        orfs = find_orfs(seq, min_length=30)
        reverse_orfs = [o for o in orfs if o.strand == "-"]
        self.assertGreater(len(reverse_orfs), 0)


if __name__ == "__main__":
    unittest.main()
    
