"""
test_alignment.py - Unit tests for Needleman-Wunsch and Smith-Waterman alignment algorithms.
"""

import unittest
from alignment import needleman_wunsch, smith_waterman, format_alignment


class TestAlignmentModule(unittest.TestCase):
    """Test suite for global and local sequence alignment functions."""

    def test_needleman_wunsch_identical(self):
        """Test Needleman-Wunsch with identical sequences."""
        res = needleman_wunsch("ATGC", "ATGC")
        self.assertEqual(res.score, 8.0)
        self.assertEqual(res.aligned_reference, "ATGC")
        self.assertEqual(res.aligned_query, "ATGC")
        self.assertEqual(res.identity, 100.0)

    def test_needleman_wunsch_gattaca(self):
        """Test Needleman-Wunsch classic educational example: GATTACA vs GATACA."""
        res = needleman_wunsch("GATTACA", "GATACA")
        # Should introduce a gap rather than multiple substitutions
        self.assertIn("-", res.aligned_reference + res.aligned_query)
        self.assertGreater(res.matches, 0)

    def test_needleman_wunsch_invalid_dna(self):
        """Test that invalid DNA raises ValueError."""
        with self.assertRaises(ValueError):
            needleman_wunsch("ATGX", "ATGC")

    def test_smith_waterman_local(self):
        """Test Smith-Waterman local alignment on sequences with strong local matching regions."""
        # Flanking unrelated sequences with a matching core "ATGCATGC"
        ref = "NNNNATGCATGCNNNN"
        query = "XXXXATGCATGCXXXX"
        res = smith_waterman(ref, query)
        self.assertGreater(res.score, 0.0)
        self.assertEqual(res.aligned_reference, "ATGCATGC")
        self.assertEqual(res.aligned_query, "ATGCATGC")

    def test_format_alignment(self):
        """Test alignment text formatting output."""
        formatted = format_alignment("ATGC", "AT-C", width=60)
        self.assertIn("Reference: ATGC", formatted)
        self.assertIn("Query:     AT-C", formatted)
        self.assertIn("|", formatted)


if __name__ == "__main__":
    unittest.main()
      
