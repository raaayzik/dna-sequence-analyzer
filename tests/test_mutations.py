"""
test_mutations.py - Unit tests for alignment-based mutation detection and summarization.
"""

import unittest
from mutations import call_mutations, summarize_mutations
from alignment import needleman_wunsch


class TestMutationsModule(unittest.TestCase):
    """Test suite for mutation calling and summary generation."""

    def test_substitution_detection(self):
        """Test detection of a single substitution."""
        ref = "ATGC"
        query = "ATCC"
        # Pre-aligned strings
        mutations = call_mutations(ref, query)
        self.assertEqual(len(mutations), 1)
        self.assertEqual(mutations[0].mutation_type, "Substitution")
        self.assertEqual(mutations[0].reference_base, "G")
        self.assertEqual(mutations[0].query_base, "C")

    def test_insertion_detection(self):
        """Test detection of an insertion."""
        aligned_ref = "AT-GC"
        aligned_query = "ATGGC"
        mutations = call_mutations(aligned_ref, aligned_query)
        self.assertEqual(len(mutations), 1)
        self.assertEqual(mutations[0].mutation_type, "Insertion")
        self.assertEqual(mutations[0].reference_base, "-")
        self.assertEqual(mutations[0].query_base, "G")

    def test_deletion_detection(self):
        """Test detection of a deletion."""
        aligned_ref = "ATGGC"
        aligned_query = "AT-GC"
        mutations = call_mutations(aligned_ref, aligned_query)
        self.assertEqual(len(mutations), 1)
        self.assertEqual(mutations[0].mutation_type, "Deletion")
        self.assertEqual(mutations[0].reference_base, "G")
        self.assertEqual(mutations[0].query_base, "-")

    def test_integration_alignment_and_mutation(self):
        """Integration test: NW alignment followed by mutation calling."""
        ref = "GATTACA"
        query = "GATACA"
        res = needleman_wunsch(ref, query)
        mutations = call_mutations(res.aligned_reference, res.aligned_query)
        summary = summarize_mutations(mutations)
        
        # GATTACA vs GATACA has a deletion of T
        self.assertGreaterEqual(summary["deletions"], 1)

    def test_multiple_mutations_summary(self):
        """Test summary counting across multiple mutation types."""
        # Aligned strings with substitution, insertion, and deletion
        aligned_ref = "A-TGCAG"
        aligned_query = "ACTG-AG"
        mutations = call_mutations(aligned_ref, aligned_query)
        summary = summarize_mutations(mutations)

        self.assertGreater(summary["total"], 0)
        self.assertIn("substitutions", summary)
        self.assertIn("insertions", summary)
        self.assertIn("deletions", summary)

    def test_empty_mutations(self):
        """Test identical sequences yield zero mutations."""
        mutations = call_mutations("ATGC", "ATGC")
        self.assertEqual(len(mutations), 0)
        summary = summarize_mutations(mutations)
        self.assertEqual(summary["total"], 0)


if __name__ == "__main__":
    unittest.main()
  
