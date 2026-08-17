"""
test_sequence.py - Unit tests for sequence cleaning, validation, statistics,
complement, transcription, and translation.
"""

import unittest
from sequence import (
    GENETIC_CODE,
    clean_sequence,
    validate_dna,
    sequence_statistics,
    complement,
    reverse_complement,
    transcribe,
    translate_dna,
)


class TestSequenceModule(unittest.TestCase):
    """Test suite for core DNA sequence functions."""

    def test_clean_sequence(self):
        """Test sequence cleaning with various whitespace and case variations."""
        self.assertEqual(clean_sequence("atgc"), "ATGC")
        self.assertEqual(clean_sequence("AT GC"), "ATGC")
        self.assertEqual(clean_sequence("atg\nc"), "ATGC")
        self.assertEqual(clean_sequence("\t  ATGC  \n"), "ATGC")
        self.assertEqual(clean_sequence("aTtGcC"), "ATTGCC")

    def test_validate_dna(self):
        """Test DNA validation for valid bases, empty inputs, and invalid characters."""
        is_valid, msg = validate_dna("ATGC")
        self.assertTrue(is_valid)
        
        is_valid, msg = validate_dna("atgc")
        self.assertTrue(is_valid)

        # Empty sequence
        is_valid, msg = validate_dna("")
        self.assertFalse(is_valid)

        # Invalid characters (X, numbers, punctuation)
        is_valid, msg = validate_dna("ATGX")
        self.assertFalse(is_valid)
        self.assertIn("X", msg)

        is_valid, msg = validate_dna("AT123GC")
        self.assertFalse(is_valid)

        is_valid, msg = validate_dna("AT.GC")
        self.assertFalse(is_valid)

    def test_sequence_statistics_balanced(self):
        """Test nucleotide statistics for a balanced sequence (AATTGGCC)."""
        stats = sequence_statistics("AATTGGCC")
        self.assertEqual(stats.length, 8)
        self.assertEqual(stats.a_count, 2)
        self.assertEqual(stats.t_count, 2)
        self.assertEqual(stats.g_count, 2)
        self.assertEqual(stats.c_count, 2)
        self.assertEqual(stats.gc_percent, 50.0)
        self.assertEqual(stats.at_percent, 50.0)

    def test_sequence_statistics_unbalanced(self):
        """Test nucleotide statistics for an unbalanced sequence (AAATGC)."""
        stats = sequence_statistics("AAATGC")
        self.assertEqual(stats.length, 6)
        self.assertEqual(stats.a_count, 3)
        self.assertEqual(stats.t_count, 1)
        self.assertEqual(stats.g_count, 1)
        self.assertEqual(stats.c_count, 1)
        self.assertEqual(stats.gc_percent, 33.33)
        self.assertEqual(stats.at_percent, 66.67)

    def test_complement(self):
        """Test DNA complement generation."""
        self.assertEqual(complement("ATGC"), "TACG")
        self.assertEqual(complement("AAAA"), "TTTT")

    def test_reverse_complement(self):
        """Test DNA reverse complement generation."""
        self.assertEqual(reverse_complement("ATGC"), "GCAT")
        self.assertEqual(reverse_complement("AAAC"), "GTTT")

    def test_transcribe(self):
        """Test DNA to RNA transcription."""
        self.assertEqual(transcribe("ATGC"), "AUGC")
        self.assertEqual(transcribe("ATTA"), "AUUA")

    def test_translate_dna(self):
        """Test DNA translation with stop codons and frames."""
        # Standard translation with stop codon stopping execution
        self.assertEqual(translate_dna("ATGAAATAG", frame=0, stop_at_stop=True), "MK")

        # Translation without stopping at stop codon
        self.assertEqual(translate_dna("ATGAAATAG", frame=0, stop_at_stop=False), "MK*")

        # Reading frames 1 and 2
        # Frame 1: TGAAATAG -> ...
        self.assertEqual(translate_dna("ATGAAATAG", frame=1, stop_at_stop=True), "E")

        # Short sequences and incomplete trailing codons
        self.assertEqual(translate_dna("AT", frame=0), "")
        self.assertEqual(translate_dna("ATGAA", frame=0), "M")

    def test_codon_table(self):
        """Verify the genetic code table properties."""
        self.assertEqual(len(GENETIC_CODE), 64)
        self.assertIn("TAA", GENETIC_CODE)
        self.assertIn("TAG", GENETIC_CODE)
        self.assertIn("TGA", GENETIC_CODE)
        self.assertEqual(GENETIC_CODE["TAA"], "*")
        self.assertEqual(GENETIC_CODE["TAG"], "*")
        self.assertEqual(GENETIC_CODE["TGA"], "*")
        self.assertEqual(GENETIC_CODE["ATG"], "M")


if __name__ == "__main__":
    unittest.main()
  
