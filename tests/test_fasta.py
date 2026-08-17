"""
test_fasta.py - Unit tests for FASTA parsing and writing.
"""

import unittest
import os
from fasta import parse_fasta, read_fasta, write_fasta, FastaRecord


class TestFastaModule(unittest.TestCase):
    """Test suite for FASTA parsing and writing functions."""

    def setUp(self):
        self.test_output_file = "data/test_output.fasta"

    def tearDown(self):
        if os.path.exists(self.test_output_file):
            os.remove(self.test_output_file)

    def test_parse_fasta_single(self):
        """Test parsing a single FASTA record."""
        text = ">gene_1 Description text\nATGCATGC"
        records = parse_fasta(text)
        self.assertEqual(len(records), 1)
        self.assertEqual(records[0].name, "gene_1")
        self.assertEqual(records[0].description, "Description text")
        self.assertEqual(records[0].sequence, "ATGCATGC")

    def test_parse_fasta_multiple_and_multiline(self):
        """Test parsing multiple multiline FASTA records."""
        text = (
            ">seq1\n"
            "ATGC\n"
            "ATGC\n"
            ">seq2 Description 2\n"
            "GGCC\n"
            "TTAA\n"
        )
        records = parse_fasta(text)
        self.assertEqual(len(records), 2)
        self.assertEqual(records[0].name, "seq1")
        self.assertEqual(records[0].sequence, "ATGCATGC")
        self.assertEqual(records[1].name, "seq2")
        self.assertEqual(records[1].description, "Description 2")
        self.assertEqual(records[1].sequence, "GGCCTTAA")

    def test_parse_fasta_edge_cases(self):
        """Test empty text, blank lines, and sequence before header."""
        self.assertEqual(parse_fasta(""), [])
        self.assertEqual(parse_fasta("   \n  "), [])

        # Sequence before header
        text = "ATGCATGC\n>header_after\nGGCC"
        records = parse_fasta(text)
        self.assertEqual(len(records), 2)
        self.assertEqual(records[0].name, "unnamed_sequence")
        self.assertEqual(records[0].sequence, "ATGCATGC")
        self.assertEqual(records[1].name, "header_after")

    def test_write_and_read_fasta_roundtrip(self):
        """Test writing FASTA records to file and reading them back."""
        records = [
            FastaRecord(name="rec1", sequence="ATGCATGCATGC", description="desc1"),
            FastaRecord(name="rec2", sequence="GGCCGGCC", description="")
        ]
        write_fasta(records, self.test_output_file, line_width=4)
        
        read_records = read_fasta(self.test_output_file)
        self.assertEqual(len(read_records), 2)
        self.assertEqual(read_records[0].name, "rec1")
        self.assertEqual(read_records[0].sequence, "ATGCATGCATGC")
        self.assertEqual(read_records[1].name, "rec2")
        self.assertEqual(read_records[1].sequence, "GGCCGGCC")


if __name__ == "__main__":
    unittest.main()
      
