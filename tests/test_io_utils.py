"""
test_io_utils.py - Unit tests for file export and I/O utility functions.
"""

import unittest
import os
from io_utils import write_text, write_json, write_csv


class TestIOUtilsModule(unittest.TestCase):
    """Test suite for I/O utilities."""

    def setUp(self):
        self.txt_path = "data/test_out.txt"
        self.json_path = "data/test_out.json"
        self.csv_path = "data/test_out.csv"

    def tearDown(self):
        for path in [self.txt_path, self.json_path, self.csv_path]:
            if os.path.exists(path):
                os.remove(path)

    def test_write_text(self):
        """Test writing plain text."""
        write_text(self.txt_path, "Hello Bioinformatics")
        self.assertTrue(os.path.exists(self.txt_path))
        with open(self.txt_path, "r", encoding="utf-8") as f:
            self.assertEqual(f.read(), "Hello Bioinformatics")

    def test_write_json(self):
        """Test writing structured data to JSON."""
        data = {"sequence": "ATGC", "length": 4, "gc": 50.0}
        write_json(self.json_path, data)
        self.assertTrue(os.path.exists(self.json_path))

    def test_write_csv(self):
        """Test writing tabular data to CSV."""
        rows = [{"pos": 1, "type": "Substitution", "ref": "A", "query": "G"}]
        write_csv(self.csv_path, rows, fieldnames=["pos", "type", "ref", "query"])
        self.assertTrue(os.path.exists(self.csv_path))


if __name__ == "__main__":
    unittest.main()
  
