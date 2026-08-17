"""
io_utils.py - Reusable file input/output utilities for exporting results
in TXT, JSON, and CSV formats.
"""

import json
import csv
import os
from typing import Any, Dict, List


def write_text(filepath: str, content: str) -> None:
    """Write plain text content to a file."""
    with open(filepath, "w", encoding="utf-8") as f:
        f.write(content)


def write_json(filepath: str, data: Any) -> None:
    """
    Write structured data to a JSON file. Automatically converts dataclasses
    to dictionaries.
    """
    with open(filepath, "w", encoding="utf-8") as f:
        json.dump(data, f, indent=4, default=lambda o: o.__dict__ if hasattr(o, '__dict__') else str(o))


def write_csv(filepath: str, rows: List[Dict[str, Any]], fieldnames: List[str]) -> None:
    """
    Write tabular data to a CSV file.
    """
    os.makedirs(os.path.dirname(os.path.abspath(filepath)) or '.', exist_ok=True)
    with open(filepath, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)
          
