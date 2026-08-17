"""
fasta.py - FASTA file and string parser, data model, and writer
for the Bioinformatics Toolkit.
"""

from dataclasses import dataclass
from typing import List
import os


@dataclass
class FastaRecord:
    """Dataclass representing a single FASTA record (name, description, sequence)."""
    name: str
    sequence: str
    description: str = ""


def parse_fasta(text: str) -> List[FastaRecord]:
    """
    Parse FASTA formatted text into a list of FastaRecord objects.

    Handles:
    - Single and multiple records
    - Multiline sequences spanning multiple lines
    - Blank lines and whitespace trimming
    - Sequence data before the first header (assigned a default name)

    Args:
        text (str): FASTA formatted string content.

    Returns:
        List[FastaRecord]: Parsed list of FASTA records.
    """
    records: List[FastaRecord] = []
    if not text or not text.strip():
        return records

    lines = text.splitlines()
    current_name = None
    current_desc = ""
    current_seq_lines: List[str] = []

    for line in lines:
        line = line.strip()
        if not line:
            continue

        if line.startswith(">"):
            # Save previous record if one exists
            if current_name is not None:
                records.append(FastaRecord(
                    name=current_name,
                    sequence="".join(current_seq_lines),
                    description=current_desc
                ))
                current_seq_lines = []

            # Parse header line (>name description)
            header_content = line[1:].strip()
            parts = header_content.split(maxsplit=1)
            current_name = parts[0] if parts else "unnamed_sequence"
            current_desc = parts[1] if len(parts) > 1 else ""
        else:
            if current_name is None:
                # Handle sequence data appearing before the first header
                current_name = "unnamed_sequence"
                current_desc = "Sequence found before first header"
            current_seq_lines.append(line)

    # Save the final record in the file/text
    if current_name is not None:
        records.append(FastaRecord(
            name=current_name,
            sequence="".join(current_seq_lines),
            description=current_desc
        ))

    return records


def read_fasta(filename: str) -> List[FastaRecord]:
    """
    Read and parse a FASTA file from disk.

    Args:
        filename (str): Path to the FASTA file.

    Returns:
        List[FastaRecord]: Parsed list of FASTA records.
    """
    if not os.path.exists(filename):
        raise FileNotFoundError(f"FASTA file not found: {filename}")

    with open(filename, "r", encoding="utf-8") as f:
        content = f.read()

    return parse_fasta(content)


def write_fasta(records: List[FastaRecord], filename: str, line_width: int = 60) -> None:
    """
    Write a list of FastaRecord objects to a file in valid FASTA format.

    Args:
        records (List[FastaRecord]): Records to write.
        filename (str): Output file path.
        line_width (int): Maximum character width for sequence lines (default: 60).
    """
    with open(filename, "w", encoding="utf-8") as f:
        for record in records:
            header = f">{record.name}"
            if record.description:
                header += f" {record.description}"
            f.write(header + "\n")

            seq = record.sequence
            for i in range(0, len(seq), line_width):
                f.write(seq[i:i+line_width] + "\n")
      
