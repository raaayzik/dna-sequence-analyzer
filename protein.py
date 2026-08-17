"""
protein.py - Protein sequence analysis module for the Bioinformatics Toolkit.
"""

from collections import Counter
from dataclasses import dataclass
from typing import Dict, List, Tuple


@dataclass
class ProteinStatistics:
    """Dataclass holding amino acid composition and frequency statistics."""
    length: int
    composition: Dict[str, int]
    percentages: Dict[str, float]
    most_common: List[Tuple[str, int]]


def analyze_protein(sequence: str) -> ProteinStatistics:
    """
    Analyze an amino acid sequence to calculate length, composition, percentages,
    and most common residues.

    Args:
        sequence (str): Input amino acid sequence string.

    Returns:
        ProteinStatistics: Dataclass containing protein compositional metrics.
    """
    cleaned = "".join(sequence.split()).upper()
    length = len(cleaned)

    if length == 0:
        return ProteinStatistics(length=0, composition={}, percentages={}, most_common=[])

    counts = Counter(cleaned)
    composition = dict(counts)
    percentages = {aa: round((count / length) * 100, 2) for aa, count in counts.items()}
    most_common = counts.most_common()

    return ProteinStatistics(
        length=length,
        composition=composition,
        percentages=percentages,
        most_common=most_common
    )
  
