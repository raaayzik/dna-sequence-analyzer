"""
orf.py - Open Reading Frame (ORF) detection module for the Bioinformatics Toolkit.
"""

from dataclasses import dataclass
from typing import List
from sequence import clean_sequence, reverse_complement, translate_dna


@dataclass
class ORFResult:
    """Dataclass holding detailed information about a detected ORF."""
    strand: str            # '+' or '-'
    frame: int             # Reading frame (1, 2, or 3)
    start: int             # 1-based inclusive start coordinate on forward reference
    end: int               # 1-based inclusive end coordinate on forward reference
    length: int            # Length in nucleotides
    dna_sequence: str      # Nucleotide sequence of the ORF
    protein_sequence: str  # Translated amino acid sequence


def find_orfs(sequence: str, min_length: int = 30) -> List[ORFResult]:
    """
    Find Open Reading Frames (ORFs) across all three forward and three reverse-complement
    reading frames.

    ORF Definition:
    - Starts with start codon: ATG
    - Ends at the first in-frame stop codon: TAA, TAG, TGA
    - Meets or exceeds the specified minimum length (in nucleotides).

    Coordinates Convention:
    - User-facing coordinates are 1-based inclusive relative to the original forward sequence.

    Args:
        sequence (str): Input DNA sequence string.
        min_length (int): Minimum length of ORF in nucleotides (default: 30).

    Returns:
        List[ORFResult]: List of detected ORF results.
    """
    cleaned = clean_sequence(sequence)
    seq_len = len(cleaned)
    orfs: List[ORFResult] = []

    # 1. Search forward strands (+1, +2, +3 -> frame indices 0, 1, 2)
    for frame_idx in range(3):
        frame_num = frame_idx + 1
        i = frame_idx
        while i < seq_len - 2:
            codon = cleaned[i:i+3]
            if codon == "ATG":
                # Search for in-frame stop codon
                found_stop = False
                stop_idx = -1
                for j in range(i + 3, seq_len - 2, 3):
                    stop_codon = cleaned[j:j+3]
                    if stop_codon in ("TAA", "TAG", "TGA"):
                        stop_idx = j
                        found_stop = True
                        break

                if found_stop:
                    orf_end_idx = stop_idx + 3  # Include stop codon in DNA sequence
                    orf_dna = cleaned[i:orf_end_idx]
                    orf_len = len(orf_dna)

                    if orf_len >= min_length:
                        # Translate without stop codon appearing in protein output
                        protein = translate_dna(orf_dna, frame=0, stop_at_stop=True)
                        orfs.append(ORFResult(
                            strand="+",
                            frame=frame_num,
                            start=i + 1,          # 1-based start
                            end=orf_end_idx,      # 1-based inclusive end
                            length=orf_len,
                            dna_sequence=orf_dna,
                            protein_sequence=protein
                        ))
            i += 3

    # 2. Search reverse-complement strands (-1, -2, -3)
    rev_seq = reverse_complement(cleaned)
    rev_len = len(rev_seq)

    for frame_idx in range(3):
        frame_num = frame_idx + 1
        i = frame_idx
        while i < rev_len - 2:
            codon = rev_seq[i:i+3]
            if codon == "ATG":
                found_stop = False
                stop_idx = -1
                for j in range(i + 3, rev_len - 2, 3):
                    stop_codon = rev_seq[j:j+3]
                    if stop_codon in ("TAA", "TAG", "TGA"):
                        stop_idx = j
                        found_stop = True
                        break

                if found_stop:
                    orf_end_idx = stop_idx + 3
                    orf_dna_rev = rev_seq[i:orf_end_idx]
                    orf_len = len(orf_dna_rev)

                    if orf_len >= min_length:
                        protein = translate_dna(orf_dna_rev, frame=0, stop_at_stop=True)

                        # Map reverse-complement coordinates back to original forward strand coordinates:
                        # orig_end = seq_len - i
                        # orig_start = seq_len - orf_end_idx + 1
                        orig_end = seq_len - i
                        orig_start = seq_len - orf_end_idx + 1

                        orfs.append(ORFResult(
                            strand="-",
                            frame=frame_num,
                            start=orig_start,
                            end=orig_end,
                            length=orf_len,
                            dna_sequence=orf_dna_rev,
                            protein_sequence=protein
                        ))
            i += 3

    return orfs
  
