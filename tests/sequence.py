"""
sequence.py - Core DNA sequence manipulation, validation, statistics,
transcription, and translation module for the Bioinformatics Toolkit.
"""

from dataclasses import dataclass
from typing import Dict, Tuple

# Standard Genetic Code mapping DNA codons to amino acids (or '*' for stop codons).
# Contains all 64 possible DNA codons.
GENETIC_CODE: Dict[str, str] = {
    # Phenylalanine (F)
    'TTT': 'F', 'TTC': 'F',
    # Leucine (L)
    'TTA': 'L', 'TTG': 'L', 'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
    # Isoleucine (I)
    'ATT': 'I', 'ATC': 'I', 'ATA': 'I',
    # Methionine / Start (M)
    'ATG': 'M',
    # Valine (V)
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
    # Serine (S)
    'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S', 'AGT': 'S', 'AGC': 'S',
    # Proline (P)
    'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    # Threonine (T)
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    # Alanine (A)
    'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    # Tyrosine (Y)
    'TAT': 'Y', 'TAC': 'Y',
    # Stop Codons (*)
    'TAA': '*', 'TAG': '*', 'TGA': '*',
    # Histidine (H)
    'CAT': 'H', 'CAC': 'H',
    # Glutamine (Q)
    'CAA': 'Q', 'CAG': 'Q',
    # Asparagine (N)
    'AAT': 'N', 'AAC': 'N',
    # Lysine (K)
    'AAA': 'K', 'AAG': 'K',
    # Aspartic Acid (D)
    'GAT': 'D', 'GAC': 'D',
    # Glutamic Acid (E)
    'GAA': 'E', 'GAG': 'E',
    # Cysteine (C)
    'TGT': 'C', 'TGC': 'C',
    # Tryptophan (W)
    'TGG': 'W',
    # Arginine (R)
    'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R', 'AGA': 'R', 'AGG': 'R',
}


@dataclass
class SequenceStatistics:
    """Dataclass holding nucleotide counts and compositional percentages."""
    length: int
    a_count: int
    t_count: int
    g_count: int
    c_count: int
    a_percent: float
    t_percent: float
    g_percent: float
    c_percent: float
    gc_percent: float
    at_percent: float


def clean_sequence(sequence: str) -> str:
    """
    Clean a DNA sequence string by converting lowercase letters to uppercase
    and removing spaces, tabs, and newline characters.

    Args:
        sequence (str): Raw input sequence string.

    Returns:
        str: Cleaned uppercase sequence without whitespace.
    """
    if not isinstance(sequence, str):
        raise TypeError("Sequence must be a string.")
    
    # Remove all whitespace characters (spaces, tabs, newlines) and uppercase
    cleaned = "".join(sequence.split()).upper()
    return cleaned


def validate_dna(sequence: str) -> Tuple[bool, str]:
    """
    Validate whether a sequence contains only standard DNA bases (A, T, G, C).
    Rejects empty sequences and invalid characters.

    Args:
        sequence (str): Input DNA sequence string.

    Returns:
        Tuple[bool, str]: (True, "Valid DNA sequence.") if valid,
                          otherwise (False, error_message).
    """
    if not sequence or not sequence.strip():
        return False, "Sequence is empty."

    cleaned = clean_sequence(sequence)
    if not cleaned:
        return False, "Sequence is empty after cleaning."

    valid_bases = {'A', 'T', 'G', 'C'}
    invalid_chars = set(cleaned) - valid_bases

    if invalid_chars:
        invalid_str = ", ".join(sorted(invalid_chars))
        return False, f"Invalid nucleotide(s): {invalid_str}"

    return True, "Valid DNA sequence."


def sequence_statistics(sequence: str) -> SequenceStatistics:
    """
    Calculate length, nucleotide counts, percentages, GC content, and AT content.

    Args:
        sequence (str): Cleaned and validated DNA sequence.

    Returns:
        SequenceStatistics: Dataclass containing compositional statistics.
    """
    cleaned = clean_sequence(sequence)
    length = len(cleaned)

    if length == 0:
        return SequenceStatistics(0, 0, 0, 0, 0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)

    a_count = cleaned.count('A')
    t_count = cleaned.count('T')
    g_count = cleaned.count('G')
    c_count = cleaned.count('C')

    a_percent = (a_count / length) * 100
    t_percent = (t_count / length) * 100
    g_percent = (g_count / length) * 100
    c_percent = (c_count / length) * 100

    gc_percent = g_percent + c_percent
    at_percent = a_percent + t_percent

    return SequenceStatistics(
        length=length,
        a_count=a_count,
        t_count=t_count,
        g_count=g_count,
        c_count=c_count,
        a_percent=round(a_percent, 2),
        t_percent=round(t_percent, 2),
        g_percent=round(g_percent, 2),
        c_percent=round(c_percent, 2),
        gc_percent=round(gc_percent, 2),
        at_percent=round(at_percent, 2)
    )


def complement(sequence: str) -> str:
    """
    Generate the DNA complement sequence using standard base pairing (A <-> T, G <-> C).

    Args:
        sequence (str): Input DNA sequence.

    Returns:
        str: Complementary DNA sequence.
    """
    cleaned = clean_sequence(sequence)
    trans_table = str.maketrans("ATGC", "TACG")
    return cleaned.translate(trans_table)


def reverse_complement(sequence: str) -> str:
    """
    Generate the reverse complement of a DNA sequence by reversing its complement.

    Args:
        sequence (str): Input DNA sequence.

    Returns:
        str: Reverse complement DNA sequence.
    """
    comp = complement(sequence)
    return comp[::-1]


def transcribe(sequence: str) -> str:
    """
    Transcribe a DNA coding strand sequence into RNA by replacing thymine (T) with uracil (U).

    Biological Convention Assumption:
    - Input represents the 5' -> 3' coding (sense) strand.
    - Thymine (T) in DNA corresponds to Uracil (U) in RNA.

    Args:
        sequence (str): Input DNA sequence.

    Returns:
        str: Transcribed RNA sequence.
    """
    cleaned = clean_sequence(sequence)
    trans_table = str.maketrans("T", "U")
    return cleaned.translate(trans_table)


def translate_dna(sequence: str, frame: int = 0, stop_at_stop: bool = True) -> str:
    """
    Translate a DNA sequence into a protein (amino acid sequence) using the standard
    genetic code across reading frames 0, 1, or 2.

    Args:
        sequence (str): Input DNA sequence.
        frame (int): Reading frame (0, 1, or 2). Default is 0.
        stop_at_stop (bool): If True, translation halts at the first stop codon (*).
                             If False, translation continues through stop codons.

    Returns:
        str: Resulting amino acid sequence.
    """
    if frame not in (0, 1, 2):
        raise ValueError("Reading frame must be 0, 1, or 2.")

    cleaned = clean_sequence(sequence)
    protein = []

    # Iterate through codons starting from the specified frame
    for i in range(frame, len(cleaned) - 2, 3):
        codon = cleaned[i:i+3]
        amino_acid = GENETIC_CODE.get(codon, 'X')  # 'X' for unknown/invalid codon

        if amino_acid == '*' and stop_at_stop:
            break
        
        protein.append(amino_acid)

    return "".join(protein)
  
