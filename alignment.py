"""
alignment.py - Needleman-Wunsch global alignment and Smith-Waterman local alignment
modules from scratch for the Bioinformatics Toolkit.

Complexity:
    Time: O(n * m) where n and m are sequence lengths.
    Space: O(n * m) for storing the dynamic programming matrices.
"""

from dataclasses import dataclass
from typing import Dict, List, Tuple
from sequence import clean_sequence, validate_dna


@dataclass
class AlignmentResult:
    """Dataclass holding detailed metrics and aligned sequences for pairwise alignment."""
    score: float
    aligned_reference: str
    aligned_query: str
    start_reference: int
    end_reference: int
    start_query: int
    end_query: int
    matches: int
    mismatches: int
    gaps: int
    alignment_length: int
    identity: float


def calculate_alignment_stats(ref_aligned: str, query_aligned: str) -> Dict[str, any]:
    """
    Calculate alignment statistics including matches, mismatches, gaps, length,
    and percent identity.

    Percent Identity Formula:
        identity = (matches / aligned non-gap positions) * 100
    """
    matches = 0
    mismatches = 0
    gaps = 0
    non_gap_positions = 0

    for r, q in zip(ref_aligned, query_aligned):
        if r == '-' or q == '-':
            gaps += 1
        else:
            non_gap_positions += 1
            if r == q:
                matches += 1
            else:
                mismatches += 1

    alignment_length = len(ref_aligned)
    identity = (matches / non_gap_positions * 100) if non_gap_positions > 0 else 0.0

    return {
        "matches": matches,
        "mismatches": mismatches,
        "gaps": gaps,
        "alignment_length": alignment_length,
        "identity": round(identity, 2)
    }


def needleman_wunsch(
    reference: str,
    query: str,
    match_score: int = 2,
    mismatch_penalty: int = -1,
    gap_penalty: int = -2
) -> AlignmentResult:
    """
    Perform Needleman-Wunsch global sequence alignment from scratch using dynamic programming.

    Recurrence relation:
        diagonal = previous_diagonal + (match_score or mismatch_penalty)
        up       = previous_up + gap_penalty
        left     = previous_left + gap_penalty
        F(i, j)  = max(diagonal, up, left)

    Tie-breaking rule:
        Diagonal > Up > Left (deterministic).

    Args:
        reference (str): Reference sequence.
        query (str): Query sequence.
        match_score (int): Score for matching bases (default: 2).
        mismatch_penalty (int): Penalty for mismatching bases (default: -1).
        gap_penalty (int): Penalty for introducing or extending a gap (default: -2).

    Returns:
        AlignmentResult: Alignment score, aligned sequences, coordinates, and statistics.
    """
    ref_clean = clean_sequence(reference)
    query_clean = clean_sequence(query)

    valid_r, msg_r = validate_dna(ref_clean)
    if not valid_r:
        raise ValueError(f"Invalid reference DNA for alignment: {msg_r}")
    valid_q, msg_q = validate_dna(query_clean)
    if not valid_q:
        raise ValueError(f"Invalid query DNA for alignment: {msg_q}")

    n = len(ref_clean)
    m = len(query_clean)

    # DP matrix: (n + 1) x (m + 1)
    dp = [[0.0] * (m + 1) for _ in range(n + 1)]
    tb = [[''] * (m + 1) for _ in range(n + 1)]

    # Initialize first row and column with gap penalties
    for i in range(n + 1):
        dp[i][0] = i * gap_penalty
        tb[i][0] = 'U'
    for j in range(m + 1):
        dp[0][j] = j * gap_penalty
        tb[0][j] = 'L'
    tb[0][0] = 'NONE'

    # Fill DP matrix
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            match_mismatch = match_score if ref_clean[i-1] == query_clean[j-1] else mismatch_penalty
            diag_score = dp[i-1][j-1] + match_mismatch
            up_score = dp[i-1][j] + gap_penalty
            left_score = dp[i][j-1] + gap_penalty

            best_score = max(diag_score, up_score, left_score)
            dp[i][j] = best_score

            # Deterministic tie-breaking: Diagonal > Up > Left
            if best_score == diag_score:
                tb[i][j] = 'D'
            elif best_score == up_score:
                tb[i][j] = 'U'
            else:
                tb[i][j] = 'L'

    # Traceback from bottom-right cell (n, m)
    i, j = n, m
    aligned_ref = []
    aligned_query = []

    while i > 0 or j > 0:
        direction = tb[i][j]
        if direction == 'D':
            aligned_ref.append(ref_clean[i-1])
            aligned_query.append(query_clean[j-1])
            i -= 1
            j -= 1
        elif direction == 'U':
            aligned_ref.append(ref_clean[i-1])
            aligned_query.append('-')
            i -= 1
        elif direction == 'L':
            aligned_ref.append('-')
            aligned_query.append(query_clean[j-1])
            j -= 1
        else:
            if i > 0:
                aligned_ref.append(ref_clean[i-1])
                aligned_query.append('-')
                i -= 1
            elif j > 0:
                aligned_ref.append('-')
                aligned_query.append(query_clean[j-1])
                j -= 1

    aligned_ref.reverse()
    aligned_query.reverse()

    ref_str = "".join(aligned_ref)
    query_str = "".join(aligned_query)

    stats = calculate_alignment_stats(ref_str, query_str)

    return AlignmentResult(
        score=dp[n][m],
        aligned_reference=ref_str,
        aligned_query=query_str,
        start_reference=1,
        end_reference=n,
        start_query=1,
        end_query=m,
        matches=stats["matches"],
        mismatches=stats["mismatches"],
        gaps=stats["gaps"],
        alignment_length=stats["alignment_length"],
        identity=stats["identity"]
    )


def smith_waterman(
    reference: str,
    query: str,
    match_score: int = 2,
    mismatch_penalty: int = -1,
    gap_penalty: int = -2
) -> AlignmentResult:
    """
    Perform Smith-Waterman local sequence alignment from scratch using dynamic programming.

    Recurrence relation (with zero reset):
        F(i, j) = max(0, diagonal, up, left)

    Why zero is used:
        Local alignment searches for regions of high similarity. Resetting negative
        scores to zero allows the algorithm to start fresh, ignoring poorly matching
        flanking regions.

    Args:
        reference (str): Reference sequence.
        query (str): Query sequence.
        match_score (int): Score for matching bases (default: 2).
        mismatch_penalty (int): Penalty for mismatching bases (default: -1).
        gap_penalty (int): Penalty for introducing or extending a gap (default: -2).

    Returns:
        AlignmentResult: Maximum score, local alignment, start/end coordinates, and statistics.
    """
    ref_clean = clean_sequence(reference)
    query_clean = clean_sequence(query)

    valid_r, msg_r = validate_dna(ref_clean)
    if not valid_r:
        raise ValueError(f"Invalid reference DNA for alignment: {msg_r}")
    valid_q, msg_q = validate_dna(query_clean)
    if not valid_q:
        raise ValueError(f"Invalid query DNA for alignment: {msg_q}")

    n = len(ref_clean)
    m = len(query_clean)

    dp = [[0.0] * (m + 1) for _ in range(n + 1)]
    tb = [[''] * (m + 1) for _ in range(n + 1)]

    max_score = 0.0
    max_i, max_j = 0, 0

    # Fill DP matrix with zero floor
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            match_mismatch = match_score if ref_clean[i-1] == query_clean[j-1] else mismatch_penalty
            diag_score = dp[i-1][j-1] + match_mismatch
            up_score = dp[i-1][j] + gap_penalty
            left_score = dp[i][j-1] + gap_penalty

            best_score = max(0.0, diag_score, up_score, left_score)
            dp[i][j] = best_score

            if best_score > max_score:
                max_score = best_score
                max_i, max_j = i, j

            if best_score == 0.0:
                tb[i][j] = 'STOP'
            elif best_score == diag_score:
                tb[i][j] = 'D'
            elif best_score == up_score:
                tb[i][j] = 'U'
            else:
                tb[i][j] = 'L'

    # Traceback from the highest scoring cell until zero is reached
    i, j = max_i, max_j
    end_ref = i
    end_query = j

    aligned_ref = []
    aligned_query = []

    while i > 0 and j > 0 and dp[i][j] > 0 and tb[i][j] != 'STOP':
        direction = tb[i][j]
        if direction == 'D':
            aligned_ref.append(ref_clean[i-1])
            aligned_query.append(query_clean[j-1])
            i -= 1
            j -= 1
        elif direction == 'U':
            aligned_ref.append(ref_clean[i-1])
            aligned_query.append('-')
            i -= 1
        elif direction == 'L':
            aligned_ref.append('-')
            aligned_query.append(query_clean[j-1])
            j -= 1
        else:
            break

    start_ref = i + 1
    start_query = j + 1

    aligned_ref.reverse()
    aligned_query.reverse()

    ref_str = "".join(aligned_ref)
    query_str = "".join(aligned_query)

    stats = calculate_alignment_stats(ref_str, query_str)

    return AlignmentResult(
        score=max_score,
        aligned_reference=ref_str,
        aligned_query=query_str,
        start_reference=start_ref if ref_str else 0,
        end_reference=end_ref if ref_str else 0,
        start_query=start_query if query_str else 0,
        end_query=end_query if query_str else 0,
        matches=stats["matches"],
        mismatches=stats["mismatches"],
        gaps=stats["gaps"],
        alignment_length=stats["alignment_length"],
        identity=stats["identity"]
    )


def format_alignment(aligned_reference: str, aligned_query: str, width: int = 60) -> str:
    """
    Format aligned reference and query sequences into a readable multi-line text representation
    with match indicator lines (|).

    Args:
        aligned_reference (str): Aligned reference sequence.
        aligned_query (str): Aligned query sequence.
        width (int): Line wrapping character width (default: 60).

    Returns:
        str: Formatted text representation.
    """
    lines = []
    n = len(aligned_reference)

    for i in range(0, n, width):
        ref_chunk = aligned_reference[i:i+width]
        query_chunk = aligned_query[i:i+width]

        match_line = []
        for r, q in zip(ref_chunk, query_chunk):
            if r == q and r != '-':
                match_line.append('|')
            else:
                match_line.append(' ')
        match_str = "".join(match_line)

        lines.append(f"Reference: {ref_chunk}")
        lines.append(f"           {match_str}")
        lines.append(f"Query:     {query_chunk}\n")

    return "\n".join(lines).strip()
  
