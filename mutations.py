"""
mutations.py - Alignment-based mutation detection and summarization module
for the Bioinformatics Toolkit.

SCIENTIFIC LIMITATION NOTE:
This is an educational sequence difference caller. It is NOT clinical software,
a diagnostic system, a validated clinical variant caller, or a production genomic
pipeline. It does not perform read alignment, sequencing quality score filtering,
or biological interpretation of variant pathogenicity.
"""

from dataclasses import dataclass
from typing import Dict, List, Any


@dataclass
class MutationRecord:
    """Dataclass holding detailed information about a detected mutation."""
    alignment_position: int  # 1-based index in the aligned string
    reference_position: int  # 1-based biological coordinate on the reference sequence
    query_position: int      # 1-based biological coordinate on the query sequence
    mutation_type: str       # 'Substitution', 'Insertion', or 'Deletion'
    reference_base: str      # Nucleotide or '-'
    query_base: str          # Nucleotide or '-'


def call_mutations(aligned_reference: str, aligned_query: str) -> List[MutationRecord]:
    """
    Scan aligned reference and query sequences to detect differences (substitutions,
    insertions, and deletions) without relying on raw unaligned indices.

    Coordinate Tracking Rules:
    - Reference coordinate increments only when the reference base is a real nucleotide (!= '-').
    - Query coordinate increments only when the query base is a real nucleotide (!= '-').

    Args:
        aligned_reference (str): Aligned reference sequence string.
        aligned_query (str): Aligned query sequence string.

    Returns:
        List[MutationRecord]: List of detected mutation records.
    """
    if len(aligned_reference) != len(aligned_query):
        raise ValueError("Aligned sequences must have equal length.")

    mutations: List[MutationRecord] = []
    ref_pos = 0
    query_pos = 0

    for idx, (r_base, q_base) in enumerate(zip(aligned_reference, aligned_query), start=1):
        has_ref = (r_base != '-')
        has_query = (q_base != '-')

        if has_ref:
            ref_pos += 1
        if has_query:
            query_pos += 1

        if r_base == q_base:
            continue

        # Determine mutation type
        if has_ref and has_query:
            m_type = "Substitution"
        elif not has_ref and has_query:
            m_type = "Insertion"
        elif has_ref and not has_query:
            m_type = "Deletion"
        else:
            # Both are gaps (edge case handling)
            continue

        mutations.append(MutationRecord(
            alignment_position=idx,
            reference_position=ref_pos,
            query_position=query_pos,
            mutation_type=m_type,
            reference_base=r_base,
            query_base=q_base
        ))

    return mutations


def summarize_mutations(mutations: List[MutationRecord]) -> Dict[str, Any]:
    """
    Produce a machine-readable summary dictionary of mutation counts.

    Args:
        mutations (List[MutationRecord]): List of detected mutation records.

    Returns:
        Dict[str, Any]: Summary containing total differences and counts by type.
    """
    substitutions = sum(1 for m in mutations if m.mutation_type == "Substitution")
    insertions = sum(1 for m in mutations if m.mutation_type == "Insertion")
    deletions = sum(1 for m in mutations if m.mutation_type == "Deletion")

    return {
        "total": len(mutations),
        "substitutions": substitutions,
        "insertions": insertions,
        "deletions": deletions
    }


def format_mutation_report(mutations: List[MutationRecord]) -> str:
    """
    Format mutation records into a human-readable text report.

    Args:
        mutations (List[MutationRecord]): List of mutation records.

    Returns:
        str: Formatted report string.
    """
    if not mutations:
        return "No differences detected (sequences are identical)."

    lines = ["==================================================", "MUTATION DETECTION REPORT", "=================================================="]
    for m in mutations:
        lines.append(f"Alignment Position {m.alignment_position} (Ref Pos: {m.reference_position}, Query Pos: {m.query_position}):")
        lines.append(f"  Type:      {m.mutation_type}")
        lines.append(f"  Reference: {m.reference_base}")
        lines.append(f"  Query:     {m.query_base}")
        lines.append("-" * 50)

    summary = summarize_mutations(mutations)
    lines.append("\nSUMMARY:")
    lines.append(f"  Total differences: {summary['total']}")
    lines.append(f"  Substitutions:     {summary['substitutions']}")
    lines.append(f"  Insertions:        {summary['insertions']}")
    lines.append(f"  Deletions:         {summary['deletions']}")
    lines.append("==================================================")

    return "\n".join(lines)
  
