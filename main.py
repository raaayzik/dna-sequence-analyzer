"""
main.py - Command-Line Interface (CLI) and interactive menu controller
for the Bioinformatics DNA Sequence Analyzer and Alignment Toolkit.
"""

import argparse
import sys
import os

from sequence import clean_sequence, validate_dna, sequence_statistics, reverse_complement, translate_dna
from fasta import read_fasta, FastaRecord
from orf import find_orfs
from alignment import needleman_wunsch, smith_waterman, format_alignment
from mutations import call_mutations, summarize_mutations, format_mutation_report
from protein import analyze_protein
from io_utils import write_text, write_json, write_csv


def interactive_menu():
    """Run the interactive command-line menu."""
    while True:
        print("\n" + "=" * 50)
        print("BIOINFORMATICS DNA SEQUENCE ANALYZER")
        print("=" * 50)
        print("1. Analyze DNA sequence")
        print("2. Global alignment")
        print("3. Local alignment")
        print("4. Mutation analysis")
        print("5. Find ORFs")
        print("6. Analyze FASTA file")
        print("7. Protein analysis")
        print("8. Export results")
        print("9. Exit")
        print("=" * 50)

        choice = input("Select an option (1-9): ").strip()

        try:
            if choice == "1":
                seq = input("Enter DNA sequence: ").strip()
                valid, msg = validate_dna(seq)
                if not valid:
                    print(f"Error: {msg}")
                    continue
                stats = sequence_statistics(seq)
                print("\n--- SEQUENCE STATISTICS ---")
                for k, v in stats.__dict__.items():
                    print(f"  {k}: {v}")

            elif choice == "2":
                ref = input("Enter reference DNA sequence: ").strip()
                query = input("Enter query DNA sequence: ").strip()
                res = needleman_wunsch(ref, query)
                print(f"\nScore: {res.score} | Identity: {res.identity}% | Matches: {res.matches} | Gaps: {res.gaps}")
                print("\nAlignment:\n" + format_alignment(res.aligned_reference, res.aligned_query))

            elif choice == "3":
                ref = input("Enter reference DNA sequence: ").strip()
                query = input("Enter query DNA sequence: ").strip()
                res = smith_waterman(ref, query)
                print(f"\nScore: {res.score} | Identity: {res.identity}%")
                print(f"Ref Coordinates: {res.start_reference}-{res.end_reference}")
                print(f"Query Coordinates: {res.start_query}-{res.end_query}")
                print("\nLocal Alignment:\n" + format_alignment(res.aligned_reference, res.aligned_query))

            elif choice == "4":
                ref = input("Enter reference DNA sequence: ").strip()
                query = input("Enter query DNA sequence: ").strip()
                res = needleman_wunsch(ref, query)
                muts = call_mutations(res.aligned_reference, res.aligned_query)
                print("\n" + format_mutation_report(muts))

            elif choice == "5":
                seq = input("Enter DNA sequence: ").strip()
                valid, msg = validate_dna(seq)
                if not valid:
                    print(f"Error: {msg}")
                    continue
                min_len = int(input("Enter minimum ORF length in nucleotides [default 30]: ") or "30")
                orfs = find_orfs(seq, min_length=min_len)
                print(f"\nFound {len(orfs)} ORF(s):")
                for idx, o in enumerate(orfs, 1):
                    print(f"[{idx}] Strand: {o.strand}, Frame: {o.frame}, Coords: {o.start}-{o.end}, Len: {o.length}nt")
                    print(f"    Protein: {o.protein_sequence[:50]}...")

            elif choice == "6":
                filepath = input("Enter path to FASTA file: ").strip()
                records = read_fasta(filepath)
                print(f"\nLoaded {len(records)} record(s) from {filepath}:")
                for rec in records:
                    stats = sequence_statistics(rec.sequence)
                    orfs = find_orfs(rec.sequence)
                    print(f"\n  Name: {rec.name} | Len: {stats.length} | GC%: {stats.gc_percent}% | ORFs: {len(orfs)}")

            elif choice == "7":
                p_seq = input("Enter amino acid sequence: ").strip()
                p_stats = analyze_protein(p_seq)
                print("\n--- PROTEIN STATISTICS ---")
                print(f"  Length: {p_stats.length}")
                print(f"  Composition: {p_stats.composition}")
                print(f"  Percentages: {p_stats.percentages}")
                print(f"  Most Common: {p_stats.most_common[:3]}")

            elif choice == "8":
                print("\nExport feature available via CLI commands or direct script usage.")

            elif choice == "9":
                print("Exiting. Goodbye!")
                sys.exit(0)
            else:
                print("Invalid option. Please enter a number between 1 and 9.")
        except Exception as e:
            print(f"\nAn error occurred: {e}")


def main():
    parser = argparse.ArgumentParser(
        description="Bioinformatics DNA Sequence Analyzer and Alignment Toolkit"
    )
    subparsers = parser.add_subparsers(dest="command", help="Sub-commands")

    # Analyze command
    p_analyze = subparsers.add_parser("analyze", help="Batch analyze sequences in a FASTA file")
    p_analyze.add_argument("fasta_file", type=str, help="Path to input FASTA file")

    # Align command
    p_align = subparsers.add_parser("align", help="Perform Needleman-Wunsch global alignment")
    p_align.add_argument("ref_file", type=str, help="Reference FASTA file")
    p_align.add_argument("query_file", type=str, help="Query FASTA file")

    # Local align command
    p_local = subparsers.add_parser("local-align", help="Perform Smith-Waterman local alignment")
    p_local.add_argument("ref_file", type=str, help="Reference FASTA file")
    p_local.add_argument("query_file", type=str, help="Query FASTA file")

    # ORF command
    p_orf = subparsers.add_parser("orf", help="Find ORFs in a FASTA file")
    p_orf.add_argument("fasta_file", type=str, help="Path to input FASTA file")
    p_orf.add_argument("--min-length", type=int, default=30, help="Minimum ORF length in nucleotides")

    # Mutations command
    p_mut = subparsers.add_parser("mutations", help="Perform alignment and call mutations")
    p_mut.add_argument("ref_file", type=str, help="Reference FASTA file")
    p_mut.add_argument("query_file", type=str, help="Query FASTA file")

    args = parser.parse_args()

    if not args.command:
        interactive_menu()
        return

    try:
        if args.command == "analyze":
            records = read_fasta(args.fasta_file)
            print(f"Batch Analysis Report for: {args.fasta_file}\n" + "=" * 60)
            for rec in records:
                stats = sequence_statistics(rec.sequence)
                orfs = find_orfs(rec.sequence)
                longest_orf = max((o.length for o in orfs), default=0)
                print(f"Sequence: {rec.name}")
                print(f"  Length: {stats.length}")
                print(f"  A/T/G/C: A={stats.a_count}, T={stats.t_count}, G={stats.g_count}, C={stats.c_count}")
                print(f"  GC%: {stats.gc_percent}% | AT%: {stats.at_percent}%")
                print(f"  Reverse Complement: {reverse_complement(rec.sequence)[:30]}...")
                print(f"  Number of ORFs: {len(orfs)} | Longest ORF: {longest_orf}nt")
                print("-" * 60)

        elif args.command == "align":
            ref_recs = read_fasta(args.ref_file)
            query_recs = read_fasta(args.query_file)
            ref_seq = ref_recs[0].sequence
            query_seq = query_recs[0].sequence

            res = needleman_wunsch(ref_seq, query_seq)
            print("GLOBAL ALIGNMENT REPORT (Needleman-Wunsch)")
            print("=" * 60)
            print(f"Reference Name: {ref_recs[0].name}")
            print(f"Query Name:     {query_recs[0].name}")
            print(f"Score:          {res.score}")
            print(f"Identity:       {res.identity}%")
            print(f"Matches:        {res.matches} | Mismatches: {res.mismatches} | Gaps: {res.gaps}")
            print("\nAlignment:\n" + format_alignment(res.aligned_reference, res.aligned_query))

        elif args.command == "local-align":
            ref_recs = read_fasta(args.ref_file)
            query_recs = read_fasta(args.query_file)
            ref_seq = ref_recs[0].sequence
            query_seq = query_recs[0].sequence

            res = smith_waterman(ref_seq, query_seq)
            print("LOCAL ALIGNMENT REPORT (Smith-Waterman)")
            print("=" * 60)
            print(f"Score:              {res.score}")
            print(f"Identity:           {res.identity}%")
            print(f"Reference Coords:   {res.start_reference} - {res.end_reference}")
            print(f"Query Coords:       {res.start_query} - {res.end_query}")
            print("\nLocal Alignment:\n" + format_alignment(res.aligned_reference, res.aligned_query))

        elif args.command == "orf":
            records = read_fasta(args.fasta_file)
            for rec in records:
                orfs = find_orfs(rec.sequence, min_length=args.min_length)
                print(f"ORF Report for: {rec.name} (Total: {len(orfs)})\n" + "=" * 60)
                for o in orfs:
                    print(f"Strand: {o.strand} | Frame: {o.frame} | Coords: {o.start}-{o.end} | Len: {o.length}")
                    print(f"DNA:     {o.dna_sequence}")
                    print(f"Protein: {o.protein_sequence}")
                    print("-" * 60)

        elif args.command == "mutations":
            ref_recs = read_fasta(args.ref_file)
            query_recs = read_fasta(args.query_file)
            res = needleman_wunsch(ref_recs[0].sequence, query_recs[0].sequence)
            muts = call_mutations(res.aligned_reference, res.aligned_query)
            print(format_mutation_report(muts))

    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
              
