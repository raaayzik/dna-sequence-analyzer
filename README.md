# Bioinformatics DNA Sequence Analyzer and Alignment Toolkit

An advanced, modular educational bioinformatics toolkit implemented in Python 3.10+. It performs DNA validation, sequence analysis, transcription, translation, open reading frame (ORF) detection, global and local pairwise alignments, and alignment-based mutation detection.

---

## Features

1. **Sequence Foundation**: Robust DNA validation, cleaning, statistics (GC/AT content), complement, reverse complement, transcription, and translation (with 64-codon table across frames 0, 1, 2).
2. **FASTA Parser & Writer**: Handles multiline sequences, blank lines, metadata descriptions, and file round-trips.
3. **Open Reading Frame (ORF) Finder**: Detects ORFs across all 3 forward and 3 reverse-complement reading frames with configurable length thresholds.
4. **Pairwise Sequence Alignment**: Implements **Needleman-Wunsch** (global) and **Smith-Waterman** (local) dynamic programming algorithms from scratch.
5. **Mutation Detection**: Identifies substitutions, insertions, and deletions directly from alignment columns.
6. **Protein Analysis**: Computes amino acid composition and frequency metrics.
7. **Robust CLI & Interactive Menu**: Supports both an interactive terminal menu and command-line arguments via `argparse`.
8. **Export Utilities**: Export results to TXT, JSON, and CSV formats.

---

## Project Structure

```text
bioinformatics_sequence_analyzer/
├── main.py
├── sequence.py
├── alignment.py
├── mutations.py
├── fasta.py
├── orf.py
├── protein.py
├── io_utils.py
├── requirements.txt
├── README.md
├── data/
│   ├── reference.fasta
│   ├── sample_sequences.fasta
│   ├── orf_examples.fasta
│   └── mutation_examples.fasta
└── tests/
    ├── test_sequence.py
    ├── test_alignment.py
    ├── test_mutations.py
    ├── test_fasta.py
    ├── test_orf.py
    ├── test_protein.py
    └── test_io_utils.py
    
