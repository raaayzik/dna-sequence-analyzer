============================================================
STAGE 5 — COMMAND-LINE APPLICATION AND REPORTING
============================================================

Stages 1–4 are complete and tested.

Now convert the underlying modules into a usable
bioinformatics command-line application.

Do NOT rewrite the existing algorithms unnecessarily.

Reuse the existing modules.


============================================================
1. PROJECT FILES
============================================================

Create:

    main.py
    io_utils.py

Update:

    README.md
    requirements.txt

Create:

    tests/test_io_utils.py


============================================================
2. INTERACTIVE MENU
============================================================

When the program is launched without arguments:

    python main.py

display:

==================================================
BIOINFORMATICS DNA SEQUENCE ANALYZER
==================================================

1. Analyze DNA sequence
2. Global alignment
3. Local alignment
4. Mutation analysis
5. Find ORFs
6. Analyze FASTA file
7. Protein analysis
8. Export results
9. Exit


Make the menu robust.


============================================================
3. COMMAND-LINE ARGUMENTS
============================================================

Use argparse.

Support commands such as:

    python main.py analyze data/sample_sequences.fasta

    python main.py align reference.fasta query.fasta

    python main.py local-align reference.fasta query.fasta

    python main.py orf data/orf_examples.fasta

    python main.py mutations reference.fasta query.fasta


Provide:

    python main.py --help

with useful documentation.


============================================================
4. FASTA BATCH ANALYSIS
============================================================

For:

    analyze file.fasta

report for each sequence:

- name
- length
- A/T/G/C counts
- GC%
- AT%
- reverse complement
- number of ORFs
- longest ORF


============================================================
5. GLOBAL ALIGNMENT COMMAND
============================================================

For:

    align reference.fasta query.fasta

perform Needleman-Wunsch alignment.

Report:

- reference name
- query name
- score
- alignment
- identity
- matches
- mismatches
- gaps


============================================================
6. LOCAL ALIGNMENT COMMAND
============================================================

For:

    local-align reference.fasta query.fasta

perform Smith-Waterman.

Report:

- score
- reference coordinates
- query coordinates
- local alignment
- identity


============================================================
7. MUTATION COMMAND
============================================================

For:

    mutations reference.fasta query.fasta

Perform:

1. Read reference.
2. Read query.
3. Validate DNA.
4. Perform global alignment.
5. Call mutations.
6. Print mutation report.
7. Print mutation summary.


============================================================
8. ORF COMMAND
============================================================

For:

    orf sequence.fasta

report all ORFs.

Include:

- strand
- frame
- coordinates
- nucleotide length
- DNA
- protein


============================================================
9. EXPORT
============================================================

Implement export functionality.

Support:

    TXT
    JSON
    CSV

Use standard-library modules where possible.

JSON should preserve structured information.

CSV should be appropriate for tabular data such as
mutation lists or sequence summaries.


============================================================
10. IO UTILITIES
============================================================

Create reusable functions such as:

    write_json(...)
    write_csv(...)
    write_text(...)

Do not put file-export logic inside biological algorithms.


============================================================
11. ERROR HANDLING
============================================================

The CLI must handle:

- missing files
- malformed FASTA
- invalid DNA
- missing arguments
- invalid commands
- empty FASTA files
- incompatible input files

Errors should be understandable to a normal user.


============================================================
12. TESTS
============================================================

Create:

    tests/test_io_utils.py

Test:

- JSON writing
- CSV writing
- text writing
- valid paths
- basic data structures

Also ensure existing tests continue to pass.


============================================================
13. NO UNNECESSARY DEPENDENCIES
============================================================

Prefer standard-library Python.

Do not introduce external dependencies unless genuinely
necessary.


============================================================
14. CODE QUALITY
============================================================

Keep:

- algorithm modules independent
- CLI code separate
- file IO separate
- biological logic reusable

The CLI should call functions from the existing modules
rather than reimplementing them.


============================================================
15. README UPDATE
============================================================

Update README.md with:

- installation
- interactive mode
- command-line mode
- examples
- export examples
- test instructions


============================================================
16. TESTING
============================================================

Run:

    python -m unittest discover -s tests -v

Then manually test:

    python main.py --help

and at least three real commands using the supplied test data.


============================================================
17. FINAL RESPONSE
============================================================

Provide:

1. Complete main.py
2. Complete io_utils.py
3. Complete tests/test_io_utils.py
4. Updated README.md
5. Updated requirements.txt
6. Any necessary changes to existing modules
7. Example CLI commands
8. Example output
9. Test results

Stop after Stage 5.
