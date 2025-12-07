# helix-hunter

#  DNA Pattern Matching and Alignment Tool

##  The Challenge & Motivation

DNA sequence analysis is foundational to bioinformatics, yet exact pattern matching in massive biological data can be computationally intensive. Researchers often need to identify specific gene sequences within larger genomes to study genetic markers, mutations, or conserved regions.

This project addresses that challenge by combining classic computer science algorithms with biological data processing — providing a fast, benchmarked, and extensible tool to perform DNA pattern matching.

---

#  Overview

This tool allows users to:

- Search for a specific DNA pattern in a `.fasta` or `.txt` file using:
  -  **Naive Brute Force**
  -  **Knuth-Morris-Pratt (KMP)**
  -  **Rabin-Karp (Rolling Hash)**
- Benchmark these algorithms by measuring:
  - Execution time
  - Number of matches found
- Optionally perform **local alignment** to extend matches, inspired by the **BLAST** (Basic Local Alignment Search Tool) methodology.

---

#  Execution
1. Command Line Interface
   
Run from terminal with the following options:

Basic usage
```Bash
python main.py --file <your_fasta_file.fasta or txt file> --pattern <ATCG> --algorithm <kmp>
```
Run all algorithms with plots
```Bash
python main.py --file test.txt/sample.fasta --pattern ATCG --algorithm all --plot
```
Enable local alignment (BLAST-inspired)
``` Bash
python main.py --file test_data/sample.fasta --pattern ATCG --extend --window 8
```

Available algorithms:

* `naive`
* `kmp`
* `rk`
* `all`

2. Gradio Web 


Launch the UI by simply running:

```bash
python main.py
```
This opens an interactive browser-based interface to upload FASTA files, enter a pattern, and select the matching algorithm.






