# Concurrent Protein Sequence Analyzer

A tool designed to analyze the amino acid composition of protein sequences from FASTA files. This project serves as a performance benchmark comparing **single-threaded**, **multi-threaded (low-level `threading` with `Queue` and `Lock`)**, and **multi-processed (`multiprocessing.Pool`)** execution modes under Python's Global Interpreter Lock (GIL) constraints.

## 📌 Project Overview

In bioinformatics, parsing large genomic or proteomic datasets (like FASTA files) can quickly become a computational bottleneck. This project explores how different concurrency models in Python handle CPU-bound tasks:

1. **Single-threaded:** Processes sequences sequentially.
2. **Multi-threaded:** Spawns low-level worker threads using a thread-safe `queue.Queue` for task distribution and a `threading.Lock` to aggregate metrics without race conditions.
3. **Multi-processed:** Utilizes a pool of independent processes to bypass the Python GIL.

Here, it automatically tracks execution benchmarks across variable data scales (e.g., 100, 1.000, and 10.000 sequences) and exports a chart (`performance_comparison.png`). The terminal output provides a comprehensive summary for each dataset, detailing the precise distribution of amino acid properties alongside total residue counts. Additionally, it logs the execution speed of each model, offering clear insights into hardware utilization.

## 🚀 Usage

To run the performance benchmark, execute the script from the terminal and pass the paths to your FASTA data files as args:

```bash
python .\protein-analyzer.py test1.100Reads test2.1000Reads test3.10kReads
```

## 🧬 Biological Context

Proteins are translated from FASTA sequences and classified dynamically based on the core physicochemical properties of their residues. Classification is taken from Pevsner J. Bioinformatics and Functional Genomics, 2009:

* **Hydrophobic:** A, I, L, M, F, W, V
* **Neutral:** N, Q, S, T, Y, C
* **Positive (Basic):** K, R, H
* **Negative (Acidic):** D, E

<img width="1071" height="1266" alt="image" src="https://github.com/user-attachments/assets/806d8a89-275a-48d8-866b-cde75a8c7b63" />
