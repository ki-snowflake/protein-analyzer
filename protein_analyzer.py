#!/usr/bin/env python3
"""
    A tool for analyzing the amino acid composition of proteins using
    different concurrency approaches (single-threaded, multi-threaded,
    and multi-processed) to benchmark Python's performance under GIL constraints.

    Usage:
        python protein_analyzer.py <fasta_10> <fasta_1000> <fasta_10k> 
"""

import multiprocessing as mp
import os
import sys
import threading
import time
from collections import defaultdict
from queue import Empty, Queue
from typing import Dict, List
import matplotlib.pyplot as plt
from Bio import SeqIO


class ProteinAnalyzer:

    HYDROPHOBIC = set("AILMFWV")
    NEUTRAL = set("NQSTYC")
    POSITIVE = set("KRH")
    NEGATIVE = set("DE")

    def __init__(self) -> None:
        """Initializes the thread-safe lock for multi-threaded aggregation."""
        self.lock = threading.Lock()

    @staticmethod
    def analyze_sequence(sequence: str) -> Dict[str, int]:
        """Analyzes the amino acid composition of a single protein sequence.

        Args:
            sequence (str): The amino acid sequence to be analyzed.

        Returns:
            Dict[str, int]: A dictionary containing counts for 'hydrophobic',
                            'neutral', 'positive', and 'negative' residues.
        """
        local_stats = {
            "hydrophobic": 0,
            "neutral": 0,
            "positive": 0,
            "negative": 0,
        }

        for aa in sequence.upper():
            if aa in ProteinAnalyzer.HYDROPHOBIC:
                local_stats["hydrophobic"] += 1
            elif aa in ProteinAnalyzer.NEUTRAL:
                local_stats["neutral"] += 1
            elif aa in ProteinAnalyzer.POSITIVE:
                local_stats["positive"] += 1
            elif aa in ProteinAnalyzer.NEGATIVE:
                local_stats["negative"] += 1

        return local_stats

    def process_file_single(self, fasta_file: str) -> Dict[str, int]:
        """Processes a FASTA file sequentially using a single thread.

        Args:
            fasta_file (str): The path to the input FASTA file.

        Returns:
            Dict[str, int]: Aggregated amino acid statistics for the entire file.
        """
        aggregated_stats: defaultdict[str, int] = defaultdict(int)
        for record in SeqIO.parse(fasta_file, "fasta"):
            seq_stats = ProteinAnalyzer.analyze_sequence(str(record.seq))
            for key, value in seq_stats.items():
                aggregated_stats[key] += value
        return dict(aggregated_stats)

    def _worker_thread(self, queue: Queue, global_stats: defaultdict[str, int]) -> None:
        """Target function for worker threads that consumes sequences from a queue.

        Args:
            queue (Queue): Thread-safe queue containing protein sequences.
            global_stats (defaultdict[str, int]): Shared dictionary to aggregate results.
        """
        local_stats: defaultdict[str, int] = defaultdict(int)

        while True:
            try:
                seq = queue.get_nowait()
            except Empty:
                break

            seq_stats = ProteinAnalyzer.analyze_sequence(seq)
            for key, value in seq_stats.items():
                local_stats[key] += value
            queue.task_done()

        # Thread-safe merge using the instance lock
        with self.lock:
            for key, value in local_stats.items():
                global_stats[key] += value

    def process_file_threaded(
        self, fasta_file: str, num_threads: int = 4
    ) -> Dict[str, int]:
        """Processes a FASTA file in parallel using multiple threads.

        Args:
            fasta_file (str): The path to the input FASTA file.
            num_threads (int): The number of working threads to spawn.

        Returns:
            Dict[str, int]: Aggregated amino acid statistics.
        """
        sequence_queue: Queue = Queue()
        global_stats: defaultdict[str, int] = defaultdict(int)

        for record in SeqIO.parse(fasta_file, "fasta"):
            sequence_queue.put(str(record.seq))

        threads = []
        for _ in range(num_threads):
            t = threading.Thread(
                target=self._worker_thread,
                args=(sequence_queue, global_stats),
            )
            threads.append(t)
            t.start()

        for t in threads:
            t.join()

        return dict(global_stats)

    def process_file_multiprocess(
        self, fasta_file: str, num_processes: int = 4
    ) -> Dict[str, int]:
        """Processes a FASTA file in parallel across multiple CPU cores
        using a multiprocessing Pool to bypass GIL.

        Args:
            fasta_file (str): The path to the input FASTA file.
            num_processes (int): The number of independent processes to spawn.

        Returns:
            Dict[str, int]: Aggregated amino acid statistics.
        """
        sequences = [
            str(record.seq) for record in SeqIO.parse(fasta_file, "fasta")
        ]

        # Running the pure static method prevents object-level serialization errors
        with mp.Pool(processes=num_processes) as pool:
            results = pool.map(ProteinAnalyzer.analyze_sequence, sequences)

        global_stats: defaultdict[str, int] = defaultdict(int)
        for seq_stats in results:
            for key, value in seq_stats.items():
                global_stats[key] += value

        return dict(global_stats)

    @staticmethod
    def plot_performance(
        reads: List[int],
        times_single: List[float],
        times_threaded: List[float],
        times_mp: List[float],
        num_workers: int,
    ) -> None:
        """Generates and saves a line plot comparing the execution times of
        single-threaded, multi-threaded, and multi-processed execution modes.
        """
        plt.figure(figsize=(10, 6))
        plt.plot(reads, times_single, "o-", label="Single-threaded")
        plt.plot(
            reads,
            times_threaded,
            "s-",
            label=f"Multi-threaded ({num_workers} threads)",
        )
        plt.plot(
            reads, times_mp, "^-", label=f"Multi-process ({num_workers} cores)"
        )

        plt.xlabel("Number of Reads (Sequences)")
        plt.ylabel("Execution Time (seconds)")
        plt.title("Performance Comparison of Concurrent Approaches")
        plt.legend()
        plt.grid(True, linestyle="--", alpha=0.6)
        plt.savefig("performance_comparison.png", dpi=300)
        plt.close()
        print(
            "\n[Success] Performance benchmark chart saved as 'performance_comparison.png'"
        )


if __name__ == "__main__":
    # Validate command line arguments
    if len(sys.argv) < 4:
        print(
            "Error: Insufficient arguments.\n"
            "Usage: python protein_analyzer.py <fasta_10> <fasta_1000> <fasta_10k>"
        )
        sys.exit(1)

    fasta_files = sys.argv[1:4]
    reads_labels = [100, 1000, 10000]

    # Verification of file existence
    for file_path in fasta_files:
        if not os.path.exists(file_path):
            print(f"Error: File '{file_path}' does not exist.")
            sys.exit(1)

    analyzer = ProteinAnalyzer()

    # Lists to store benchmarking metrics
    times_single = []
    times_threaded = []
    times_mp = []

    workers_count = 4

    # Execute benchmarks
    for fasta_file, read_count in zip(fasta_files, reads_labels):
        print(f"\n" + "=" * 50)
        print(f"Benchmarking dataset with {read_count} reads:")
        print("=" * 50)

        # Single-threaded
        start = time.time()
        results_single = analyzer.process_file_single(fasta_file)
        elapsed_single = time.time() - start
        times_single.append(elapsed_single)
        print(f"-> Single-threaded finished in: {elapsed_single:.4f} seconds")

        print("   Amino acid distribution in this file:")
        print(f"     - Hydrophobic: {results_single.get('hydrophobic', 0)}")
        print(f"     - Neutral:     {results_single.get('neutral', 0)}")
        print(f"     - Positive:    {results_single.get('positive', 0)}")
        print(f"     - Negative:    {results_single.get('negative', 0)}")
        print(f"     - Total typed: {sum(results_single.values())}")

        # Multi-threaded
        start = time.time()
        results_threaded = analyzer.process_file_threaded(
            fasta_file, num_threads=workers_count
        )
        elapsed_threaded = time.time() - start
        times_threaded.append(elapsed_threaded)
        print(f"-> Multi-threaded finished in:  {elapsed_threaded:.4f} seconds")

        # Multi-process
        start = time.time()
        results_mp = analyzer.process_file_multiprocess(
            fasta_file, num_processes=workers_count
        )
        elapsed_mp = time.time() - start
        times_mp.append(elapsed_mp)
        print(f"-> Multi-processed finished in: {elapsed_mp:.4f} seconds")

    # To save the chart
    ProteinAnalyzer.plot_performance(
        reads_labels, times_single, times_threaded, times_mp, workers_count
    )
