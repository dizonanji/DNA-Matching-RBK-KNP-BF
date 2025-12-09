"""
DNA Pattern Matching Benchmark using Rabin-Karp only
Empirical analysis on three genomes: DNA sequence, chimpanzee, and dog.
This version cleans the input sequence so matcher.py only gets A/C/G/T uppercase.
"""

from typing import List, Dict, Tuple
from benchmark import Benchmark
from matcher import PatternMatcher
import statistics
import matplotlib.pyplot as plt
import sys

ALLOWED_BASES = set("ACGT")


def load_and_clean_sequence(file_path: str) -> str:
    """
    Read a DNA file, remove FASTA headers and whitespace, uppercase,
    and filter out any character that is not A/C/G/T.
    Returns the cleaned sequence.
    """
    raw_chunks = []
    try:
        with open(file_path, "r") as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                if line.startswith(">"):   # skip FASTA headers
                    continue
                raw_chunks.append(line)
    except FileNotFoundError:
        print(f"Error: file not found: {file_path}", file=sys.stderr)
        raise

    raw_text = "".join(raw_chunks).upper()

    # Count and filter non-ACGT characters
    total_len = len(raw_text)
    filtered_chars = [ch for ch in raw_text if ch in ALLOWED_BASES]
    filtered_len = len(filtered_chars)
    removed = total_len - filtered_len

    cleaned = "".join(filtered_chars)

    if removed > 0:
        print(
            f"Warning ({file_path}): removed {removed:,} / {total_len:,} characters "
            f"that were not A/C/G/T (e.g. N, <, >, digits, whitespace)."
        )
    if filtered_len == 0:
        raise ValueError(
            f"No valid A/C/G/T bases remain after cleaning the input file: {file_path}"
        )

    return cleaned


def get_valid_pattern() -> str:
    """
    Prompt user for a pattern and require it contains only A/C/G/T.
    Returns the uppercased pattern.
    """
    while True:
        pattern = input("\nEnter pattern to search (DNA motif): ").strip().upper()
        if pattern == "":
            print("Pattern cannot be empty. Try again.")
            continue
        invalid = [ch for ch in pattern if ch not in ALLOWED_BASES]
        if invalid:
            print(
                f"Pattern contains invalid characters: {sorted(set(invalid))}. "
                f"Please enter only A, C, G, or T."
            )
            continue
        return pattern


def analyze_single_file(
    label: str, text: str, pattern: str, runs: int = 3
) -> Tuple[int, int, List[int], List[float]]:
    """
    Run Rabin-Karp on a single cleaned DNA sequence.

    Returns:
        total_length: length of the cleaned text
        num_matches: number of matches in the full sequence
        sizes: list of input sizes used in scaling experiment
        avg_times: average runtime for each size in `sizes`
    """
    print(f"\n=== {label.upper()} DNA SEQUENCE ANALYSIS ===")
    print("=" * 60)
    print(f"Sequence length (cleaned): {len(text):,} bp")
    print(f"Pattern searched: {pattern}\n")

    # ------------------ FULL TEXT BENCHMARK ------------------
    times_full = []
    match_positions: List[int] = []

    for _ in range(runs):
        t, output = Benchmark.time_function(
            PatternMatcher.rabin_karp_match, text, pattern
        )
        times_full.append(t)
        match_positions = output  # keep last match result

    avg_time_full = statistics.mean(times_full)

    print(f"Total matches found: {len(match_positions)}\n")
    if len(match_positions) > 20:
        print("Match positions (first 20):", match_positions[:20])
        print(f"... and {len(match_positions) - 20} more\n")
    else:
        print("Match positions:", match_positions, "\n")

    print("RABIN-KARP PERFORMANCE (Full Text)")
    print("=" * 60)
    print(f"| {'Algorithm':<12} | {'Avg Time (s)':<12} | {'Matches':<10} |")
    print("|" + "-" * 56 + "|")
    print(f"| {'Rabin-Karp':<12} | {avg_time_full:<12.6f} | {len(match_positions):<10} |")

    # ===================================================================
    # 📈 AVERAGE RUNTIME VS INPUT SIZE (Scaling for Rabin-Karp)
    # ===================================================================

    fractions = [0.2, 0.4, 0.6, 0.8, 1.0]
    sizes = sorted(
        set(
            max(len(pattern) + 1, int(len(text) * frac))
            for frac in fractions
            if int(len(text) * frac) >= len(pattern) + 1
        )
    )

    # Ensure at least 3 test cases even for short sequences
    if len(sizes) < 3:
        n = len(text)
        sizes = sorted(
            set(
                [
                    max(len(pattern) + 1, n // 3),
                    max(len(pattern) + 1, (2 * n) // 3),
                    n,
                ]
            )
        )

    runs_per_case = max(3, runs)

    print("\nSCALING EXPERIMENT (Rabin-Karp only)")
    print(f"(Each test case averaged over {runs_per_case} runs)")
    print("=" * 60)

    avg_times_vs_n: List[float] = []

    for n in sizes:
        subtext = text[:n]
        print(f"\nInput size n = {n} bp")
        times = []
        for _ in range(runs_per_case):
            t, _ = Benchmark.time_function(
                PatternMatcher.rabin_karp_match, subtext, pattern
            )
            times.append(t)
        avg_t = statistics.mean(times)
        avg_times_vs_n.append(avg_t)
        print(f"Average time: {avg_t:.6f} s")

    return len(text), len(match_positions), sizes, avg_times_vs_n


def analyze_dna_multiple(runs: int = 3):
    """
    Run Rabin-Karp empirical analysis on three DNA sequences:
    - dna_sequence.txt
    - chimpanzee.txt
    - dog.txt
    """

    # You can rename these to match your actual filenames
    datasets = [
        ("Human / DNA_Sequence", "dna_sequence.txt"),
        ("Chimpanzee", "chimpanzee.txt"),
        ("Dog", "dog.txt"),
    ]

    pattern = get_valid_pattern()

    results_per_dataset: Dict[str, Dict[str, object]] = {}

    for label, path in datasets:
        try:
            text = load_and_clean_sequence(path)
        except Exception as e:
            print(f"\nSkipping {label} ({path}) due to error: {e}")
            continue

        total_len, num_matches, sizes, avg_times = analyze_single_file(
            label, text, pattern, runs=runs
        )

        results_per_dataset[label] = {
            "length": total_len,
            "matches": num_matches,
            "sizes": sizes,
            "times": avg_times,
        }

    # ===================== PLOT: AVG TIME VS INPUT SIZE =====================

    if not results_per_dataset:
        print("\nNo datasets were successfully processed. Exiting.")
        return

    plt.figure(figsize=(10, 6))

    # Assign distinct colors/markers per dataset
    style_cycle = [
        {"color": "blue", "marker": "o"},
        {"color": "green", "marker": "^"},
        {"color": "red", "marker": "s"},
    ]

    for (label, _), style in zip(datasets, style_cycle):
        if label not in results_per_dataset:
            continue
        sizes = results_per_dataset[label]["sizes"]
        times = results_per_dataset[label]["times"]
        plt.plot(
            sizes,
            times,
            label=f"{label} (Rabin-Karp)",
            linewidth=1.5,
            marker=style["marker"],
            markersize=6,
            color=style["color"],
        )

    plt.title("Average Runtime vs Input Size (Rabin-Karp, DNA Pattern Matching)")
    plt.xlabel("Input size n (length of DNA text in bp)")
    plt.ylabel("Average running time (seconds)")
    plt.grid(linestyle="--", alpha=0.3)
    plt.legend(title="Genome")

    info_lines = []
    for label, info in results_per_dataset.items():
        info_lines.append(f"{label}: {info['length']:,} bp, matches={info['matches']}")

    info_text = (
        f"Pattern: \"{pattern}\" ({len(pattern)} bp)\n"
        + "\n".join(info_lines)
        + f"\nRuns per test case: {max(3, runs)}"
    )

    plt.text(
        0.02,
        0.95,
        info_text,
        transform=plt.gca().transAxes,
        fontsize=9,
        verticalalignment="top",
        bbox=dict(
            boxstyle="round,pad=0.3",
            facecolor="white",
            alpha=0.8,
            edgecolor="gray",
        ),
    )

    plt.tight_layout()
    plt.show()


# ======================== RUN ===========================
if __name__ == "__main__":
    analyze_dna_multiple()
