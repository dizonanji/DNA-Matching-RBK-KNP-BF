"""
DNA Pattern Matching Benchmark using Naive, KMP & Rabin-Karp
Plots runtime results and uses benchmark.py for timing.
This version cleans the input sequence so matcher.py only gets A/C/G/T uppercase.
"""

from typing import List
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
            f"Warning: removed {removed:,} / {total_len:,} characters "
            f"that were not A/C/G/T (e.g. N, <, >, digits, whitespace)."
        )
    if filtered_len == 0:
        raise ValueError(
            "No valid A/C/G/T bases remain after cleaning the input file."
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


def analyze_dna(file_path="human-dna.txt", runs=3):
    # Load & clean
    text = load_and_clean_sequence(file_path)
    pattern = get_valid_pattern()

    print("\nDNA SEQUENCE ANALYSIS RESULTS")
    print("=" * 60)
    print(f"Sequence length (cleaned): {len(text):,} bp")
    print(f"Pattern searched: {pattern}\n")

    algos = {
        "Naive": PatternMatcher.naive_match,
        "KMP": PatternMatcher.kmp_match,
        "Rabin-Karp": PatternMatcher.rabin_karp_match,
    }

    results = {}
    match_positions: List[int] = []

    # ------------------ MAIN runs-RUN BENCHMARK ON FULL TEXT ------------------
    for name, func in algos.items():
        times = []
        for _ in range(runs):
            t, output = Benchmark.time_function(func, text, pattern)
            times.append(t)
            match_positions = output  # keep last match result

        results[name] = (statistics.mean(times), len(match_positions))

    print(f"Total matches found: {len(match_positions)}\n")
    if len(match_positions) > 20:
        print("Match positions (first 20):", match_positions[:20])
        print(f"... and {len(match_positions) - 20} more\n")
    else:
        print("Match positions:", match_positions, "\n")

    # ================= PRINT TIMING TABLE (FULL TEXT) =================
    print("ALGORITHM PERFORMANCE COMPARISON (Full Text)")
    print("=" * 60)
    print(f"| {'Algorithm':<12} | {'Avg Time (s)':<12} | {'Matches':<10} |")
    print("|" + "-" * 56 + "|")

    for name, (t, m) in results.items():
        print(f"| {name:<12} | {t:<12.6f} | {m:<10} |")

    fastest = min(results, key=lambda x: results[x][0])
    print("\nFastest Algorithm (full text):", fastest)

    # ===================================================================
    # 📈 AVERAGE RUNTIME VS INPUT SIZE (each test case has ≥3 runs)
    # ===================================================================

    # Define different input sizes as fractions of the full cleaned sequence
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

    runs_per_case = max(3, runs)  # each test case must have at least 3 runs

    print("\nSCALING EXPERIMENT: Average runtime vs. input size")
    print(f"(Each test case averaged over {runs_per_case} runs)")
    print("=" * 60)

    # algo -> list of avg times for each input size
    avg_times_vs_n = {name: [] for name in algos.keys()}

    # Tabulated results per test size
    for n in sizes:
        subtext = text[:n]
        print(f"\nInput size n = {n} bp")
        print(f"| {'Algorithm':<12} | {'Avg Time (s)':<12} |")
        print("|" + "-" * 32 + "|")

        for name, func in algos.items():
            times = []
            for _ in range(runs_per_case):
                t, _ = Benchmark.time_function(func, subtext, pattern)
                times.append(t)
            avg_t = statistics.mean(times)
            avg_times_vs_n[name].append(avg_t)
            print(f"| {name:<12} | {avg_t:<12.6f} |")

    # ===================== PLOT: AVG TIME VS INPUT SIZE =====================

    plt.figure(figsize=(10, 6))

    # Retain the original color & marker assignment
    line_styles = {
        "Naive": {"color": "blue", "marker": "o"},
        "KMP": {"color": "green", "marker": "^"},
        "Rabin-Karp": {"color": "red", "marker": "D"},
    }

    for name in algos.keys():
        style = line_styles.get(name, {"color": "black", "marker": "o"})
        plt.plot(
            sizes,
            avg_times_vs_n[name],
            label=name,
            linewidth=1.5,
            marker=style["marker"],
            markersize=6,
            color=style["color"],
        )

    plt.title("Average Runtime vs Input Size (DNA Pattern Matching)")
    plt.xlabel("Input size n (length of DNA text in bp)")
    plt.ylabel("Average running time (seconds)")
    plt.grid(linestyle="--", alpha=0.3)
    plt.legend(title="Algorithm")

    info_text = (
        f"Full Text Length: {len(text):,} bp\n"
        f"Pattern: \"{pattern}\" ({len(pattern)} bp)\n"
        f"Runs per test case: {runs_per_case}"
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
    analyze_dna()
