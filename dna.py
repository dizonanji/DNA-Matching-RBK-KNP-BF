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

    # Helpful diagnostics printed to stderr so they don't mix with normal output if desired
    if removed > 0:
        print(f"Warning: removed {removed:,} / {total_len:,} characters "
              f"that were not A/C/G/T (e.g. N, <, >, digits, whitespace).")
    if filtered_len == 0:
        raise ValueError("No valid A/C/G/T bases remain after cleaning the input file.")

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
            print(f"Pattern contains invalid characters: {sorted(set(invalid))}. "
                  f"Please enter only A, C, G, or T.")
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
        "Rabin-Karp": PatternMatcher.rabin_karp_match
    }

    results = {}
    match_positions = []

    # ------------------ MAIN 3-RUN BENCHMARK ------------------
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

    # ================= PRINT TIMING TABLE =================
    print("ALGORITHM PERFORMANCE COMPARISON")
    print("=" * 60)
    print(f"| {'Algorithm':<12} | {'Avg Time (s)':<12} | {'Matches':<10} |")
    print("|" + "-" * 56 + "|")

    for name, (t, m) in results.items():
        print(f"| {name:<12} | {t:<12.6f} | {m:<10} |")

    fastest = min(results, key=lambda x: results[x][0])
    print("\nFastest Algorithm:", fastest)

    # ===================================================================
    # 📈 5 RUN BENCHMARK — MATPLOTLIB PLOT
    # ===================================================================
    runs_plot = 5
    raw_times = {}

    print("\nBenchmarking 5 runs per algorithm...")
    for name, func in algos.items():
        times = []
        for _ in range(runs_plot):
            t, output = Benchmark.time_function(func, text, pattern)
            times.append(t)
        raw_times[name] = times

    x_runs = list(range(1, runs_plot + 1))

    plt.figure(figsize=(10, 6))
    line_styles = {
        "Naive":      {"color": "blue",  "marker": "o"},
        "KMP":        {"color": "green", "marker": "^"},
        "Rabin-Karp": {"color": "red",   "marker": "D"}
    }

    for algo, times in raw_times.items():
        style = line_styles.get(algo, {"marker": "o", "color": "black"})
        plt.plot(
            x_runs, times,
            label=algo,
            linewidth=1,
            marker=style["marker"],
            markersize=5,
            color=style["color"]
        )
        for x, t in zip(x_runs, times):
            plt.text(
                x, t, f"{t:.6f}",
                ha='center', va='bottom',
                fontsize=8, color=style["color"]
            )

    plt.title("Algorithm Runtime Comparison (DNA Sequence)")
    plt.xlabel("Run Number")
    plt.ylabel("Execution Time (seconds)")
    plt.xticks(x_runs)
    plt.grid(linestyle='--', alpha=0.3)
    plt.legend(title="Algorithm")

    info_text = (
        f"Text Length: {len(text):,} bp\n"
        f"Pattern: \"{pattern}\" ({len(pattern)} bp)"
    )
    plt.text(
        0.02, 0.95,
        info_text,
        transform=plt.gca().transAxes,
        fontsize=9,
        verticalalignment="top",
        bbox=dict(boxstyle="round,pad=0.3",
                  facecolor="white", alpha=0.8, edgecolor="gray")
    )

    plt.tight_layout()
    plt.show()


# ======================== RUN ===========================
if __name__ == "__main__":
    analyze_dna()
