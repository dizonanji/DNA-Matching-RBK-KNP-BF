"""
DNA Pattern Matching Benchmark using Naive, KMP & Rabin-Karp
Plots runtime results and uses benchmark.py for timing.
"""

from typing import List
from benchmark import Benchmark
import statistics
import matplotlib.pyplot as plt


class PatternMatcher:

    @staticmethod
    def naive_match(text: str, pattern: str) -> List[int]:
        matches = []
        t_len, p_len = len(text), len(pattern)

        for i in range(t_len - p_len + 1):
            if text[i:i+p_len] == pattern:
                matches.append(i)
        return matches

    @staticmethod
    def compute_lps(pattern: str) -> List[int]:
        lps = [0] * len(pattern)
        length = 0
        i = 1

        while i < len(pattern):
            if pattern[i] == pattern[length]:
                length += 1
                lps[i] = length
                i += 1
            else:
                if length != 0:
                    length = lps[length - 1]
                else:
                    lps[i] = 0
                    i += 1
        return lps

    @staticmethod
    def kmp_match(text: str, pattern: str) -> List[int]:
        matches = []
        lps = PatternMatcher.compute_lps(pattern)
        i = j = 0

        while i < len(text):
            if text[i] == pattern[j]:
                i += 1
                j += 1

            if j == len(pattern):
                matches.append(i-j)
                j = lps[j-1]

            elif i < len(text) and pattern[j] != text[i]:
                if j != 0:
                    j = lps[j-1]
                else:
                    i += 1

        return matches

    @staticmethod
    def rabin_karp_match(text: str, pattern: str, prime=101) -> List[int]:
        matches = []
        t_len, p_len = len(text), len(pattern)
        if p_len > t_len: return matches
        
        bases = sorted(set(text + pattern))  # DNA alphabet detection
        map_val = {b: i+1 for i, b in enumerate(bases)}
        base = len(bases)

        p_hash = t_hash = 0
        h = 1

        for _ in range(p_len-1):
            h = (h * base) % prime

        for i in range(p_len):
            p_hash = (base*p_hash + map_val[pattern[i]]) % prime
            t_hash = (base*t_hash + map_val[text[i]]) % prime

        for i in range(t_len - p_len + 1):
            if p_hash == t_hash and text[i:i+p_len] == pattern:
                matches.append(i)

            if i < t_len - p_len:
                t_hash = (base*(t_hash - map_val[text[i]]*h) + map_val[text[i+p_len]]) % prime
                if t_hash < 0: t_hash += prime

        return matches



# =======================================================================
#                        MAIN BENCHMARK + PLOTTER
# =======================================================================

def analyze_dna(file_path="human-dna.txt", runs=3):

    with open(file_path, "r") as f:
        text = f.read().replace("\n", "").strip()

    pattern = input("\nEnter pattern to search (DNA motif): ").strip()

    print("\nDNA SEQUENCE ANALYSIS RESULTS")
    print("="*60)
    print(f"Sequence length: {len(text):,} bp")
    print(f"Pattern searched: {pattern}\n")

    algos = {
        "Naive": PatternMatcher.naive_match,
        "KMP": PatternMatcher.kmp_match,
        "Rabin-Karp": PatternMatcher.rabin_karp_match
    }

    results = {}
    match_positions = None

    for name, func in algos.items():
        times = []
        for _ in range(runs):
            t, output = Benchmark.time_function(func, text, pattern)
            times.append(t)
            match_positions = output

        results[name] = (statistics.mean(times), len(match_positions))

    print(f"Total matches found: {len(match_positions)}\n")
    if len(match_positions) > 20:
        print("Match positions (first 20):", match_positions[:20])
        print(f"... and {len(match_positions)-20} more\n")
    else:
        print("Match positions:", match_positions, "\n")

    # ================= PRINT TIMING TABLE =================
    print("ALGORITHM PERFORMANCE COMPARISON")
    print("="*60)
    print(f"| {'Algorithm':<10} | {'Avg Time (s)':<12} | {'Matches':<10} |")
    print("|" + "-"*56 + "|")

    for name, (t, m) in results.items():
        print(f"| {name:<10} | {t:<12.6f} | {m:<10} |")

    fastest = min(results, key=lambda x: results[x][0])
    print("\nFastest Algorithm:", fastest)

    # ===================================================================
    # 📈 5 RUN BENCHMARK — MATPLOTLIB PLOT
    # ===================================================================
    runs = 5
    raw_times = {}

    print("\nBenchmarking 5 runs per algorithm...")
    for name, func in algos.items():
        times = []
        for _ in range(runs):
            t, output = Benchmark.time_function(func, text, pattern)
            times.append(t)
        raw_times[name] = times

    x_runs = list(range(1, runs + 1))

    plt.figure(figsize=(10, 6))
    line_styles = {
        "Naive":      {"color": "blue",  "marker": "o"},
        "KMP":        {"color": "green", "marker": "^"},
        "Rabin-Karp": {"color": "red",   "marker": "D"}
    }

    for algo, times in raw_times.items():
        style = line_styles.get(algo)
        plt.plot(
            x_runs, times,
            label=algo,
            linewidth=1,
            marker=style["marker"],
            markersize=5,
            color=style["color"]
        )
        for x, t in zip(x_runs, times):
            plt.text(x, t, f"{t:.6f}",
                     ha='center', va='bottom',
                     fontsize=8, color=style["color"])

    plt.title("Algorithm Runtime Comparison (DNA Sequence)")
    plt.xlabel("Run Number")
    plt.ylabel("Execution Time (seconds)")
    plt.xticks(x_runs)
    plt.grid(linestyle='--', alpha=0.3)
    plt.legend(title="Algorithm")

    # DNA version annotation
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
        bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.8, edgecolor="gray")
    )

    plt.tight_layout()
    plt.show()



# ======================== RUN ===========================
if __name__ == "__main__":
    analyze_dna()
