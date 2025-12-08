"""
Protein Pattern Matching Benchmark using Naive, KMP & Rabin-Karp
Supports MULTIPLE protein motifs (patterns) entered by the user.
Plots runtime results and uses benchmark.py for timing.
"""

from typing import List, Dict
from benchmark import Benchmark
import statistics
import matplotlib.pyplot as plt


class PatternMatcher:
    """Collection of pattern matching algorithms for protein sequences."""

    @staticmethod
    def naive_match(text: str, pattern: str) -> List[int]:
        """
        Naive (brute force) pattern matching algorithm.

        Args:
            text: Protein sequence to search in
            pattern: Pattern/motif to search for

        Returns:
            List of starting positions where pattern matches
        """
        matches: List[int] = []
        text_len = len(text)
        pattern_len = len(pattern)

        for i in range(text_len - pattern_len + 1):
            if text[i:i + pattern_len] == pattern:
                matches.append(i)

        return matches

    @staticmethod
    def compute_lps_array(pattern: str) -> List[int]:
        """
        Compute Longest Proper Prefix which is also Suffix (LPS) array for KMP.

        Args:
            pattern: Pattern string

        Returns:
            LPS array
        """
        length = 0  # Length of the previous longest prefix suffix
        lps = [0] * len(pattern)
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
        """
        Knuth-Morris-Pratt pattern matching algorithm.

        Args:
            text: Protein sequence to search in
            pattern: Pattern/motif to search for

        Returns:
            List of starting positions where pattern matches
        """
        matches: List[int] = []
        text_len = len(text)
        pattern_len = len(pattern)

        if pattern_len == 0 or text_len == 0 or pattern_len > text_len:
            return matches

        # Compute LPS array
        lps = PatternMatcher.compute_lps_array(pattern)

        i = 0  # Index for text
        j = 0  # Index for pattern

        while i < text_len:
            if pattern[j] == text[i]:
                i += 1
                j += 1

            if j == pattern_len:
                matches.append(i - j)
                j = lps[j - 1]
            elif i < text_len and pattern[j] != text[i]:
                if j != 0:
                    j = lps[j - 1]
                else:
                    i += 1

        return matches

    @staticmethod
    def rabin_karp_match(text: str, pattern: str, prime: int = 101) -> List[int]:
        """
        Rabin-Karp pattern matching algorithm using rolling hash.

        Args:
            text: Protein sequence to search in
            pattern: Pattern/motif to search for
            prime: Prime number for hashing

        Returns:
            List of starting positions where pattern matches
        """
        matches: List[int] = []
        text_len = len(text)
        pattern_len = len(pattern)

        if pattern_len == 0 or text_len == 0 or pattern_len > text_len:
            return matches

        # Auto-detect alphabet from protein sequence and pattern
        alphabet = sorted(set(text + pattern))
        char_map = {ch: i + 1 for i, ch in enumerate(alphabet)}
        base = len(alphabet)

        # Calculate hash values
        pattern_hash = 0
        text_hash = 0
        h = 1

        # Calculate h = pow(base, pattern_len-1) % prime
        for _ in range(pattern_len - 1):
            h = (h * base) % prime

        # Calculate hash value of pattern and first window of text
        for i in range(pattern_len):
            pattern_hash = (base * pattern_hash + char_map[pattern[i]]) % prime
            text_hash = (base * text_hash + char_map[text[i]]) % prime

        # Slide the pattern over text one by one
        for i in range(text_len - pattern_len + 1):
            # Check if hash values match
            if pattern_hash == text_hash:
                # Check characters one by one for exact match
                if text[i:i + pattern_len] == pattern:
                    matches.append(i)

            # Calculate hash value for next window
            if i < text_len - pattern_len:
                text_hash = (
                    base * (text_hash - char_map[text[i]] * h)
                    + char_map[text[i + pattern_len]]
                ) % prime

                # Convert negative hash to positive
                if text_hash < 0:
                    text_hash += prime

        return matches


# =======================================================================
#                     MULTI-PATTERN MATCHING CLASS
# =======================================================================

class MultiPatternMatcher:
    """
    Wrapper class to search for MULTIPLE patterns in the same text.
    It reuses the single-pattern algorithms from PatternMatcher.

    All methods return:
        Dict[str, List[int]]
    mapping each pattern -> list of starting indices where it appears.
    """

    @staticmethod
    def naive_multi(text: str, patterns: List[str]) -> Dict[str, List[int]]:
        results: Dict[str, List[int]] = {}
        for pat in patterns:
            results[pat] = PatternMatcher.naive_match(text, pat)
        return results

    @staticmethod
    def kmp_multi(text: str, patterns: List[str]) -> Dict[str, List[int]]:
        results: Dict[str, List[int]] = {}
        for pat in patterns:
            results[pat] = PatternMatcher.kmp_match(text, pat)
        return results

    @staticmethod
    def rabin_karp_multi(text: str, patterns: List[str], prime: int = 101) -> Dict[str, List[int]]:
        results: Dict[str, List[int]] = {}
        for pat in patterns:
            results[pat] = PatternMatcher.rabin_karp_match(text, pat, prime=prime)
        return results


# =======================================================================
#                  MAIN BENCHMARK + PLOTS (MULTIPLE PATTERNS)
# =======================================================================

def analyze_proteins(file_path="protein.txt", runs=3):

    with open(file_path, "r") as f:
        text = f.read().replace("\n", "").strip()

    raw = input(
        "\nEnter protein motifs to search (comma-separated, e.g. KDEL,GGH,HLH): "
    ).strip()

    patterns = [p.strip() for p in raw.split(",") if p.strip()]

    if not patterns:
        print("No valid patterns entered. Exiting.")
        return

    print("\nPROTEIN SEQUENCE ANALYSIS RESULTS")
    print("=" * 60)
    print(f"Sequence length: {len(text):,} amino acids")
    print(f"Patterns searched ({len(patterns)}): {', '.join(patterns)}\n")

    algos = {
        "Naive": MultiPatternMatcher.naive_multi,
        "KMP": MultiPatternMatcher.kmp_multi,
        "Rabin-Karp": MultiPatternMatcher.rabin_karp_multi
    }

    # ================= MULTI-PATTERN FULL-TEXT BENCHMARK =================

    results = {}                 # algorithm -> (avg_time, total_matches_across_patterns)
    matches_per_algo = {}        # algorithm -> {pattern: [positions]}

    for name, func in algos.items():
        times = []
        last_output: Dict[str, List[int]] = {}

        for _ in range(runs):
            t, output = Benchmark.time_function(func, text, patterns)
            times.append(t)
            last_output = output   # dict: pattern -> positions

        total_matches = sum(len(v) for v in last_output.values())
        results[name] = (statistics.mean(times), total_matches)
        matches_per_algo[name] = last_output

    # ================= PRINT MATCH SUMMARY (using Naive as reference) =================

    ref_algo = "Naive" if "Naive" in matches_per_algo else list(matches_per_algo.keys())[0]
    ref_matches = matches_per_algo[ref_algo]

    total_ref_matches = sum(len(v) for v in ref_matches.values())
    print(f"Total matches found across all patterns ({ref_algo}): {total_ref_matches}\n")

    # Print per-pattern positions (capped)
    max_positions_to_show = 20
    for pat in patterns:
        positions = ref_matches.get(pat, [])
        print(f"Pattern: {pat}")
        if not positions:
            print("  No matches found.\n")
        elif len(positions) > max_positions_to_show:
            print("  Match positions (first 20):", positions[:max_positions_to_show])
            print(f"  ... and {len(positions) - max_positions_to_show} more\n")
        else:
            print("  Match positions:", positions, "\n")

    # ================= PRINT TIMING TABLE (FULL TEXT) =================

    print("ALGORITHM PERFORMANCE COMPARISON (Full Text, all patterns)")
    print("=" * 60)
    print(f"| {'Algorithm':<10} | {'Avg Time (s)':<12} | {'Total Matches':<14} |")
    print("|" + "-" * 62 + "|")

    for name, (t, m) in results.items():
        print(f"| {name:<10} | {t:<12.6f} | {m:<14} |")

    fastest = min(results, key=lambda x: results[x][0])
    print("\nFastest Algorithm (full text, all patterns):", fastest)

    # ===================================================================
    # 📈 AVERAGE RUNTIME VS INPUT SIZE (at least 3 runs per test case)
    # ===================================================================

    max_pattern_len = max(len(p) for p in patterns)

    fractions = [0.2, 0.4, 0.6, 0.8, 1.0]
    sizes = sorted(set(
        max(max_pattern_len + 1, int(len(text) * frac))
        for frac in fractions
        if int(len(text) * frac) >= max_pattern_len + 1
    ))

    # Guarantee at least 3 test cases
    if len(sizes) < 3:
        n = len(text)
        sizes = sorted(set([
            max(max_pattern_len + 1, n // 3),
            max(max_pattern_len + 1, (2 * n) // 3),
            n
        ]))

    runs_per_case = max(3, runs)  # ensure ≥ 3 runs per test case

    print("\nSCALING EXPERIMENT: Average runtime vs. input size")
    print(f"(Each test case averaged over {runs_per_case} runs; multiple patterns)")
    print("=" * 60)

    # algo -> list of avg times for each size
    avg_times_vs_n = {name: [] for name in algos.keys()}

    # Tabulated results for each test case
    for n in sizes:
        subtext = text[:n]
        print(f"\nInput size n = {n} characters")
        print(f"| {'Algorithm':<10} | {'Avg Time (s)':<12} |")
        print("|" + "-" * 30 + "|")

        for name, func in algos.items():
            times = []
            for _ in range(runs_per_case):
                t, _ = Benchmark.time_function(func, subtext, patterns)
                times.append(t)
            avg_t = statistics.mean(times)
            avg_times_vs_n[name].append(avg_t)
            print(f"| {name:<10} | {avg_t:<12.6f} |")

    # ===================== PLOT: AVG TIME VS INPUT SIZE =====================

    plt.figure(figsize=(10, 6))

    # Original color + marker assignments (UNCHANGED)
    line_styles = {
        "Naive":      {"color": "blue",  "marker": "o"},
        "KMP":        {"color": "green", "marker": "^"},
        "Rabin-Karp": {"color": "red",   "marker": "D"}
    }

    for name in algos.keys():
        style = line_styles.get(name, {"color": "black", "marker": "o"})
        plt.plot(
            sizes,
            avg_times_vs_n[name],
            marker=style["marker"],
            color=style["color"],
            linewidth=1.5,
            markersize=6,
            label=name
        )

    plt.title("Average Runtime vs Input Size (Protein Multi-Pattern Matching)")
    plt.xlabel("Input size n (length of text in amino acids)")
    plt.ylabel("Average running time (seconds)")
    plt.grid(linestyle="--", alpha=0.3)
    plt.legend(title="Algorithm")

    info_text = (
        f"Full Text Length: {len(text):,} aa\n"
        f"Patterns: {len(patterns)} motifs\n"
        f"Runs per test case: {runs_per_case}"
    )

    plt.text(
        0.02, 0.95,
        info_text,
        transform=plt.gca().transAxes,
        fontsize=9,
        verticalalignment="top",
        bbox=dict(
            boxstyle="round,pad=0.3",
            facecolor="white",
            alpha=0.8,
            edgecolor="gray"
        )
    )

    plt.tight_layout()
    plt.show()


# ======================== RUN ===========================
if __name__ == "__main__":
    analyze_proteins()
