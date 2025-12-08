"""
Protein Pattern Matching Benchmark using Naive, KMP & Rabin-Karp
Uses benchmark.py for time measurement.
"""

from typing import List
from benchmark import Benchmark   # <====== IMPORT HERE
import statistics


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
        
        # auto-map amino acids to numeric values
        amino = sorted(set(text + pattern))
        map_val = {aa: i+1 for i, aa in enumerate(amino)}
        base = len(amino)

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
                if t_hash < 0:
                    t_hash += prime

        return matches



# ======================== RUN ANALYSIS ===========================

def analyze_proteins(file_path="protein.txt", runs=3):

    with open(file_path, "r") as f:
        text = f.read().replace("\n", "").strip()

    pattern = input("\nEnter pattern to search (protein motif): ").strip()

    print("\nPROTEIN SEQUENCE ANALYSIS RESULTS")
    print("="*60)
    print(f"Sequence length: {len(text):,} amino acids")
    print(f"Pattern searched: {pattern}\n")


    algos = {
        "NAIVE": PatternMatcher.naive_match,
        "KMP": PatternMatcher.kmp_match,
        "RABIN-KARP": PatternMatcher.rabin_karp_match
    }

    results = {}

    for name, func in algos.items():
        times = []
        matches = None

        for _ in range(runs):
            t, output = Benchmark.time_function(func, text, pattern)
            times.append(t)
            matches = len(output)

        results[name] = (statistics.mean(times), matches)
        match_positions = output


    # Print matches summary
    print(f"Total matches found: {len(match_positions)}\n")

    if len(match_positions) > 20:
        print("Match positions (first 20):", match_positions[:20])
        print(f"... and {len(match_positions)-20} more\n")
    else:
        print("Match positions:", match_positions, "\n")


    # Print timing table
    print("ALGORITHM PERFORMANCE COMPARISON")
    print("="*60)
    print(f"| {'Algorithm':<10} | {'Time (sec)':<12} | {'Matches':<10} |")
    print("|" + "-"*56 + "|")

    for name, (t, m) in results.items():
        print(f"| {name:<10} | {t:<12.6f} | {m:<10} |")

    fastest = min(results, key=lambda x: results[x][0])
    print("\nFastest Algorithm:", fastest, "\n")


# ======================== RUN ===========================
if __name__ == "__main__":
    analyze_proteins()
