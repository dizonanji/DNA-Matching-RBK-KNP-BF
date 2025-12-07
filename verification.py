import random
import csv
import json
import platform
import re
import plotly.graph_objects as go
from matcher import PatternMatcher
from benchmark import Benchmark

# ------------------------------
# FLAGS AND CONFIGURATION
# ------------------------------
# Control which sequences to use
USE_PREDEFINED = False       # Include predefined test cases
USE_MANUAL = False          # Include manual DNA + pattern
USE_RANDOM = False           # Include randomly generated sequences
USE_FILE = True            # Include DNA from file

# DNA / Pattern files
DNA_INPUT_FILE = "dog.txt"  # File containing DNA sequence
DNA_PATTERN_FILE = None     # Optional: file containing pattern
DNA_PATTERN_FOR_FILE = "GCTA" # Manual pattern if no pattern file

# Manual sequence
MANUAL_DNA_TEXT = "ATCGATCGGCTAATCGGCTAGCTAATCG"
MANUAL_DNA_PATTERN = "GCTA"

# Random settings
RANDOM_SIZES = [4, 10, 25, 50, 80, 100, 250, 500, 1000, 5000, 10000, 25000, 50000]
PATTERN_LENGTH = 4
REPETITIONS = 3
EXPANSION_FACTORS = [1, 2, 5]

# Predefined test cases
PREDEFINED_TESTS = [
    {"desc": "Simple match", "text": "ATCGATCGA", "pattern": "ATC", "expected": [0, 4]},
    {"desc": "No match", "text": "GATTACA", "pattern": "CCC", "expected": []},
    {"desc": "Pattern longer than text", "text": "ATC", "pattern": "ATCG", "expected": []},
    {"desc": "Multiple overlapping matches", "text": "ATATATAT", "pattern": "ATAT", "expected": [0, 2, 4]},
    {"desc": "Highly repetitive DNA", "text": "A"*17, "pattern": "AAAA", "expected": list(range(14))}
]

# ------------------------------
# PREDEFINED TEST CASES
# ------------------------------
PREDEFINED_TESTS = [
    {"desc": "Simple match", "text": "ATCGATCGA", "pattern": "ATC", "expected": [0, 4]},
    {"desc": "No match", "text": "GATTACA", "pattern": "CCC", "expected": []},
    {"desc": "Pattern longer than text", "text": "ATC", "pattern": "ATCG", "expected": []},
    {"desc": "Multiple overlapping matches", "text": "ATATATAT", "pattern": "ATAT", "expected": [0, 2, 4]},
    {"desc": "Highly repetitive DNA", "text": "A" * 17, "pattern": "AAAA", "expected": list(range(14))}
]

# ------------------------------
# DNA GENERATION & LOADERS
# ------------------------------
def generate_dna_sequence(length):
    return ''.join(random.choice("ATCG") for _ in range(length))

def generate_dna_pattern(length):
    return ''.join(random.choice("ATCG") for _ in range(length))

def load_dna_from_file(file_path):
    """
    Load a DNA sequence from a file.
    - Accepts FASTA (lines starting with '>') or plain text/multiline files.
    - Removes any characters other than A/T/C/G.
    - Returns uppercase cleaned sequence.
    """
    with open(file_path, 'r') as f:
        raw = f.read().upper()

    # Remove FASTA header lines (starting with '>')
    lines = raw.splitlines()
    seq_lines = [line for line in lines if not line.startswith(">")]

    # Join and filter only valid bases
    joined = "".join(seq_lines)
    cleaned = re.sub(r"[^ATCG]", "", joined)

    removed = len(joined) - len(cleaned)
    print(f"Loaded DNA file: {len(cleaned)} bases (removed {removed} invalid characters)")

    if len(cleaned) == 0:
        raise ValueError("ERROR: File does not contain valid A/T/C/G bases.")

    return cleaned

def load_pattern_from_file(file_path):
    """Load a simple pattern string from a text file (single line or multiline joined)."""
    with open(file_path, 'r') as f:
        raw = f.read().upper()
    # Keep only A/T/C/G for safety
    cleaned = re.sub(r"[^ATCG]", "", raw)
    return cleaned.strip()

# ------------------------------
# RUN TESTS
# ------------------------------
# ------------------------------
# RUN TESTS
# ------------------------------
def run_tests():
    results = []
    sequences = []

    # 1) USER FILE
    if USE_FILE and DNA_INPUT_FILE:
        text = load_dna_from_file(DNA_INPUT_FILE)
        if DNA_PATTERN_FILE:
            pattern = load_pattern_from_file(DNA_PATTERN_FILE)
        else:
            pattern = DNA_PATTERN_FOR_FILE
        # Repeat the file test REPETITIONS times
        for rep in range(REPETITIONS):
            sequences.append((text, pattern, None, f"User DNA File (Run {rep+1})"))

    # 2) MANUAL SEQUENCE
    if USE_MANUAL:
        sequences.append((MANUAL_DNA_TEXT, MANUAL_DNA_PATTERN, None, "Manual DNA"))

    # 3) PREDEFINED TEST CASES
    if USE_PREDEFINED:
        for test in PREDEFINED_TESTS:
            for factor in EXPANSION_FACTORS:
                text = test["text"] * factor
                pattern = test["pattern"]
                expected = [i*factor for i in test["expected"]]
                desc = f"{test['desc']} x{factor}"
                sequences.append((text, pattern, expected, desc))

    # 4) RANDOM SEQUENCES
    if USE_RANDOM:
        for size in RANDOM_SIZES:
            text = generate_dna_sequence(size)
            pattern = generate_dna_pattern(PATTERN_LENGTH)
            sequences.append((text, pattern, None, "Random sequence"))

    # ------------------------------
    # BENCHMARK FUNCTIONS
    # ------------------------------
    for idx, (text, pattern, expected, desc) in enumerate(sequences):
        size = len(text)
        print(f"\nRunning test {idx + 1} | {desc} | Text size: {size}")

        naive_times, kmp_times, rk_times = [], [], []

        for run in range(REPETITIONS):
            print(f"  Run {run + 1}/{REPETITIONS}")
            n_time, naive_matches = Benchmark.time_function(PatternMatcher.naive_match, text, pattern)
            k_time, kmp_matches = Benchmark.time_function(PatternMatcher.kmp_match, text, pattern)
            r_time, rk_matches = Benchmark.time_function(PatternMatcher.rabin_karp_match, text, pattern)

            naive_times.append(n_time)
            kmp_times.append(k_time)
            rk_times.append(r_time)

            if run == 0:
                correct_kmp = (kmp_matches == naive_matches)
                correct_rk  = (rk_matches == naive_matches)
                num_naive = len(naive_matches)
                num_kmp = len(kmp_matches)
                num_rk = len(rk_matches)
                kmp_diff = len(set(kmp_matches) ^ set(naive_matches))
                rk_diff = len(set(rk_matches) ^ set(naive_matches))
                if expected is not None:
                    rk_correct_vs_expected = (rk_matches == expected)
                else:
                    rk_correct_vs_expected = None

        avg_naive = sum(naive_times)/REPETITIONS
        avg_kmp   = sum(kmp_times)/REPETITIONS
        avg_rk    = sum(rk_times)/REPETITIONS

        print(f"Naive avg time: {avg_naive:.6f}s | Matches: {num_naive}")
        print(f"KMP avg time:   {avg_kmp:.6f}s | Matches: {num_kmp} | Correct: {correct_kmp}")
        print(f"Rabin–Karp avg: {avg_rk:.6f}s | Matches: {num_rk} | Correct vs naive: {correct_rk} | Correct vs expected: {rk_correct_vs_expected}")

        results.append([
            size,
            text,
            pattern,
            desc,
            json.dumps(naive_matches),
            json.dumps(kmp_matches),
            json.dumps(rk_matches),
            num_naive,
            num_kmp,
            num_rk,
            correct_kmp,
            correct_rk,
            kmp_diff,
            rk_diff,
            avg_naive,
            avg_kmp,
            avg_rk,
            json.dumps(expected) if expected else None
        ])
    return results


# ------------------------------
# WRITE CSV
# ------------------------------
def save_results(results, filename="benchmark_results.csv"):
    with open(filename, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow([
            "Text Size",
            "DNA Text",
            "Pattern",
            "Description",
            "Naive Matches",
            "KMP Matches",
            "RK Matches",
            "#Naive",
            "#KMP",
            "#RK",
            "KMP Correct",
            "RK Correct",
            "KMP Diff",
            "RK Diff",
            "Naive Avg Time",
            "KMP Avg Time",
            "RK Avg Time",
            "Expected Matches"
        ])
        writer.writerows(results)
    print(f"\nResults saved to {filename}")

# ------------------------------
# PLOT RESULTS
# ------------------------------
def plot_results(results):
    # Separate results by type
    predefined = [r for r in results if "Random" not in r[3] and "User DNA File" not in r[3] and "Manual DNA" not in r[3]]
    manual_seq = [r for r in results if "Manual DNA" in r[3]]
    random_seq = [r for r in results if "Random" in r[3]]
    user_file  = [r for r in results if "User DNA File" in r[3]]

    def plot_group(group, title):
        if not group:
            return

        group_sizes = [r[0] for r in group]
        naive_times = [r[14] for r in group]
        kmp_times   = [r[15] for r in group]
        rk_times    = [r[16] for r in group]

        fig = go.Figure()

        fig.add_trace(go.Scatter(
            x=group_sizes, y=naive_times, mode='lines+markers', name='Naive',
            marker=dict(symbol='circle', size=8), line=dict(color='blue')
        ))
        fig.add_trace(go.Scatter(
            x=group_sizes, y=kmp_times, mode='lines+markers', name='KMP',
            marker=dict(symbol='triangle-up', size=8), line=dict(color='green')
        ))
        fig.add_trace(go.Scatter(
            x=group_sizes, y=rk_times, mode='lines+markers', name='Rabin-Karp',
            marker=dict(symbol='diamond', size=8), line=dict(color='red')
        ))

        fig.update_layout(
            title=title,
            xaxis_title='DNA Sequence Length',
            yaxis_title='Average Execution Time (seconds)',
            legend_title='Algorithm',
            template='plotly_white',
            hovermode='x unified'
        )

        fig.show()

    # Plot according to flags
    if USE_PREDEFINED:
        plot_group(predefined, "Predefined Sequence Tests")
    if USE_MANUAL:
        plot_group(manual_seq, "Manual Sequence Tests")
    if USE_RANDOM:
        plot_group(random_seq, "Randomly Generated DNA Sequence Tests")
    if USE_FILE:
        plot_group(user_file, f"User DNA File: {DNA_INPUT_FILE}")


# ------------------------------
# MACHINE & PARAMETER INFO
# ------------------------------
def print_environment_info():
    print("\n--- Environment & Parameters ---")
    print(f"Python version: {platform.python_version()}")
    print(f"Platform: {platform.platform()}")
    print(f"Manual sequence used: {USE_MANUAL}")
    print(f"Random sequence sizes: {RANDOM_SIZES}")
    print(f"Pattern length: {PATTERN_LENGTH}")
    print(f"Repetitions per test: {REPETITIONS}")
    print(f"Predefined test cases: {len(PREDEFINED_TESTS)}")
    print(f"Expansion factors: {EXPANSION_FACTORS}")
    if DNA_INPUT_FILE:
        print(f"User DNA file: {DNA_INPUT_FILE}")
    if USE_FILE and DNA_PATTERN_FILE:
        print(f"User pattern file: {DNA_PATTERN_FILE}")

# ------------------------------
# MAIN
# ------------------------------
if __name__ == "__main__":
    print_environment_info()
    results = run_tests()
    save_results(results)
    plot_results(results)
