"""
Command-line interface for DNA pattern alignment tool with matplotlib visualization.
Updated to support both FASTA and TXT files.
"""
import argparse
import matplotlib.pyplot as plt
import numpy as np
from typing import List, Tuple, Dict, Any
from parsers import DNAParser  # Updated import
from matcher import PatternMatcher
from aligner import LocalAligner
from benchmark import Benchmark

class DNACli:
    """Command-line interface for DNA pattern alignment tool with visualization."""
    
    def __init__(self):
        self.parser = DNAParser()  # Updated to use DNAParser
        self.matcher = PatternMatcher()
        self.aligner = LocalAligner()
        self.benchmark = Benchmark()
    
    def run_algorithm(self, algorithm: str, text: str, pattern: str) -> Tuple[float, List[int]]:
        """Run specified algorithm and return timing and results."""
        if algorithm == 'naive':
            return self.benchmark.time_function(self.matcher.naive_match, text, pattern)
        elif algorithm == 'kmp':
            return self.benchmark.time_function(self.matcher.kmp_match, text, pattern)
        elif algorithm == 'rk':
            return self.benchmark.time_function(self.matcher.rabin_karp_match, text, pattern)
        else:
            raise ValueError(f"Unknown algorithm: {algorithm}")
    
    def visualize_results(self, results: Dict[str, Dict[str, Any]], 
                         dna_sequence: str, pattern: str, save_plots: bool = False):
        """Create comprehensive matplotlib visualizations of the results."""
        
        # Set up the plotting style
        plt.style.use('default')
        fig = plt.figure(figsize=(16, 12))
        
        # 1. Performance Comparison Bar Chart
        ax1 = plt.subplot(2, 3, 1)
        algorithms = list(results.keys())
        times = [results[algo]['time'] for algo in algorithms]
        colors = ['#FF6B6B', '#4ECDC4', '#45B7D1']
        
        bars = ax1.bar(algorithms, times, color=colors[:len(algorithms)])
        ax1.set_title('Algorithm Performance Comparison', fontsize=14, fontweight='bold')
        ax1.set_ylabel('Execution Time (seconds)')
        ax1.set_xlabel('Algorithm')
        
        # Add value labels on bars
        for bar, time in zip(bars, times):
            height = bar.get_height()
            ax1.text(bar.get_x() + bar.get_width()/2., height,
                    f'{time:.6f}s', ha='center', va='bottom', fontsize=10)
        
        # 2. Match Count Comparison
        ax2 = plt.subplot(2, 3, 2)
        match_counts = [len(results[algo]['matches']) for algo in algorithms]
        
        bars2 = ax2.bar(algorithms, match_counts, color=colors[:len(algorithms)])
        ax2.set_title('Matches Found by Algorithm', fontsize=14, fontweight='bold')
        ax2.set_ylabel('Number of Matches')
        ax2.set_xlabel('Algorithm')
        
        # Add value labels on bars
        for bar, count in zip(bars2, match_counts):
            height = bar.get_height()
            ax2.text(bar.get_x() + bar.get_width()/2., height,
                    f'{count}', ha='center', va='bottom', fontsize=10)
        
        # 3. Performance vs Matches Scatter Plot
        ax3 = plt.subplot(2, 3, 3)
        scatter = ax3.scatter(times, match_counts, c=colors[:len(algorithms)], 
                             s=100, alpha=0.7, edgecolors='black')
        
        for i, algo in enumerate(algorithms):
            ax3.annotate(algo.upper(), (times[i], match_counts[i]), 
                        xytext=(5, 5), textcoords='offset points', fontsize=10)
        
        ax3.set_title('Performance vs Accuracy', fontsize=14, fontweight='bold')
        ax3.set_xlabel('Execution Time (seconds)')
        ax3.set_ylabel('Number of Matches')
        ax3.grid(True, alpha=0.3)
        
        # 4. Match Distribution Along Sequence
        ax4 = plt.subplot(2, 3, 4)
        if results:
            # Use the first algorithm's results for match distribution
            first_algo = list(results.keys())[0]
            matches = results[first_algo]['matches']
            
            if matches:
                # Create histogram of match positions
                bins = min(50, len(dna_sequence) // 1000 + 1)
                ax4.hist(matches, bins=bins, alpha=0.7, color=colors[0], edgecolor='black')
                ax4.set_title(f'Match Distribution Along Sequence\n({first_algo.upper()})', 
                             fontsize=14, fontweight='bold')
                ax4.set_xlabel('Sequence Position')
                ax4.set_ylabel('Number of Matches')
                ax4.grid(True, alpha=0.3)
            else:
                ax4.text(0.5, 0.5, 'No matches found', ha='center', va='center', 
                        transform=ax4.transAxes, fontsize=12)
                ax4.set_title('Match Distribution', fontsize=14, fontweight='bold')
        
        # 5. Relative Performance (if multiple algorithms)
        ax5 = plt.subplot(2, 3, 5)
        if len(algorithms) > 1:
            # Normalize times relative to fastest algorithm
            min_time = min(times)
            relative_times = [t / min_time for t in times]
            
            bars3 = ax5.bar(algorithms, relative_times, color=colors[:len(algorithms)])
            ax5.set_title('Relative Performance\n(vs Fastest Algorithm)', 
                         fontsize=14, fontweight='bold')
            ax5.set_ylabel('Relative Speed Factor')
            ax5.set_xlabel('Algorithm')
            ax5.axhline(y=1, color='red', linestyle='--', alpha=0.7, label='Baseline')
            ax5.legend()
            
            # Add value labels
            for bar, rel_time in zip(bars3, relative_times):
                height = bar.get_height()
                ax5.text(bar.get_x() + bar.get_width()/2., height,
                        f'{rel_time:.2f}x', ha='center', va='bottom', fontsize=10)
        else:
            ax5.text(0.5, 0.5, 'Single Algorithm\nNo Comparison Available', 
                    ha='center', va='center', transform=ax5.transAxes, fontsize=12)
            ax5.set_title('Relative Performance', fontsize=14, fontweight='bold')
        
        # 6. Summary Statistics Table
        ax6 = plt.subplot(2, 3, 6)
        ax6.axis('tight')
        ax6.axis('off')
        
        # Create summary table data
        table_data = []
        headers = ['Algorithm', 'Time (s)', 'Matches', 'Speed Rank']
        
        # Sort algorithms by performance for ranking
        sorted_algos = sorted(algorithms, key=lambda x: results[x]['time'])
        
        for i, algo in enumerate(algorithms):
            rank = sorted_algos.index(algo) + 1
            table_data.append([
                algo.upper(),
                f"{results[algo]['time']:.6f}",
                str(len(results[algo]['matches'])),
                str(rank)
            ])
        
        table = ax6.table(cellText=table_data, colLabels=headers,
                         cellLoc='center', loc='center',
                         colColours=['lightgray'] * len(headers))
        table.auto_set_font_size(False)
        table.set_fontsize(10)
        table.scale(1, 2)
        
        ax6.set_title('Performance Summary', fontsize=14, fontweight='bold', pad=20)
        
        # Add overall title and metadata
        fig.suptitle(f'DNA Pattern Alignment Analysis\nPattern: {pattern} | Sequence Length: {len(dna_sequence):,}', 
                     fontsize=16, fontweight='bold', y=0.98)
        
        plt.tight_layout()
        plt.subplots_adjust(top=0.92)
        
        if save_plots:
            plt.savefig('dna_alignment_analysis.png', dpi=300, bbox_inches='tight')
            print("Plots saved as 'dna_alignment_analysis.png'")
        
        plt.show()
    
    def create_simple_performance_chart(self, results: Dict[str, Dict[str, Any]]):
        """Create a simple performance comparison chart."""
        algorithms = list(results.keys())
        times = [results[algo]['time'] for algo in algorithms]
        matches = [len(results[algo]['matches']) for algo in algorithms]
        
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
        
        # Performance chart
        bars1 = ax1.bar(algorithms, times, color=['#FF6B6B', '#4ECDC4', '#45B7D1'][:len(algorithms)])
        ax1.set_title('Execution Time Comparison')
        ax1.set_ylabel('Time (seconds)')
        
        for bar, time in zip(bars1, times):
            height = bar.get_height()
            ax1.text(bar.get_x() + bar.get_width()/2., height,
                    f'{time:.4f}s', ha='center', va='bottom')
        
        # Matches chart
        bars2 = ax2.bar(algorithms, matches, color=['#FF6B6B', '#4ECDC4', '#45B7D1'][:len(algorithms)])
        ax2.set_title('Matches Found')
        ax2.set_ylabel('Number of Matches')
        
        for bar, count in zip(bars2, matches):
            height = bar.get_height()
            ax2.text(bar.get_x() + bar.get_width()/2., height,
                    f'{count}', ha='center', va='bottom')
        
        plt.tight_layout()
        plt.show()
    
    def detect_file_type(self, file_path: str) -> str:
        """Detect if file is FASTA or TXT format."""
        import os
        
        # Check by extension first
        file_ext = os.path.splitext(file_path)[1].lower()
        if file_ext in ['.fasta', '.fa', '.fas']:
            return 'FASTA'
        elif file_ext in ['.txt', '.text']:
            return 'TXT'
        
        # If extension is unclear, check content
        try:
            with open(file_path, 'r') as f:
                first_line = f.readline().strip()
                if first_line.startswith('>'):
                    return 'FASTA'
                else:
                    return 'TXT'
        except:
            return 'Unknown'
    
    def run_cli(self):
        """Main CLI execution function."""
        parser = argparse.ArgumentParser(description='DNA Pattern Alignment Tool (Supports FASTA and TXT files)')
        parser.add_argument('--file', required=True, help='Path to DNA sequence file (FASTA or TXT)')
        parser.add_argument('--pattern', required=True, help='DNA pattern to search for')
        parser.add_argument('--algorithm', choices=['naive', 'kmp', 'rk', 'all'], 
                          default='all', help='Algorithm to use')
        parser.add_argument('--extend', action='store_true', 
                          help='Enable local extension alignment')
        parser.add_argument('--window', type=int, default=5, 
                          help='Window size for extension (default: 5)')
        parser.add_argument('--plot', action='store_true', 
                          help='Generate matplotlib visualizations')
        parser.add_argument('--save-plots', action='store_true',
                          help='Save plots to file')
        parser.add_argument('--simple-plot', action='store_true',
                          help='Generate simple performance plots only')
        
        args = parser.parse_args()
        
        try:
            # Detect and parse file
            file_type = self.detect_file_type(args.file)
            print(f"Detected file type: {file_type}")
            print(f"Parsing DNA sequence file: {args.file}")
            
            # Use the universal parser that handles both formats
            dna_sequence = self.parser.parse_file(args.file)
            print(f"Loaded sequence of length: {len(dna_sequence)}")
            
            # Validate pattern
            pattern = args.pattern.upper()
            if not self.parser.validate_dna_sequence(pattern):
                raise ValueError("Invalid DNA pattern. Use only A, T, C, G characters.")
            
            print(f"Searching for pattern: {pattern}")
            print("-" * 50)
            
            # Run algorithms
            algorithms = ['naive', 'kmp', 'rk'] if args.algorithm == 'all' else [args.algorithm]
            results = {}
            
            for algo in algorithms:
                exec_time, matches = self.run_algorithm(algo, dna_sequence, pattern)
                results[algo] = {'time': exec_time, 'matches': matches}
                
                print(f"{algo.upper()} Algorithm:")
                print(f"  Time: {exec_time:.6f} seconds")
                print(f"  Matches found: {len(matches)}")
                if matches:
                    print(f"  Positions: {matches[:10]}{'...' if len(matches) > 10 else ''}")
                
                # Perform local extension if requested
                if args.extend and matches:
                    print("  Local alignments:")
                    for i, pos in enumerate(matches[:3]):  # Show first 3 alignments
                        alignment = self.aligner.extend_alignment(
                            dna_sequence, pattern, pos, args.window
                        )
                        print(f"    Position {pos}:")
                        print(f"      Text:    {alignment.extended_text}")
                        print(f"      Pattern: {alignment.extended_pattern}")
                        print(f"      Score: {alignment.score}")
                
                print()
            
            # Text-based performance summary
            if len(results) > 1:
                print("Performance Summary:")
                print("| Algorithm   | Time (s) | Matches |")
                print("|-------------|----------|---------|")
                for algo, data in results.items():
                    print(f"| {algo.ljust(11)} | {data['time']:.6f} | {len(data['matches'])} |")
            
            # Generate plots if requested
            if args.plot or args.simple_plot:
                print("\nGenerating visualizations...")
                if args.simple_plot:
                    self.create_simple_performance_chart(results)
                else:
                    self.visualize_results(results, dna_sequence, pattern, args.save_plots)
        
        except Exception as e:
            print(f"Error: {e}")
            return 1
        
        return 0

if __name__ == '__main__':
    cli = DNACli()
    exit(cli.run_cli())