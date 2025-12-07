"""
Gradio-based web interface for DNA pattern alignment tool with matplotlib visualizations.
"""
try:
    import gradio as gr
    GRADIO_AVAILABLE = True
except ImportError:
    GRADIO_AVAILABLE = False

import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg')  # Use non-interactive backend for web
import numpy as np
import io
import base64
from typing import Dict, Any, List, Tuple
from parsers import DNAParser
from matcher import PatternMatcher
from aligner import LocalAligner
from benchmark import Benchmark

class DNAGradioUI:
    """Gradio-based web interface for DNA pattern alignment tool with visualizations."""
    
    def __init__(self):
        self.parser = DNAParser()
        self.matcher = PatternMatcher()
        self.aligner = LocalAligner()
        self.benchmark = Benchmark()
    
    def create_performance_chart(self, results: Dict[str, Dict[str, Any]]) -> plt.Figure:
        """Create performance comparison chart."""
        algorithms = list(results.keys())
        times = [results[algo]['time'] for algo in algorithms]
        matches = [len(results[algo]['matches']) for algo in algorithms]
        
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
        
        # Performance chart
        colors = ['#FF6B6B', '#4ECDC4', '#45B7D1']
        bars1 = ax1.bar(algorithms, times, color=colors[:len(algorithms)])
        ax1.set_title('Algorithm Performance Comparison', fontweight='bold', fontsize=14)
        ax1.set_ylabel('Execution Time (seconds)')
        ax1.set_xlabel('Algorithm')
        
        # Add value labels on bars
        for bar, time in zip(bars1, times):
            height = bar.get_height()
            ax1.text(bar.get_x() + bar.get_width()/2., height,
                    f'{time:.6f}s', ha='center', va='bottom', fontsize=10)
        
        # Matches chart
        bars2 = ax2.bar(algorithms, matches, color=colors[:len(algorithms)])
        ax2.set_title('Matches Found by Algorithm', fontweight='bold', fontsize=14)
        ax2.set_ylabel('Number of Matches')
        ax2.set_xlabel('Algorithm')
        
        for bar, count in zip(bars2, matches):
            height = bar.get_height()
            ax2.text(bar.get_x() + bar.get_width()/2., height,
                    f'{count}', ha='center', va='bottom', fontsize=10)
        
        plt.tight_layout()
        return fig
    
    def create_comprehensive_analysis(self, results: Dict[str, Dict[str, Any]], 
                                    dna_sequence: str, pattern: str) -> plt.Figure:
        """Create comprehensive analysis dashboard."""
        fig = plt.figure(figsize=(16, 12))
        algorithms = list(results.keys())
        times = [results[algo]['time'] for algo in algorithms]
        matches = [len(results[algo]['matches']) for algo in algorithms]
        colors = ['#FF6B6B', '#4ECDC4', '#45B7D1', '#FFD93D']
        
        # 1. Performance Bar Chart
        ax1 = plt.subplot(2, 3, 1)
        bars1 = ax1.bar(algorithms, times, color=colors[:len(algorithms)])
        ax1.set_title('Execution Time Comparison', fontweight='bold')
        ax1.set_ylabel('Time (seconds)')
        
        for bar, time in zip(bars1, times):
            height = bar.get_height()
            ax1.text(bar.get_x() + bar.get_width()/2., height,
                    f'{time:.4f}s', ha='center', va='bottom', fontsize=9)
        
        # 2. Matches Bar Chart
        ax2 = plt.subplot(2, 3, 2)
        bars2 = ax2.bar(algorithms, matches, color=colors[:len(algorithms)])
        ax2.set_title('Matches Found', fontweight='bold')
        ax2.set_ylabel('Number of Matches')
        
        for bar, count in zip(bars2, matches):
            height = bar.get_height()
            ax2.text(bar.get_x() + bar.get_width()/2., height,
                    f'{count}', ha='center', va='bottom', fontsize=9)
        
        # 3. Performance vs Matches Scatter
        ax3 = plt.subplot(2, 3, 3)
        scatter = ax3.scatter(times, matches, c=colors[:len(algorithms)], 
                             s=100, alpha=0.7, edgecolors='black')
        
        for i, algo in enumerate(algorithms):
            ax3.annotate(algo.upper(), (times[i], matches[i]), 
                        xytext=(5, 5), textcoords='offset points', fontsize=10)
        
        ax3.set_title('Performance vs Accuracy', fontweight='bold')
        ax3.set_xlabel('Execution Time (seconds)')
        ax3.set_ylabel('Matches Found')
        ax3.grid(True, alpha=0.3)
        
        # 4. Match Distribution Histogram
        ax4 = plt.subplot(2, 3, 4)
        if results and algorithms:
            first_algo = algorithms[0]
            algo_matches = results[first_algo]['matches']
            
            if algo_matches:
                bins = min(30, len(dna_sequence) // 1000 + 1)
                ax4.hist(algo_matches, bins=bins, alpha=0.7, color=colors[0], edgecolor='black')
                ax4.set_title(f'Match Distribution\n({first_algo.upper()})', fontweight='bold')
                ax4.set_xlabel('Sequence Position')
                ax4.set_ylabel('Frequency')
            else:
                ax4.text(0.5, 0.5, 'No matches found', ha='center', va='center', 
                        transform=ax4.transAxes, fontsize=12)
                ax4.set_title('Match Distribution', fontweight='bold')
        
        # 5. Relative Performance
        ax5 = plt.subplot(2, 3, 5)
        if len(algorithms) > 1:
            min_time = min(times)
            relative_times = [t / min_time for t in times]
            
            bars3 = ax5.bar(algorithms, relative_times, color=colors[:len(algorithms)])
            ax5.set_title('Relative Performance\n(vs Fastest)', fontweight='bold')
            ax5.set_ylabel('Speed Factor')
            ax5.axhline(y=1, color='red', linestyle='--', alpha=0.7)
            
            for bar, rel_time in zip(bars3, relative_times):
                height = bar.get_height()
                ax5.text(bar.get_x() + bar.get_width()/2., height,
                        f'{rel_time:.2f}x', ha='center', va='bottom', fontsize=9)
        else:
            ax5.text(0.5, 0.5, 'Single Algorithm', ha='center', va='center', 
                    transform=ax5.transAxes, fontsize=12)
            ax5.set_title('Relative Performance', fontweight='bold')
        
        # 6. Sequence Statistics Pie Chart
        ax6 = plt.subplot(2, 3, 6)
        if dna_sequence:
            base_counts = {'A': dna_sequence.count('A'), 'T': dna_sequence.count('T'),
                          'C': dna_sequence.count('C'), 'G': dna_sequence.count('G')}
            
            labels = list(base_counts.keys())
            sizes = list(base_counts.values())
            colors_pie = ['#FF9999', '#66B2FF', '#99FF99', '#FFCC99']
            
            wedges, texts, autotexts = ax6.pie(sizes, labels=labels, colors=colors_pie,
                                              autopct='%1.1f%%', startangle=90)
            ax6.set_title('DNA Base Composition', fontweight='bold')
        
        plt.suptitle(f'DNA Pattern Analysis Dashboard\nPattern: {pattern} | Sequence Length: {len(dna_sequence):,}', 
                     fontsize=16, fontweight='bold', y=0.98)
        plt.tight_layout()
        plt.subplots_adjust(top=0.92)
        
        return fig
    
    def create_match_timeline(self, results: Dict[str, Dict[str, Any]], 
                             dna_sequence: str, pattern: str) -> plt.Figure:
        """Create match timeline visualization."""
        fig, ax = plt.subplots(figsize=(14, 6))
        
        if not results:
            ax.text(0.5, 0.5, 'No results to display', ha='center', va='center', 
                   transform=ax.transAxes, fontsize=14)
            return fig
        
        # Use first algorithm's matches for timeline
        first_algo = list(results.keys())[0]
        matches = results[first_algo]['matches']
        
        if matches:
            # Create timeline plot
            y_pos = [1] * len(matches)
            colors = plt.cm.viridis(np.linspace(0, 1, len(matches)))
            
            scatter = ax.scatter(matches, y_pos, c=colors, s=50, alpha=0.7)
            
            # Add some statistical lines
            if len(matches) > 1:
                ax.axvline(x=np.mean(matches), color='red', linestyle='--', 
                          label=f'Mean: {np.mean(matches):.0f}', alpha=0.7)
                ax.axvline(x=np.median(matches), color='orange', linestyle='--', 
                          label=f'Median: {np.median(matches):.0f}', alpha=0.7)
            
            ax.set_xlim(0, len(dna_sequence))
            ax.set_ylim(0.5, 1.5)
            ax.set_xlabel('Sequence Position')
            ax.set_title(f'Match Timeline - {pattern} ({first_algo.upper()})\n{len(matches)} matches found', 
                        fontweight='bold')
            ax.legend()
            ax.grid(True, alpha=0.3)
            
            # Remove y-axis ticks as they're not meaningful
            ax.set_yticks([])
        else:
            ax.text(0.5, 0.5, 'No matches found', ha='center', va='center', 
                   transform=ax.transAxes, fontsize=14)
            ax.set_title(f'Match Timeline - {pattern}', fontweight='bold')
        
        plt.tight_layout()
        return fig
    
    def process_file_and_pattern(self, file_obj, pattern: str, algorithm: str, 
                                enable_extension: bool, window_size: int):
        """Main processing function for Gradio interface with visualizations."""
        try:
            if file_obj is None:
                return "Please upload a FASTA file.", "", "", None, None, None
            
            # Read uploaded file
            dna_sequence = self.parser.parse_fasta(file_obj.name)
            
            # Validate pattern
            pattern = pattern.upper().strip()
            if not pattern:
                return "Please enter a DNA pattern.", "", "", None, None, None
            
            if not self.parser.validate_dna_sequence(pattern):
                return "Invalid DNA pattern. Use only A, T, C, G characters.", "", "", None, None, None
            
            # Run algorithms
            algorithms = ['naive', 'kmp', 'rk'] if algorithm == 'all' else [algorithm]
            results = {}
            all_matches = []
            
            for algo in algorithms:
                if algo == 'naive':
                    exec_time, matches = self.benchmark.time_function(
                        self.matcher.naive_match, dna_sequence, pattern)
                elif algo == 'kmp':
                    exec_time, matches = self.benchmark.time_function(
                        self.matcher.kmp_match, dna_sequence, pattern)
                elif algo == 'rk':
                    exec_time, matches = self.benchmark.time_function(
                        self.matcher.rabin_karp_match, dna_sequence, pattern)
                
                results[algo] = {'time': exec_time, 'matches': matches}
                if not all_matches:  # Use first algorithm's matches for display
                    all_matches = matches
            
            # Generate results summary
            summary = f"🧬 DNA SEQUENCE ANALYSIS RESULTS\n"
            summary += f"{'='*50}\n"
            summary += f"Sequence length: {len(dna_sequence):,} bases\n"
            summary += f"Pattern searched: {pattern}\n"
            summary += f"Total matches found: {len(all_matches)}\n\n"
            
            if all_matches:
                summary += f"📍 Match positions (first 20): {all_matches[:20]}"
                if len(all_matches) > 20:
                    summary += f"\n   ... and {len(all_matches) - 20} more positions"
                summary += "\n\n"
            
            # Performance table
            table_data = []
            for algo, data in results.items():
                table_data.append([algo.upper(), f"{data['time']:.6f}", len(data['matches'])])
            
            performance_table = "⚡ ALGORITHM PERFORMANCE COMPARISON\n"
            performance_table += "="*50 + "\n"
            performance_table += "| Algorithm | Time (seconds) | Matches Found |\n"
            performance_table += "|-----------|----------------|---------------|\n"
            for row in table_data:
                performance_table += f"| {row[0]:9} | {row[1]:14} | {row[2]:13} |\n"
            
            # Find fastest algorithm
            if len(results) > 1:
                fastest = min(results.keys(), key=lambda x: results[x]['time'])
                performance_table += f"\n🏆 Fastest Algorithm: {fastest.upper()}\n"
            
            # Alignment results
            alignment_results = ""
            if enable_extension and all_matches:
                alignment_results = "🔍 LOCAL ALIGNMENT EXTENSIONS\n"
                alignment_results += "="*50 + "\n"
                alignment_results += f"Showing first 5 matches with {window_size}-base extensions:\n\n"
                
                for i, pos in enumerate(all_matches[:5]):
                    alignment = self.aligner.extend_alignment(
                        dna_sequence, pattern, pos, window_size
                    )
                    alignment_results += f"Match #{i+1} at position {pos}:\n"
                    alignment_results += f"  Extended Text:    {alignment.extended_text}\n"
                    alignment_results += f"  Extended Pattern: {alignment.extended_pattern}\n"
                    alignment_results += f"  Alignment Score:  {alignment.score}\n"
                    alignment_results += "-" * 40 + "\n"
            
            # Create visualizations
            performance_chart = self.create_performance_chart(results)
            comprehensive_chart = self.create_comprehensive_analysis(results, dna_sequence, pattern)
            timeline_chart = self.create_match_timeline(results, dna_sequence, pattern)
            
            return (summary, performance_table, alignment_results, 
                   performance_chart, comprehensive_chart, timeline_chart)
        
        except Exception as e:
            error_msg = f"❌ ERROR: {str(e)}"
            return error_msg, "", "", None, None, None
    
    def create_interface(self):
        """Create and return Gradio interface with visualizations."""
        if not GRADIO_AVAILABLE:
            raise ImportError("Gradio is not available. Install with 'pip install gradio'")
        
        # Custom CSS for better styling
        custom_css = """
        .gradio-container {
            max-width: 1200px !important;
        }
        .plot-container {
            border: 2px solid #e1e5e9;
            border-radius: 8px;
            padding: 10px;
            margin: 10px 0;
        }
        """
        
        with gr.Blocks(title="DNA Pattern Alignment Tool") as interface:
            gr.Markdown("""
            # 🧬 DNA Pattern Alignment Analysis Tool
            
            Upload a FASTA file and search for DNA patterns using different algorithms. 
            Get comprehensive performance analysis with interactive visualizations!
            """)
            
            with gr.Row():
                with gr.Column(scale=1):
                    gr.Markdown("### 📁 Input Configuration")
                    
                    file_input = gr.File(
                        label="Upload FASTA File",
                        file_types=[".fasta", ".fa", ".fas", ".txt"],
                        file_count="single"
                    )
                    
                    pattern_input = gr.Textbox(
                        label="DNA Pattern to Search",
                        placeholder="Enter DNA sequence (e.g., ATCG, GCTA)",
                        lines=1,
                        info="Use only A, T, C, G characters"
                    )
                    
                    algorithm_dropdown = gr.Dropdown(
                        choices=["all", "naive", "kmp", "rk"],
                        label="Algorithm Selection",
                        value="all",
                        info="Choose specific algorithm or run all for comparison"
                    )
                    
                    with gr.Row():
                        extend_checkbox = gr.Checkbox(
                            label="Enable Local Extension Analysis",
                            value=False,
                            info="Extend matches for local alignment"
                        )
                        window_slider = gr.Slider(
                            minimum=1,
                            maximum=20,
                            value=5,
                            step=1,
                            label="Extension Window Size",
                            info="Number of bases to extend on each side"
                        )
                    
                    submit_btn = gr.Button("🔍 Analyze DNA Pattern", 
                                         variant="primary", 
                                         size="lg")
                
                with gr.Column(scale=2):
                    gr.Markdown("### 📊 Analysis Results")
                    
                    with gr.Tabs():
                        with gr.TabItem("📈 Summary"):
                            results_output = gr.Textbox(
                                label="Results Summary",
                                lines=12,
                                max_lines=15,
                                
                            )
                        
                        with gr.TabItem("⚡ Performance"):
                            performance_output = gr.Textbox(
                                label="Performance Comparison",
                                lines=10,
                                
                            )
                        
                        with gr.TabItem("🔍 Alignments"):
                            alignment_output = gr.Textbox(
                                label="Local Alignment Extensions",
                                lines=15,
                                max_lines=20,
                                
                            )
            
            # Visualization section
            gr.Markdown("### 📊 Performance Visualizations")
            
            with gr.Tabs():
                with gr.TabItem("📊 Performance Charts"):
                    performance_plot = gr.Plot(
                        label="Algorithm Performance Comparison",
                        elem_classes=["plot-container"]
                    )
                
                with gr.TabItem("📈 Comprehensive Analysis"):
                    comprehensive_plot = gr.Plot(
                        label="Complete Analysis Dashboard",
                        elem_classes=["plot-container"]
                    )
                
                with gr.TabItem("📍 Match Timeline"):
                    timeline_plot = gr.Plot(
                        label="Match Distribution Timeline",
                        elem_classes=["plot-container"]
                    )
            
            # Examples section
            gr.Markdown("### 💡 Example Patterns to Try")
            gr.Examples(
                examples=[
                    ["ATCG", "all", False, 5],
                    ["GCTA", "kmp", True, 8],
                    ["AAAA", "rk", False, 3],
                    ["CGCG", "all", True, 10]
                ],
                inputs=[pattern_input, algorithm_dropdown, extend_checkbox, window_slider],
                label="Click to load example configurations"
            )
            
            # Connect the submit button to the processing function
            submit_btn.click(
                fn=self.process_file_and_pattern,
                inputs=[
                    file_input,
                    pattern_input,
                    algorithm_dropdown,
                    extend_checkbox,
                    window_slider
                ],
                outputs=[
                    results_output, 
                    performance_output, 
                    alignment_output,
                    performance_plot,
                    comprehensive_plot,
                    timeline_plot
                ]
            )
            
            gr.Markdown("""
            ---
            ### ℹ️ Algorithm Information
            - **Naive**: Simple brute-force pattern matching
            - **KMP**: Knuth-Morris-Pratt algorithm with preprocessing
            - **RK**: Rabin-Karp rolling hash algorithm
            """)
        
        return interface

# Usage example
if __name__ == "__main__":
    if GRADIO_AVAILABLE:
        ui = DNAGradioUI()
        interface = ui.create_interface()
        interface.launch(share=True, debug=True)
    else:
        print("Gradio is not installed. Please install with: pip install gradio")