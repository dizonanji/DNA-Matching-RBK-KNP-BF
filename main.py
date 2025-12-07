# ===== main.py =====
"""
Main entry point for the DNA Pattern Alignment Tool.
Supports both FASTA and TXT files in CLI and UI modes.
"""
import sys
import os
from c import DNACli
from ui import DNAGradioUI, GRADIO_AVAILABLE


def print_usage():
    """Print usage instructions."""
    print("\n" + "="*60)
    print("DNA Pattern Alignment Tool")
    print("="*60)
    print("Supports both FASTA (.fasta, .fa, .fas) and TXT (.txt) files")
    print("\nCLI Usage:")
    print("python main.py --file <sequence_file> --pattern <DNA_pattern> [options]")
    print("\nExamples:")
    print("  # Basic usage with TXT file")
    print("  python main.py --file sequence.txt --pattern ATCG")
    print("")
    print("  # Use specific algorithm with FASTA file")
    print("  python main.py --file data.fasta --pattern GCTA --algorithm kmp")
    print("")
    print("  # Run all algorithms with visualization")
    print("  python main.py --file dna.txt --pattern AATTCC --algorithm all --plot")
    print("")
    print("  # With local alignment extension")
    print("  python main.py --file genome.fasta --pattern ATCG --extend --window 10")
    print("")
    print("Available options:")
    print("  --file          Path to DNA sequence file (FASTA or TXT)")
    print("  --pattern       DNA pattern to search for (A, T, C, G only)")
    print("  --algorithm     Algorithm: naive, kmp, rk, or all (default: all)")
    print("  --extend        Enable local extension alignment")
    print("  --window        Window size for extension (default: 5)")
    print("  --plot          Generate comprehensive visualizations")
    print("  --simple-plot   Generate simple performance plots")
    print("  --save-plots    Save plots to file")
    print("")
    print("UI Mode:")
    print("  python main.py    (no arguments - launches web interface)")
    print("="*60)


def check_dependencies():
    """Check if required modules are available."""
    missing_deps = []
    
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        missing_deps.append("matplotlib")
    
    try:
        import numpy as np
    except ImportError:
        missing_deps.append("numpy")
    
    # Check for optional dependencies
    optional_missing = []
    if not GRADIO_AVAILABLE:
        optional_missing.append("gradio")
    
    if missing_deps:
        print("ERROR: Missing required dependencies:")
        for dep in missing_deps:
            print(f"  - {dep}")
        print("\nInstall with: pip install " + " ".join(missing_deps))
        return False
    
    if optional_missing:
        print("NOTE: Optional dependencies not available:")
        for dep in optional_missing:
            print(f"  - {dep} (for web UI)")
        print("Install with: pip install " + " ".join(optional_missing))
    
    return True


def main():
    """Main entry point for the application."""
    try:
        if len(sys.argv) > 1:
            # CLI mode
            print("DNA Pattern Alignment Tool - CLI Mode")
            print("Supports FASTA and TXT files")
            print("-" * 40)
            
            # Check for help flags
            if '--help' in sys.argv or '-h' in sys.argv:
                print_usage()
                return 0
            
            # Check dependencies
            if not check_dependencies():
                return 1
            
            cli = DNACli()
            return cli.run_cli()
        else:
            # UI mode
            print("DNA Pattern Alignment Tool - UI Mode")
            print("-" * 40)
            
            if GRADIO_AVAILABLE:
                print("Launching web interface...")
                print("The interface supports both FASTA and TXT files.")
                print("You can upload files or paste sequences directly.")
                
                ui = DNAGradioUI()
                interface = ui.create_interface()
                interface.launch(share=False, inbrowser=True)
            else:
                print("Gradio not available for web UI.")
                print("You can still use the CLI mode.")
                print_usage()
                print("\nTo install Gradio for web interface:")
                print("pip install gradio")
                return 1
        
        return 0
    
    except KeyboardInterrupt:
        print("\n\nOperation cancelled by user.")
        return 0
    except ImportError as e:
        print(f"Import error: {e}")
        print("Please ensure all required modules are installed.")
        return 1
    except Exception as e:
        print(f"Unexpected error: {e}")
        return 1


if __name__ == "__main__":
    exit(main())