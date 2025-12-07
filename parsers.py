"""
DNA sequence parser for both FASTA and TXT files.
"""
import os
import re

class DNAParser:
    """Parser for DNA sequence files (FASTA and TXT formats)."""
    
    @staticmethod
    def parse_fasta(file_path: str) -> str:
        """
        Parse FASTA file and return concatenated DNA sequence.
        
        Args:
            file_path: Path to FASTA file
            
        Returns:
            Clean DNA sequence string (uppercase, no whitespace)
        """
        if not os.path.exists(file_path):
            raise FileNotFoundError(f"FASTA file not found: {file_path}")
        
        sequence_lines = []
        
        with open(file_path, 'r') as file:
            for line in file:
                line = line.strip()
                # Skip header lines starting with '>'
                if not line.startswith('>') and line:
                    # Remove any non-DNA characters and convert to uppercase
                    clean_line = re.sub(r'[^ATCG]', '', line.upper())
                    sequence_lines.append(clean_line)
        
        if not sequence_lines:
            raise ValueError("No valid DNA sequence found in FASTA file")
        
        return ''.join(sequence_lines)
    
    @staticmethod
    def parse_txt(file_path: str) -> str:
        """
        Parse TXT file and return concatenated DNA sequence.
        
        Args:
            file_path: Path to TXT file containing DNA sequence
            
        Returns:
            Clean DNA sequence string (uppercase, no whitespace)
        """
        if not os.path.exists(file_path):
            raise FileNotFoundError(f"TXT file not found: {file_path}")
        
        sequence_lines = []
        
        with open(file_path, 'r') as file:
            for line in file:
                line = line.strip()
                if line:  # Skip empty lines
                    # Remove any non-DNA characters and convert to uppercase
                    clean_line = re.sub(r'[^ATCG]', '', line.upper())
                    if clean_line:  # Only add if there are valid DNA characters
                        sequence_lines.append(clean_line)
        
        if not sequence_lines:
            raise ValueError("No valid DNA sequence found in TXT file")
        
        return ''.join(sequence_lines)
    
    @staticmethod
    def parse_file(file_path: str) -> str:
        """
        Parse DNA sequence file (auto-detects FASTA or TXT format).
        
        Args:
            file_path: Path to DNA sequence file
            
        Returns:
            Clean DNA sequence string (uppercase, no whitespace)
        """
        if not os.path.exists(file_path):
            raise FileNotFoundError(f"File not found: {file_path}")
        
        # Determine file type by extension
        file_ext = os.path.splitext(file_path)[1].lower()
        
        if file_ext in ['.fasta', '.fa', '.fas']:
            return DNAParser.parse_fasta(file_path)
        elif file_ext in ['.txt', '.text']:
            return DNAParser.parse_txt(file_path)
        else:
            # Try to auto-detect by checking first line
            with open(file_path, 'r') as file:
                first_line = file.readline().strip()
                if first_line.startswith('>'):
                    return DNAParser.parse_fasta(file_path)
                else:
                    return DNAParser.parse_txt(file_path)
    
    @staticmethod
    def validate_dna_sequence(sequence: str) -> bool:
        """Validate that sequence contains only valid DNA bases."""
        return bool(re.match(r'^[ATCG]+$', sequence.upper()))
    
    @staticmethod
    def get_sequence_stats(sequence: str) -> dict:
        """
        Get basic statistics about the DNA sequence.
        
        Args:
            sequence: DNA sequence string
            
        Returns:
            Dictionary with sequence statistics
        """
        if not DNAParser.validate_dna_sequence(sequence):
            raise ValueError("Invalid DNA sequence")
        
        sequence = sequence.upper()
        length = len(sequence)
        
        return {
            'length': length,
            'A_count': sequence.count('A'),
            'T_count': sequence.count('T'),
            'C_count': sequence.count('C'),
            'G_count': sequence.count('G'),
            'GC_content': (sequence.count('G') + sequence.count('C')) / length * 100 if length > 0 else 0
        }

# Example usage
if __name__ == "__main__":
    # Example for TXT file
    try:
        # Parse a TXT file
        txt_sequence = DNAParser.parse_txt("dna_sequence.txt")
        print(f"TXT Sequence: {txt_sequence[:50]}...")
        print(f"Length: {len(txt_sequence)}")
        
        # Get statistics
        stats = DNAParser.get_sequence_stats(txt_sequence)
        print(f"Stats: {stats}")
        
        # Auto-detect file format
        sequence = DNAParser.parse_file("dna_sequence.txt")
        print(f"Auto-detected sequence: {sequence[:50]}...")
        
    except FileNotFoundError as e:
        print(f"Error: {e}")
    except ValueError as e:
        print(f"Error: {e}")