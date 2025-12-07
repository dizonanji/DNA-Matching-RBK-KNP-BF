"""
Local sequence alignment with extension and scoring.
"""
from dataclasses import dataclass

@dataclass
class AlignmentResult:
    """Data class for alignment results."""
    position: int
    extended_text: str
    extended_pattern: str
    score: int
    alignment_length: int


class LocalAligner:
    """Local sequence alignment with simple scoring."""
    
    def __init__(self, match_score: int = 1, mismatch_score: int = -1, gap_score: int = -2):
        self.match_score = match_score
        self.mismatch_score = mismatch_score
        self.gap_score = gap_score
    
    def extend_alignment(self, text: str, pattern: str, match_pos: int, 
                        window_size: int) -> AlignmentResult:
        """
        Perform local extension around a match position.
        
        Args:
            text: Full DNA sequence
            pattern: Original pattern
            match_pos: Position where pattern was found
            window_size: Number of bases to extend on each side
            
        Returns:
            AlignmentResult with extended alignment and score
        """
        text_len = len(text)
        pattern_len = len(pattern)
        
        # Calculate extension boundaries
        left_start = max(0, match_pos - window_size)
        right_end = min(text_len, match_pos + pattern_len + window_size)
        
        # Extract extended regions
        extended_text = text[left_start:right_end]
        
        # Create extended pattern with gaps for unmatched regions
        pattern_start = match_pos - left_start
        pattern_end = pattern_start + pattern_len
        
        extended_pattern = '-' * pattern_start + pattern + '-' * (len(extended_text) - pattern_end)
        
        # Calculate alignment score
        score = self._calculate_score(extended_text, extended_pattern)
        
        return AlignmentResult(
            position=match_pos,
            extended_text=extended_text,
            extended_pattern=extended_pattern,
            score=score,
            alignment_length=len(extended_text)
        )
    
    def _calculate_score(self, seq1: str, seq2: str) -> int:
        """Calculate alignment score between two sequences."""
        score = 0
        
        for i in range(len(seq1)):
            if i < len(seq2):
                if seq2[i] == '-':
                    score += self.gap_score
                elif seq1[i] == seq2[i]:
                    score += self.match_score
                else:
                    score += self.mismatch_score
            else:
                score += self.gap_score
        
        return score