"""
Pattern matching algorithms for DNA sequences.
"""
from typing import List

class PatternMatcher:
    """Collection of pattern matching algorithms for DNA sequences."""
    
    @staticmethod
    def naive_match(text: str, pattern: str) -> List[int]:
        """
        Naive (brute force) pattern matching algorithm.
        
        Args:
            text: DNA sequence to search in
            pattern: Pattern to search for
            
        Returns:
            List of starting positions where pattern matches
        """
        matches = []
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
            text: DNA sequence to search in
            pattern: Pattern to search for
            
        Returns:
            List of starting positions where pattern matches
        """
        matches = []
        text_len = len(text)
        pattern_len = len(pattern)
        
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
            text: DNA sequence to search in
            pattern: Pattern to search for
            prime: Prime number for hashing
            
        Returns:
            List of starting positions where pattern matches
        """
        matches = []
        text_len = len(text)
        pattern_len = len(pattern)
        
        if pattern_len > text_len:
            return matches
        
        # Hash mapping for DNA bases
        base_hash = {'A': 1, 'T': 2, 'C': 3, 'G': 4}
        base = 4  # Number of characters in DNA alphabet
        
        # Calculate hash values
        pattern_hash = 0
        text_hash = 0
        h = 1
        
        # Calculate h = pow(base, pattern_len-1) % prime
        for i in range(pattern_len - 1):
            h = (h * base) % prime
        
        # Calculate hash value of pattern and first window of text
        for i in range(pattern_len):
            pattern_hash = (base * pattern_hash + base_hash[pattern[i]]) % prime
            text_hash = (base * text_hash + base_hash[text[i]]) % prime
        
        # Slide the pattern over text one by one
        for i in range(text_len - pattern_len + 1):
            # Check if hash values match
            if pattern_hash == text_hash:
                # Check characters one by one for exact match
                if text[i:i + pattern_len] == pattern:
                    matches.append(i)
            
            # Calculate hash value for next window
            if i < text_len - pattern_len:
                text_hash = (base * (text_hash - base_hash[text[i]] * h) + 
                           base_hash[text[i + pattern_len]]) % prime
                
                # Convert negative hash to positive
                if text_hash < 0:
                    text_hash += prime
        
        return matches

