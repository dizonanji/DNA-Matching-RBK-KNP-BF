# ===== benchmark.py =====
"""
Timing and performance measurement utilities.
"""
import time
from typing import Tuple, Callable, Any

class Benchmark:
    """Timing and performance measurement utilities."""
    
    @staticmethod
    def time_function(func: Callable, *args, **kwargs) -> Tuple[float, Any]:
        """
        Time a function execution.
        
        Args:
            func: Function to time
            *args: Function arguments
            **kwargs: Function keyword arguments
            
        Returns:
            Tuple of (execution_time, result)
        """
        start_time = time.time()
        result = func(*args, **kwargs)
        end_time = time.time()
        
        return end_time - start_time, result

