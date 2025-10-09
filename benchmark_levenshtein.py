"""
Benchmark different Levenshtein distance implementations for performance comparison.
"""
import time
import random
import string
from typing import List, Tuple, Callable


def generate_test_data(num_pairs: int = 1000, min_length: int = 5, max_length: int = 50) -> List[Tuple[str, str]]:
    """Generate random string pairs for testing."""
    test_pairs = []
    for _ in range(num_pairs):
        len1 = random.randint(min_length, max_length)
        len2 = random.randint(min_length, max_length)
        
        str1 = ''.join(random.choices(string.ascii_lowercase, k=len1))
        str2 = ''.join(random.choices(string.ascii_lowercase, k=len2))
        
        test_pairs.append((str1, str2))
    
    return test_pairs


def benchmark_implementation(name: str, func: Callable, test_data: List[Tuple[str, str]]) -> Tuple[float, bool]:
    """Benchmark a Levenshtein implementation."""
    try:
        start_time = time.time()
        
        for str1, str2 in test_data:
            _ = func(str1, str2)
        
        end_time = time.time()
        elapsed = end_time - start_time
        
        print(f"{name:20}: {elapsed:.4f} seconds ({len(test_data)/elapsed:.0f} ops/sec)")
        return elapsed, True
        
    except ImportError as e:
        print(f"{name:20}: Not available ({e})")
        return float('inf'), False
    except Exception as e:
        print(f"{name:20}: Error ({e})")
        return float('inf'), False


def pure_python_levenshtein(s1: str, s2: str) -> int:
    """Pure Python implementation (our original algorithm)."""
    if len(s1) < len(s2):
        return pure_python_levenshtein(s2, s1)

    if len(s2) == 0:
        return len(s1)

    previous_row = list(range(len(s2) + 1))
    for i, c1 in enumerate(s1):
        current_row = [i + 1]
        for j, c2 in enumerate(s2):
            insertions = previous_row[j + 1] + 1
            deletions = current_row[j] + 1
            substitutions = previous_row[j] + (c1 != c2)
            current_row.append(min(insertions, deletions, substitutions))
        previous_row = current_row
    
    return previous_row[-1]


def main():
    print("Levenshtein Distance Implementation Benchmark")
    print("=" * 50)
    
    # Generate test data
    print("Generating test data...")
    test_data = generate_test_data(num_pairs=1000, min_length=10, max_length=30)
    print(f"Created {len(test_data)} string pairs for testing\n")
    
    results = {}
    
    # Test Pure Python implementation
    results["Pure Python"] = benchmark_implementation(
        "Pure Python", 
        pure_python_levenshtein, 
        test_data
    )
    
    # Test python-Levenshtein (current)
    def test_python_levenshtein(s1, s2):
        import Levenshtein
        return Levenshtein.distance(s1, s2)
    
    results["python-Levenshtein"] = benchmark_implementation(
        "python-Levenshtein",
        test_python_levenshtein,
        test_data
    )
    
    # Test RapidFuzz
    def test_rapidfuzz(s1, s2):
        from rapidfuzz.distance import Levenshtein as RFLevenshtein
        return RFLevenshtein.distance(s1, s2)
    
    results["RapidFuzz"] = benchmark_implementation(
        "RapidFuzz",
        test_rapidfuzz,
        test_data
    )
    
    # Test editdistance
    def test_editdistance(s1, s2):
        import editdistance
        return editdistance.eval(s1, s2)
    
    results["editdistance"] = benchmark_implementation(
        "editdistance",
        test_editdistance,
        test_data
    )
    
    # Test difflib (standard library)
    def test_difflib(s1, s2):
        # Note: difflib doesn't directly compute Levenshtein, this is approximate
        import difflib
        return len([x for x in difflib.unified_diff(s1, s2) if x.startswith('+')])
    
    # results["difflib (approx)"] = benchmark_implementation(
    #     "difflib (approx)",
    #     test_difflib,
    #     test_data
    # )
    
    print("\nPerformance Summary:")
    print("-" * 30)
    
    # Sort by performance (excluding failed imports)
    successful_results = [(name, time_taken) for name, (time_taken, success) in results.items() if success]
    successful_results.sort(key=lambda x: x[1])
    
    if successful_results:
        fastest_name, fastest_time = successful_results[0]
        print(f"Fastest: {fastest_name}")
        
        for name, time_taken in successful_results:
            if name == fastest_name:
                print(f"  {name}: baseline")
            else:
                slowdown = time_taken / fastest_time if fastest_time > 0 else 1
                print(f"  {name}: {slowdown:.1f}x slower")


if __name__ == "__main__":
    main()