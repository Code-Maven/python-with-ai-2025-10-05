import Levenshtein


def calculate_word_distance(word1: str, word2: str) -> int:
    """
    Calculate the Levenshtein distance between two words.
    
    The Levenshtein distance is the minimum number of single-character edits
    (insertions, deletions, or substitutions) required to change one word into another.
    
    Uses the optimized python-Levenshtein library for fast computation.
    
    Args:
        word1 (str): First word
        word2 (str): Second word
    
    Returns:
        int: The Levenshtein distance between the two words
    """
    # Convert to lowercase for case-insensitive comparison
    word1 = word1.lower().strip()
    word2 = word2.lower().strip()
    
    # Use the optimized Levenshtein library
    return Levenshtein.distance(word1, word2)


def main() -> None:
    try:
        # Get input from user
        print("Word Distance Calculator")
        print("-" * 25)
        print("This program calculates the Levenshtein distance between two words.")
        print("The distance is the minimum number of edits needed to transform one word into another.\n")
        
        word1 = input("Enter the first word: ").strip()
        word2 = input("Enter the second word: ").strip()
        
        # Validate input
        if not word1 or not word2:
            print("Error: Both words must be non-empty.")
            return
            
        if not word1.replace(' ', '').isalpha() or not word2.replace(' ', '').isalpha():
            print("Error: Words should contain only alphabetic characters and spaces.")
            return
        
        # Calculate distance
        distance = calculate_word_distance(word1, word2)
        
        # Print results
        print(f"\nWord Distance Analysis:")
        print(f"First word: '{word1}'")
        print(f"Second word: '{word2}'")
        print(f"Levenshtein distance: {distance}")
        
        # Provide interpretation
        if distance == 0:
            print("Interpretation: The words are identical.")
        elif distance == 1:
            print("Interpretation: The words differ by exactly one character.")
        elif distance <= 3:
            print("Interpretation: The words are quite similar.")
        else:
            print("Interpretation: The words are quite different.")
        
    except KeyboardInterrupt:
        print("\nProgram interrupted by user.")


if __name__ == "__main__":
    main()