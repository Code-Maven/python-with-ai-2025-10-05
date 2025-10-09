import pytest
from unittest.mock import patch
import sys
import io
from word_distance import calculate_word_distance, main


# Tests for calculate_word_distance function

def test_calculate_word_distance_identical_words():
    """Test distance between identical words."""
    distance = calculate_word_distance("hello", "hello")
    assert distance == 0


def test_calculate_word_distance_case_insensitive():
    """Test that distance calculation is case-insensitive."""
    distance = calculate_word_distance("Hello", "HELLO")
    assert distance == 0


def test_calculate_word_distance_single_substitution():
    """Test distance with single character substitution."""
    distance = calculate_word_distance("cat", "bat")
    assert distance == 1


def test_calculate_word_distance_single_insertion():
    """Test distance with single character insertion."""
    distance = calculate_word_distance("cat", "cats")
    assert distance == 1


def test_calculate_word_distance_single_deletion():
    """Test distance with single character deletion."""
    distance = calculate_word_distance("cats", "cat")
    assert distance == 1


def test_calculate_word_distance_multiple_operations():
    """Test distance with multiple edit operations."""
    distance = calculate_word_distance("kitten", "sitting")
    assert distance == 3


def test_calculate_word_distance_completely_different():
    """Test distance between completely different words."""
    distance = calculate_word_distance("abc", "xyz")
    assert distance == 3


def test_calculate_word_distance_empty_to_word():
    """Test distance from empty string to word."""
    distance = calculate_word_distance("", "hello")
    assert distance == 5


def test_calculate_word_distance_word_to_empty():
    """Test distance from word to empty string."""
    distance = calculate_word_distance("hello", "")
    assert distance == 5


def test_calculate_word_distance_both_empty():
    """Test distance between two empty strings."""
    distance = calculate_word_distance("", "")
    assert distance == 0


def test_calculate_word_distance_with_spaces():
    """Test distance calculation with words containing spaces."""
    distance = calculate_word_distance("hello world", "hello world")
    assert distance == 0


def test_calculate_word_distance_whitespace_handling():
    """Test that leading/trailing whitespace is handled."""
    distance = calculate_word_distance("  hello  ", "hello")
    assert distance == 0


# Tests for main function

def test_main_valid_input():
    """Test main function with valid input."""
    user_inputs = ["hello", "world"]
    expected_output = [
        "Word Distance Calculator",
        "Word Distance Analysis:",
        "First word: 'hello'",
        "Second word: 'world'",
        "Levenshtein distance: 4"
    ]
    
    with patch('builtins.input', side_effect=user_inputs):
        with patch('sys.stdout', new_callable=io.StringIO) as mock_stdout:
            main()
            output = mock_stdout.getvalue()
            
            for expected_line in expected_output:
                assert expected_line in output


def test_main_identical_words():
    """Test main function with identical words."""
    user_inputs = ["test", "test"]
    
    with patch('builtins.input', side_effect=user_inputs):
        with patch('sys.stdout', new_callable=io.StringIO) as mock_stdout:
            main()
            output = mock_stdout.getvalue()
            
            assert "Levenshtein distance: 0" in output
            assert "The words are identical" in output


def test_main_empty_first_word():
    """Test main function with empty first word."""
    user_inputs = ["", "hello"]
    
    with patch('builtins.input', side_effect=user_inputs):
        with patch('sys.stdout', new_callable=io.StringIO) as mock_stdout:
            main()
            output = mock_stdout.getvalue()
            
            assert "Error: Both words must be non-empty" in output


def test_main_empty_second_word():
    """Test main function with empty second word."""
    user_inputs = ["hello", ""]
    
    with patch('builtins.input', side_effect=user_inputs):
        with patch('sys.stdout', new_callable=io.StringIO) as mock_stdout:
            main()
            output = mock_stdout.getvalue()
            
            assert "Error: Both words must be non-empty" in output


def test_main_whitespace_only_words():
    """Test main function with whitespace-only words."""
    user_inputs = ["   ", "hello"]
    
    with patch('builtins.input', side_effect=user_inputs):
        with patch('sys.stdout', new_callable=io.StringIO) as mock_stdout:
            main()
            output = mock_stdout.getvalue()
            
            assert "Error: Both words must be non-empty" in output


def test_main_invalid_characters():
    """Test main function with invalid characters."""
    user_inputs = ["hello123", "world"]
    
    with patch('builtins.input', side_effect=user_inputs):
        with patch('sys.stdout', new_callable=io.StringIO) as mock_stdout:
            main()
            output = mock_stdout.getvalue()
            
            assert "Error: Words should contain only alphabetic characters and spaces" in output


def test_main_keyboard_interrupt():
    """Test main function handles KeyboardInterrupt gracefully."""
    with patch('builtins.input', side_effect=KeyboardInterrupt):
        with patch('sys.stdout', new_callable=io.StringIO) as mock_stdout:
            main()
            output = mock_stdout.getvalue()
            
            assert "Program interrupted by user" in output


def test_main_similar_words():
    """Test main function interpretation for similar words."""
    user_inputs = ["cat", "bat"]
    
    with patch('builtins.input', side_effect=user_inputs):
        with patch('sys.stdout', new_callable=io.StringIO) as mock_stdout:
            main()
            output = mock_stdout.getvalue()
            
            assert "Levenshtein distance: 1" in output
            assert "differ by exactly one character" in output


def test_main_quite_similar_words():
    """Test main function interpretation for quite similar words."""
    user_inputs = ["hellx", "hallo"]
    
    with patch('builtins.input', side_effect=user_inputs):
        with patch('sys.stdout', new_callable=io.StringIO) as mock_stdout:
            main()
            output = mock_stdout.getvalue()
            
            assert "quite similar" in output


def test_main_quite_different_words():
    """Test main function interpretation for quite different words."""
    user_inputs = ["python", "java"]
    
    with patch('builtins.input', side_effect=user_inputs):
        with patch('sys.stdout', new_callable=io.StringIO) as mock_stdout:
            main()
            output = mock_stdout.getvalue()
            
            assert "quite different" in output