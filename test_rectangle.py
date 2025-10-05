import pytest
from unittest.mock import patch
import sys
import io
from rectangle import calculate_rectangle_properties, main


# Tests for calculate_rectangle_properties function

def test_calculate_rectangle_properties_positive_integers():
    """Test calculation with positive integer values."""
    area, perimeter = calculate_rectangle_properties(5, 3)
    assert area == 15
    assert perimeter == 16


def test_calculate_rectangle_properties_positive_floats():
    """Test calculation with positive float values."""
    area, perimeter = calculate_rectangle_properties(2.5, 4.0)
    assert area == 10.0
    assert perimeter == 13.0


def test_calculate_rectangle_properties_square():
    """Test calculation with a square (equal width and height)."""
    area, perimeter = calculate_rectangle_properties(4, 4)
    assert area == 16
    assert perimeter == 16


def test_calculate_rectangle_properties_decimal_precision():
    """Test calculation with decimal values requiring precision."""
    area, perimeter = calculate_rectangle_properties(1.5, 2.3)
    assert area == pytest.approx(3.45, abs=1e-2)
    assert perimeter == pytest.approx(7.6, abs=1e-1)


def test_calculate_rectangle_properties_large_numbers():
    """Test calculation with large numbers."""
    area, perimeter = calculate_rectangle_properties(1000, 2000)
    assert area == 2000000
    assert perimeter == 6000


def test_calculate_rectangle_properties_small_numbers():
    """Test calculation with very small numbers."""
    area, perimeter = calculate_rectangle_properties(0.1, 0.2)
    assert area == pytest.approx(0.02, abs=1e-3)
    assert perimeter == pytest.approx(0.6, abs=1e-1)


# Tests for main function

@patch('builtins.input')
@patch('sys.stdout', new_callable=io.StringIO)
def test_main_valid_input(mock_stdout, mock_input):
    """Test main function with valid input."""
    mock_input.side_effect = ['5', '3']
    
    main()
    
    output = mock_stdout.getvalue()
    assert "Rectangle Calculator" in output
    assert "Width: 5.0" in output
    assert "Height: 3.0" in output
    assert "Area: 15.0" in output
    assert "Perimeter: 16.0" in output


@patch('builtins.input')
@patch('sys.stdout', new_callable=io.StringIO)
def test_main_negative_width(mock_stdout, mock_input):
    """Test main function with negative width."""
    mock_input.side_effect = ['-5', '3']
    
    main()
    
    output = mock_stdout.getvalue()
    assert "Error: Width and height must be positive numbers." in output


@patch('builtins.input')
@patch('sys.stdout', new_callable=io.StringIO)
def test_main_negative_height(mock_stdout, mock_input):
    """Test main function with negative height."""
    mock_input.side_effect = ['5', '-3']
    
    main()
    
    output = mock_stdout.getvalue()
    assert "Error: Width and height must be positive numbers." in output


@patch('builtins.input')
@patch('sys.stdout', new_callable=io.StringIO)
def test_main_zero_width(mock_stdout, mock_input):
    """Test main function with zero width."""
    mock_input.side_effect = ['0', '3']
    
    main()
    
    output = mock_stdout.getvalue()
    assert "Error: Width and height must be positive numbers." in output


@patch('builtins.input')
@patch('sys.stdout', new_callable=io.StringIO)
def test_main_zero_height(mock_stdout, mock_input):
    """Test main function with zero height."""
    mock_input.side_effect = ['5', '0']
    
    main()
    
    output = mock_stdout.getvalue()
    assert "Error: Width and height must be positive numbers." in output


@patch('builtins.input')
@patch('sys.stdout', new_callable=io.StringIO)
def test_main_invalid_width_input(mock_stdout, mock_input):
    """Test main function with invalid (non-numeric) width input."""
    mock_input.side_effect = ['abc', '3']
    
    main()
    
    output = mock_stdout.getvalue()
    assert "Error: Please enter valid numbers for width and height." in output


@patch('builtins.input')
@patch('sys.stdout', new_callable=io.StringIO)
def test_main_invalid_height_input(mock_stdout, mock_input):
    """Test main function with invalid (non-numeric) height input."""
    mock_input.side_effect = ['5', 'xyz']
    
    main()
    
    output = mock_stdout.getvalue()
    assert "Error: Please enter valid numbers for width and height." in output


@patch('builtins.input')
@patch('sys.stdout', new_callable=io.StringIO)
def test_main_decimal_input(mock_stdout, mock_input):
    """Test main function with decimal input."""
    mock_input.side_effect = ['2.5', '4.2']
    
    main()
    
    output = mock_stdout.getvalue()
    assert "Width: 2.5" in output
    assert "Height: 4.2" in output
    assert "Area: 10.5" in output
    assert "Perimeter: 13.4" in output


@patch('builtins.input')
@patch('sys.stdout', new_callable=io.StringIO)
def test_main_keyboard_interrupt(mock_stdout, mock_input):
    """Test main function with keyboard interrupt."""
    mock_input.side_effect = KeyboardInterrupt()
    
    main()
    
    output = mock_stdout.getvalue()
    assert "Program interrupted by user." in output


# Edge case tests

def test_very_large_numbers():
    """Test with very large numbers."""
    area, perimeter = calculate_rectangle_properties(1e6, 1e6)
    assert area == 1e12
    assert perimeter == 4e6


def test_very_small_numbers():
    """Test with very small numbers."""
    area, perimeter = calculate_rectangle_properties(1e-6, 1e-6)
    assert area == pytest.approx(1e-12, abs=1e-15)
    assert perimeter == pytest.approx(4e-6, abs=1e-10)


def test_return_type():
    """Test that function returns correct data types."""
    result = calculate_rectangle_properties(5, 3)
    assert isinstance(result, tuple)
    assert len(result) == 2
    area, perimeter = result
    assert isinstance(area, (int, float))
    assert isinstance(perimeter, (int, float))