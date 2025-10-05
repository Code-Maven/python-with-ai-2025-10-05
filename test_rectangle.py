import unittest
from unittest.mock import patch, call
import sys
import io
from rectangle import calculate_rectangle_properties, main


class TestRectangleProperties(unittest.TestCase):
    """Test cases for the calculate_rectangle_properties function."""
    
    def test_calculate_rectangle_properties_positive_integers(self):
        """Test calculation with positive integer values."""
        area, perimeter = calculate_rectangle_properties(5, 3)
        self.assertEqual(area, 15)
        self.assertEqual(perimeter, 16)
    
    def test_calculate_rectangle_properties_positive_floats(self):
        """Test calculation with positive float values."""
        area, perimeter = calculate_rectangle_properties(2.5, 4.0)
        self.assertEqual(area, 10.0)
        self.assertEqual(perimeter, 13.0)
    
    def test_calculate_rectangle_properties_square(self):
        """Test calculation with a square (equal width and height)."""
        area, perimeter = calculate_rectangle_properties(4, 4)
        self.assertEqual(area, 16)
        self.assertEqual(perimeter, 16)
    
    def test_calculate_rectangle_properties_decimal_precision(self):
        """Test calculation with decimal values requiring precision."""
        area, perimeter = calculate_rectangle_properties(1.5, 2.3)
        self.assertAlmostEqual(area, 3.45, places=2)
        self.assertAlmostEqual(perimeter, 7.6, places=1)
    
    def test_calculate_rectangle_properties_large_numbers(self):
        """Test calculation with large numbers."""
        area, perimeter = calculate_rectangle_properties(1000, 2000)
        self.assertEqual(area, 2000000)
        self.assertEqual(perimeter, 6000)
    
    def test_calculate_rectangle_properties_small_numbers(self):
        """Test calculation with very small numbers."""
        area, perimeter = calculate_rectangle_properties(0.1, 0.2)
        self.assertAlmostEqual(area, 0.02, places=3)
        self.assertAlmostEqual(perimeter, 0.6, places=1)


class TestMainFunction(unittest.TestCase):
    """Test cases for the main function."""
    
    @patch('builtins.input')
    @patch('sys.stdout', new_callable=io.StringIO)
    def test_main_valid_input(self, mock_stdout, mock_input):
        """Test main function with valid input."""
        mock_input.side_effect = ['5', '3']
        
        main()
        
        output = mock_stdout.getvalue()
        self.assertIn("Rectangle Calculator", output)
        self.assertIn("Width: 5.0", output)
        self.assertIn("Height: 3.0", output)
        self.assertIn("Area: 15.0", output)
        self.assertIn("Perimeter: 16.0", output)
    
    @patch('builtins.input')
    @patch('sys.stdout', new_callable=io.StringIO)
    def test_main_negative_width(self, mock_stdout, mock_input):
        """Test main function with negative width."""
        mock_input.side_effect = ['-5', '3']
        
        main()
        
        output = mock_stdout.getvalue()
        self.assertIn("Error: Width and height must be positive numbers.", output)
    
    @patch('builtins.input')
    @patch('sys.stdout', new_callable=io.StringIO)
    def test_main_negative_height(self, mock_stdout, mock_input):
        """Test main function with negative height."""
        mock_input.side_effect = ['5', '-3']
        
        main()
        
        output = mock_stdout.getvalue()
        self.assertIn("Error: Width and height must be positive numbers.", output)
    
    @patch('builtins.input')
    @patch('sys.stdout', new_callable=io.StringIO)
    def test_main_zero_width(self, mock_stdout, mock_input):
        """Test main function with zero width."""
        mock_input.side_effect = ['0', '3']
        
        main()
        
        output = mock_stdout.getvalue()
        self.assertIn("Error: Width and height must be positive numbers.", output)
    
    @patch('builtins.input')
    @patch('sys.stdout', new_callable=io.StringIO)
    def test_main_zero_height(self, mock_stdout, mock_input):
        """Test main function with zero height."""
        mock_input.side_effect = ['5', '0']
        
        main()
        
        output = mock_stdout.getvalue()
        self.assertIn("Error: Width and height must be positive numbers.", output)
    
    @patch('builtins.input')
    @patch('sys.stdout', new_callable=io.StringIO)
    def test_main_invalid_width_input(self, mock_stdout, mock_input):
        """Test main function with invalid (non-numeric) width input."""
        mock_input.side_effect = ['abc', '3']
        
        main()
        
        output = mock_stdout.getvalue()
        self.assertIn("Error: Please enter valid numbers for width and height.", output)
    
    @patch('builtins.input')
    @patch('sys.stdout', new_callable=io.StringIO)
    def test_main_invalid_height_input(self, mock_stdout, mock_input):
        """Test main function with invalid (non-numeric) height input."""
        mock_input.side_effect = ['5', 'xyz']
        
        main()
        
        output = mock_stdout.getvalue()
        self.assertIn("Error: Please enter valid numbers for width and height.", output)
    
    @patch('builtins.input')
    @patch('sys.stdout', new_callable=io.StringIO)
    def test_main_decimal_input(self, mock_stdout, mock_input):
        """Test main function with decimal input."""
        mock_input.side_effect = ['2.5', '4.2']
        
        main()
        
        output = mock_stdout.getvalue()
        self.assertIn("Width: 2.5", output)
        self.assertIn("Height: 4.2", output)
        self.assertIn("Area: 10.5", output)
        self.assertIn("Perimeter: 13.4", output)
    
    @patch('builtins.input')
    @patch('sys.stdout', new_callable=io.StringIO)
    def test_main_keyboard_interrupt(self, mock_stdout, mock_input):
        """Test main function with keyboard interrupt."""
        mock_input.side_effect = KeyboardInterrupt()
        
        main()
        
        output = mock_stdout.getvalue()
        self.assertIn("Program interrupted by user.", output)


class TestEdgeCases(unittest.TestCase):
    """Test cases for edge cases and boundary conditions."""
    
    def test_very_large_numbers(self):
        """Test with very large numbers."""
        area, perimeter = calculate_rectangle_properties(1e6, 1e6)
        self.assertEqual(area, 1e12)
        self.assertEqual(perimeter, 4e6)
    
    def test_very_small_numbers(self):
        """Test with very small numbers."""
        area, perimeter = calculate_rectangle_properties(1e-6, 1e-6)
        self.assertAlmostEqual(area, 1e-12, places=15)
        self.assertAlmostEqual(perimeter, 4e-6, places=10)
    
    def test_return_type(self):
        """Test that function returns correct data types."""
        result = calculate_rectangle_properties(5, 3)
        self.assertIsInstance(result, tuple)
        self.assertEqual(len(result), 2)
        area, perimeter = result
        self.assertIsInstance(area, (int, float))
        self.assertIsInstance(perimeter, (int, float))


if __name__ == '__main__':
    # Run the tests
    unittest.main(verbosity=2)