import pytest
import math
from unittest.mock import patch, mock_open
import sys
import io
from pie_chart import calculate_pie_chart_segments, generate_svg_path, create_pie_chart_svg, main


# Tests for calculate_pie_chart_segments function

def test_calculate_pie_chart_segments_basic():
    """Test basic segment calculation."""
    data = [("A", 50), ("B", 50)]
    segments = calculate_pie_chart_segments(data)
    
    assert len(segments) == 2
    assert segments[0]['label'] == "A"
    assert segments[0]['value'] == 50
    assert segments[0]['percentage'] == 50.0
    assert segments[0]['start_angle'] == 0
    assert segments[0]['end_angle'] == 180
    assert segments[0]['angle'] == 180
    
    assert segments[1]['label'] == "B"
    assert segments[1]['value'] == 50
    assert segments[1]['percentage'] == 50.0
    assert segments[1]['start_angle'] == 180
    assert segments[1]['end_angle'] == 360
    assert segments[1]['angle'] == 180


def test_calculate_pie_chart_segments_unequal():
    """Test segment calculation with unequal values."""
    data = [("Small", 10), ("Large", 90)]
    segments = calculate_pie_chart_segments(data)
    
    assert len(segments) == 2
    assert segments[0]['percentage'] == 10.0
    assert segments[0]['angle'] == 36.0  # 10% of 360
    assert segments[1]['percentage'] == 90.0
    assert segments[1]['angle'] == 324.0  # 90% of 360


def test_calculate_pie_chart_segments_multiple():
    """Test segment calculation with multiple segments."""
    data = [("A", 25), ("B", 25), ("C", 25), ("D", 25)]
    segments = calculate_pie_chart_segments(data)
    
    assert len(segments) == 4
    for i, segment in enumerate(segments):
        assert segment['percentage'] == 25.0
        assert segment['angle'] == 90.0
        assert segment['start_angle'] == i * 90
        assert segment['end_angle'] == (i + 1) * 90


def test_calculate_pie_chart_segments_zero_total():
    """Test that zero total raises ValueError."""
    data = [("A", 0), ("B", 0)]
    with pytest.raises(ValueError, match="Total value cannot be zero"):
        calculate_pie_chart_segments(data)


def test_calculate_pie_chart_segments_empty():
    """Test that empty data raises ValueError."""
    data = []
    with pytest.raises(ValueError, match="Total value cannot be zero"):
        calculate_pie_chart_segments(data)


def test_calculate_pie_chart_segments_single():
    """Test segment calculation with single segment."""
    data = [("Only", 100)]
    segments = calculate_pie_chart_segments(data)
    
    assert len(segments) == 1
    assert segments[0]['percentage'] == 100.0
    assert segments[0]['angle'] == 360.0
    assert segments[0]['start_angle'] == 0
    assert segments[0]['end_angle'] == 360


def test_calculate_pie_chart_segments_decimal_values():
    """Test segment calculation with decimal values."""
    data = [("A", 33.33), ("B", 66.67)]
    segments = calculate_pie_chart_segments(data)
    
    assert len(segments) == 2
    assert abs(segments[0]['percentage'] - 33.33) < 0.01
    assert abs(segments[1]['percentage'] - 66.67) < 0.01


# Tests for generate_svg_path function

def test_generate_svg_path_quarter_circle():
    """Test SVG path generation for quarter circle."""
    path = generate_svg_path(100, 100, 50, 0, 90)
    
    # Should start from center, go to top, arc to right
    assert "M 100 100" in path  # Move to center
    assert "L 100.0 50.0" in path   # Line to top (radius 50 up from center)
    assert "A 50 50" in path    # Arc with radius 50
    assert "Z" in path          # Close path


def test_generate_svg_path_half_circle():
    """Test SVG path generation for half circle."""
    path = generate_svg_path(100, 100, 50, 0, 180)
    
    assert "M 100 100" in path
    assert "A 50 50 0 0 1" in path  # Small arc (< 180 degrees logic should be 0)
    assert "Z" in path


def test_generate_svg_path_large_arc():
    """Test SVG path generation for large arc (> 180 degrees)."""
    path = generate_svg_path(100, 100, 50, 0, 270)
    
    assert "M 100 100" in path
    assert "A 50 50 0 1 1" in path  # Large arc flag should be 1
    assert "Z" in path


def test_generate_svg_path_full_circle():
    """Test SVG path generation for full circle."""
    path = generate_svg_path(100, 100, 50, 0, 360)
    
    assert "M 100 100" in path
    assert "A 50 50 0 1 1" in path  # Large arc
    assert "Z" in path


def test_generate_svg_path_different_center():
    """Test SVG path generation with different center coordinates."""
    path = generate_svg_path(200, 150, 75, 0, 90)
    
    assert "M 200 150" in path  # Should use specified center
    assert "A 75 75" in path    # Should use specified radius


# Tests for create_pie_chart_svg function

def test_create_pie_chart_svg_basic():
    """Test basic SVG creation."""
    data = [("A", 50), ("B", 50)]
    svg = create_pie_chart_svg(data)
    
    # Check SVG structure
    assert '<?xml version="1.0" encoding="UTF-8"?>' in svg
    assert '<svg' in svg
    assert '</svg>' in svg
    assert 'xmlns="http://www.w3.org/2000/svg"' in svg
    
    # Check for pie chart elements
    assert '<path' in svg
    assert 'fill=' in svg
    assert 'stroke=' in svg


def test_create_pie_chart_svg_with_title():
    """Test SVG creation with custom title."""
    data = [("Test", 100)]
    svg = create_pie_chart_svg(data, title="Custom Title")
    
    assert "Custom Title" in svg


def test_create_pie_chart_svg_with_dimensions():
    """Test SVG creation with custom dimensions."""
    data = [("Test", 100)]
    svg = create_pie_chart_svg(data, width=800, height=600)
    
    assert 'width="800"' in svg
    assert 'height="600"' in svg
    assert 'viewBox="0 0 800 600"' in svg


def test_create_pie_chart_svg_multiple_segments():
    """Test SVG creation with multiple segments."""
    data = [("A", 25), ("B", 25), ("C", 25), ("D", 25)]
    svg = create_pie_chart_svg(data)
    
    # Should have 4 path elements
    path_count = svg.count('<path')
    assert path_count == 4
    
    # Should have legend entries
    assert "A: 25" in svg
    assert "B: 25" in svg
    assert "C: 25" in svg
    assert "D: 25" in svg


def test_create_pie_chart_svg_legend():
    """Test that legend is created correctly."""
    data = [("Alpha", 60), ("Beta", 40)]
    svg = create_pie_chart_svg(data)
    
    # Check for legend elements
    assert '<rect' in svg  # Legend color boxes
    assert 'Alpha: 60 (60.0%)' in svg
    assert 'Beta: 40 (40.0%)' in svg


def test_create_pie_chart_svg_tooltips():
    """Test that tooltips are included."""
    data = [("Item", 100)]
    svg = create_pie_chart_svg(data)
    
    assert '<title>' in svg
    assert 'Item: 100' in svg


# Tests for main function

def test_main_valid_input():
    """Test main function with valid input."""
    user_inputs = ["Test Chart", "Item A", "50", "Item B", "30", "Item C", "20", "", "test_output.svg"]
    
    with patch('builtins.input', side_effect=user_inputs):
        with patch('builtins.open', mock_open()) as mock_file:
            with patch('sys.stdout', new_callable=io.StringIO) as mock_stdout:
                main()
                
                # Check that file was opened for writing
                mock_file.assert_called_once_with("test_output.svg", 'w', encoding='utf-8')
                
                # Check output messages
                output = mock_stdout.getvalue()
                assert "SVG Pie Chart Generator" in output
                assert "Total segments: 3" in output
                assert "Total value: 100" in output


def test_main_default_title():
    """Test main function with default title."""
    user_inputs = ["", "Item", "100", "", "chart.svg"]
    
    with patch('builtins.input', side_effect=user_inputs):
        with patch('builtins.open', mock_open()) as mock_file:
            with patch('sys.stdout', new_callable=io.StringIO):
                main()
                
                # Should use default title
                handle = mock_file.return_value.__enter__.return_value
                written_content = ''.join(call.args[0] for call in handle.write.call_args_list)
                assert "Pie Chart" in written_content


def test_main_default_filename():
    """Test main function with default filename."""
    user_inputs = ["Test", "Item", "50", "", ""]
    
    with patch('builtins.input', side_effect=user_inputs):
        with patch('builtins.open', mock_open()) as mock_file:
            with patch('sys.stdout', new_callable=io.StringIO):
                main()
                
                # Should use default filename
                mock_file.assert_called_once_with("pie_chart.svg", 'w', encoding='utf-8')


def test_main_filename_auto_extension():
    """Test main function adds .svg extension if missing."""
    user_inputs = ["Test", "Item", "50", "", "myfile"]
    
    with patch('builtins.input', side_effect=user_inputs):
        with patch('builtins.open', mock_open()) as mock_file:
            with patch('sys.stdout', new_callable=io.StringIO):
                main()
                
                # Should add .svg extension
                mock_file.assert_called_once_with("myfile.svg", 'w', encoding='utf-8')


def test_main_negative_values():
    """Test main function handles negative values."""
    user_inputs = ["Test", "Item", "-10", "", "test.svg"]
    
    with patch('builtins.input', side_effect=user_inputs):
        with patch('builtins.open', mock_open()) as mock_file:
            with patch('sys.stdout', new_callable=io.StringIO) as mock_stdout:
                main()
                
                output = mock_stdout.getvalue()
                assert "Warning: Negative values will be treated as zero" in output


def test_main_invalid_value():
    """Test main function handles invalid numeric input."""
    user_inputs = ["Test", "Item", "abc", "Item", "50", "", "test.svg"]
    
    with patch('builtins.input', side_effect=user_inputs):
        with patch('builtins.open', mock_open()) as mock_file:
            with patch('sys.stdout', new_callable=io.StringIO) as mock_stdout:
                main()
                
                output = mock_stdout.getvalue()
                assert "Error: Please enter a valid number" in output


def test_main_no_data():
    """Test main function handles no data input."""
    user_inputs = ["Test", "", "test.svg"]
    
    with patch('builtins.input', side_effect=user_inputs):
        with patch('sys.stdout', new_callable=io.StringIO) as mock_stdout:
            main()
            
            output = mock_stdout.getvalue()
            assert "Error: No data provided" in output


def test_main_all_zero_values():
    """Test main function handles all zero values."""
    user_inputs = ["Test", "Item A", "0", "Item B", "0", "", "test.svg"]
    
    with patch('builtins.input', side_effect=user_inputs):
        with patch('sys.stdout', new_callable=io.StringIO) as mock_stdout:
            main()
            
            output = mock_stdout.getvalue()
            assert "Error: All values are zero" in output


def test_main_keyboard_interrupt():
    """Test main function handles KeyboardInterrupt gracefully."""
    with patch('builtins.input', side_effect=KeyboardInterrupt):
        with patch('sys.stdout', new_callable=io.StringIO) as mock_stdout:
            main()
            
            output = mock_stdout.getvalue()
            assert "Program interrupted by user" in output


def test_main_data_summary():
    """Test main function shows data summary."""
    user_inputs = ["Summary Test", "A", "30", "B", "70", "", "summary.svg"]
    
    with patch('builtins.input', side_effect=user_inputs):
        with patch('builtins.open', mock_open()):
            with patch('sys.stdout', new_callable=io.StringIO) as mock_stdout:
                main()
                
                output = mock_stdout.getvalue()
                assert "Data summary:" in output
                assert "A: 30.0 (30.0%)" in output
                assert "B: 70.0 (70.0%)" in output


# Edge case and integration tests

def test_pie_chart_segments_precision():
    """Test that segment calculations maintain precision."""
    data = [("A", 1), ("B", 2), ("C", 3)]  # Should be 16.67%, 33.33%, 50%
    segments = calculate_pie_chart_segments(data)
    
    total_percentage = sum(seg['percentage'] for seg in segments)
    assert abs(total_percentage - 100.0) < 0.001  # Should sum to 100%
    
    total_angle = sum(seg['angle'] for seg in segments)
    assert abs(total_angle - 360.0) < 0.001  # Should sum to 360°


def test_svg_path_coordinates_reasonable():
    """Test that generated coordinates are reasonable."""
    path = generate_svg_path(200, 200, 100, 0, 90)
    
    # Extract coordinates and check they're reasonable
    assert "200" in path  # Center coordinates should appear
    assert "100" in path or "300" in path  # Should have points at radius distance


def test_full_workflow():
    """Test complete workflow from data to SVG."""
    data = [("Test A", 40), ("Test B", 60)]
    svg = create_pie_chart_svg(data, title="Integration Test")
    
    # Verify complete SVG structure
    assert svg.startswith('<?xml version="1.0" encoding="UTF-8"?>')
    assert svg.endswith('</svg>')
    assert "Integration Test" in svg
    assert svg.count('<path') == 2  # Two segments
    assert 'fill="#FF6B6B"' in svg  # First color
    assert 'fill="#4ECDC4"' in svg  # Second color