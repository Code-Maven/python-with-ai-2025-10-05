import math
from typing import List, Tuple, Dict


def calculate_pie_chart_segments(data: List[Tuple[str, float]]) -> List[Dict]:
    """
    Calculate pie chart segments with angles and positions.
    
    Args:
        data: List of tuples (label, value) for pie chart segments
    
    Returns:
        List of dictionaries with segment information
    """
    total = sum(value for _, value in data)
    if total == 0:
        raise ValueError("Total value cannot be zero")
    
    segments = []
    current_angle = 0
    
    for label, value in data:
        # Calculate percentage and angle
        percentage = (value / total) * 100
        angle = (value / total) * 360
        
        # Calculate end angle
        end_angle = current_angle + angle
        
        # Store segment information
        segments.append({
            'label': label,
            'value': value,
            'percentage': percentage,
            'start_angle': current_angle,
            'end_angle': end_angle,
            'angle': angle
        })
        
        current_angle = end_angle
    
    return segments


def generate_svg_path(center_x: float, center_y: float, radius: float, 
                     start_angle: float, end_angle: float) -> str:
    """
    Generate SVG path for a pie chart segment.
    
    Args:
        center_x, center_y: Center coordinates of the pie chart
        radius: Radius of the pie chart
        start_angle, end_angle: Angles in degrees
    
    Returns:
        SVG path string
    """
    # Convert degrees to radians
    start_rad = math.radians(start_angle - 90)  # -90 to start from top
    end_rad = math.radians(end_angle - 90)
    
    # Calculate start and end points
    start_x = center_x + radius * math.cos(start_rad)
    start_y = center_y + radius * math.sin(start_rad)
    end_x = center_x + radius * math.cos(end_rad)
    end_y = center_y + radius * math.sin(end_rad)
    
    # Determine if arc should be large (> 180 degrees)
    large_arc = 1 if (end_angle - start_angle) > 180 else 0
    
    # Create SVG path
    path = f"M {center_x} {center_y} "  # Move to center
    path += f"L {start_x} {start_y} "   # Line to start point
    path += f"A {radius} {radius} 0 {large_arc} 1 {end_x} {end_y} "  # Arc
    path += "Z"  # Close path
    
    return path


def create_pie_chart_svg(data: List[Tuple[str, float]], 
                        width: int = 400, 
                        height: int = 400,
                        title: str = "Pie Chart") -> str:
    """
    Create an SVG pie chart.
    
    Args:
        data: List of tuples (label, value) for pie chart segments
        width, height: SVG dimensions
        title: Chart title
    
    Returns:
        Complete SVG string
    """
    # Calculate segments
    segments = calculate_pie_chart_segments(data)
    
    # Chart settings
    center_x = width // 2
    center_y = height // 2
    radius = min(width, height) // 3
    
    # Color palette
    colors = [
        "#FF6B6B", "#4ECDC4", "#45B7D1", "#96CEB4", "#FFEAA7",
        "#DDA0DD", "#98D8C8", "#F7DC6F", "#BB8FCE", "#85C1E9"
    ]
    
    # Start SVG
    svg = f'''<?xml version="1.0" encoding="UTF-8"?>
<svg width="{width}" height="{height}" viewBox="0 0 {width} {height}" 
     xmlns="http://www.w3.org/2000/svg">
  
  <!-- Title -->
  <text x="{center_x}" y="30" text-anchor="middle" 
        font-family="Arial, sans-serif" font-size="18" font-weight="bold" 
        fill="#333">{title}</text>
  
  <!-- Pie chart segments -->
'''
    
    # Add each segment
    for i, segment in enumerate(segments):
        color = colors[i % len(colors)]
        path = generate_svg_path(center_x, center_y, radius,
                               segment['start_angle'], segment['end_angle'])
        
        svg += f'''  <path d="{path}" fill="{color}" stroke="white" stroke-width="2">
    <title>{segment['label']}: {segment['value']} ({segment['percentage']:.1f}%)</title>
  </path>
'''
    
    # Add legend
    legend_x = center_x + radius + 40
    legend_y = center_y - (len(segments) * 12)
    
    svg += "\n  <!-- Legend -->\n"
    for i, segment in enumerate(segments):
        color = colors[i % len(colors)]
        y_pos = legend_y + i * 25
        
        # Legend color box
        svg += f'''  <rect x="{legend_x}" y="{y_pos}" width="15" height="15" 
        fill="{color}" stroke="#333" stroke-width="1"/>
'''
        
        # Legend text
        svg += f'''  <text x="{legend_x + 20}" y="{y_pos + 12}" 
        font-family="Arial, sans-serif" font-size="12" fill="#333">
    {segment['label']}: {segment['value']} ({segment['percentage']:.1f}%)
  </text>
'''
    
    # Close SVG
    svg += "\n</svg>"
    
    return svg


def main() -> None:
    try:
        print("SVG Pie Chart Generator")
        print("-" * 22)
        print("This program creates an SVG pie chart from your data.\n")
        
        # Get chart title
        DEFAULT_TITLE = "Pie Chart"

        title = input(f"Enter chart title (or press Enter for '{DEFAULT_TITLE}'): ").strip()
        if not title:
            title = DEFAULT_TITLE
        
        # Get data
        data = []
        print("\nEnter your data (label and value). Press Enter with empty label to finish:")
        
        while True:
            label = input("Label: ").strip()
            if not label:
                break
                
            try:
                value = float(input("Value: ").strip())
                if value < 0:
                    print("Warning: Negative values will be treated as zero.")
                    value = 0
                data.append((label, value))
            except ValueError:
                print("Error: Please enter a valid number for value.")
                continue
        
        # Validate data
        if not data:
            print("Error: No data provided.")
            return
            
        if sum(value for _, value in data) == 0:
            print("Error: All values are zero.")
            return
        
        # Generate SVG
        svg_content = create_pie_chart_svg(data, title=title)
        
        # Save to file
        filename = input(f"\nEnter filename (or press Enter for 'pie_chart.svg'): ").strip()
        if not filename:
            filename = "pie_chart.svg"
        elif not filename.endswith('.svg'):
            filename += '.svg'
        
        with open(filename, 'w', encoding='utf-8') as f:
            f.write(svg_content)
        
        print(f"\nPie chart saved as '{filename}'")
        print(f"Total segments: {len(data)}")
        print(f"Total value: {sum(value for _, value in data)}")
        
        # Show data summary
        print("\nData summary:")
        for label, value in data:
            percentage = (value / sum(v for _, v in data)) * 100
            print(f"  {label}: {value} ({percentage:.1f}%)")
        
    except KeyboardInterrupt:
        print("\nProgram interrupted by user.")
    except Exception as e:
        print(f"Error: {e}")


if __name__ == "__main__":
    main()