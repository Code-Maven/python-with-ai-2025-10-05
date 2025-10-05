def calculate_rectangle_properties(width, height):
    """
    Calculate the area and perimeter of a rectangle.
    
    Args:
        width (float): Width of the rectangle
        height (float): Height of the rectangle
    
    Returns:
        tuple: (area, perimeter)
    """
    area = width * height
    perimeter = 2 * (width + height)
    return area, perimeter


def main():
    try:
        # Get input from user
        print("Rectangle Calculator")
        print("-" * 20)
        
        width = float(input("Enter the width of the rectangle: "))
        height = float(input("Enter the height of the rectangle: "))
        
        # Validate input
        if width <= 0 or height <= 0:
            print("Error: Width and height must be positive numbers.")
            return
        
        # Calculate properties
        area, perimeter = calculate_rectangle_properties(width, height)
        
        # Print results
        print(f"\nRectangle Properties:")
        print(f"Width: {width}")
        print(f"Height: {height}")
        print(f"Area: {area}")
        print(f"Perimeter: {perimeter}")
        
    except ValueError:
        print("Error: Please enter valid numbers for width and height.")
    except KeyboardInterrupt:
        print("\nProgram interrupted by user.")


if __name__ == "__main__":
    main()