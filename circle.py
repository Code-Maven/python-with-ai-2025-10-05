import math


def calculate_circle_properties(radius: float) -> tuple[float, float]:
    """
    Calculate the area and circumference of a circle.
    
    Args:
        radius (float): Radius of the circle
    
    Returns:
        tuple: (area, circumference)
    """
    area = math.pi * radius * radius
    circumference = 2 * math.pi * radius
    return area, circumference


def main() -> None:
    # area, circumference = calculate_circle_properties("abc")

    try:
        # Get input from user
        print("Circle Calculator")
        print("-" * 16)
        
        radius = float(input("Enter the radius of the circle: "))
        
        # Validate input
        if radius <= 0:
            print("Error: Radius must be a positive number.")
            return
        
        # Calculate properties
        area, circumference = calculate_circle_properties(radius)
        
        # Print results
        print(f"\nCircle Properties:")
        print(f"Radius: {radius}")
        print(f"Area: {area:.2f}")
        print(f"Circumference: {circumference:.2f}")
        
    except ValueError:
        print("Error: Please enter a valid number for radius.")
    except KeyboardInterrupt:
        print("\nProgram interrupted by user.")


if __name__ == "__main__":
    main()