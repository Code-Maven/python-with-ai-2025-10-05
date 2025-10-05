import sys
import argparse


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


def parse_arguments():
    """
    Parse command line arguments for rectangle dimensions.
    
    Returns:
        argparse.Namespace: Parsed arguments containing width and height
    """
    parser = argparse.ArgumentParser(
        description='Calculate the area and perimeter of a rectangle',
        epilog='Example: python rectangle_argv.py 5.0 3.0'
    )
    
    parser.add_argument(
        'width',
        type=float,
        help='Width of the rectangle (must be positive)'
    )
    
    parser.add_argument(
        'height', 
        type=float,
        help='Height of the rectangle (must be positive)'
    )
    
    parser.add_argument(
        '-v', '--verbose',
        action='store_true',
        help='Enable verbose output'
    )
    
    return parser.parse_args()


def main():
    try:
        # Parse command line arguments
        args = parse_arguments()
        
        width = args.width
        height = args.height
        verbose = args.verbose
        
        # Validate input
        if width <= 0:
            print(f"Error: Width must be positive, got {width}", file=sys.stderr)
            sys.exit(1)
        
        if height <= 0:
            print(f"Error: Height must be positive, got {height}", file=sys.stderr)
            sys.exit(1)
        
        # Calculate properties
        area, perimeter = calculate_rectangle_properties(width, height)
        
        # Print results
        if verbose:
            print("Rectangle Calculator")
            print("-" * 20)
            print(f"Width: {width}")
            print(f"Height: {height}")
            print(f"Area: {area}")
            print(f"Perimeter: {perimeter}")
        else:
            print(f"Area: {area}")
            print(f"Perimeter: {perimeter}")
            
    except ValueError as e:
        print(f"Error: Invalid number format - {e}", file=sys.stderr)
        sys.exit(1)
    except SystemExit:
        # Re-raise SystemExit to allow argparse to handle help/error messages
        raise
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()