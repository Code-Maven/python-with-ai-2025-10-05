import flet as ft


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


def main(page: ft.Page):
    # Configure the page
    page.title = "Rectangle Calculator - Modern Flet UI"
    page.theme_mode = ft.ThemeMode.LIGHT
    page.padding = 20
    
    # Create input fields
    width_field = ft.TextField(
        label="Width",
        hint_text="Enter rectangle width",
        width=200,
    )
    
    height_field = ft.TextField(
        label="Height", 
        hint_text="Enter rectangle height",
        width=200,
    )
    
    # Result displays
    area_result = ft.Text("Area: --", size=18, weight=ft.FontWeight.BOLD)
    perimeter_result = ft.Text("Perimeter: --", size=18, weight=ft.FontWeight.BOLD)
    
    # Status text
    status_text = ft.Text(
        "Enter width and height, then calculate",
        size=14,
        color=ft.Colors.GREY_600,
    )
    
    def calculate_clicked(e):
        try:
            # Get and validate input
            width_str = width_field.value.strip() if width_field.value else ""
            height_str = height_field.value.strip() if height_field.value else ""
            
            if not width_str:
                show_error("Please enter a width value")
                return
            
            if not height_str:
                show_error("Please enter a height value")
                return
            
            # Convert to float
            try:
                width = float(width_str)
                height = float(height_str)
            except ValueError:
                show_error("Please enter valid numbers")
                return
            
            # Validate positive numbers
            if width <= 0:
                show_error("Width must be positive")
                return
            
            if height <= 0:
                show_error("Height must be positive")
                return
            
            # Calculate results
            area, perimeter = calculate_rectangle_properties(width, height)
            
            # Update results
            area_result.value = f"Area: {area:.2f}"
            area_result.color = ft.Colors.BLUE
            perimeter_result.value = f"Perimeter: {perimeter:.2f}"
            perimeter_result.color = ft.Colors.GREEN
            
            # Update status
            status_text.value = f"✅ Calculated for rectangle {width} × {height}"
            status_text.color = ft.Colors.GREEN_600
            
            page.update()
            
        except Exception as ex:
            show_error(f"An error occurred: {str(ex)}")
    
    def clear_clicked(e):
        width_field.value = ""
        height_field.value = ""
        area_result.value = "Area: --"
        area_result.color = ft.Colors.BLACK
        perimeter_result.value = "Perimeter: --"
        perimeter_result.color = ft.Colors.BLACK
        status_text.value = "Enter width and height, then calculate"
        status_text.color = ft.Colors.GREY_600
        page.update()
    
    def show_error(message):
        status_text.value = f"❌ Error: {message}"
        status_text.color = ft.Colors.RED_600
        page.update()
    
    # Create buttons
    calculate_btn = ft.ElevatedButton(
        text="Calculate",
        on_click=calculate_clicked,
        bgcolor=ft.Colors.BLUE,
        color=ft.Colors.WHITE,
    )
    
    clear_btn = ft.OutlinedButton(
        text="Clear",
        on_click=clear_clicked,
    )
    
    # Theme toggle
    def toggle_theme(e):
        page.theme_mode = ft.ThemeMode.DARK if page.theme_mode == ft.ThemeMode.LIGHT else ft.ThemeMode.LIGHT
        page.update()
    
    theme_btn = ft.IconButton(
        icon=ft.Icons.DARK_MODE,
        tooltip="Toggle Dark/Light Mode",
        on_click=toggle_theme,
    )
    
    # Build the page layout
    page.add(
        ft.Column([
            # Header
            ft.Row([
                ft.Text(
                    "Rectangle Calculator",
                    size=24,
                    weight=ft.FontWeight.BOLD,
                ),
                theme_btn,
            ], alignment=ft.MainAxisAlignment.SPACE_BETWEEN),
            
            ft.Divider(),
            
            # Input section
            ft.Text("Enter Dimensions:", size=16, weight=ft.FontWeight.W_500),
            width_field,
            height_field,
            
            ft.Divider(),
            
            # Buttons
            ft.Row([
                calculate_btn,
                clear_btn,
            ], alignment=ft.MainAxisAlignment.CENTER),
            
            ft.Divider(),
            
            # Results
            ft.Text("Results:", size=16, weight=ft.FontWeight.W_500),
            area_result,
            perimeter_result,
            
            ft.Divider(),
            
            # Status
            status_text,
            
        ], spacing=10, scroll=ft.ScrollMode.AUTO)
    )


if __name__ == "__main__":
    ft.app(target=main)