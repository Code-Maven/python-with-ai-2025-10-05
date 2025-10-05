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
        color=ft.colors.GREY_600,
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
            area_result.color = ft.colors.BLUE
            perimeter_result.value = f"Perimeter: {perimeter:.2f}"
            perimeter_result.color = ft.colors.GREEN
            
            # Update status
            status_text.value = f"✅ Calculated for rectangle {width} × {height}"
            status_text.color = ft.colors.GREEN_600
            
            page.update()
            
        except Exception as ex:
            show_error(f"An error occurred: {str(ex)}")
    
    def clear_clicked(e):
        width_field.value = ""
        height_field.value = ""
        area_result.value = "Area: --"
        area_result.color = ft.colors.BLACK
        perimeter_result.value = "Perimeter: --"
        perimeter_result.color = ft.colors.BLACK
        status_text.value = "Enter width and height, then calculate"
        status_text.color = ft.colors.GREY_600
        page.update()
    
    def show_error(message):
        status_text.value = f"❌ Error: {message}"
        status_text.color = ft.colors.RED_600
        page.update()
    
    # Create buttons
    calculate_btn = ft.ElevatedButton(
        text="Calculate",
        on_click=calculate_clicked,
        bgcolor=ft.colors.BLUE,
        color=ft.colors.WHITE,
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
        icon=ft.icons.DARK_MODE,
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
    
    # Result display cards
    area_card = ft.Card(
        content=ft.Container(
            content=ft.Column([
                ft.Text("Area", size=14, color=ft.colors.GREY_600),
                ft.Text("--", size=24, weight=ft.FontWeight.BOLD, color=ft.colors.BLUE),
            ], tight=True),
            padding=15,
            width=180,
        ),
        elevation=3,
    )
    
    perimeter_card = ft.Card(
        content=ft.Container(
            content=ft.Column([
                ft.Text("Perimeter", size=14, color=ft.colors.GREY_600),
                ft.Text("--", size=24, weight=ft.FontWeight.BOLD, color=ft.colors.GREEN),
            ], tight=True),
            padding=15,
            width=180,
        ),
        elevation=3,
    )
    
    # Status text
    status_text = ft.Text(
        "Enter width and height, then calculate",
        size=12,
        color=ft.colors.GREY_600,
        text_align=ft.TextAlign.CENTER,
    )
    
    def calculate_clicked(e):
        try:
            # Get and validate input
            width_str = width_field.value.strip() if width_field.value else ""
            height_str = height_field.value.strip() if height_field.value else ""
            
            if not width_str:
                show_error("Please enter a width value")
                width_field.focus()
                return
            
            if not height_str:
                show_error("Please enter a height value")
                height_field.focus()
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
                width_field.focus()
                width_field.select_all()
                return
            
            if height <= 0:
                show_error("Height must be positive")
                height_field.focus() 
                height_field.select_all()
                return
            
            # Calculate results
            area, perimeter = calculate_rectangle_properties(width, height)
            
            # Update result cards
            area_card.content.content.controls[1].value = f"{area:.2f}"
            perimeter_card.content.content.controls[1].value = f"{perimeter:.2f}"
            
            # Update status
            status_text.value = f"✅ Calculated for rectangle {width} × {height}"
            status_text.color = ft.colors.GREEN_600
            
            page.update()
            
        except Exception as ex:
            show_error(f"An error occurred: {str(ex)}")
    
    def clear_clicked(e):
        width_field.value = ""
        height_field.value = ""
        area_card.content.content.controls[1].value = "--"
        perimeter_card.content.content.controls[1].value = "--"
        status_text.value = "Enter width and height, then calculate"
        status_text.color = ft.colors.GREY_600
        width_field.focus()
        page.update()
    
    def show_error(message):
        status_text.value = f"❌ Error: {message}"
        status_text.color = ft.colors.RED_600
        page.update()
        
        # Also show a snack bar for better UX
        page.show_snack_bar(
            ft.SnackBar(
                content=ft.Text(message),
                bgcolor=ft.colors.RED_400,
            )
        )
    
    def on_enter_key(e):
        if e.key == "Enter":
            calculate_clicked(e)
    
    # Add enter key support
    width_field.on_submit = lambda e: calculate_clicked(e)
    height_field.on_submit = lambda e: calculate_clicked(e)
    
    # Create buttons
    calculate_btn = ft.ElevatedButton(
        text="Calculate",
        icon=ft.icons.CALCULATE,
        on_click=calculate_clicked,
        style=ft.ButtonStyle(
            color=ft.colors.WHITE,
            bgcolor=ft.colors.BLUE,
            elevation=5,
        ),
        width=120,
        height=45,
    )
    
    clear_btn = ft.OutlinedButton(
        text="Clear",
        icon=ft.icons.CLEAR,
        on_click=clear_clicked,
        width=120,
        height=45,
    )
    
    # Theme toggle button
    def toggle_theme(e):
        page.theme_mode = "dark" if page.theme_mode == "light" else "light"
        page.update()
    
    theme_btn = ft.IconButton(
        icon=ft.icons.DARK_MODE,
        tooltip="Toggle Dark/Light Mode",
        on_click=toggle_theme,
    )
    
    # Build the page layout
    page.add(
        ft.Column([
            # Header
            ft.Container(
                content=ft.Row([
                    ft.Text(
                        "Rectangle Calculator",
                        size=28,
                        weight=ft.FontWeight.BOLD,
                        color=ft.colors.PRIMARY,
                    ),
                    theme_btn,
                ], alignment=ft.MainAxisAlignment.SPACE_BETWEEN),
                padding=ft.padding.only(bottom=20),
            ),
            
            # Input section
            ft.Card(
                content=ft.Container(
                    content=ft.Column([
                        ft.Text("Dimensions", size=18, weight=ft.FontWeight.W_500),
                        ft.Divider(height=10),
                        ft.Row([
                            width_field,
                            height_field,
                        ], alignment=ft.MainAxisAlignment.SPACE_BETWEEN),
                    ]),
                    padding=20,
                ),
                elevation=2,
            ),
            
            # Buttons
            ft.Container(
                content=ft.Row([
                    calculate_btn,
                    clear_btn,
                ], alignment=ft.MainAxisAlignment.CENTER),
                padding=ft.padding.symmetric(vertical=15),
            ),
            
            # Results section
            ft.Card(
                content=ft.Container(
                    content=ft.Column([
                        ft.Text("Results", size=18, weight=ft.FontWeight.W_500),
                        ft.Divider(height=10),
                        ft.Row([
                            area_card,
                            perimeter_card,
                        ], alignment=ft.MainAxisAlignment.SPACE_BETWEEN),
                    ]),
                    padding=20,
                ),
                elevation=2,
            ),
            
            # Status
            ft.Container(
                content=status_text,
                padding=ft.padding.only(top=10),
                alignment=ft.alignment.center,
            ),
            
        ], scroll=ft.ScrollMode.AUTO)
    )
    
    # Focus on first field
    width_field.focus()


if __name__ == "__main__":
    ft.app(target=main)