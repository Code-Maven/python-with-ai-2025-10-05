import tkinter as tk
from tkinter import ttk, messagebox
import sys


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


class RectangleCalculatorGUI:
    """GUI application for rectangle calculations."""
    
    def __init__(self, root):
        self.root = root
        self.root.title("Rectangle Calculator")
        self.root.geometry("400x300")
        self.root.resizable(True, True)
        
        # Configure style
        style = ttk.Style()
        style.theme_use('clam')
        
        self.setup_widgets()
        self.setup_layout()
        
        # Focus on width entry initially
        self.width_entry.focus()
    
    def setup_widgets(self):
        """Create and configure all GUI widgets."""
        # Main frame
        self.main_frame = ttk.Frame(self.root, padding="20")
        
        # Title
        self.title_label = ttk.Label(
            self.main_frame, 
            text="Rectangle Calculator",
            font=("Arial", 16, "bold")
        )
        
        # Input section
        self.input_frame = ttk.LabelFrame(self.main_frame, text="Dimensions", padding="10")
        
        # Width input
        self.width_label = ttk.Label(self.input_frame, text="Width:")
        self.width_var = tk.StringVar()
        self.width_entry = ttk.Entry(
            self.input_frame, 
            textvariable=self.width_var,
            width=15
        )
        
        # Height input
        self.height_label = ttk.Label(self.input_frame, text="Height:")
        self.height_var = tk.StringVar()
        self.height_entry = ttk.Entry(
            self.input_frame, 
            textvariable=self.height_var,
            width=15
        )
        
        # Buttons
        self.button_frame = ttk.Frame(self.main_frame)
        
        self.calculate_btn = ttk.Button(
            self.button_frame,
            text="Calculate",
            command=self.calculate,
            style="Accent.TButton"
        )
        
        self.clear_btn = ttk.Button(
            self.button_frame,
            text="Clear",
            command=self.clear_all
        )
        
        # Results section
        self.results_frame = ttk.LabelFrame(self.main_frame, text="Results", padding="10")
        
        # Area result
        self.area_label = ttk.Label(self.results_frame, text="Area:")
        self.area_result = ttk.Label(
            self.results_frame, 
            text="--",
            font=("Arial", 10, "bold"),
            foreground="blue"
        )
        
        # Perimeter result
        self.perimeter_label = ttk.Label(self.results_frame, text="Perimeter:")
        self.perimeter_result = ttk.Label(
            self.results_frame, 
            text="--",
            font=("Arial", 10, "bold"),
            foreground="blue"
        )
        
        # Status bar
        self.status_var = tk.StringVar()
        self.status_var.set("Enter width and height, then click Calculate")
        self.status_label = ttk.Label(
            self.main_frame, 
            textvariable=self.status_var,
            relief="sunken",
            padding="5"
        )
        
        # Bind Enter key to calculate
        self.root.bind('<Return>', lambda event: self.calculate())
        self.root.bind('<KP_Enter>', lambda event: self.calculate())
    
    def setup_layout(self):
        """Arrange widgets using grid layout."""
        # Main frame
        self.main_frame.pack(fill="both", expand=True)
        
        # Title
        self.title_label.grid(row=0, column=0, columnspan=2, pady=(0, 20))
        
        # Input section
        self.input_frame.grid(row=1, column=0, columnspan=2, sticky="ew", pady=(0, 10))
        
        # Width input
        self.width_label.grid(row=0, column=0, sticky="w", padx=(0, 10))
        self.width_entry.grid(row=0, column=1, sticky="ew")
        
        # Height input
        self.height_label.grid(row=1, column=0, sticky="w", padx=(0, 10), pady=(10, 0))
        self.height_entry.grid(row=1, column=1, sticky="ew", pady=(10, 0))
        
        # Configure input frame columns
        self.input_frame.columnconfigure(1, weight=1)
        
        # Buttons
        self.button_frame.grid(row=2, column=0, columnspan=2, pady=(0, 10))
        self.calculate_btn.pack(side="left", padx=(0, 10))
        self.clear_btn.pack(side="left")
        
        # Results section
        self.results_frame.grid(row=3, column=0, columnspan=2, sticky="ew", pady=(0, 10))
        
        # Area result
        self.area_label.grid(row=0, column=0, sticky="w", padx=(0, 10))
        self.area_result.grid(row=0, column=1, sticky="w")
        
        # Perimeter result
        self.perimeter_label.grid(row=1, column=0, sticky="w", padx=(0, 10), pady=(10, 0))
        self.perimeter_result.grid(row=1, column=1, sticky="w", pady=(10, 0))
        
        # Status bar
        self.status_label.grid(row=4, column=0, columnspan=2, sticky="ew", pady=(10, 0))
        
        # Configure main frame columns
        self.main_frame.columnconfigure(0, weight=1)
    
    def calculate(self):
        """Calculate and display rectangle properties."""
        try:
            # Get input values
            width_str = self.width_var.get().strip()
            height_str = self.height_var.get().strip()
            
            # Validate input
            if not width_str:
                self.show_error("Please enter a width value")
                self.width_entry.focus()
                return
            
            if not height_str:
                self.show_error("Please enter a height value")
                self.height_entry.focus()
                return
            
            # Convert to float
            try:
                width = float(width_str)
                height = float(height_str)
            except ValueError:
                self.show_error("Please enter valid numbers for width and height")
                return
            
            # Validate positive numbers
            if width <= 0:
                self.show_error("Width must be a positive number")
                self.width_entry.focus()
                self.width_entry.select_range(0, tk.END)
                return
            
            if height <= 0:
                self.show_error("Height must be a positive number")
                self.height_entry.focus()
                self.height_entry.select_range(0, tk.END)
                return
            
            # Calculate properties
            area, perimeter = calculate_rectangle_properties(width, height)
            
            # Display results
            self.area_result.config(text=f"{area}")
            self.perimeter_result.config(text=f"{perimeter}")
            
            # Update status
            self.status_var.set(f"Calculated for rectangle {width} × {height}")
            
        except Exception as e:
            self.show_error(f"An error occurred: {str(e)}")
    
    def clear_all(self):
        """Clear all input fields and results."""
        self.width_var.set("")
        self.height_var.set("")
        self.area_result.config(text="--")
        self.perimeter_result.config(text="--")
        self.status_var.set("Enter width and height, then click Calculate")
        self.width_entry.focus()
    
    def show_error(self, message):
        """Display error message and update status."""
        messagebox.showerror("Input Error", message)
        self.status_var.set(f"Error: {message}")


def main():
    """Main function to run the GUI application."""
    try:
        # Create main window
        root = tk.Tk()
        
        # Create application
        app = RectangleCalculatorGUI(root)
        
        # Set minimum window size
        root.minsize(350, 280)
        
        # Center window on screen
        root.update_idletasks()
        width = root.winfo_width()
        height = root.winfo_height()
        x = (root.winfo_screenwidth() // 2) - (width // 2)
        y = (root.winfo_screenheight() // 2) - (height // 2)
        root.geometry(f"{width}x{height}+{x}+{y}")
        
        # Start the GUI event loop
        root.mainloop()
        
    except Exception as e:
        print(f"Error starting GUI application: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()