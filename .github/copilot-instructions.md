# Python Learning Project with AI Tools

## Project Context
This is a Python educational project focused on learning Python through AI-assisted development. It demonstrates progressive complexity from basic CLI programs to GUI applications, following a tutorial-based learning path (in Hebrew).

## Architecture & Core Patterns

### Single-Function Core Pattern
All rectangle calculators share the same core function:
```python
def calculate_rectangle_properties(width, height):
    area = width * height
    perimeter = 2 * (width + height)
    return area, perimeter
```
**Key**: Always reuse this exact function signature across all rectangle implementations to maintain consistency.

### Progressive UI Implementations
The project follows a specific progression pattern:
1. `rectangle.py` - Interactive CLI with input validation
2. `rectangle_argv.py` - Command-line arguments using argparse
3. `rectangle_gui.py` - Tkinter GUI with ttk styling
4. `rectangle_flet.py` - Modern Flet-based GUI (multi-implementation file)
5. `rectangle_flet_fixed.py` - Refined Flet implementation

## Development Workflow

### Environment Management
- **Always use `uv` for package management** (not pip/venv)
- Python version: 3.12+ (specified in `.python-version`)
- Dependencies managed in `pyproject.toml`
- Run commands with: `uv run python <file.py>`

### Testing Approach
- Uses **pytest with plain functions** (not unittest classes)
- Test file: `test_rectangle.py` 
- Pattern: Test both the core function and main() with mocking
- Example test structure:
```python
def test_calculate_rectangle_properties_positive_integers():
    area, perimeter = calculate_rectangle_properties(5, 3)
    assert area == 15
    assert perimeter == 16
```

### GUI Framework Progression
1. **Tkinter** (`rectangle_gui.py`) - Basic GUI, class-based structure
2. **Flet** (`rectangle_flet*.py`) - Modern web-based GUI, function-based structure

## Project-Specific Conventions

### File Naming Pattern
- Base functionality: `rectangle.py`
- Variants: `rectangle_<interface_type>.py` (argv, gui, flet)
- Fixed versions: `rectangle_<type>_fixed.py`

### Error Handling Standards
- CLI: Try-catch with user-friendly error messages
- GUI: Visual feedback with status text and message boxes
- Always validate positive numbers for dimensions

### Documentation Style
- Function docstrings with Args/Returns sections
- Inline comments for GUI setup and validation logic
- Class docstrings for GUI classes: `"""GUI application for rectangle calculations."""`

## Key Integration Points

### Shared Core Logic
Every rectangle calculator must import and use the same `calculate_rectangle_properties()` function to ensure consistent calculations across all interfaces.

### GUI State Management
- Flet implementations use nested functions for event handlers
- Tkinter uses class methods for state management
- Both patterns include clear/reset functionality

### Command-Line Interface
`rectangle_argv.py` uses argparse with:
- Positional arguments for width/height
- Optional verbose flag (`-v, --verbose`)
- Built-in help and validation

## Development Guidelines

1. **Reuse the core calculation function** - never reimplement the math
2. **Follow the progressive complexity pattern** when adding new interfaces
3. **Use `uv run` for all Python execution** 
4. **Write pytest functions, not unittest classes**
5. **Include comprehensive input validation** in all implementations
6. **Maintain consistent error messaging** across all interfaces

## Testing Commands
```bash
uv run pytest                    # Run all tests
uv run python rectangle.py       # Interactive CLI
uv run python rectangle_argv.py 5 3  # Command-line usage
uv run python rectangle_gui.py   # Tkinter GUI
uv run python rectangle_flet.py  # Modern Flet GUI
```