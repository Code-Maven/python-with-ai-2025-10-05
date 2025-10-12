# Python Learning Project with AI Tools

## Project Context
This is a Python educational project focused on learning Python through AI-assisted development. It demonstrates progressive complexity from basic CLI programs to GUI applications, following a tutorial-based learning path (in Hebrew).

## Architecture & Core Patterns

### Single-Function Core Pattern

## Development Workflow

### Environment Management
- **Always use `uv` for package management** (not pip/venv)
- Python version: 3.12+ (specified in `.python-version`)
- Dependencies managed in `pyproject.toml`
- Run commands with: `uv run python <file.py>`

### Testing Approach
- Uses **pytest with plain functions** (not unittest classes)
- Pattern: Test both the core function and main() with mocking
- Example test structure:
```python
def test_calculate_rectangle_properties_positive_integers():
    area, perimeter = calculate_rectangle_properties(5, 3)
    assert area == 15
    assert perimeter == 16
```

## Project-Specific Conventions

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

### Command-Line Interface
- use argparse

## Development Guidelines

1. **Reuse the core calculation function** - never reimplement the math
1. **Use `uv run` for all Python execution**
1. **Write pytest functions, not unittest classes**
1. **Include comprehensive input validation** in all implementations
1. **Maintain consistent error messaging** across all interfaces
1. Add type annotation to every function.
1. Make sure none of the function names are defined more than once.
1. Write tests for every function.

## Testing Commands

```bash
uv run pytest                    # Run a
```
