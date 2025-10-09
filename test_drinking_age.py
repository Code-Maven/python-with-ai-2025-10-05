import pytest
from unittest.mock import patch
import sys
import io
from drinking_age import check_drinking_age, main


# Tests for check_drinking_age function

def test_check_drinking_age_us_legal():
    """Test legal drinking age in US (21)."""
    can_drink, explanation = check_drinking_age(21, "US")
    assert can_drink == True
    assert "Yes, you can legally drink" in explanation
    assert "legal age: 21" in explanation


def test_check_drinking_age_us_underage():
    """Test underage person in US."""
    can_drink, explanation = check_drinking_age(18, "US")
    assert can_drink == False
    assert "you must wait 3 more years" in explanation
    assert "legal age: 21" in explanation


def test_check_drinking_age_uk_legal():
    """Test legal drinking age in UK (18)."""
    can_drink, explanation = check_drinking_age(18, "UK")
    assert can_drink == True
    assert "Yes, you can legally drink" in explanation
    assert "legal age: 18" in explanation


def test_check_drinking_age_germany_beer_wine_legal():
    """Test legal drinking age for beer/wine in Germany (16)."""
    can_drink, explanation = check_drinking_age(16, "DE", "beer_wine")
    assert can_drink == True
    assert "Yes, you can legally drink beer and wine" in explanation
    assert "legal age: 16" in explanation


def test_check_drinking_age_germany_spirits_legal():
    """Test legal drinking age for spirits in Germany (18)."""
    can_drink, explanation = check_drinking_age(18, "DE", "spirits")
    assert can_drink == True
    assert "Yes, you can legally drink spirits" in explanation
    assert "legal age: 18" in explanation


def test_check_drinking_age_germany_spirits_underage():
    """Test underage for spirits in Germany."""
    can_drink, explanation = check_drinking_age(17, "DE", "spirits")
    assert can_drink == False
    assert "you must wait 1 more year" in explanation
    assert "spirits" in explanation


def test_check_drinking_age_germany_all_alcohol_17():
    """Test 17-year-old with all alcohol types in Germany."""
    can_drink, explanation = check_drinking_age(17, "DE", "all")
    assert can_drink == True
    assert "beer and wine" in explanation
    assert "must wait 1 more year for spirits" in explanation


def test_check_drinking_age_germany_all_alcohol_18():
    """Test 18-year-old with all alcohol types in Germany."""
    can_drink, explanation = check_drinking_age(18, "DE", "all")
    assert can_drink == True
    assert "all types of alcohol" in explanation
    assert "beer/wine: 16+, spirits: 18+" in explanation


def test_check_drinking_age_germany_all_alcohol_underage():
    """Test underage for all alcohol in Germany."""
    can_drink, explanation = check_drinking_age(15, "DE", "all")
    assert can_drink == False
    assert "you must wait 1 more year" in explanation
    assert "beer/wine" in explanation


def test_check_drinking_age_case_insensitive():
    """Test that country codes are case-insensitive."""
    can_drink1, _ = check_drinking_age(21, "us")
    can_drink2, _ = check_drinking_age(21, "US")
    can_drink3, _ = check_drinking_age(21, "Us")
    
    assert can_drink1 == can_drink2 == can_drink3 == True


def test_check_drinking_age_invalid_country():
    """Test with invalid country code."""
    can_drink, explanation = check_drinking_age(25, "XX")
    assert can_drink == False
    assert "not in our database" in explanation


def test_check_drinking_age_edge_case_exact_age():
    """Test exact legal age boundary."""
    can_drink, explanation = check_drinking_age(21, "US")
    assert can_drink == True
    
    can_drink, explanation = check_drinking_age(20, "US")
    assert can_drink == False


def test_check_drinking_age_single_year_wait():
    """Test singular 'year' vs plural 'years' in message."""
    can_drink, explanation = check_drinking_age(20, "US")
    assert "1 more year to" in explanation  # singular
    
    can_drink, explanation = check_drinking_age(19, "US")
    assert "2 more years to" in explanation  # plural


def test_check_drinking_age_multiple_countries():
    """Test different countries with same age."""
    # 18-year-old in different countries
    can_drink_uk, _ = check_drinking_age(18, "UK")
    can_drink_fr, _ = check_drinking_age(18, "FR")
    can_drink_au, _ = check_drinking_age(18, "AU")
    can_drink_us, _ = check_drinking_age(18, "US")
    
    assert can_drink_uk == can_drink_fr == can_drink_au == True
    assert can_drink_us == False


def test_check_drinking_age_very_old():
    """Test with very old age."""
    can_drink, explanation = check_drinking_age(99, "US")
    assert can_drink == True
    assert "Yes, you can legally drink" in explanation


# Tests for main function

def test_main_legal_age_us():
    """Test main function with legal age in US."""
    user_inputs = ["21", "US"]
    
    with patch('builtins.input', side_effect=user_inputs):
        with patch('sys.stdout', new_callable=io.StringIO) as mock_stdout:
            main()
            output = mock_stdout.getvalue()
            
            assert "Legal Drinking Age Checker" in output
            assert "Age: 21" in output
            assert "Country: US" in output
            assert "Yes, you can legally drink" in output
            assert "drink responsibly" in output


def test_main_underage_us():
    """Test main function with underage person in US."""
    user_inputs = ["18", "US"]
    
    with patch('builtins.input', side_effect=user_inputs):
        with patch('sys.stdout', new_callable=io.StringIO) as mock_stdout:
            main()
            output = mock_stdout.getvalue()
            
            assert "Age: 18" in output
            assert "Country: US" in output
            assert "you must wait 3 more years" in output


def test_main_default_country():
    """Test main function with default country (US)."""
    user_inputs = ["25", ""]  # Empty string for country
    
    with patch('builtins.input', side_effect=user_inputs):
        with patch('sys.stdout', new_callable=io.StringIO) as mock_stdout:
            main()
            output = mock_stdout.getvalue()
            
            assert "Country: US" in output
            assert "Yes, you can legally drink" in output


def test_main_young_person_advisory():
    """Test main function with very young person."""
    user_inputs = ["14", "US"]
    
    with patch('builtins.input', side_effect=user_inputs):
        with patch('sys.stdout', new_callable=io.StringIO) as mock_stdout:
            main()
            output = mock_stdout.getvalue()
            
            assert "not recommended for anyone under 16" in output


def test_main_teen_advisory():
    """Test main function with teenager."""
    user_inputs = ["17", "DE", "1"]  # Legal for beer/wine in Germany but still gets advisory
    
    with patch('builtins.input', side_effect=user_inputs):
        with patch('sys.stdout', new_callable=io.StringIO) as mock_stdout:
            main()
            output = mock_stdout.getvalue()
            
            assert "responsible alcohol consumption is important" in output


def test_main_negative_age():
    """Test main function with negative age."""
    user_inputs = ["-5", "US"]
    
    with patch('builtins.input', side_effect=user_inputs):
        with patch('sys.stdout', new_callable=io.StringIO) as mock_stdout:
            main()
            output = mock_stdout.getvalue()
            
            assert "Error: Age cannot be negative" in output


def test_main_unrealistic_age():
    """Test main function with unrealistic age."""
    user_inputs = ["200", "US"]
    
    with patch('builtins.input', side_effect=user_inputs):
        with patch('sys.stdout', new_callable=io.StringIO) as mock_stdout:
            main()
            output = mock_stdout.getvalue()
            
            assert "Error: Please enter a realistic age" in output


def test_main_invalid_age_format():
    """Test main function with invalid age format."""
    user_inputs = ["twenty", "US"]
    
    with patch('builtins.input', side_effect=user_inputs):
        with patch('sys.stdout', new_callable=io.StringIO) as mock_stdout:
            main()
            output = mock_stdout.getvalue()
            
            assert "Error: Please enter a valid number for age" in output


def test_main_keyboard_interrupt():
    """Test main function handles KeyboardInterrupt gracefully."""
    with patch('builtins.input', side_effect=KeyboardInterrupt):
        with patch('sys.stdout', new_callable=io.StringIO) as mock_stdout:
            main()
            output = mock_stdout.getvalue()
            
            assert "Program interrupted by user" in output


def test_main_invalid_country():
    """Test main function with invalid country code."""
    user_inputs = ["25", "XYZ"]
    
    with patch('builtins.input', side_effect=user_inputs):
        with patch('sys.stdout', new_callable=io.StringIO) as mock_stdout:
            main()
            output = mock_stdout.getvalue()
            
            assert "not in our database" in output


def test_main_various_countries():
    """Test main function with different countries."""
    # Test Germany (lower drinking age for beer/wine)
    user_inputs = ["17", "DE", "1"]  # Choose beer/wine option
    
    with patch('builtins.input', side_effect=user_inputs):
        with patch('sys.stdout', new_callable=io.StringIO) as mock_stdout:
            main()
            output = mock_stdout.getvalue()
            
            assert "Country: DE" in output
            assert "Yes, you can legally drink" in output


def test_main_edge_case_zero_age():
    """Test main function with zero age."""
    user_inputs = ["0", "US"]
    
    with patch('builtins.input', side_effect=user_inputs):
        with patch('sys.stdout', new_callable=io.StringIO) as mock_stdout:
            main()
            output = mock_stdout.getvalue()
            
            assert "you must wait 21 more years" in output
            assert "not recommended for anyone under 16" in output


def test_main_boundary_ages():
    """Test main function with boundary ages for advisories."""
    # Test age 16 (boundary for first advisory)
    user_inputs = ["16", "US"]
    
    with patch('builtins.input', side_effect=user_inputs):
        with patch('sys.stdout', new_callable=io.StringIO) as mock_stdout:
            main()
            output = mock_stdout.getvalue()
            
            assert "responsible alcohol consumption is important" in output
            assert "not recommended for anyone under 16" not in output


def test_main_whitespace_handling():
    """Test main function handles whitespace in input."""
    user_inputs = ["  21  ", "  US  "]
    
    with patch('builtins.input', side_effect=user_inputs):
        with patch('sys.stdout', new_callable=io.StringIO) as mock_stdout:
            main()
            output = mock_stdout.getvalue()
            
            assert "Age: 21" in output
            assert "Country: US" in output
            assert "Yes, you can legally drink" in output


def test_main_germany_spirits_choice():
    """Test main function with Germany spirits choice."""
    user_inputs = ["17", "DE", "2"]  # Choose spirits option
    
    with patch('builtins.input', side_effect=user_inputs):
        with patch('sys.stdout', new_callable=io.StringIO) as mock_stdout:
            main()
            output = mock_stdout.getvalue()
            
            assert "Country: DE" in output
            assert "you must wait 1 more year" in output
            assert "spirits" in output


def test_main_germany_all_alcohol_choice():
    """Test main function with Germany all alcohol choice."""
    user_inputs = ["17", "DE", "3"]  # Choose all alcohol types
    
    with patch('builtins.input', side_effect=user_inputs):
        with patch('sys.stdout', new_callable=io.StringIO) as mock_stdout:
            main()
            output = mock_stdout.getvalue()
            
            assert "Country: DE" in output
            assert "beer and wine" in output
            assert "must wait 1 more year for spirits" in output


def test_main_germany_default_choice():
    """Test main function with Germany default choice (empty input)."""
    user_inputs = ["18", "DE", ""]  # Empty choice defaults to all
    
    with patch('builtins.input', side_effect=user_inputs):
        with patch('sys.stdout', new_callable=io.StringIO) as mock_stdout:
            main()
            output = mock_stdout.getvalue()
            
            assert "Country: DE" in output
            assert "all types of alcohol" in output