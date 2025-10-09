def check_drinking_age(age: int, country: str = "US") -> tuple[bool, str]:
    """
    Check if a person can legally drink alcohol based on age and country.
    
    Args:
        age (int): Person's age in years
        country (str): Country code (US, UK, DE, etc.) - defaults to US
    
    Returns:
        tuple: (can_drink: bool, explanation: str)
    """
    # Define legal drinking ages by country
    drinking_ages = {
        "US": 21,    # United States
        "UK": 18,    # United Kingdom
        "CA": 19,    # Canada (varies by province, using majority)
        "DE": 16,    # Germany (beer/wine), 18 for spirits
        "FR": 18,    # France
        "IT": 18,    # Italy
        "JP": 20,    # Japan
        "AU": 18,    # Australia
        "BR": 18,    # Brazil
        "IN": 21,    # India (varies by state)
    }
    
    country = country.upper()
    
    if country not in drinking_ages:
        return False, f"Legal drinking age for {country} is not in our database."
    
    legal_age = drinking_ages[country]
    
    if age >= legal_age:
        return True, f"Yes, you can legally drink alcohol in {country} (legal age: {legal_age})."
    else:
        years_to_wait = legal_age - age
        return False, f"No, you must wait {years_to_wait} more year{'s' if years_to_wait != 1 else ''} to legally drink in {country} (legal age: {legal_age})."


def main() -> None:
    try:
        # Get input from user
        print("Legal Drinking Age Checker")
        print("-" * 26)
        print("This program checks if you can legally drink alcohol based on your age and country.\n")
        
        # Get age
        age_input = input("Enter your age: ").strip()
        age = int(age_input)
        
        # Validate age
        if age < 0:
            print("Error: Age cannot be negative.")
            return
        elif age > 150:
            print("Error: Please enter a realistic age.")
            return
        
        # Get country (optional)
        country_input = input("Enter your country code (US, UK, CA, DE, etc.) or press Enter for US: ").strip()
        country = country_input if country_input else "US"
        
        # Check drinking eligibility
        can_drink, explanation = check_drinking_age(age, country)
        
        # Print results
        print(f"\nDrinking Age Analysis:")
        print(f"Age: {age}")
        print(f"Country: {country.upper()}")
        print(f"Result: {explanation}")
        
        # Additional context
        if age < 16:
            print("\nNote: Regardless of local laws, alcohol consumption is not recommended for anyone under 16.")
        elif age < 18:
            print("\nNote: Even where legal, responsible alcohol consumption is important.")
        elif can_drink:
            print("\nNote: Please drink responsibly if you choose to consume alcohol.")
        
    except ValueError:
        print("Error: Please enter a valid number for age.")
    except KeyboardInterrupt:
        print("\nProgram interrupted by user.")


if __name__ == "__main__":
    main()