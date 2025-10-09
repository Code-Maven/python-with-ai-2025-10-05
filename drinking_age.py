def check_drinking_age(age: int, country: str = "US", alcohol_type: str = "all") -> tuple[bool, str]:
    """
    Check if a person can legally drink alcohol based on age and country.
    
    Args:
        age (int): Person's age in years
        country (str): Country code (US, UK, DE, etc.) - defaults to US
        alcohol_type (str): Type of alcohol ("all", "beer_wine", "spirits") - defaults to "all"
    
    Returns:
        tuple: (can_drink: bool, explanation: str)
    """
    # Define legal drinking ages by country
    drinking_ages = {
        "US": 21,    # United States
        "UK": 18,    # United Kingdom
        "CA": 19,    # Canada (varies by province, using majority)
        "DE": {"beer_wine": 16, "spirits": 18},  # Germany has different ages
        "FR": 18,    # France
        "IT": 18,    # Italy
        "JP": 20,    # Japan
        "AU": 18,    # Australia
        "BR": 18,    # Brazil
        "IN": 21,    # India (varies by state)
    }
    
    country = country.upper()
    alcohol_type = alcohol_type.lower()
    
    if country not in drinking_ages:
        return False, f"Legal drinking age for {country} is not in our database."
    
    # Handle Germany's special case
    if country == "DE" and isinstance(drinking_ages[country], dict):
        if alcohol_type == "beer_wine":
            legal_age = drinking_ages[country]["beer_wine"]
            alcohol_desc = "beer and wine"
        elif alcohol_type == "spirits":
            legal_age = drinking_ages[country]["spirits"]
            alcohol_desc = "spirits (hard liquor)"
        else:  # "all" or any other value
            # For "all" alcohol types, use the higher age requirement
            beer_wine_age = drinking_ages[country]["beer_wine"]
            spirits_age = drinking_ages[country]["spirits"]
            
            if age >= spirits_age:
                return True, f"Yes, you can legally drink all types of alcohol in {country} (beer/wine: {beer_wine_age}+, spirits: {spirits_age}+)."
            elif age >= beer_wine_age:
                years_to_spirits = spirits_age - age
                return True, f"Yes, you can legally drink beer and wine in {country} (age {beer_wine_age}+), but must wait {years_to_spirits} more year{'s' if years_to_spirits != 1 else ''} for spirits (age {spirits_age}+)."
            else:
                years_to_beer_wine = beer_wine_age - age
                return False, f"No, you must wait {years_to_beer_wine} more year{'s' if years_to_beer_wine != 1 else ''} to legally drink beer/wine in {country} (age {beer_wine_age}+)."
    else:
        # Handle all other countries with single age requirement
        legal_age = drinking_ages[country]
        alcohol_desc = "alcohol"
    
    if age >= legal_age:
        return True, f"Yes, you can legally drink {alcohol_desc} in {country} (legal age: {legal_age})."
    else:
        years_to_wait = legal_age - age
        return False, f"No, you must wait {years_to_wait} more year{'s' if years_to_wait != 1 else ''} to legally drink {alcohol_desc} in {country} (legal age: {legal_age})."


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
        
        # Get alcohol type for Germany
        alcohol_type = "all"
        if country.upper() == "DE":
            print("\nGermany has different legal ages for different types of alcohol:")
            print("1. Beer and wine (age 16+)")
            print("2. Spirits/hard liquor (age 18+)")
            print("3. All alcohol types")
            
            alcohol_choice = input("What type are you asking about? (1/2/3 or press Enter for all): ").strip()
            if alcohol_choice == "1":
                alcohol_type = "beer_wine"
            elif alcohol_choice == "2":
                alcohol_type = "spirits"
            else:
                alcohol_type = "all"
        
        # Check drinking eligibility
        can_drink, explanation = check_drinking_age(age, country, alcohol_type)
        
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