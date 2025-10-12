#!/usr/bin/env python3
"""
Example script demonstrating how to use the GenBankAnalyzer class programmatically.
"""

from genbank_analyzer import GenBankAnalyzer


def main():
    """Demonstrate programmatic usage of GenBankAnalyzer."""
    
    # Example file path - replace with your own GenBank file
    file_path = "biodata/ON205948.1.gb"
    
    try:
        # Create analyzer instance
        print("Creating GenBank analyzer...")
        analyzer = GenBankAnalyzer(file_path)
        
        # Basic file information
        print("\n=== Basic Information ===")
        info = analyzer.get_basic_info()
        print(f"Accession: {info.get('accession', 'N/A')}")
        print(f"Organism: {info.get('source', 'N/A')}")
        print(f"Definition: {info.get('definition', 'N/A')}")
        
        # Sequence statistics
        print(f"\n=== Sequence Statistics ===")
        print(f"Length: {analyzer.length:,} nucleotides")
        print(f"GC Content: {analyzer.get_gc_content():.1f}%")
        
        # Nucleotide composition
        print(f"\n=== Nucleotide Composition ===")
        composition = analyzer.get_nucleotide_composition()
        for nucleotide, count in sorted(composition.items()):
            percentage = (count / analyzer.length) * 100
            print(f"{nucleotide.upper()}: {count:,} ({percentage:.1f}%)")
        
        # Features summary
        print(f"\n=== Features Summary ===")
        feature_summary = analyzer.get_feature_summary()
        for feat_type, count in sorted(feature_summary.items()):
            print(f"{feat_type}: {count}")
        
        # Validation
        print(f"\n=== Validation ===")
        is_valid_extension = analyzer.validate_file_extension()
        print(f"Valid file extension: {is_valid_extension}")
        
        length_matches, declared_length = analyzer.validate_sequence_length()
        if declared_length:
            print(f"Length validation: {'✓ PASS' if length_matches else '✗ FAIL'}")
            print(f"Declared: {declared_length}, Actual: {analyzer.length}")
        
        # Sequence preview
        print(f"\n=== Sequence Preview ===")
        preview = analyzer.get_sequence_preview(100)
        print(f"First 100 nucleotides:")
        print(f"{preview}")
        
        # Individual feature details (first 3 CDS features)
        print(f"\n=== CDS Features (first 3) ===")
        features = analyzer.get_features()
        cds_count = 0
        for feature in features:
            if feature['type'] == 'CDS' and cds_count < 3:
                print(f"CDS {cds_count + 1}:")
                print(f"  Location: {feature['location']}")
                if 'product' in feature:
                    print(f"  Product: {feature['product']}")
                if 'gene' in feature:
                    print(f"  Gene: {feature['gene']}")
                cds_count += 1
        
        print(f"\n=== Analysis Complete ===")
        print(f"Successfully analyzed {analyzer.length:,} nucleotides from {len(features)} features")
        
    except FileNotFoundError:
        print(f"Error: File '{file_path}' not found.")
        print("Please make sure you have a GenBank file at the specified path.")
    except Exception as e:
        print(f"Error analyzing file: {e}")


if __name__ == "__main__":
    main()