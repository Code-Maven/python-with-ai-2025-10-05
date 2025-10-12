import sys
from pathlib import Path
from typing import Optional
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


def parse_genbank_sequence(file_path: str) -> tuple[str, int]:
    """
    Parse a GenBank format file and extract the sequence using BioPython.
    
    Args:
        file_path (str): Path to the GenBank format file
    
    Returns:
        tuple: (sequence, length) where sequence is the DNA/RNA sequence
               and length is the number of nucleotides
    """
    try:
        record = SeqIO.read(file_path, "genbank")
        sequence = str(record.seq).lower()
        return sequence, len(sequence)
    
    except FileNotFoundError:
        raise FileNotFoundError(f"File not found: {file_path}")
    except Exception as e:
        raise Exception(f"Error reading GenBank file {file_path}: {e}")


def get_genbank_info(file_path: str) -> dict:
    """
    Extract basic information from GenBank file header using BioPython.
    
    Args:
        file_path (str): Path to the GenBank format file
    
    Returns:
        dict: Dictionary with basic file information
    """
    info = {}
    
    try:
        record: SeqRecord = SeqIO.read(file_path, "genbank")
        
        # Basic identifiers
        info['locus'] = record.name if record.name else record.id
        info['accession'] = record.id
        info['definition'] = record.description
        info['declared_length'] = len(record.seq)
        
        # Extract source organism from annotations
        if 'source' in record.annotations:
            info['source'] = record.annotations['source']
        elif 'organism' in record.annotations:
            info['source'] = record.annotations['organism']
        
        # Extract additional useful information
        if 'date' in record.annotations:
            info['date'] = record.annotations['date']
        
        if 'molecule_type' in record.annotations:
            info['molecule_type'] = record.annotations['molecule_type']
        
        if 'topology' in record.annotations:
            info['topology'] = record.annotations['topology']
    
    except Exception as e:
        print(f"Warning: Could not extract header information: {e}")
    
    return info


def get_genbank_features(file_path: str) -> list:
    """
    Extract feature information from GenBank file using BioPython.
    
    Args:
        file_path (str): Path to the GenBank format file
    
    Returns:
        list: List of feature dictionaries with type, location, and qualifiers
    """
    features = []
    
    try:
        record: SeqRecord = SeqIO.read(file_path, "genbank")
        
        for feature in record.features:
            feature_info = {
                'type': feature.type,
                'location': str(feature.location),
                'strand': getattr(feature.location, 'strand', None)
            }
            
            # Extract common qualifiers
            if feature.qualifiers:
                for key in ['gene', 'product', 'protein_id', 'translation', 'note']:
                    if key in feature.qualifiers:
                        feature_info[key] = feature.qualifiers[key][0]  # Take first value
            
            features.append(feature_info)
    
    except Exception as e:
        print(f"Warning: Could not extract feature information: {e}")
    
    return features


def main() -> None:
    try:
        print("GenBank Sequence Length Analyzer")
        print("-" * 32)
        print("This program reads GenBank format files and analyzes the sequence length.\n")
        
        # Get file path from user or command line argument
        if len(sys.argv) > 1:
            file_path = sys.argv[1]
        else:
            file_path = input("Enter the path to the GenBank file: ").strip()
        
        # Validate file exists
        if not Path(file_path).exists():
            print(f"Error: File '{file_path}' does not exist.")
            return
        
        # Validate file extension (optional warning)
        if not file_path.lower().endswith(('.gb', '.gbk', '.genbank')):
            print(f"Warning: File '{file_path}' doesn't have a typical GenBank extension (.gb, .gbk, .genbank)")
        
        print(f"Analyzing file: {file_path}")
        print("-" * 50)
        
        # Extract header information
        info = get_genbank_info(file_path)
        
        if info:
            print("File Information:")
            if 'locus' in info:
                print(f"  Locus: {info['locus']}")
            if 'accession' in info:
                print(f"  Accession: {info['accession']}")
            if 'source' in info:
                print(f"  Source: {info['source']}")
            if 'definition' in info:
                print(f"  Definition: {info['definition']}")
            if 'declared_length' in info:
                print(f"  Declared length: {info['declared_length']} bp")
            print()
        
        # Parse the sequence
        sequence, actual_length = parse_genbank_sequence(file_path)
        
        if not sequence:
            print("Warning: No sequence data found in the file.")
            return
        
        # Display results
        print("Sequence Analysis:")
        print(f"  Actual sequence length: {actual_length} nucleotides")
        
        # Verify against declared length if available
        if info.get('declared_length'):
            declared_length = info['declared_length']
            if actual_length == declared_length:
                print(f"  ✓ Length matches declared length ({declared_length} bp)")
            else:
                print(f"  ⚠ Length mismatch! Declared: {declared_length} bp, Actual: {actual_length} bp")
        
        # Show sequence composition
        if sequence:
            # Count all possible nucleotides including ambiguous ones
            nucleotides = ['a', 't', 'g', 'c', 'u', 'n']  # Include U for RNA
            composition = {}
            
            for nuc in nucleotides:
                count = sequence.count(nuc)
                if count > 0:  # Only show nucleotides that are present
                    composition[nuc] = count
            
            # Count any other characters (ambiguous nucleotides)
            total_counted = sum(composition.values())
            if total_counted < actual_length:
                other_count = actual_length - total_counted
                composition['other'] = other_count
            
            print(f"\nNucleotide Composition:")
            for nucleotide, count in composition.items():
                percentage = (count / actual_length) * 100 if actual_length > 0 else 0
                print(f"  {nucleotide.upper()}: {count:,} ({percentage:.1f}%)")
            
            # Calculate GC content (works for both DNA and RNA)
            gc_count = composition.get('g', 0) + composition.get('c', 0)
            gc_content = (gc_count / actual_length) * 100 if actual_length > 0 else 0
            print(f"\nGC Content: {gc_content:.1f}%")
        
        # Show features information
        features = get_genbank_features(file_path)
        if features:
            print(f"\nFeatures Found ({len(features)} total):")
            
            # Group features by type for better display
            feature_types = {}
            for feature in features:
                feat_type = feature['type']
                if feat_type not in feature_types:
                    feature_types[feat_type] = []
                feature_types[feat_type].append(feature)
            
            # Display feature summary
            for feat_type, feat_list in feature_types.items():
                print(f"  {feat_type}: {len(feat_list)}")
            
            # Show detailed information for important features
            important_features = ['CDS', 'gene', 'mRNA', 'tRNA', 'rRNA']
            shown_features = 0
            max_features_to_show = 5
            
            print(f"\nDetailed Feature Information (showing up to {max_features_to_show}):")
            for feature in features:
                if shown_features >= max_features_to_show:
                    remaining = len(features) - shown_features
                    if remaining > 0:
                        print(f"  ... and {remaining} more features")
                    break
                
                if feature['type'] in important_features or shown_features < 3:
                    print(f"  {feature['type']} at {feature['location']}")
                    if 'gene' in feature:
                        print(f"    Gene: {feature['gene']}")
                    if 'product' in feature:
                        print(f"    Product: {feature['product']}")
                    shown_features += 1
        
        # Show sequence preview
        if len(sequence) > 0:
            preview_length = min(60, len(sequence))
            print(f"\nSequence Preview (first {preview_length} nucleotides):")
            print(f"  {sequence[:preview_length].upper()}")
            if len(sequence) > preview_length:
                print("  ...")
        
    except FileNotFoundError as e:
        print(f"Error: {e}")
    except Exception as e:
        print(f"Error: {e}")
    except KeyboardInterrupt:
        print("\nProgram interrupted by user.")


if __name__ == "__main__":
    main()