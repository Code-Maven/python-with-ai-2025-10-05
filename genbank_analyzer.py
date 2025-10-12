import sys
import re
from pathlib import Path
from typing import Optional


def parse_genbank_sequence(file_path: str) -> tuple[str, int]:
    """
    Parse a GenBank format file and extract the sequence.
    
    Args:
        file_path (str): Path to the GenBank format file
    
    Returns:
        tuple: (sequence, length) where sequence is the DNA/RNA sequence
               and length is the number of nucleotides
    """
    sequence = ""
    in_sequence_section = False
    
    try:
        with open(file_path, 'r', encoding='utf-8') as file:
            for line in file:
                line = line.strip()
                
                # Check if we've reached the sequence section
                if line.startswith('ORIGIN'):
                    in_sequence_section = True
                    continue
                
                # Check if we've reached the end of the file
                if line.startswith('//'):
                    break
                
                # If we're in the sequence section, extract the sequence
                if in_sequence_section:
                    # Remove line numbers and spaces, keep only nucleotides
                    # GenBank sequence lines look like: "     1 atgcgtacgt acgtacgtac gtacgtacgt"
                    cleaned_line = re.sub(r'^\s*\d+\s*', '', line)  # Remove line numbers
                    cleaned_line = re.sub(r'\s+', '', cleaned_line)  # Remove all spaces
                    sequence += cleaned_line.lower()
    
    except FileNotFoundError:
        raise FileNotFoundError(f"File not found: {file_path}")
    except Exception as e:
        raise Exception(f"Error reading file {file_path}: {e}")
    
    return sequence, len(sequence)


def get_genbank_info(file_path: str) -> dict:
    """
    Extract basic information from GenBank file header.
    
    Args:
        file_path (str): Path to the GenBank format file
    
    Returns:
        dict: Dictionary with basic file information
    """
    info = {}
    
    try:
        with open(file_path, 'r', encoding='utf-8') as file:
            for line in file:
                line = line.strip()
                
                if line.startswith('LOCUS'):
                    # Extract locus name and sequence length from LOCUS line
                    # Format: LOCUS       ON205948                6076 bp    RNA     linear   VRL 16-JUN-2022
                    parts = line.split()
                    if len(parts) >= 1:
                        info['locus'] = parts[1] if len(parts) > 1 else ''
                        # Look for the pattern: number followed by 'bp'
                        for i, part in enumerate(parts):
                            if part == 'bp' and i > 0:
                                # The number should be in the previous part
                                if parts[i-1].isdigit():
                                    info['declared_length'] = int(parts[i-1])
                                    break
                
                elif line.startswith('DEFINITION'):
                    info['definition'] = line[10:].strip()
                
                elif line.startswith('ACCESSION'):
                    info['accession'] = line[9:].strip()
                
                elif line.startswith('SOURCE'):
                    info['source'] = line[6:].strip()
                
                elif line.startswith('ORIGIN'):
                    break
    
    except Exception as e:
        print(f"Warning: Could not extract header information: {e}")
    
    return info


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
            composition = {
                'a': sequence.count('a'),
                't': sequence.count('t'),
                'g': sequence.count('g'),
                'c': sequence.count('c'),
                'n': sequence.count('n')  # Unknown nucleotides
            }
            
            print(f"\nNucleotide Composition:")
            for nucleotide, count in composition.items():
                percentage = (count / actual_length) * 100 if actual_length > 0 else 0
                print(f"  {nucleotide.upper()}: {count:,} ({percentage:.1f}%)")
            
            # Calculate GC content
            gc_count = composition['g'] + composition['c']
            gc_content = (gc_count / actual_length) * 100 if actual_length > 0 else 0
            print(f"\nGC Content: {gc_content:.1f}%")
        
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