import sys
from pathlib import Path
from typing import Optional, Dict, List, Tuple
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


class GenBankAnalyzer:
    """
    A class for analyzing GenBank format files using BioPython.
    
    This class provides methods to parse GenBank files, extract sequence information,
    analyze features, and calculate sequence composition statistics.
    """
    
    def __init__(self, file_path: str) -> None:
        """
        Initialize the GenBank analyzer with a file path.
        
        Args:
            file_path (str): Path to the GenBank format file
        
        Raises:
            FileNotFoundError: If the file doesn't exist
            ValueError: If the file path is invalid
        """
        self.file_path = Path(file_path)
        if not self.file_path.exists():
            raise FileNotFoundError(f"File not found: {file_path}")
        
        self._record: Optional[SeqRecord] = None
        self._sequence: Optional[str] = None
        self._info: Optional[Dict] = None
        self._features: Optional[List[Dict]] = None
        self._composition: Optional[Dict] = None
    
    @property
    def record(self) -> SeqRecord:
        """
        Get the BioPython SeqRecord object (lazy loading).
        
        Returns:
            SeqRecord: The parsed GenBank record
        """
        if self._record is None:
            try:
                self._record = SeqIO.read(str(self.file_path), "genbank")
            except Exception as e:
                raise Exception(f"Error reading GenBank file {self.file_path}: {e}")
        return self._record
    
    @property
    def sequence(self) -> str:
        """
        Get the DNA/RNA sequence as a lowercase string.
        
        Returns:
            str: The nucleotide sequence
        """
        if self._sequence is None:
            self._sequence = str(self.record.seq).lower()
        return self._sequence
    
    @property
    def length(self) -> int:
        """
        Get the length of the sequence.
        
        Returns:
            int: Number of nucleotides in the sequence
        """
        return len(self.sequence)
    
    def get_basic_info(self) -> Dict[str, str]:
        """
        Extract basic information from GenBank file header.
        
        Returns:
            Dict[str, str]: Dictionary with basic file information
        """
        if self._info is None:
            info = {}
            try:
                record = self.record
                
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
            
            self._info = info
        return self._info
    
    def get_features(self) -> List[Dict]:
        """
        Extract feature information from GenBank file.
        
        Returns:
            List[Dict]: List of feature dictionaries with type, location, and qualifiers
        """
        if self._features is None:
            features = []
            try:
                record = self.record
                
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
            
            self._features = features
        return self._features
    
    def get_nucleotide_composition(self) -> Dict[str, int]:
        """
        Calculate nucleotide composition of the sequence.
        
        Returns:
            Dict[str, int]: Dictionary with nucleotide counts
        """
        if self._composition is None:
            sequence = self.sequence
            # Count all possible nucleotides including ambiguous ones
            nucleotides = ['a', 't', 'g', 'c', 'u', 'n']  # Include U for RNA
            composition = {}
            
            for nuc in nucleotides:
                count = sequence.count(nuc)
                if count > 0:  # Only include nucleotides that are present
                    composition[nuc] = count
            
            # Count any other characters (ambiguous nucleotides)
            total_counted = sum(composition.values())
            if total_counted < len(sequence):
                other_count = len(sequence) - total_counted
                composition['other'] = other_count
            
            self._composition = composition
        return self._composition
    
    def get_gc_content(self) -> float:
        """
        Calculate GC content as a percentage.
        
        Returns:
            float: GC content percentage (0-100)
        """
        composition = self.get_nucleotide_composition()
        gc_count = composition.get('g', 0) + composition.get('c', 0)
        return (gc_count / self.length) * 100 if self.length > 0 else 0.0
    
    def validate_file_extension(self) -> bool:
        """
        Check if the file has a typical GenBank extension.
        
        Returns:
            bool: True if extension is typical GenBank format
        """
        return self.file_path.suffix.lower() in ['.gb', '.gbk', '.genbank']
    
    def validate_sequence_length(self) -> Tuple[bool, Optional[int]]:
        """
        Validate that actual sequence length matches declared length.
        
        Returns:
            Tuple[bool, Optional[int]]: (matches, declared_length)
        """
        info = self.get_basic_info()
        declared_length = info.get('declared_length')
        if declared_length is not None:
            return self.length == declared_length, declared_length
        return True, None  # No declared length to compare
    
    def get_feature_summary(self) -> Dict[str, int]:
        """
        Get a summary of features grouped by type.
        
        Returns:
            Dict[str, int]: Dictionary with feature types and their counts
        """
        features = self.get_features()
        feature_types = {}
        for feature in features:
            feat_type = feature['type']
            feature_types[feat_type] = feature_types.get(feat_type, 0) + 1
        return feature_types
    
    def get_sequence_preview(self, length: int = 60) -> str:
        """
        Get a preview of the sequence.
        
        Args:
            length (int): Number of nucleotides to include in preview
        
        Returns:
            str: Uppercase sequence preview
        """
        sequence = self.sequence
        preview_length = min(length, len(sequence))
        return sequence[:preview_length].upper()


class GenBankAnalyzerCLI:
    """
    Command-line interface for the GenBank analyzer.
    """
    
    def __init__(self) -> None:
        """Initialize the CLI."""
        self.analyzer: Optional[GenBankAnalyzer] = None
    
    def get_file_path(self) -> str:
        """
        Get file path from command line argument or user input.
        
        Returns:
            str: Path to the GenBank file
        """
        if len(sys.argv) > 1:
            return sys.argv[1]
        else:
            return input("Enter the path to the GenBank file: ").strip()
    
    def display_header(self) -> None:
        """Display the program header."""
        print("GenBank Sequence Length Analyzer")
        print("-" * 32)
        print("This program reads GenBank format files and analyzes the sequence length.\n")
    
    def display_file_info(self, analyzer: GenBankAnalyzer) -> None:
        """
        Display basic file information.
        
        Args:
            analyzer (GenBankAnalyzer): The analyzer instance
        """
        info = analyzer.get_basic_info()
        
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
    
    def display_sequence_analysis(self, analyzer: GenBankAnalyzer) -> None:
        """
        Display sequence analysis results.
        
        Args:
            analyzer (GenBankAnalyzer): The analyzer instance
        """
        print("Sequence Analysis:")
        print(f"  Actual sequence length: {analyzer.length} nucleotides")
        
        # Verify against declared length if available
        length_matches, declared_length = analyzer.validate_sequence_length()
        if declared_length is not None:
            if length_matches:
                print(f"  ✓ Length matches declared length ({declared_length} bp)")
            else:
                print(f"  ⚠ Length mismatch! Declared: {declared_length} bp, Actual: {analyzer.length} bp")
    
    def display_nucleotide_composition(self, analyzer: GenBankAnalyzer) -> None:
        """
        Display nucleotide composition analysis.
        
        Args:
            analyzer (GenBankAnalyzer): The analyzer instance
        """
        composition = analyzer.get_nucleotide_composition()
        
        if composition:
            print(f"\nNucleotide Composition:")
            for nucleotide, count in composition.items():
                percentage = (count / analyzer.length) * 100 if analyzer.length > 0 else 0
                print(f"  {nucleotide.upper()}: {count:,} ({percentage:.1f}%)")
            
            gc_content = analyzer.get_gc_content()
            print(f"\nGC Content: {gc_content:.1f}%")
    
    def display_features_analysis(self, analyzer: GenBankAnalyzer) -> None:
        """
        Display features analysis.
        
        Args:
            analyzer (GenBankAnalyzer): The analyzer instance
        """
        features = analyzer.get_features()
        if features:
            print(f"\nFeatures Found ({len(features)} total):")
            
            # Display feature summary
            feature_summary = analyzer.get_feature_summary()
            for feat_type, count in feature_summary.items():
                print(f"  {feat_type}: {count}")
            
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
    
    def display_sequence_preview(self, analyzer: GenBankAnalyzer) -> None:
        """
        Display sequence preview.
        
        Args:
            analyzer (GenBankAnalyzer): The analyzer instance
        """
        if analyzer.length > 0:
            preview = analyzer.get_sequence_preview(60)
            print(f"\nSequence Preview (first {len(preview)} nucleotides):")
            print(f"  {preview}")
            if analyzer.length > len(preview):
                print("  ...")
    
    def run(self) -> None:
        """Run the command-line interface."""
        try:
            self.display_header()
            
            # Get file path
            file_path = self.get_file_path()
            
            # Create analyzer
            self.analyzer = GenBankAnalyzer(file_path)
            
            # Validate file extension (optional warning)
            if not self.analyzer.validate_file_extension():
                print(f"Warning: File '{file_path}' doesn't have a typical GenBank extension (.gb, .gbk, .genbank)")
            
            print(f"Analyzing file: {file_path}")
            print("-" * 50)
            
            # Display analysis results
            self.display_file_info(self.analyzer)
            self.display_sequence_analysis(self.analyzer)
            
            if not self.analyzer.sequence:
                print("Warning: No sequence data found in the file.")
                return
            
            self.display_nucleotide_composition(self.analyzer)
            self.display_features_analysis(self.analyzer)
            self.display_sequence_preview(self.analyzer)
            
        except FileNotFoundError as e:
            print(f"Error: {e}")
        except Exception as e:
            print(f"Error: {e}")
        except KeyboardInterrupt:
            print("\nProgram interrupted by user.")


def main() -> None:
    """Main entry point for the application."""
    cli = GenBankAnalyzerCLI()
    cli.run()


if __name__ == "__main__":
    main()