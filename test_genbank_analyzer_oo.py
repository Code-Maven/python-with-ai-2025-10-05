import pytest
from unittest.mock import patch, MagicMock
from pathlib import Path
import tempfile
import os
from genbank_analyzer import GenBankAnalyzer, GenBankAnalyzerCLI


class TestGenBankAnalyzer:
    """Test cases for the GenBankAnalyzer class."""
    
    def test_init_file_not_found(self):
        """Test initialization with non-existent file."""
        with pytest.raises(FileNotFoundError):
            GenBankAnalyzer("nonexistent_file.gb")
    
    def test_init_valid_file(self):
        """Test initialization with valid file."""
        with tempfile.NamedTemporaryFile(suffix=".gb", delete=False) as tmp:
            tmp.write(b"test content")
            tmp_path = tmp.name
        
        try:
            analyzer = GenBankAnalyzer(tmp_path)
            assert analyzer.file_path == Path(tmp_path)
        finally:
            os.unlink(tmp_path)
    
    @patch('genbank_analyzer.SeqIO.read')
    def test_record_property(self, mock_seqio_read):
        """Test the record property with lazy loading."""
        with tempfile.NamedTemporaryFile(suffix=".gb", delete=False) as tmp:
            tmp.write(b"test content")
            tmp_path = tmp.name
        
        try:
            mock_record = MagicMock()
            mock_seqio_read.return_value = mock_record
            
            analyzer = GenBankAnalyzer(tmp_path)
            
            # First access should call SeqIO.read
            record1 = analyzer.record
            assert record1 == mock_record
            mock_seqio_read.assert_called_once()
            
            # Second access should use cached value
            record2 = analyzer.record
            assert record2 == mock_record
            assert mock_seqio_read.call_count == 1  # Still only called once
        finally:
            os.unlink(tmp_path)
    
    @patch('genbank_analyzer.SeqIO.read')
    def test_sequence_property(self, mock_seqio_read):
        """Test the sequence property."""
        with tempfile.NamedTemporaryFile(suffix=".gb", delete=False) as tmp:
            tmp.write(b"test content")
            tmp_path = tmp.name
        
        try:
            mock_record = MagicMock()
            mock_record.seq = "ATCGATCG"
            mock_seqio_read.return_value = mock_record
            
            analyzer = GenBankAnalyzer(tmp_path)
            assert analyzer.sequence == "atcgatcg"  # Should be lowercase
        finally:
            os.unlink(tmp_path)
    
    @patch('genbank_analyzer.SeqIO.read')
    def test_length_property(self, mock_seqio_read):
        """Test the length property."""
        with tempfile.NamedTemporaryFile(suffix=".gb", delete=False) as tmp:
            tmp.write(b"test content")
            tmp_path = tmp.name
        
        try:
            mock_record = MagicMock()
            mock_record.seq = "ATCGATCG"  # 8 nucleotides
            mock_seqio_read.return_value = mock_record
            
            analyzer = GenBankAnalyzer(tmp_path)
            assert analyzer.length == 8
        finally:
            os.unlink(tmp_path)
    
    @patch('genbank_analyzer.SeqIO.read')
    def test_get_basic_info(self, mock_seqio_read):
        """Test basic info extraction."""
        with tempfile.NamedTemporaryFile(suffix=".gb", delete=False) as tmp:
            tmp.write(b"test content")
            tmp_path = tmp.name
        
        try:
            mock_record = MagicMock()
            mock_record.name = "TEST123"
            mock_record.id = "TEST123.1"
            mock_record.description = "Test sequence"
            mock_record.seq = "ATCGATCG"
            mock_record.annotations = {
                'source': 'Test organism',
                'date': '01-JAN-2023',
                'molecule_type': 'DNA'
            }
            mock_seqio_read.return_value = mock_record
            
            analyzer = GenBankAnalyzer(tmp_path)
            info = analyzer.get_basic_info()
            
            assert info['locus'] == "TEST123"
            assert info['accession'] == "TEST123.1"
            assert info['definition'] == "Test sequence"
            assert info['declared_length'] == 8
            assert info['source'] == "Test organism"
            assert info['date'] == "01-JAN-2023"
            assert info['molecule_type'] == "DNA"
        finally:
            os.unlink(tmp_path)
    
    @patch('genbank_analyzer.SeqIO.read')
    def test_get_features(self, mock_seqio_read):
        """Test feature extraction."""
        with tempfile.NamedTemporaryFile(suffix=".gb", delete=False) as tmp:
            tmp.write(b"test content")
            tmp_path = tmp.name
        
        try:
            # Mock feature
            mock_feature = MagicMock()
            mock_feature.type = "CDS"
            mock_feature.location = "[100:200](+)"
            mock_feature.location.strand = 1
            mock_feature.qualifiers = {
                'gene': ['test_gene'],
                'product': ['test_protein']
            }
            
            mock_record = MagicMock()
            mock_record.features = [mock_feature]
            mock_seqio_read.return_value = mock_record
            
            analyzer = GenBankAnalyzer(tmp_path)
            features = analyzer.get_features()
            
            assert len(features) == 1
            assert features[0]['type'] == "CDS"
            assert features[0]['gene'] == "test_gene"
            assert features[0]['product'] == "test_protein"
        finally:
            os.unlink(tmp_path)
    
    @patch('genbank_analyzer.SeqIO.read')
    def test_get_nucleotide_composition(self, mock_seqio_read):
        """Test nucleotide composition calculation."""
        with tempfile.NamedTemporaryFile(suffix=".gb", delete=False) as tmp:
            tmp.write(b"test content")
            tmp_path = tmp.name
        
        try:
            mock_record = MagicMock()
            mock_record.seq = "ATCGATCG"  # 2A, 2T, 2G, 2C
            mock_seqio_read.return_value = mock_record
            
            analyzer = GenBankAnalyzer(tmp_path)
            composition = analyzer.get_nucleotide_composition()
            
            assert composition['a'] == 2
            assert composition['t'] == 2
            assert composition['g'] == 2
            assert composition['c'] == 2
        finally:
            os.unlink(tmp_path)
    
    @patch('genbank_analyzer.SeqIO.read')
    def test_get_gc_content(self, mock_seqio_read):
        """Test GC content calculation."""
        with tempfile.NamedTemporaryFile(suffix=".gb", delete=False) as tmp:
            tmp.write(b"test content")
            tmp_path = tmp.name
        
        try:
            mock_record = MagicMock()
            mock_record.seq = "GCGC"  # 100% GC
            mock_seqio_read.return_value = mock_record
            
            analyzer = GenBankAnalyzer(tmp_path)
            gc_content = analyzer.get_gc_content()
            
            assert gc_content == 100.0
        finally:
            os.unlink(tmp_path)
    
    def test_validate_file_extension(self):
        """Test file extension validation."""
        with tempfile.NamedTemporaryFile(suffix=".gb", delete=False) as tmp:
            tmp.write(b"test content")
            tmp_path = tmp.name
        
        try:
            analyzer = GenBankAnalyzer(tmp_path)
            assert analyzer.validate_file_extension() is True
        finally:
            os.unlink(tmp_path)
        
        with tempfile.NamedTemporaryFile(suffix=".txt", delete=False) as tmp:
            tmp.write(b"test content")
            tmp_path = tmp.name
        
        try:
            analyzer = GenBankAnalyzer(tmp_path)
            assert analyzer.validate_file_extension() is False
        finally:
            os.unlink(tmp_path)
    
    @patch('genbank_analyzer.SeqIO.read')
    def test_validate_sequence_length(self, mock_seqio_read):
        """Test sequence length validation."""
        with tempfile.NamedTemporaryFile(suffix=".gb", delete=False) as tmp:
            tmp.write(b"test content")
            tmp_path = tmp.name
        
        try:
            mock_record = MagicMock()
            mock_record.seq = "ATCG"  # 4 nucleotides
            mock_record.name = "TEST"
            mock_record.id = "TEST"
            mock_record.description = "Test"
            mock_record.annotations = {}
            mock_seqio_read.return_value = mock_record
            
            analyzer = GenBankAnalyzer(tmp_path)
            matches, declared = analyzer.validate_sequence_length()
            
            # Should match because declared_length will be 4 (len(seq))
            assert matches is True
            assert declared == 4
        finally:
            os.unlink(tmp_path)
    
    @patch('genbank_analyzer.SeqIO.read')
    def test_get_feature_summary(self, mock_seqio_read):
        """Test feature summary generation."""
        with tempfile.NamedTemporaryFile(suffix=".gb", delete=False) as tmp:
            tmp.write(b"test content")
            tmp_path = tmp.name
        
        try:
            # Mock multiple features
            mock_feature1 = MagicMock()
            mock_feature1.type = "CDS"
            mock_feature1.location = "[100:200](+)"
            mock_feature1.qualifiers = {}
            
            mock_feature2 = MagicMock()
            mock_feature2.type = "CDS"
            mock_feature2.location = "[300:400](+)"
            mock_feature2.qualifiers = {}
            
            mock_feature3 = MagicMock()
            mock_feature3.type = "gene"
            mock_feature3.location = "[100:400](+)"
            mock_feature3.qualifiers = {}
            
            mock_record = MagicMock()
            mock_record.features = [mock_feature1, mock_feature2, mock_feature3]
            mock_seqio_read.return_value = mock_record
            
            analyzer = GenBankAnalyzer(tmp_path)
            summary = analyzer.get_feature_summary()
            
            assert summary['CDS'] == 2
            assert summary['gene'] == 1
        finally:
            os.unlink(tmp_path)
    
    @patch('genbank_analyzer.SeqIO.read')
    def test_get_sequence_preview(self, mock_seqio_read):
        """Test sequence preview generation."""
        with tempfile.NamedTemporaryFile(suffix=".gb", delete=False) as tmp:
            tmp.write(b"test content")
            tmp_path = tmp.name
        
        try:
            mock_record = MagicMock()
            mock_record.seq = "atcgatcgatcg"
            mock_seqio_read.return_value = mock_record
            
            analyzer = GenBankAnalyzer(tmp_path)
            
            # Test with default length
            preview = analyzer.get_sequence_preview()
            assert preview == "ATCGATCGATCG"  # Should be uppercase
            
            # Test with custom length
            preview_short = analyzer.get_sequence_preview(5)
            assert preview_short == "ATCGA"
        finally:
            os.unlink(tmp_path)


class TestGenBankAnalyzerCLI:
    """Test cases for the GenBankAnalyzerCLI class."""
    
    def test_init(self):
        """Test CLI initialization."""
        cli = GenBankAnalyzerCLI()
        assert cli.analyzer is None
    
    @patch('sys.argv', ['genbank_analyzer.py', 'test.gb'])
    def test_get_file_path_from_argv(self):
        """Test getting file path from command line arguments."""
        cli = GenBankAnalyzerCLI()
        file_path = cli.get_file_path()
        assert file_path == 'test.gb'
    
    @patch('sys.argv', ['genbank_analyzer.py'])
    @patch('builtins.input', return_value='input_file.gb')
    def test_get_file_path_from_input(self, mock_input):
        """Test getting file path from user input."""
        cli = GenBankAnalyzerCLI()
        file_path = cli.get_file_path()
        assert file_path == 'input_file.gb'
        mock_input.assert_called_once()
    
    def test_display_header(self, capsys):
        """Test header display."""
        cli = GenBankAnalyzerCLI()
        cli.display_header()
        
        captured = capsys.readouterr()
        assert "GenBank Sequence Length Analyzer" in captured.out
        assert "analyzes the sequence length" in captured.out
    
    @patch('sys.argv', ['genbank_analyzer.py'])
    @patch('builtins.input', return_value='nonexistent.gb')
    def test_run_file_not_found(self, mock_input, capsys):
        """Test CLI run with non-existent file."""
        cli = GenBankAnalyzerCLI()
        cli.run()
        
        captured = capsys.readouterr()
        assert "Error:" in captured.out


if __name__ == "__main__":
    pytest.main([__file__, "-v"])