import pytest
from unittest.mock import patch, mock_open, MagicMock
from pathlib import Path
import tempfile
import os
from ncbi_genbank_downloader import (
    setup_entrez_email,
    search_ncbi_nucleotides,
    download_genbank_file,
    get_sequence_info,
    main
)


def test_setup_entrez_email():
    """Test that Entrez email is set up correctly."""
    with patch('ncbi_genbank_downloader.Entrez') as mock_entrez:
        setup_entrez_email()
        assert hasattr(mock_entrez, 'email')


@patch('ncbi_genbank_downloader.Entrez.esearch')
@patch('ncbi_genbank_downloader.Entrez.read')
def test_search_ncbi_nucleotides_success(mock_read, mock_esearch):
    """Test successful NCBI search."""
    # Mock the search handle and results
    mock_handle = MagicMock()
    mock_esearch.return_value = mock_handle
    mock_read.return_value = {"IdList": ["123456", "789012"]}
    
    result = search_ncbi_nucleotides("test virus", 2)
    
    assert result == ["123456", "789012"]
    mock_esearch.assert_called_once_with(
        db="nucleotide",
        term="test virus",
        retmax=2,
        sort="relevance"
    )
    mock_handle.close.assert_called_once()


@patch('ncbi_genbank_downloader.Entrez.esearch')
@patch('ncbi_genbank_downloader.Entrez.read')
def test_search_ncbi_nucleotides_no_results(mock_read, mock_esearch):
    """Test NCBI search with no results."""
    mock_handle = MagicMock()
    mock_esearch.return_value = mock_handle
    mock_read.return_value = {"IdList": []}
    
    result = search_ncbi_nucleotides("nonexistent virus", 5)
    
    assert result == []


@patch('ncbi_genbank_downloader.Entrez.esearch')
def test_search_ncbi_nucleotides_error(mock_esearch):
    """Test NCBI search with error."""
    mock_esearch.side_effect = Exception("Network error")
    
    result = search_ncbi_nucleotides("test virus", 2)
    
    assert result == []


def test_download_genbank_file_success():
    """Test successful GenBank file download."""
    with tempfile.TemporaryDirectory() as temp_dir:
        # Mock Entrez.efetch
        with patch('ncbi_genbank_downloader.Entrez.efetch') as mock_efetch:
            mock_handle = MagicMock()
            mock_handle.read.return_value = "LOCUS       TEST123                 1000 bp    DNA     linear   VRL 01-JAN-2023\nORIGIN\n        1 atcgatcgat\n//"
            mock_efetch.return_value = mock_handle
            
            result = download_genbank_file("TEST123", temp_dir)
            
            expected_path = Path(temp_dir) / "TEST123.gb"
            assert result == str(expected_path)
            assert expected_path.exists()
            
            # Check file content
            with open(expected_path, 'r') as f:
                content = f.read()
                assert "LOCUS       TEST123" in content


def test_download_genbank_file_no_data():
    """Test GenBank file download with no data."""
    with tempfile.TemporaryDirectory() as temp_dir:
        with patch('ncbi_genbank_downloader.Entrez.efetch') as mock_efetch:
            mock_handle = MagicMock()
            mock_handle.read.return_value = ""  # No data
            mock_efetch.return_value = mock_handle
            
            result = download_genbank_file("TEST123", temp_dir)
            
            assert result is None


def test_download_genbank_file_error():
    """Test GenBank file download with error."""
    with tempfile.TemporaryDirectory() as temp_dir:
        with patch('ncbi_genbank_downloader.Entrez.efetch') as mock_efetch:
            mock_efetch.side_effect = Exception("Network error")
            
            result = download_genbank_file("TEST123", temp_dir)
            
            assert result is None


@patch('ncbi_genbank_downloader.SeqIO.read')
def test_get_sequence_info_success(mock_seqio_read):
    """Test successful sequence info extraction."""
    # Mock SeqRecord
    mock_record = MagicMock()
    mock_record.id = "TEST123"
    mock_record.description = "Test sequence"
    mock_record.seq = "ATCGATCG"  # 8 bp
    mock_record.annotations = {"organism": "Test virus"}
    mock_seqio_read.return_value = mock_record
    
    result = get_sequence_info("test.gb")
    
    expected = {
        'accession': 'TEST123',
        'description': 'Test sequence',
        'length': 8,
        'organism': 'Test virus'
    }
    assert result == expected


@patch('ncbi_genbank_downloader.SeqIO.read')
def test_get_sequence_info_error(mock_seqio_read):
    """Test sequence info extraction with error."""
    mock_seqio_read.side_effect = Exception("Parse error")
    
    result = get_sequence_info("test.gb")
    
    assert 'error' in result
    assert "Parse error" in result['error']


@patch('ncbi_genbank_downloader.setup_entrez_email')
@patch('ncbi_genbank_downloader.search_ncbi_nucleotides')
@patch('ncbi_genbank_downloader.download_genbank_file')
@patch('ncbi_genbank_downloader.get_sequence_info')
@patch('sys.argv', ['ncbi_genbank_downloader.py', 'test virus', '2'])
def test_main_command_line_success(mock_get_info, mock_download, mock_search, mock_setup):
    """Test main function with command line arguments."""
    # Mock successful execution
    mock_search.return_value = ["123456", "789012"]
    mock_download.side_effect = ["file1.gb", "file2.gb"]
    mock_get_info.side_effect = [
        {'accession': '123456', 'description': 'Test 1', 'length': 1000, 'organism': 'Virus 1'},
        {'accession': '789012', 'description': 'Test 2', 'length': 2000, 'organism': 'Virus 2'}
    ]
    
    # Capture output by patching print
    with patch('builtins.print') as mock_print:
        main()
    
    mock_setup.assert_called_once()
    mock_search.assert_called_once_with('test virus', 2)
    assert mock_download.call_count == 2


@patch('ncbi_genbank_downloader.setup_entrez_email')
@patch('ncbi_genbank_downloader.search_ncbi_nucleotides')
@patch('builtins.input')
@patch('sys.argv', ['ncbi_genbank_downloader.py'])
def test_main_interactive_success(mock_input, mock_search, mock_setup):
    """Test main function with interactive input."""
    # Mock user input
    mock_input.side_effect = ['influenza virus', '3']
    mock_search.return_value = []  # No results to keep test simple
    
    with patch('builtins.print'):
        main()
    
    mock_setup.assert_called_once()
    mock_search.assert_called_once_with('influenza virus', 3)


@patch('sys.argv', ['ncbi_genbank_downloader.py', 'test', 'invalid_number'])
def test_main_invalid_number():
    """Test main function with invalid number argument."""
    with patch('builtins.print') as mock_print:
        main()
    
    # Check that error message was printed
    calls = [str(call) for call in mock_print.call_args_list]
    assert any("must be an integer" in call for call in calls)


@patch('builtins.input')
@patch('sys.argv', ['ncbi_genbank_downloader.py'])
def test_main_empty_search_term(mock_input):
    """Test main function with empty search term."""
    mock_input.return_value = ""  # Empty search term
    
    with patch('builtins.print') as mock_print:
        main()
    
    calls = [str(call) for call in mock_print.call_args_list]
    assert any("cannot be empty" in call for call in calls)


@patch('builtins.input')
@patch('sys.argv', ['ncbi_genbank_downloader.py'])
def test_main_invalid_interactive_number(mock_input):
    """Test main function with invalid number in interactive mode."""
    mock_input.side_effect = ['test virus', 'not_a_number']
    
    with patch('builtins.print') as mock_print:
        main()
    
    calls = [str(call) for call in mock_print.call_args_list]
    assert any("valid number" in call for call in calls)


@patch('builtins.input')
@patch('sys.argv', ['ncbi_genbank_downloader.py'])
def test_main_number_out_of_range(mock_input):
    """Test main function with number out of valid range."""
    mock_input.side_effect = ['test virus', '150']  # > 100
    
    with patch('builtins.print') as mock_print:
        main()
    
    calls = [str(call) for call in mock_print.call_args_list]
    assert any("between 1 and 100" in call for call in calls)


@patch('ncbi_genbank_downloader.setup_entrez_email')
@patch('ncbi_genbank_downloader.search_ncbi_nucleotides')
@patch('sys.argv', ['ncbi_genbank_downloader.py', 'test', '2'])
def test_main_keyboard_interrupt(mock_search, mock_setup):
    """Test main function handling keyboard interrupt."""
    mock_search.side_effect = KeyboardInterrupt()
    
    with patch('builtins.print') as mock_print:
        main()
    
    calls = [str(call) for call in mock_print.call_args_list]
    assert any("interrupted by user" in call for call in calls)


if __name__ == "__main__":
    pytest.main([__file__])