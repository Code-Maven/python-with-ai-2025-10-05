import sys
import os
from pathlib import Path
from typing import List, Optional
from Bio import Entrez
from Bio import SeqIO
import time


def setup_entrez_email() -> None:
    """
    Set up Entrez email for NCBI API access.
    NCBI requires an email address for API access.
    """
    # You should set your email address here for NCBI API access
    # This is required by NCBI's usage guidelines
    Entrez.email = "gabor@szabgba.com"  # Replace with your actual email


def search_ncbi_nucleotides(search_term: str, max_results: int = 10) -> List[str]:
    """
    Search NCBI nucleotide database for sequences matching the search term.
    
    Args:
        search_term (str): Search term to look for in NCBI nucleotide database
        max_results (int): Maximum number of results to return
    
    Returns:
        List[str]: List of GenBank accession IDs
    """
    try:
        print(f"Searching NCBI nucleotide database for: '{search_term}'")
        print(f"Maximum results requested: {max_results}")
        
        # Search the nucleotide database
        search_handle = Entrez.esearch(
            db="nucleotide",
            term=search_term,
            retmax=max_results,
            sort="relevance"
        )
        
        search_results = Entrez.read(search_handle)
        search_handle.close()
        
        id_list = search_results["IdList"]
        
        if not id_list:
            print("No results found for the search term.")
            return []
        
        print(f"Found {len(id_list)} sequences matching your search.")
        return id_list
    
    except Exception as e:
        print(f"Error searching NCBI database: {e}")
        return []


def download_genbank_file(accession_id: str, output_dir: str) -> Optional[str]:
    """
    Download a GenBank file for a given accession ID.
    
    Args:
        accession_id (str): NCBI accession ID
        output_dir (str): Directory to save the downloaded file
    
    Returns:
        Optional[str]: Path to the downloaded file, or None if failed
    """
    try:
        # Fetch the sequence record
        fetch_handle = Entrez.efetch(
            db="nucleotide",
            id=accession_id,
            rettype="gb",
            retmode="text"
        )
        
        # Read the GenBank data
        genbank_data = fetch_handle.read()
        fetch_handle.close()
        
        if not genbank_data:
            print(f"  Warning: No data retrieved for {accession_id}")
            return None
        
        # Create output directory if it doesn't exist
        Path(output_dir).mkdir(parents=True, exist_ok=True)
        
        # Save to file
        output_file = Path(output_dir) / f"{accession_id}.gb"
        with open(output_file, 'w', encoding='utf-8') as f:
            f.write(genbank_data)
        
        print(f"  ✓ Downloaded: {output_file}")
        return str(output_file)
    
    except Exception as e:
        print(f"  ✗ Error downloading {accession_id}: {e}")
        return None


def get_sequence_info(file_path: str) -> dict:
    """
    Extract basic information from a downloaded GenBank file.
    
    Args:
        file_path (str): Path to the GenBank file
    
    Returns:
        dict: Basic sequence information
    """
    try:
        record = SeqIO.read(file_path, "genbank")
        return {
            'accession': record.id,
            'description': record.description,
            'length': len(record.seq),
            'organism': record.annotations.get('organism', 'Unknown')
        }
    except Exception as e:
        return {'error': str(e)}


def main() -> None:
    """
    Main function to handle user input and coordinate the download process.
    """
    try:
        print("NCBI GenBank Downloader")
        print("=" * 30)
        print("This program searches NCBI's nucleotide database and downloads GenBank files.\n")
        
        # Set up Entrez email (required by NCBI)
        setup_entrez_email()
        
        # Get search parameters from user or command line
        if len(sys.argv) >= 3:
            search_term = sys.argv[1]
            try:
                num_results = int(sys.argv[2])
            except ValueError:
                print("Error: Number of results must be an integer.")
                return
        else:
            search_term = input("Enter search term: ").strip()
            if not search_term:
                print("Error: Search term cannot be empty.")
                return
            
            try:
                num_results = int(input("Enter number of results to download (1-100): ").strip())
                if num_results < 1 or num_results > 100:
                    print("Error: Number of results must be between 1 and 100.")
                    return
            except ValueError:
                print("Error: Please enter a valid number.")
                return
        
        print(f"\nSearch term: '{search_term}'")
        print(f"Number of results to download: {num_results}")
        print("-" * 50)
        
        # Search NCBI database
        accession_ids = search_ncbi_nucleotides(search_term, num_results)
        
        if not accession_ids:
            print("No sequences found. Try a different search term.")
            return
        
        # Create output directory
        output_dir = "downloaded_genbank"
        print(f"\nDownloading {len(accession_ids)} GenBank files to '{output_dir}/' directory...")
        print("-" * 50)
        
        downloaded_files = []
        
        # Download each file with a small delay to be respectful to NCBI servers
        for i, accession_id in enumerate(accession_ids, 1):
            print(f"[{i}/{len(accession_ids)}] Downloading {accession_id}...")
            
            file_path = download_genbank_file(accession_id, output_dir)
            if file_path:
                downloaded_files.append(file_path)
                
                # Get basic info about the sequence
                info = get_sequence_info(file_path)
                if 'error' not in info:
                    print(f"      {info['organism']} - {info['length']} bp")
                    print(f"      {info['description'][:80]}{'...' if len(info['description']) > 80 else ''}")
            
            # Small delay between downloads to be respectful to NCBI servers
            if i < len(accession_ids):
                time.sleep(0.5)
        
        # Summary
        print("\n" + "=" * 50)
        print(f"Download Summary:")
        print(f"  Requested: {len(accession_ids)} files")
        print(f"  Downloaded: {len(downloaded_files)} files")
        print(f"  Output directory: {output_dir}/")
        
        if downloaded_files:
            print(f"\nDownloaded files:")
            for file_path in downloaded_files:
                print(f"  - {Path(file_path).name}")
            
            print(f"\nYou can now analyze these files using:")
            print(f"  uv run python genbank_analyzer.py {output_dir}/<filename>.gb")
    
    except KeyboardInterrupt:
        print("\nDownload interrupted by user.")
    except Exception as e:
        print(f"Error: {e}")


if __name__ == "__main__":
    main()