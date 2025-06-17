# Peptide Finder Tool

This tool identifies peptide sequences between `AQ` and `AQAQ` motifs in DNA sequences stored in FASTQ files. It translates each sequence into all six reading frames, highlights matched peptide segments, and reports their frequencies.

## Features

- Parses uploaded FASTQ files
- Translates DNA into six reading frames
- Detects peptides between `AQ...AQAQ` motifs
- Highlights identified peptides in the protein sequence
- Counts and exports unique peptides to a downloadable CSV
- Displays results per record and per reading frame

## Technologies Used

- Python 3.x
- Flask
- Biopython

## Installation

1. Clone the repository:
   ```bash
   git clone https://github.com/yourusername/peptide-finder-tool.git
   cd peptide-finder-tool
   
2. (Optional) Create a virtual environment:
    ```bash
    python -m venv venv
    source venv/bin/activate  # On Windows: venv\Scripts\activate
    
3. Install dependencies:
   ```bash
   pip install -r requirements.txt

