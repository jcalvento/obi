# Obi 1

## Setup
- Download and install BLAST+ locally from https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download.
Version should be >= 2.10.1
- Install [Conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/) or [Miniconda](https://docs.conda.io/en/latest/miniconda.html) (if you want a light weight version)
- Create conda env with Python 3.7.* using requirements from environment.yml file. i.e. `conda env create --file environment.yml python=3.7`
- Download Swissprot and Uniprot-PDB mapping Databases. You can use fetch_db.py script to do so. `python obi/fetch_db.py`

## Usage
- Run main script with `python main.py`. You can check usage with `python main.py --help`