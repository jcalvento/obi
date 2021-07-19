# Obi 1

![Status](https://github.com/jcalvento/tesina/workflows/Obi%201/badge.svg)
[![Anaconda-Server Badge](https://anaconda.org/jcalvento/obi/badges/installer/conda.svg)](https://conda.anaconda.org/jcalvento)

## Setup
- Install [Conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/) or [Miniconda](https://docs.conda.io/en/latest/miniconda.html) (if you want a lightweight version)
- Create conda env with Python 3.7.* using requirements from environment.yml file. i.e. `conda env create --file environment.yml python=3.7`
- Download Swissprot and Uniprot-PDB mapping Databases. You can use fetch_db.py script to do so. `python src/scripts/fetch_db.py`

## Usage
- Run main script with `python src/scripts/main.py`. You can check usage with `python src/scripts/main.py --help`
```commandline
usage: main.py [-h] [--fasta FASTA] [--email EMAIL] [--output-path OUTPUT_PATH] [--blast BLAST] [--min-identity MIN_IDENTITY] [--max-evalue MAX_EVALUE] [--min-coverage MIN_COVERAGE] [--max-gaps MAX_GAPS] [--include-analysis]
               [--mode MODE]

optional arguments:
  -h, --help            show this help message and exit
  --fasta FASTA         Path of the Fasta file containing the problem sequence
  --email EMAIL         Email to register while using Entrez
  --output-path OUTPUT_PATH
                        Directory where the program creates output files. Default: ./results/<fasta-filename>
  --blast BLAST         Path to Blast Swissprot DB. Default: ./swissprot/swissprot
  --min-identity MIN_IDENTITY
                        Minimum identity percentage for Blast query expressed as float. Default: 0.4
  --max-evalue MAX_EVALUE
                        Max evalue permitted for Blast results. Default: 0.005
  --min-coverage MIN_COVERAGE
                        Min coverage permitted for Blast results. Default: 90
  --max-gaps MAX_GAPS   Maximum number of gaps permitted for Blast results. Default: 6
  --include-analysis    If present, runs positive selection analysis, which includes running Hyphy and getting final results. If not present, result includes until nucleotide alignment
  --mode MODE           Only used when --include-analysis is present. "local" if you want to run Hyphy locally, "remote" to use Datamonkey instead. Default: local
```
- To run only analysis (if previous step didn't include it), `python src/scripts/analysis.py`.  You can check usage with `python src/scripts/analysis.py --help`
```commandline
usage: analysis.py [-h] [--input-path INPUT_PATH] [--mode MODE] [--api-key API_KEY] [--email EMAIL]

optional arguments:
  -h, --help            show this help message and exit
  --input-path INPUT_PATH
                        Path of alignment preparation result
  --mode MODE           "local" if you want to run Hyphy locally, "remote" to use Datamonkey instead. Default: local
  --api-key API_KEY     Required when running in remote mode. To get one go to http://datamonkey.org/apiKey
  --email EMAIL         Required when running in remote mode. You will get notified once job is done
```
- To resume remote analysis `python src/scripts/resume_analysis.py`. You can check usage with `python src/scripts/resume_analysis.py --help`
```commandline
usage: resume_analysis.py [-h] [--input-path INPUT_PATH]

optional arguments:
  -h, --help            show this help message and exit
  --input-path INPUT_PATH
                        Path of alignment preparation result
```
