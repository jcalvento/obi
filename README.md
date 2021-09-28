<p align="center">
  <img src="img/logo.png" alt="Logo">
</p>

# Obi
![Status](https://github.com/jcalvento/tesina/workflows/Obi%201/badge.svg)
[![Anaconda-Server Badge](https://anaconda.org/jcalvento/obi/badges/installer/conda.svg)](https://conda.anaconda.org/jcalvento)
[![Anaconda-Server Badge](https://anaconda.org/jcalvento/obi/badges/latest_release_date.svg)](https://conda.anaconda.org/jcalvento)

## Setup
- Install [Conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/) or [Miniconda](https://docs.conda.io/en/latest/miniconda.html) (if you want a lightweight version)
- Create conda env with Python 3.7.* using requirements from environment.yml file. i.e. `conda env create --file environment.yml python=3.7`
- Download [Swissprot DB](https://ftp.ncbi.nlm.nih.gov/blast/db/swissprot.tar.gz) and [Uniprot <-> PDB mappings](http://ftp.ebi.ac.uk/pub/databases/msd/sifts/flatfiles/csv/pdb_chain_uniprot.csv.gz) from NCBI FTP server. You can use fetch_db.py script to do so (`python src/scripts/fetch_db.py`) or simply run `obi fetch-db` if you installed it as a Conda package. 

## Usage
For each script you can get arguments description and usage by using `--help` argument. i.e. `obi --help`

### Analysis entrypoint
To run the whole analysis or first part depending on the arguments, you can either use `python obi/scripts/main.py` or `obi` if you installed it program as a Conda package.

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
### Only analysis
To run only the analysis (if previous step didn't include it), `python obi/scripts/analysis.py` or `obi analysis` (if installed as Conda package).

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
### Resume remote analysis
Run `python obi/scripts/resume_analysis.py` or `obi resume-analysis` (if installed as Conda package).

```commandline
usage: resume_analysis.py [-h] [--input-path INPUT_PATH]

optional arguments:
  -h, --help            show this help message and exit
  --input-path INPUT_PATH
                        Path of alignment preparation result
```

## Conda Build
- Run `conda build conda.recipe -c bioconda -c anaconda -c conda-forge -c defaults`
- Convert build to different for distributions (OSx, Linux & Windows, remove current one): `conda convert path/to/obi/build/obi-x.x.x-x.tar.bz2 -p linux-64 -p win-64 -p osx-64 -p ...`
- Upload: `anaconda upload /path/to/build.tar.bz2`

## Conda Install
Conda package can be found in https://anaconda.org/jcalvento/obi. To install, run:
- `conda install -c jcalvento obi --channel bioconda`