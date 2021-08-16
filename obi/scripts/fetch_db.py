import gzip
import os
import pathlib
import shutil
import tarfile
from contextlib import closing
from urllib import request


def fetch_db():
    print("Fetching Swissprot DB")
    swissprot_tar_file = f"{pathlib.Path(__file__).parent.parent.absolute()}/swissprot.tar.gz"
    with closing(request.urlopen('https://ftp.ncbi.nlm.nih.gov/blast/db/swissprot.tar.gz')) as r:
        with open(swissprot_tar_file, 'wb') as f:
            shutil.copyfileobj(r, f)
    tf = tarfile.open(swissprot_tar_file)
    swissprot_db_dir = f"{pathlib.Path(__file__).parent.parent.absolute()}/swissprot/"
    tf.extractall(path=swissprot_db_dir)
    os.remove(swissprot_tar_file)
    print(f"Swissprot DB downloaded successfully. Path: {swissprot_db_dir}")
    print("Fetching PDB to Uniprot CSV")
    pdb_chain_uniprot_file = f"{pathlib.Path(__file__).parent.parent.absolute()}/pdb_chain_uniprot.csv"
    with closing(request.urlopen(
            'http://ftp.ebi.ac.uk/pub/databases/msd/sifts/flatfiles/csv/pdb_chain_uniprot.csv.gz')) as r:
        with open(pdb_chain_uniprot_file, 'wb') as f:
            with gzip.open(r, 'rb') as gzip_file:
                shutil.copyfileobj(gzip_file, f)
    print(f"PDB to Uniprot CSV downloaded successfully. Path {pdb_chain_uniprot_file}")


if __name__ == '__main__':
    fetch_db()
