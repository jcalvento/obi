import os
import pathlib
import shutil
import tarfile
from contextlib import closing
from urllib import request

if __name__ == '__main__':
    swissprot_tar_file = f"{pathlib.Path(__file__).parent.parent.absolute()}/swissprot.tar.gz"
    with closing(request.urlopen('https://ftp.ncbi.nlm.nih.gov/blast/db/swissprot.tar.gz')) as r:
        with open(swissprot_tar_file, 'wb') as f:
            shutil.copyfileobj(r, f)
    tf = tarfile.open(swissprot_tar_file)
    tf.extractall(path=f"{pathlib.Path(__file__).parent.parent.absolute()}/swissprot/")
    os.remove(swissprot_tar_file)
