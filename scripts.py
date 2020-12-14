import datetime
import os

from obi.blast import BlastResultsError, Blast
from obi.entrez import InvalidEntrezIds
from obi.obi_1 import Obi
from obi.positive_selection import HyphyError
from obi.utils import create_results_dir, root_path

if __name__ == '__main__':
    init_time = datetime.datetime.now()
    print(f"Whole process init time: {init_time}")
    email = 'juliancalvento@gmail.com'
    root_path = root_path()
    fastas_dir = f"{root_path}/fasta"
    files = os.listdir(fastas_dir)
    failed_count = 0
    blast = Blast(f'{root_path}/swissprot/swissprot')

    for file in ['P00766.fasta']:
        print(f"Init time: {datetime.datetime.now()}")
        results_dir = create_results_dir(file)
        obi_1 = Obi(
            blast=blast, email=email, uniprot_pdb_csv_path=f"{root_path}/pdb_chain_uniprot.csv", results_dir=results_dir
        )
        try:
            obi_1.run(f"{fastas_dir}/{file}")
        except (BlastResultsError, HyphyError, InvalidEntrezIds) as e:
            print(e.message)
            failed_count += 1
            continue
        print(f"Finish time: {datetime.datetime.now()}")
    finish_time = datetime.datetime.now()
    print(f"Whole process finish time: {finish_time}, took {finish_time - init_time}, failed: {failed_count}")
