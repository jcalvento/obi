import datetime
import os

from obi import utils
from obi.blast import BlastResultsError, Blast
from obi.entrez import InvalidEntrezIds
from obi.obi_1 import Obi
from obi.utils import create_results_dir

if __name__ == '__main__':
    print(f"Whole process init time: {datetime.datetime.now()}")
    email = 'juliancalvento@gmail.com'
    root_path = utils.root_path()
    fastas_dir = f"{root_path}/fasta"
    files = os.listdir(fastas_dir)
    failed_count = 0
    blast = Blast(f'{root_path}/swissprot/swissprot')

    for file in ['P02769.fasta']:
        print(f"Init time: {datetime.datetime.now()}")
        results_dir = create_results_dir(file)
        obi_1 = Obi(
            blast=blast, email=email, uniprot_pdb_csv_path=f"{root_path}/pdb_chain_uniprot.csv", results_dir=results_dir
        )
        try:
            obi_1.run(f"{fastas_dir}/{file}")
        except BlastResultsError as e:
            print(e.message)
            failed_count += 1
            continue
        except InvalidEntrezIds:
            failed_count += 1
        print(f"Finish time: {datetime.datetime.now()}")
    print(f"Whole process finish time: {datetime.datetime.now()}, failed: {failed_count}")
