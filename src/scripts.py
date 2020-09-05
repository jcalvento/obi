import datetime
import os

from src import utils
from src.alignment_preparation import AlignmentPreparation
from src.blast import BlastResultsError, Blast
from src.entrez import InvalidEntrezIds
from src.hyphy import Hyphy


def create_results_dir(file):
    filename = file.split('/')[-1].split('.')[0]
    dir_path = f'{root_path}/results/{filename}'
    # if os.path.isdir(dir_path):
    #     shutil.rmtree(dir_path)
    # os.mkdir(dir_path)
    return dir_path


if __name__ == '__main__':
    print(f"Whole process init time: {datetime.datetime.now()}")
    email = 'juliancalvento@gmail.com'
    root_path = utils.root_path()
    fastas_dir = f"{root_path}/fasta"
    files = os.listdir(fastas_dir)
    failed_count = 0
    blast = Blast(f'{root_path}/swissprot/swissprot')

    for file in [files[0]]:
        print(f"Init time: {datetime.datetime.now()}")
        results_dir = create_results_dir(file)
        try:
            fasta_file = blast.run(f"{fastas_dir}/{file}", results_dir)
            # None
        except BlastResultsError as e:
            print(e.message)
            failed_count += 1
            continue
        try:
            alignment_preparation = AlignmentPreparation(fasta_file, results_dir, email, f"{root_path}/pdb_chain_uniprot.csv")
            nucleotide_alignment = alignment_preparation.run()
            hyphy_result = Hyphy(nucleotide_alignment).run(1000)
            hyphy_result
            print(alignment_preparation.pdb_mapping)
        except InvalidEntrezIds:
            failed_count += 1
        print(f"Finish time: {datetime.datetime.now()}")
    print(f"Whole process finish time: {datetime.datetime.now()}, failed: {failed_count}")
