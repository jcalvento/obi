import datetime
import os

from obi import utils
from obi.alignment_preparation import AlignmentPreparation
from obi.blast import BlastResultsError, Blast
from obi.entrez import InvalidEntrezIds
from obi.hyphy import Hyphy
from obi.positive_selection_report import PositiveSelectionReport
from obi.sifts import Sifts


def create_results_dir(file):
    filename = file.split('/')[-1].split('.')[0]
    dir_path = f'{root_path}/results/{filename}'
    # if os.path.isdir(dir_path):
    #     shutil.rmtree(dir_path)
    if not os.path.isdir(dir_path):
        os.mkdir(dir_path)
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
        except BlastResultsError as e:
            print(e.message)
            failed_count += 1
            continue
        try:
            alignment_preparation_result = AlignmentPreparation(
                fasta_file, results_dir, email, f"{root_path}/pdb_chain_uniprot.csv"
            ).run()
            hyphy_result = Hyphy(alignment_preparation_result.nucleotide_alignment_path).run(1000)
            print(alignment_preparation_result.uniprot_pdb_mapping)
            sifts = Sifts()
            pdb_mappings = list(map(
                lambda mapping: sifts.map_to(mapping['pdb']),
                alignment_preparation_result.uniprot_pdb_mapping
            ))
            report = PositiveSelectionReport(alignment_preparation_result, hyphy_result, pdb_mappings).generate()
            print(report)

        except InvalidEntrezIds:
            failed_count += 1
        print(f"Finish time: {datetime.datetime.now()}")
    print(f"Whole process finish time: {datetime.datetime.now()}, failed: {failed_count}")
