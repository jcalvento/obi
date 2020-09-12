import datetime
import os

from obi import utils
from obi.alignment_preparation import AlignmentPreparation
from obi.blast import BlastResultsError, Blast
from obi.entrez import InvalidEntrezIds
from obi.hyphy import Hyphy
from obi.sifts import Sifts
from obi.utils import detect


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
            # None
        except BlastResultsError as e:
            print(e.message)
            failed_count += 1
            continue
        try:
            alignment_preparation = AlignmentPreparation(fasta_file, results_dir, email, f"{root_path}/pdb_chain_uniprot.csv")
            nucleotide_alignment = alignment_preparation.run()
            hyphy_result = Hyphy(nucleotide_alignment).run(1000)
            print(alignment_preparation.pdb_mapping)
            sifts = Sifts()
            pdb_mappings = list(map(lambda mapping: sifts.map_to(mapping['pdb']), alignment_preparation.pdb_mapping))
            positive_selection_rows = []
            for index, row in enumerate(hyphy_result):
                if row['p-value'] <= 0.1:
                    row['index'] = index
                    positive_selection_rows.append(row)
            result = []
            for uniprot_id, alignment in alignment_preparation.nucleotide_alignments.items():
                acc_number = detect(
                    lambda mapping: uniprot_id.startswith(mapping.from_id),
                    alignment_preparation.uniprot_entrez_mapping
                ).to_id
                rows = []
                codons = [alignment[index:index + 3] for index in range(0, len(alignment), 3)]
                for selection_row in positive_selection_rows:
                    row = {
                        "uniprot_id": uniprot_id,
                        'acc_number': acc_number,
                        'index': selection_row['index'],
                        'p-value': selection_row['p-value'],
                        'codon': codons[selection_row['index']],  # Creo que no se necesita del alineamiento si no de la cadena
                        'aa': alignment_preparation.alignments[uniprot_id][selection_row['index']]
                    }
                    rows.append(row)
                result.append({uniprot_id: rows})
            print(result)

        except InvalidEntrezIds:
            failed_count += 1
        print(f"Finish time: {datetime.datetime.now()}")
    print(f"Whole process finish time: {datetime.datetime.now()}, failed: {failed_count}")
