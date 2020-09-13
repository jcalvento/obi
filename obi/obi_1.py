from obi.alignment_preparation import AlignmentPreparation
from obi.hyphy import Hyphy
from obi.positive_selection_report import PositiveSelectionReport
from obi.sifts import Sifts


class Obi:
    def __init__(self, blast, email, uniprot_pdb_csv_path, results_dir):
        self._results_dir = results_dir
        self._email = email
        self._blast = blast
        self._uniprot_pdb_csv_path = uniprot_pdb_csv_path

    def run(self, input_fasta_file):
        fasta_file = self._blast.run(input_fasta_file, self._results_dir)
        alignment_preparation_result = AlignmentPreparation(
            fasta_file, self._results_dir, self._email, self._uniprot_pdb_csv_path
        ).run()
        hyphy_result = Hyphy(alignment_preparation_result.nucleotide_alignment_path).run(1000)
        print(alignment_preparation_result.uniprot_pdb_mapping)
        sifts = Sifts()
        pdb_mappings = {}
        for uniprot_id, pdbs in alignment_preparation_result.uniprot_pdb_mapping.items():
            pdbs_data = list(map(
                lambda pdb_id: sifts.map_to(pdb_id),
                pdbs
            ))
            pdb_mappings[uniprot_id] = pdbs_data
        report = PositiveSelectionReport(alignment_preparation_result, hyphy_result, pdb_mappings).generate()
        print(report)
        return report
