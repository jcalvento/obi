from obi.alignment_preparation import AlignmentPreparation
from obi.positive_selection import PositiveSelection


class Obi:
    def __init__(self, blast, email, uniprot_pdb_csv_path, results_dir):
        self._results_dir = results_dir
        self._email = email
        self._blast = blast
        self._uniprot_pdb_csv_path = uniprot_pdb_csv_path

    def run(self, input_fasta_file, include_analysis=True):
        fasta_file = self._blast.run(input_fasta_file, self._results_dir)
        alignment_preparation_result = AlignmentPreparation(
            fasta_file, self._results_dir, self._email, self._uniprot_pdb_csv_path
        ).run()

        if include_analysis:
            return PositiveSelection().analyse(self._results_dir, alignment_preparation_result)
        print("Positive selection analysis skipped")
