from obi.alignment_preparation import AlignmentPreparation
from obi.positive_selection import PositiveSelectionAnalyzer


class Obi:
    def __init__(self, blast, email, uniprot_pdb_csv_path, results_dir):
        self._results_dir = results_dir
        self._email = email
        self._blast = blast
        self._uniprot_pdb_csv_path = uniprot_pdb_csv_path

    def run(self, input_fasta_file, include_analysis=True, min_identity=None, max_evalue=None, min_coverage=None,
            max_gaps=None, mode=PositiveSelectionAnalyzer.LOCAL_MODE):
        fasta_file = self._blast.run(
            input_fasta_file, self._results_dir, min_identity=min_identity,
            max_evalue=max_evalue, min_coverage=min_coverage, max_gaps=max_gaps
        )
        alignment_preparation_result = AlignmentPreparation(
            fasta_file, self._results_dir, self._email, self._uniprot_pdb_csv_path
        ).run()

        if include_analysis:
            return PositiveSelectionAnalyzer.for_mode(mode).analyse(
                self._results_dir, alignment_preparation_result
            )
        print("Positive selection analysis skipped")
