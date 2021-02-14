from unittest import mock
from unittest.mock import Mock

from obi.blast import Blast
from obi.obi_1 import Obi
from tests.utils import get_resource, results_dir


class TestObi:
    @mock.patch('obi.blast.Blast.run')
    @mock.patch('obi.alignment_preparation.AlignmentPreparation.run')
    @mock.patch('obi.positive_selection.PositiveSelectionAnalyzer.analyse')
    def test_when_include_analysis_is_false_it_does_not_run_analysis_generation(self, mock_positive_selection,
                                                                                mock_alignment_preparation,
                                                                                mock_blast_run):
        blast_result = get_resource("valid_blast.fasta")
        mock_blast_run.return_value = blast_result
        blast = Blast('mocked')
        input_fasta = get_resource("valid_blast.fasta")
        email = "test@email.com"
        uniprot_pdb_csv_path = "uniprot_csv.csv"

        Obi(blast, email, uniprot_pdb_csv_path, results_dir()).run(input_fasta, include_analysis=False)

        mock_blast_run.assert_called_once_with(input_fasta, results_dir(), min_identity=None,
                                               max_evalue=None, min_coverage=None, max_gaps=None)
        mock_alignment_preparation.assert_called_once()
        mock_positive_selection.assert_not_called()

    @mock.patch('obi.blast.Blast.run')
    @mock.patch('obi.alignment_preparation.AlignmentPreparation.run')
    @mock.patch('obi.positive_selection.PositiveSelectionAnalyzer.analyse')
    def test_when_include_analysis_is_true_it_runs_analysis_generation(self, mock_positive_selection, mock_alignment_preparation, mock_blast_run):
        blast_result = get_resource("valid_blast.fasta")
        mock_blast_run.return_value = blast_result
        alignment_preparation_result = Mock()
        mock_alignment_preparation.return_value = alignment_preparation_result
        blast = Blast('mocked')
        input_fasta = get_resource("valid_blast.fasta")
        email = "test@email.com"
        uniprot_pdb_csv_path = "uniprot_csv.csv"

        Obi(blast, email, uniprot_pdb_csv_path, results_dir()).run(input_fasta, include_analysis=True)

        mock_blast_run.assert_called_once_with(input_fasta, results_dir(), min_identity=None,
                                               max_evalue=None, min_coverage=None, max_gaps=None)
        mock_alignment_preparation.assert_called_once()
        mock_positive_selection.assert_called_once_with(results_dir(), alignment_preparation_result)

    @mock.patch('obi.blast.Blast.run')
    @mock.patch('obi.alignment_preparation.AlignmentPreparation.run')
    @mock.patch('obi.positive_selection.PositiveSelectionAnalyzer.analyse')
    def test_when_blast_params_are_included_it_forwards_them_to_blast(self, mock_positive_selection, mock_alignment_preparation, mock_blast_run):
        blast_result = get_resource("valid_blast.fasta")
        mock_blast_run.return_value = blast_result
        blast = Blast('mocked')
        input_fasta = get_resource("valid_blast.fasta")
        email = "test@email.com"
        uniprot_pdb_csv_path = "uniprot_csv.csv"

        Obi(blast, email, uniprot_pdb_csv_path, results_dir()).run(
            input_fasta, include_analysis=False, min_identity=0.5, max_evalue=0.003, min_coverage=95, max_gaps=20
        )

        mock_blast_run.assert_called_once_with(
            input_fasta, results_dir(), min_identity=0.5, max_evalue=0.003, min_coverage=95, max_gaps=20
        )
        mock_alignment_preparation.assert_called_once()
        mock_positive_selection.assert_not_called()
