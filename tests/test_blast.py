import os

from pytest import fixture, fail
from pytest_mock import mock

from obi.blast import Blast, BlastResultsError
from obi.utils import root_path
from tests.utils import results_dir, get_resource


class TestBlast:
    @fixture
    def blast(self):
        return Blast(f'{root_path()}/swissprot/swissprot')

    def test_when_there_are_no_results_for_the_given_id_it_returns_and_error(self, blast):
        try:
            blast.run(get_resource("invalid_blast_file.fasta"), results_dir())
            fail()
        except BlastResultsError as e:
            assert e.message == "There is none blast result that meets filtering criteria"

    @mock.patch('obi.blast.Blast._Blast__blastp')
    def test_when_there_is_a_result_but_its_evalue_is_higher_than_required_it_is_not_included_on_the_result(self, mock_blastp, blast):
        mock_blastp.return_value = open(get_resource("blast_high_evalue.xml"))

        try:
            blast.run(input_fasta=get_resource("valid_blast.fasta"), results_dir=results_dir(), max_evalue=0.4)
            fail()
        except BlastResultsError as e:
            assert e.message == "There is none blast result that meets filtering criteria"

    @mock.patch('obi.blast.Blast._Blast__blastp')
    def test_when_there_is_a_result_but_its_coverage_is_less_than_required_it_is_not_included_on_the_result(self, mock_blastp, blast):
        mock_blastp.return_value = open(get_resource("blast_low_coverage.xml"))

        try:
            blast.run(input_fasta=get_resource("valid_blast.fasta"), results_dir=results_dir(), min_coverage=85)
            fail()
        except BlastResultsError as e:
            assert e.message == "There is none blast result that meets filtering criteria"

    @mock.patch('obi.blast.Blast._Blast__blastp')
    def test_when_there_is_a_result_but_its_identity_is_less_than_required_it_is_not_included_on_the_result(self, mock_blastp, blast):
        mock_blastp.return_value = open(get_resource("blast_low_identity.xml"))

        try:
            blast.run(input_fasta=get_resource("valid_blast.fasta"), results_dir=results_dir(), min_identity=0.3)
            fail()
        except BlastResultsError as e:
            assert e.message == "There is none blast result that meets filtering criteria"

    @mock.patch('obi.blast.Blast._Blast__blastp')
    def test_when_there_is_a_result_but_its_gap_percentage_is_higher_than_required_it_is_not_included_on_the_result(self, mock_blastp, blast):
        mock_blastp.return_value = open(get_resource("blast_with_gaps.xml"))

        try:
            blast.run(input_fasta=get_resource("valid_blast.fasta"), results_dir=results_dir(), max_gaps=5)
            fail()
        except BlastResultsError as e:
            assert e.message == "There is none blast result that meets filtering criteria"

    @mock.patch('obi.blast.Blast._Blast__blastp')
    @mock.patch('obi.blast.Blast._Blast__blastdbcmd')
    def test_when_there_are_valid_results_it_returns_the_path_of_the_resulting_fasta_with_the_sequences_that_meet_the_criteria(self, mock_blastdbcmd, mock_blastp, blast):
        mock_blastp.return_value = open(get_resource("blastp_result.xml"))
        mock_blastdbcmd.return_value = open(get_resource("blastcmd_result.fasta")).read()

        fasta_result_path = blast.run(get_resource("valid_blast.fasta"), results_dir())

        assert f"{results_dir()}/blast_result.fasta" == fasta_result_path
        with open(fasta_result_path, "r") as f:
            result_content = f.read()
            assert result_content == mock_blastdbcmd.return_value

        os.remove(fasta_result_path)
