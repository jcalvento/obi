from pytest import fixture
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
        except BlastResultsError as e:
            assert e.message == "There is none blast result that meets filtering criteria"

    @mock.patch('obi.blast.Blast._Blast__blastp')
    @mock.patch('obi.blast.Blast._Blast__blastdbcmd')
    def test_when_there_is_a_result_but_is_does_not_have_the_correct_evalue_it_is_not_included_on_the_result(self, mock_blastdbcmd, mock_blastp, blast):
        mock_blastp.return_value = open(get_resource("blast_high_evalue.xml"))
        mock_blastdbcmd.return_value = open(get_resource("blastcmd_result.fasta")).read()

        try:
            blast.run(get_resource("invalid_blast_file.fasta"), results_dir())
        except BlastResultsError as e:
            assert e.message == "There is none blast result that meets filtering criteria"
