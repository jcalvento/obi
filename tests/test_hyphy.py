import os
import subprocess

from pytest_mock import mock

from obi.hyphy import Hyphy
from tests.utils import get_resource


class TestHyphy:

    @mock.patch('os.popen')
    @mock.patch('subprocess.call')
    def test_it_generates_iq_tree_and_then_runs_hyphy_using_it_along_with_the_alignment(self, mocked_popen, mocked_subprocess):
        nucleotide_alignment_path = get_resource("nucleotide_alignment.fasta")
        boostrap = 2000

        hyphy_result = Hyphy(nucleotide_alignment_path).run(boostrap)

        expected_result = [
            {'# branches under selection': 0,
            '&beta;<sup>+</sup>': 0,
            '&beta;<sup>-</sup>': 0,
            'FEL LogL': 0,
            'LRT': 0,
            'MEME LogL': 0,
            'Total branch length': 0,
            'alpha;': 0,
            'p-value': 0.05,
            'p<sup>+</sup>': 0,
            'p<sup>-</sup>': 1},
            {'# branches under selection': 0,
            '&beta;<sup>+</sup>': 0,
            '&beta;<sup>-</sup>': 0,
            'FEL LogL': 0,
            'LRT': 0,
            'MEME LogL': 0,
            'Total branch length': 0,
            'alpha;': 0,
            'p-value': 1,
            'p<sup>+</sup>': 0,
            'p<sup>-</sup>': 1}
        ]
        os.popen.assert_called_once_with(f'iqtree -s {nucleotide_alignment_path} -bb {boostrap}')
        subprocess.call.assert_called_once_with(['hyphy', 'meme', '--alignment', nucleotide_alignment_path, '-bb', f'{nucleotide_alignment_path}.treefile'])
        assert expected_result == hyphy_result

    @mock.patch('os.popen')
    @mock.patch('subprocess.call')
    def test_it_generates_iq_tree_and_then_runs_hyphy_with_the_1000_as_the_default_boostrap_value(self, mocked_popen, mocked_subprocess):
        nucleotide_alignment_path = get_resource("nucleotide_alignment.fasta")

        Hyphy(nucleotide_alignment_path).run()

        os.popen.assert_called_once_with(f'iqtree -s {nucleotide_alignment_path} -bb 1000')
        subprocess.call.assert_called_once_with(['hyphy', 'meme', '--alignment', nucleotide_alignment_path, '-bb', f'{nucleotide_alignment_path}.treefile'])
