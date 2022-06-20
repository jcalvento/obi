import json
import os
import subprocess

import pytest
from pytest import fail
from pytest_mock import mock

from obi.src.hyphy import Hyphy, InvalidApiKeyError, HyphyAPIError, HyphyJobNotReady
from tests.utils import get_resource, results_dir, results_path


class TestHyphy:
    @mock.patch('os.popen')
    @mock.patch('subprocess.call')
    def test_it_generates_iq_tree_and_then_runs_hyphy_locally_using_it_along_with_the_alignment(self, mocked_popen,
                                                                                                mocked_subprocess):
        nucleotide_alignment_path = get_resource("nucleotide_alignment.fasta")
        bootstrap = 2000

        hyphy_result = Hyphy.local().run(nucleotide_alignment_path, bootstrap=bootstrap)

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
        os.popen.assert_called_once_with(f'iqtree -s {nucleotide_alignment_path} -bb {bootstrap}')
        subprocess.call.assert_called_once_with(
            ['hyphy', 'meme', '--alignment', nucleotide_alignment_path, '-bb', f'{nucleotide_alignment_path}.treefile'])
        assert expected_result == hyphy_result

    @mock.patch('os.popen')
    @mock.patch('subprocess.call')
    def test_it_generates_iq_tree_and_then_runs_hyphy_locally_with_the_1000_as_the_default_bootstrap_value(self,
                                                                                                          mocked_popen,
                                                                                                          mocked_subprocess):
        nucleotide_alignment_path = get_resource("nucleotide_alignment.fasta")

        Hyphy.local().run(nucleotide_alignment_path)

        os.popen.assert_called_once_with(f'iqtree -s {nucleotide_alignment_path} -bb 1000')
        subprocess.call.assert_called_once_with(
            ['hyphy', 'meme', '--alignment', nucleotide_alignment_path, '-bb', f'{nucleotide_alignment_path}.treefile'])

    @pytest.fixture(scope="module")
    def vcr_config(self):
        return {"record_mode": "once"}

    @pytest.mark.vcr()
    @pytest.mark.default_cassette("submit_job.yaml")
    def test_it_runs_hyphy_generating_a_new_job_in_datamonkey_api(self):
        nucleotide_alignment_path = get_resource("nucleotide_alignment.fasta")

        result = Hyphy.remote().run(nucleotide_alignment_path, "juliancalvento@gmail.com")

        assert result == {
            'id': '62b000b1925b0370fdb762ad',
            'status': 'queue',
            'mail': 'mail@mail.com',
            'created': '2022-06-20T05:08:01.030Z',
        }

    def test_it_fails_when_email_is_not_present_for_remote_job(self):
        try:
            Hyphy.remote().run("some/path.fasta", None)
            fail()
        except HyphyAPIError as e:
            assert e.message == "Email is required to submit a remote job"

    @pytest.mark.vcr()
    @pytest.mark.default_cassette("invalid_job.yaml")
    def test_it_runs_hyphy_generating_a_new_job_in_datamonkey_api(self):
        try:
            Hyphy.remote().run(get_resource("alignment_preparation_result.json"), "mail@mail.com")
            fail()
        except HyphyAPIError as e:
            assert e.message == 'There was an error submitting the job. Error: An unexpected error occured when ' \
                                'parsing the sequence alignment!'

    @pytest.mark.vcr()
    @pytest.mark.default_cassette("job_done.yaml")
    def test_it_returns_job_results_once_done(self):
        hyphy_remote = Hyphy.remote()

        result = hyphy_remote.job_result(results_dir(), "6028aa9c328ea41e07fe788d")

        with open(results_path("nucleotide_alignment.fasta.MEME.json"), "r") as f:
            job_data = json.load(f)
        assert result == hyphy_remote._parsed_response(job_data)

    @pytest.mark.vcr()
    @pytest.mark.default_cassette("job_queue.yaml")
    def test_it_fails_when_job_is_still_queued(self):
        job_id = "6028b28f328ea41e07fe78d9"
        try:
            Hyphy.remote().job_result(results_dir(), job_id)
            fail()
        except HyphyJobNotReady as e:
            assert e.message == f"Job {job_id} is not ready, status: queue. Try again later."
