import json
import os
import subprocess

import pytest
from pytest import fail
from pytest_mock import mock

from src.obi.hyphy import Hyphy, InvalidApiKeyError, HyphyAPIError, HyphyJobNotReady
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

        result = Hyphy.remote().run(nucleotide_alignment_path, "your-api-key", "juliancalvento@gmail.com")

        assert result == {
            'id': '6028aa9c328ea41e07fe788d',
            'status': 'queue',
            'time_stamp': '2021-02-14T04:44:12.117Z',
            'url': 'datamonkey.org/MEME/6028aa9c328ea41e07fe788d'
        }

    @mock.patch('src.obi.hyphy.RemoteHyphy._key_info')
    def test_it_raises_an_error_when_api_key_is_expired(self, mock_key_info):
        mock_key_info.return_value = {
            "_id": "Your API Key",
            "__v": 3,
            "associated_job_ids": [
                "datamonkey.org/FEL/5f0f0d95ba5dd5a981b52fbb",
                "datamonkey.org/FEL/5f0f0da7ba5dd5a981b52fd3",
                "datamonkey.org/FEL/5f0f0da8ba5dd5a981b52feb"
            ],
            "job_remaining": 0,
            "job_request_made": 100,
            "job_request_limit": 100,
            "expires": "2020-08-04T14:03:26.431Z",
            "created": "2020-07-15T14:04:16.881Z"
        }
        api_key = "6020736d27e5ee138cbb0123"

        try:
            Hyphy.remote().run("something.fasta", api_key, "email@mail.com")
            fail()
        except InvalidApiKeyError as e:
            assert e.message == f"Key {api_key} expired, to get a new one http://datamonkey.org/apiKey"

    @pytest.mark.vcr()
    @pytest.mark.default_cassette("invalid_api_key.yaml")
    def test_it_runs_hyphy_generating_a_new_job_in_datamonkey_api(self):
        api_key = "somethinginvalid"

        try:
            Hyphy.remote().run("some/path.fasta", api_key, "mail@mail.com")
            fail()
        except HyphyAPIError as e:
            assert e.message.startswith('Datamonkey API failed:')

    def test_it_fails_when_api_key_is_not_present_for_remote_job(self):
        try:
            Hyphy.remote().run("some/path.fasta", None, "mail@mail.com")
            fail()
        except InvalidApiKeyError as e:
            assert e.message == "API Key is required to submit a remote job"

    def test_it_fails_when_email_is_not_present_for_remote_job(self):
        try:
            Hyphy.remote().run("some/path.fasta", "thekey", None)
            fail()
        except HyphyAPIError as e:
            assert e.message == "Email is required to submit a remote job"

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
