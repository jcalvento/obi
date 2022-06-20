import json

import pytest
from _pytest.outcomes import fail

from obi.src.alignment_preparation import AlignmentPreparationResultSchema
from obi.src.hyphy import HyphyJobNotReady
from obi.src.positive_selection import PositiveSelectionAnalyzer
from tests.utils import results_dir, get_resource


class TestRemoteAnalyzer:
    @pytest.mark.vcr()
    @pytest.mark.default_cassette("submit_job.yaml")
    def test_creates_a_new_datamonkey_job(self, mocker):
        alignment_file = open(get_resource("nucleotide_alignment.fasta"), 'rb')
        mocker.patch('obi.src.hyphy.RemoteHyphy._RemoteHyphy__file', return_value=alignment_file)
        results_dir_path = results_dir()
        file = open(get_resource("/alignment_preparation_result.json"), 'r')
        data = json.load(file)
        alignment_preparation_result = AlignmentPreparationResultSchema().load(data)
        email = "mail@mail.com"
        analyser = PositiveSelectionAnalyzer.for_mode(PositiveSelectionAnalyzer.REMOTE_MODE)

        result = analyser.analyse(results_dir_path, alignment_preparation_result, email=email)

        assert result == {
          "id": "62b000b1925b0370fdb762ad",
          "status": "queue",
          'mail': 'mail@mail.com',
          'created': '2022-06-20T05:08:01.030Z',
        }

    @pytest.mark.vcr()
    @pytest.mark.default_cassette("job_queue.yaml")
    def test_raises_an_error_when_the_job_is_still_in_progress(self):
        results_dir_path = results_dir()
        file = open(get_resource("/alignment_preparation_result.json"), 'r')
        data = json.load(file)
        alignment_preparation_result = AlignmentPreparationResultSchema().load(data)
        analyser = PositiveSelectionAnalyzer.for_mode(PositiveSelectionAnalyzer.REMOTE_MODE)
        job_id = "62b000b1925b0370fdb762ad"

        try:
            analyser.resume(results_dir_path, alignment_preparation_result)
            fail()
        except HyphyJobNotReady as e:
            assert e.message == f"Job {job_id} is not ready, status: queue. Try again later."

    @pytest.mark.vcr()
    @pytest.mark.default_cassette("job_done.yaml")
    def test_returns_the_result_when_the_job_is_done(self, mocker):
        results_dir_path = results_dir()
        file = open(get_resource("/alignment_preparation_result.json"), 'r')
        data = json.load(file)
        alignment_preparation_result = AlignmentPreparationResultSchema().load(data)
        analyser = PositiveSelectionAnalyzer.for_mode(PositiveSelectionAnalyzer.REMOTE_MODE)
        mock_report = mocker.patch(
            'obi.src.positive_selection.PositiveSelectionAnalyzer._process_hyphy_result_and_generate_report'
        )

        result = analyser.resume(results_dir_path, alignment_preparation_result)

        mock_report.assert_called_once()
