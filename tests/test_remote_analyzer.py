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
        mocker.patch('obi.src.hyphy.upload', return_value="hyphy.com/file")
        results_dir_path = results_dir()
        file = open(get_resource("/alignment_preparation_result.json"), 'r')
        data = json.load(file)
        alignment_preparation_result = AlignmentPreparationResultSchema().load(data)
        api_key = "6020736d27e5ee138cbb0123"
        email = "mail@mail.com"
        analyser = PositiveSelectionAnalyzer.for_mode(PositiveSelectionAnalyzer.REMOTE_MODE)

        result = analyser.analyse(results_dir_path, alignment_preparation_result, api_key=api_key, email=email)

        assert result == {
          "time_stamp": "2021-02-14T04:44:12.117Z",
          "id": "6028aa9c328ea41e07fe788d",
          "status": "queue",
          "url": "datamonkey.org/MEME/6028aa9c328ea41e07fe788d"
        }

    @pytest.mark.vcr()
    @pytest.mark.default_cassette("job_queue.yaml")
    def test_raises_an_error_when_the_job_is_still_in_progress(self):
        results_dir_path = results_dir()
        file = open(get_resource("/alignment_preparation_result.json"), 'r')
        data = json.load(file)
        alignment_preparation_result = AlignmentPreparationResultSchema().load(data)
        analyser = PositiveSelectionAnalyzer.for_mode(PositiveSelectionAnalyzer.REMOTE_MODE)
        job_id = "6028aa9c328ea41e07fe788d"

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
