import json

from obi.src.alignment_preparation import AlignmentPreparationResultSchema
from obi.src.positive_selection import PositiveSelectionAnalyzer
from tests.utils import results_dir, get_resource


class TestLocalAnalyzer:
    def test_it_runs_positive_selection_analysis_localy(self, mocker):
        hyphy_result = [
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
        ]
        mock_hyphy = mocker.patch('obi.src.hyphy.LocalHyphy.run', return_value=hyphy_result)
        mock_report = mocker.patch(
            'obi.src.positive_selection.PositiveSelectionAnalyzer._process_hyphy_result_and_generate_report'
        )
        analyser = PositiveSelectionAnalyzer.for_mode(PositiveSelectionAnalyzer.LOCAL_MODE)
        results_dir_path = results_dir()
        file = open(get_resource("/alignment_preparation_result.json"), 'r')
        data = json.load(file)
        alignment_preparation_result = AlignmentPreparationResultSchema().load(data)

        analyser.analyse(results_dir_path, alignment_preparation_result)

        mock_report.assert_called_once_with(alignment_preparation_result, hyphy_result, results_dir_path)
