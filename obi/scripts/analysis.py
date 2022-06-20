import argparse
import json

from obi.src.alignment_preparation import AlignmentPreparationResultSchema
from obi.src.positive_selection import PositiveSelectionAnalyzer


def analysis():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input-path', help='Path of alignment preparation result')
    parser.add_argument(
        '--mode', help='"local" if you want to run Hyphy locally, "remote" to use Datamonkey instead. Default: local',
        default=PositiveSelectionAnalyzer.LOCAL_MODE
    )
    parser.add_argument(
        "--email", help='Required when running in remote mode. You will get notified once job is done'
    )
    args = parser.parse_args()
    with open(args.input_path + "/alignment_preparation_result.json", 'r') as file_content:
        data = json.load(file_content)
        alignment_preparation_result = AlignmentPreparationResultSchema().load(data)
        PositiveSelectionAnalyzer.for_mode(args.mode).analyse(args.input_path, alignment_preparation_result, args.email)


if __name__ == "__main__":
    analysis()
