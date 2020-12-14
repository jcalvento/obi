import argparse
import json

from obi.alignment_preparation import AlignmentPreparationResultSchema
from obi.positive_selection import PositiveSelection

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--input-path', help='Path of alignment preparation result')

    args = parser.parse_args()

    with open(args.input_path + "/alignment_preparation_result.json", 'r') as file_content:
        data = json.load(file_content)
        alignment_preparation_result = AlignmentPreparationResultSchema().load(data)
        PositiveSelection().analyse(args.input_path, alignment_preparation_result)
