import argparse
import json

from src.obi.alignment_preparation import AlignmentPreparationResultSchema
from src.obi.hyphy import HyphyJobNotReady
from src.obi.positive_selection import PositiveSelectionAnalyzer

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--input-path', help='Path of alignment preparation result')
    args = parser.parse_args()

    with open(args.input_path + "/alignment_preparation_result.json", 'r') as file_content:
        data = json.load(file_content)
        alignment_preparation_result = AlignmentPreparationResultSchema().load(data)
        try:
            PositiveSelectionAnalyzer.remote().resume(args.input_path, alignment_preparation_result)
        except HyphyJobNotReady as e:
            print(e.message)
