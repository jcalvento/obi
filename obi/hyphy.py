import json
import os
import subprocess
from functools import reduce


class Hyphy:
    def __init__(self, nucleotide_alignment_path):
        self._nucleotide_alignment_path = nucleotide_alignment_path

    def run(self, bootstrap=1000):
        tree_path = self._generate_tree(bootstrap)
        self._hyphy(tree_path)

        return self._parsed_response()

    def _generate_tree(self, boostrap):
        os.popen(f'iqtree -s {self._nucleotide_alignment_path} -bb {boostrap}').read()

        return f'{self._nucleotide_alignment_path}.treefile'

    def _hyphy(self, tree_path):
        subprocess.call(['hyphy', 'meme', '--alignment', self._nucleotide_alignment_path, '-bb', tree_path])

    def _parsed_response(self):
        with open(self._nucleotide_alignment_path + ".MEME.json", 'r') as hyphy_result:
            data = json.load(hyphy_result)
            mle_report = data['MLE']
            content = mle_report['content']['0']  # puede haber más de 1?
            headers = mle_report['headers']
            return list(map(lambda row: self._hyphy_row(row, headers), content))

    def _hyphy_row(self, row, headers):
        def map_headers_with_row(result, index_header):
            index, header = index_header
            result[header[0]] = row[index]

            return result

        return reduce(map_headers_with_row, enumerate(headers), {})
