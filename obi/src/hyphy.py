import json
import os
import subprocess
from abc import abstractmethod
from functools import reduce

import requests
from requests_toolbelt import MultipartEncoder


class Hyphy:
    @staticmethod
    def local():
        return LocalHyphy()

    @staticmethod
    def remote():
        return RemoteHyphy()

    @abstractmethod
    def run(self, *args):
        pass

    def _parsed_response(self, data):
        mle_report = data['MLE']
        content = mle_report['content']['0']
        headers = mle_report['headers']
        return list(map(lambda row: self._hyphy_row(row, headers), content))

    def _hyphy_row(self, row, headers):
        def map_headers_with_row(result, index_header):
            index, header = index_header
            result[header[0]] = row[index]

            return result

        return reduce(map_headers_with_row, enumerate(headers), {})


class LocalHyphy(Hyphy):
    def run(self, nucleotide_alignment_path, api_key=None, email=None, bootstrap=1000):
        tree_path = self._generate_tree(nucleotide_alignment_path, bootstrap)
        self._hyphy(nucleotide_alignment_path, tree_path)

        with open(nucleotide_alignment_path + ".MEME.json", 'r') as hyphy_result:
            data = json.load(hyphy_result)

            return self._parsed_response(data)

    def _generate_tree(self, nucleotide_alignment_path, bootstrap):
        os.popen(f'iqtree -s {nucleotide_alignment_path} -bb {bootstrap}').read()

        return f'{nucleotide_alignment_path}.treefile'

    def _hyphy(self, nucleotide_alignment_path, tree_path):
        subprocess.call(['hyphy', 'meme', '--alignment', nucleotide_alignment_path, '-bb', tree_path])


class RemoteHyphy(Hyphy):
    DATAMONKEY_BASE_URL = "http://datamonkey.org/api/v1"

    def run(self, nucleotide_alignment_path, api_key, email, bootstrap=None):
        self._validate_params(api_key, email)
        file_url = upload(nucleotide_alignment_path)
        url = f'{self.DATAMONKEY_BASE_URL}/submit'
        body = {
            'api_key': api_key,
            "method": "MEME",
            "fastaLoc": file_url,
            "mail": email,
            "gencodeid": "Universal",
            "fileExtension": "fasta"
        }
        response = requests.post(url, json=body)

        return response.json()

    def job_result(self, results_dir, job_id):
        response = requests.get('http://datamonkey.org/api/v1/status', json={'method': 'MEME', 'id': job_id})
        data = response.json()

        if data['status'] == 'completed':
            analysis_content = requests.get(f"http://{data['url']}/results").json()

            with open(f"{results_dir}/nucleotide_alignment.fasta.MEME.json", "w") as f:
                f.write(json.dumps(analysis_content, indent=2))

            return self._parsed_response(analysis_content)
        else:
            raise HyphyJobNotReady(f"Job {job_id} is not ready, status: {data['status']}. Try again later.")

    def _key_info(self, api_key):
        response = requests.post(f"{self.DATAMONKEY_BASE_URL}/keyinfo", {'api_key': api_key})

        if response.status_code != 200:
            raise HyphyAPIError(f"Datamonkey API failed: {response.text}")
        return response.json()

    def _validate_params(self, api_key, email):
        if not api_key:
            raise InvalidApiKeyError("API Key is required to submit a remote job")

        if not email:
            raise HyphyAPIError("Email is required to submit a remote job")

        key_info = self._key_info(api_key)
        if key_info["job_remaining"] == 0:
            raise InvalidApiKeyError(f"Key {api_key} expired, to get a new one http://datamonkey.org/apiKey")


class HyphyJobNotReady(RuntimeError):
    def __init__(self, message):
        self.message = message


class InvalidApiKeyError(RuntimeError):
    def __init__(self, message):
        self.message = message


class HyphyAPIError(RuntimeError):
    def __init__(self, message):
        self.message = message


def upload(file_path):
    file = open(file_path, 'rb')
    data = {"reqtype": "fileupload", "time": "1h", "fileToUpload": (file.name, file, 'application/fasta')}
    encoder = MultipartEncoder(fields=data)
    response = requests.post(
        'https://litterbox.catbox.moe/resources/internals/api.php',
        data=encoder,
        headers={'Content-Type': encoder.content_type}
    )
    return response.text
