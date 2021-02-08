from functools import reduce
import requests


def upload(file_path):
    # https://github.com/espebra/filebin#web-service
    with open(file_path, "r") as f:
        content = f.read()
        response = requests.post(
            url='https://filebin.net/',
            data=content,
            headers={'Content-Type': 'application/octet-stream', 'filename': f.name.split("/")[-1]}
        )
        link = response.json()['links'][0]['href']
        html = requests.get(link).text
        html_link = link + '?t='
        return html[html.find(html_link):html.find(html_link) + len(html_link) + 8]


class Hyphy:
    def __init__(self, nucleotide_alignment_path):
        self._nucleotide_alignment_path = nucleotide_alignment_path

    def run(self):
        file_url = upload(self._nucleotide_alignment_path)
        # file_url = 'https://filebin.net/cxnvt3xfabuukw3x/nucleotide_alignment.fasta?t=wsa9tj86'
        url = 'http://datamonkey.org/api/v1/submit'
        body = {
            'api_key': '6020736d27e5ee138cbb0b3c',
            "method": "MEME",
            "fastaLoc": file_url,
            "mail": "juliancalvento@gmail.com",
            "gencodeid": "Universal",
            "fileExtension": "fasta"
        }
        req = requests.post(url, json=body)
        data = req.json()
        print(data)

    def check_status(self):
        # {'time_stamp': '2021-02-08T00:08:11.876Z', 'id': '602080eb27e5ee138cbb0bb7', 'status': 'queue',
        # 'url': 'datamonkey.org/MEME/602080eb27e5ee138cbb0bb7'}
        response = requests.get('http://datamonkey.org/api/v1/status', json={'method': 'MEME', 'id': '602080eb27e5ee138cbb0bb7'})
        data = response.json()
        if data['status'] == 'completed':
            analysis_content = requests.get(f"http://{data['url']}/results").json()
            print(self.parse_response(analysis_content))

    def key_info(self):
        response = requests.post("http://datamonkey.org/api/v1/keyinfo", {'api_key': '6020736d27e5ee138cbb0b3c'})
        print(response.json())

    def parse_response(self, analysis_content):
        mle_report = analysis_content['MLE']
        content = mle_report['content']['0']  # puede haber más de 1?
        headers = mle_report['headers']
        return list(map(lambda row: self._hyphy_row(row, headers), content))

    def _hyphy_row(self, row, headers):
        def map_headers_with_row(result, index_header):
            index, header = index_header
            result[header[0]] = row[index]

            return result

        return reduce(map_headers_with_row, enumerate(headers), {})


#  Para usar este metodo, tener en cuenta servicio donde subir el archivo
# https://github.com/espebra/filebin#web-service
# https://gofile.io/welcome
# Hyphy API https://github.com/veg/datamonkey-js/wiki/API-Submit
# El token es a mano, tendria que pedirlo por param
if __name__ == '__main__':
    # Hyphy("/Users/julian/Documents/UNQ/tesis/pruebas/results/P00766/nucleotide_alignment.fasta").check_status()
    Hyphy("/Users/julian/Documents/UNQ/tesis/pruebas/results/P00766/nucleotide_alignment.fasta").run()
    # Hyphy("/Users/julian/Documents/UNQ/tesis/pruebas/results/P00784/nucleotide_alignment.fasta").key_info()
