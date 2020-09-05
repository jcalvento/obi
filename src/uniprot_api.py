import requests


def parse_response(uniprot_response):
    return list(map(
        lambda result_line: UniprotIdMapping(result_line), uniprot_response.text.split("\n")[1:-1]
    ))


class UniprotIdMapping:
    def __init__(self, tab_separated_ids):
        self.from_id, self.to_id = tab_separated_ids.split("\t")


class UniprotAPIClient:
    UNIPROT_UPLOADLISTS_URL = 'https://www.uniprot.org/uploadlists/'

    def refseq_ids(self, uniprot_ids):
        uniprot_response = requests.post(
            url=self.UNIPROT_UPLOADLISTS_URL,
            data={
                'from': 'ACC+ID',
                'to': 'REFSEQ_NT_ID',
                'format': 'tab',
                'query': " ".join(uniprot_ids)
            }
        )

        assert uniprot_response.status_code == requests.status_codes.codes.ok

        return parse_response(uniprot_response)
