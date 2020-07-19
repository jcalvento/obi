import requests


class UniprotIdMapping:
    def __init__(self, tab_separated_ids):
        self.from_id, self.to_id = tab_separated_ids.split("\t")


class UniprotAPIClient:
    UNIPROT_UPLOADLISTS_URL = 'https://www.uniprot.org/uploadlists/'

    def refseq_ids(self, fasta_file):
        uniprot_ids = self.__parse_uniprot_ids(fasta_file)

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

        return self.__parse_response(uniprot_response)

    def __parse_response(self, uniprot_response):
        return list(map(
            lambda result_line: UniprotIdMapping(result_line), uniprot_response.text.split("\n")[1:-1]
        ))

    def __parse_uniprot_ids(self, fasta_file):
        uniprot_ids = []
        with open(fasta_file, "r") as f:
            for line in f.readlines():
                if line.startswith(">"):
                    uniprot_ids.append(line.split(" ")[0].replace(">", ""))
        return uniprot_ids
