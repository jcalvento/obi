import requests
from marshmallow import Schema, fields, post_load


def parse_response(uniprot_response):
    return list(map(
        lambda result_line: UniprotIdMapping.for_tab_separated_ids(result_line), uniprot_response.text.split("\n")[1:-1]
    ))


class UniprotIdMappingSchema(Schema):
    from_id = fields.Str()
    to_id = fields.Str()

    @post_load
    def deserialize(self, data, **kwargs):
        return UniprotIdMapping(**data)


class UniprotIdMapping:
    @classmethod
    def for_tab_separated_ids(cls, tab_separated_ids):
        return cls(*tab_separated_ids.split("\t"))

    def __init__(self, from_id, to_id):
        self.from_id = from_id
        self.to_id = to_id


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
