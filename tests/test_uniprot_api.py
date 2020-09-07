from src.uniprot_api import UniprotAPIClient


class TestUniprotApiClient:
    def test_given_a_list_of_uniprot_ids_it_returns_the_mapping_to_refseq_nucleotide_id(self):
        api = UniprotAPIClient()
        uniprot_ids = ['P14080.2', 'O65493.1', 'Q9LM66.2']

        ids_mapping = api.refseq_ids(uniprot_ids)

        mapping = ids_mapping[0]
        another_mapping = ids_mapping[1]
        assert len(ids_mapping) == 2
        assert 'O65493' == mapping.from_id
        assert 'NM_119701.4' == mapping.to_id
        assert 'Q9LM66' == another_mapping.from_id
        assert 'NM_101938.5' == another_mapping.to_id

    def test_given_an_id_with_no_refseq_mapping_it_returns_an_empty_list(self):
        api = UniprotAPIClient()
        uniprot_ids = ['P14080.2']

        ids_mapping = api.refseq_ids(uniprot_ids)

        assert len(ids_mapping) == 0

    def test_given_an_empty_list_of_ids_it_returns_another_empty_list(self):
        api = UniprotAPIClient()

        ids_mapping = api.refseq_ids([])

        assert len(ids_mapping) == 0
