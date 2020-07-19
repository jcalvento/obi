import urllib
from datetime import datetime

from Bio import Entrez

from src.utils import detect


class InvalidEntrezIds(RuntimeError):
    def __init__(self, message):
        self.message = message


class EntrezApiClient:
    def __init__(self, email):
        self.__email = email
        self.__num_calls = 0
        self.__last_reset = datetime.now()
        Entrez.email = self.__email

    def efetch(self, ids):
        if (datetime.now() - self.__last_reset).seconds < 1 and self.__num_calls >= 3:
            raise RuntimeError("Entrez api calls rate limit exceeded. Try again in a few seconds")
        elif (datetime.now() - self.__last_reset).seconds < 1 and self.__num_calls < 3:
            self.__num_calls += 1
            return self.__efetch(ids)
        else:
            self.__num_calls = 0
            self.__last_reset = datetime.now()
            return self.__efetch(ids)

    def __efetch(self, ids):
        fetch_response = Entrez.efetch(db="nucleotide", id=ids, rettype="gb", retmode="xml")
        return Entrez.read(fetch_response)


class EntrezElement:
    def __init__(self, uniprot_id, location, translation, sequence):
        self.sequence = sequence
        self.translation = translation
        self.location = location
        self.uniprot_id = uniprot_id


class EntrezDB:
    def __init__(self):
        self._entrez_api_client = EntrezApiClient('juliancalvento@gmail.com')

    def fetch_cds(self, ids_mapping):
        if not ids_mapping:
            raise InvalidEntrezIds("Empty ids lists")
        entrez_ids = list(map(lambda id_mapping: id_mapping.to_id, ids_mapping))
        try:
            entrez_response = self._entrez_api_client.efetch(entrez_ids)
        except urllib.error.HTTPError:
            print(f"There's been an error with Entrez, ids: {entrez_ids}")
            return
        return self.__parse_response(entrez_response, ids_mapping)

    def __parse_response(self, entrez_response, ids_mapping):
        return list(map(lambda element: EntrezElement(
            uniprot_id=self.__uniprot_id(element['GBSeq_other-seqids'], ids_mapping),
            location=self.__cds_location(element['GBSeq_feature-table']),
            translation=self.__entrez_translation(element['GBSeq_feature-table']),
            sequence=element['GBSeq_sequence']
        ), entrez_response))

    def __uniprot_id(self, sequence_ids, ids_mapping):
        entrez_id = detect(
            lambda seq_id: seq_id.startswith('ref|'),
            sequence_ids,
            if_not_none=lambda seq_id: seq_id.replace('ref|', '').replace('|', '')
        )
        return detect(
            lambda id_mapping: entrez_id in id_mapping.to_id,
            ids_mapping, if_not_none=lambda id_mapping: id_mapping.from_id
        )

    def __cds_location(self, features):
        return detect(
            self.is_cds_key,
            features, if_not_none=lambda feature: feature['GBFeature_location']
        )

    def is_cds_key(self, feature):
        return feature['GBFeature_key'] == 'CDS'

    def __entrez_translation(self, features):
        feature_qualifiers = detect(
            self.is_cds_key,
            features, if_not_none=lambda feature: feature['GBFeature_quals']
        )
        return detect(
            lambda qualifier: qualifier['GBQualifier_name'] == 'translation',
            feature_qualifiers, if_not_none=lambda qualifier: qualifier['GBQualifier_value']
        )
