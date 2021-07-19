import urllib
from datetime import datetime

from Bio import Entrez

from src.obi.utils import detect


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
    def __init__(self, uniprot_id, location, translation, sequence, locus_version):
        self.locus_version = locus_version
        self.sequence = sequence
        self.translation = translation
        self.location = location
        self.uniprot_id = uniprot_id


class EntrezDB:
    def __init__(self, email):
        self._entrez_api_client = EntrezApiClient(email)

    def fetch_cds(self, ids_mapping):
        if not ids_mapping:
            raise InvalidEntrezIds("Empty ids lists")
        entrez_ids = list(map(lambda id_mapping: id_mapping.to_id.split(".")[0], ids_mapping))
        try:
            entrez_response = self._entrez_api_client.efetch(entrez_ids)
        except urllib.error.HTTPError:
            print(f"There's been an error with Entrez, ids: {entrez_ids}")
            return
        return self._parse_response(entrez_response, ids_mapping)

    def _parse_response(self, entrez_response, ids_mapping):
        parsed_response = []
        for element in entrez_response:
            try:
                if not element.get('GBSeq_sequence'):
                    raise InvalidEntrezIds("Sequence not found")  # TODO: Test validations
                parsed_response.append(
                    EntrezElement(
                        uniprot_id=self._uniprot_id(element['GBSeq_locus'], ids_mapping),
                        location=self._cds_location(element['GBSeq_feature-table']),
                        translation=self._entrez_translation(element['GBSeq_feature-table']),
                        sequence=element['GBSeq_sequence'],
                        locus_version=element['GBSeq_accession-version']
                    )
                )
            except InvalidEntrezIds as e:
                print(f"{self._uniprot_id(element['GBSeq_locus'], ids_mapping)} {e}")

        return parsed_response

    def _uniprot_id(self, primary_accession, ids_mapping):
        return detect(
            lambda id_mapping: primary_accession in id_mapping.to_id.split(".")[0],
            ids_mapping, if_not_none=lambda id_mapping: id_mapping.from_id
        )

    def _cds_location(self, features):
        return detect(
            self._is_cds_key,
            features, if_not_none=lambda feature: feature['GBFeature_location'],
            if_none=lambda: self.__throw(InvalidEntrezIds("CDS not found"))
        )

    def _is_cds_key(self, feature):
        return feature['GBFeature_key'] == 'CDS'

    def _entrez_translation(self, features):
        feature_qualifiers = detect(
            self._is_cds_key,
            features, if_not_none=lambda feature: feature['GBFeature_quals']
        )
        return detect(
            lambda qualifier: qualifier['GBQualifier_name'] == 'translation',
            feature_qualifiers, if_not_none=lambda qualifier: qualifier['GBQualifier_value']
        )

    def __throw(self, exception):
        raise exception
