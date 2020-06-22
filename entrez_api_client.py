from datetime import datetime

from Bio import Entrez


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
