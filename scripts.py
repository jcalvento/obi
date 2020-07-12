import csv
import os
import subprocess
import datetime
import urllib

import requests

from blast import Blast, EmptyBlastResultError, BlastResultsError
from entrez_api_client import EntrezApiClient


class InvalidEntrezIds(RuntimeError):
    def __init__(self, message):
        self.message = message


def map_uniprot_ids_to_entrez_ids(fasta_file):
    # Mappeo de Uniprot a Entrez
    uniprot_ids = []
    with open(fasta_file, "r") as f:
        for line in f.readlines():
            if line.startswith(">"):
                uniprot_ids.append(line.split(" ")[0].replace(">", ""))

    uniprot_response = requests.post(
        url='https://www.uniprot.org/uploadlists/',
        data={
            'from': 'ACC+ID',
            'to': 'P_ENTREZGENEID',
            'format': 'tab',
            'query': " ".join(uniprot_ids)
        }
    )

    assert uniprot_response.status_code == requests.status_codes.codes.ok

    return list(map(
        lambda result_line: {"uniprot_id": result_line.split("\t")[0], "entrez_id": result_line.split("\t")[1]},
        uniprot_response.text.split("\n")[1:-1]
    ))


def cd_hit(input_fasta):
    cd_hit_output_file = 'results/clust100_' + input_fasta.split('/')[-1].split('.')[0]
    subprocess.call(['./../cdhit/cd-hit', '-i', input_fasta, '-o', cd_hit_output_file, '-c', '1'])

    return cd_hit_output_file


def clustal(input_fasta):
    clustal_output = f"results/clustalo_aligned{input_fasta.split('/')[-1].split('.')[0]}.fasta"
    subprocess.call(['clustalo', '-i', input_fasta, '-o', clustal_output, "--force"])

    return clustal_output


def fetch_cds_from_entrez(ids):
    if not ids:
        raise InvalidEntrezIds("Empty ids lists")
    entrez_ids = list(map(lambda entrez_dict: entrez_dict['entrez_id'], ids))
    entrez_api_client = EntrezApiClient('juliancalvento@gmail.com')
    try:
        entez_response = entrez_api_client.efetch(entrez_ids)
    except urllib.error.HTTPError:
        print(f"There's been an error with Entrez, ids: {entrez_ids}")
        return
    # # del fetch response[0] sacar 'GBSeq_feature-table', y de la parte de CDS 'GBFeature_location', después cortar la secuencia desde inicio - 1 al fin
    # # La secuencia está en 'GBSeq_sequence'
    return entez_response


def alignment_preparation(fasta_file):
    cd_hit_output_file = cd_hit(fasta_file)
    uniprot_to_entrez = map_uniprot_ids_to_entrez_ids(cd_hit_output_file)
    fetch_cds_from_entrez(uniprot_to_entrez)
    # clustal_output = clustal(cd_hit_output_file)


if __name__ == '__main__':
    files = os.listdir('./fasta')
    print(f"Whole process init time: {datetime.datetime.now()}")
    failed_count = 0
    for file in files:
        print(f"Init time: {datetime.datetime.now()}")
        try:
            fasta_file = Blast().run(f"fasta/{file}", '../swissprot/swissprot')
        except (EmptyBlastResultError, BlastResultsError) as e:
            print(e.message)
            failed_count += 1
            continue
        # Checkpoint
        try:
            alignment_preparation(fasta_file)
        except InvalidEntrezIds:
            failed_count += 1
        print(f"Finish time: {datetime.datetime.now()}")
    print(f"Whole process finish time: {datetime.datetime.now()}, failed: {failed_count}")
