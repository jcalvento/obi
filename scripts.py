import subprocess
import requests

from blast import Blast
from entrez_api_client import EntrezApiClient


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
    cd_hit_output_file = 'clust100_' + input_fasta
    subprocess.call(['./../cdhit/cd-hit', '-i', input_fasta, '-o', cd_hit_output_file, '-c', '1'])

    return cd_hit_output_file


def clustal(input_fasta):
    clustal_output = "clustalo_aligned.fasta"
    subprocess.call(['clustalo', '-i', input_fasta, '-o', clustal_output, "--force"])

    return clustal_output


def fetch_cds_from_entrez(ids):
    entrez_ids = list(map(lambda entrez_dict: entrez_dict['entrez_id'], ids))
    entrez_api_client = EntrezApiClient('juliancalvento@gmail.com')
    entez_response = entrez_api_client.efetch(entrez_ids)
    # # del fetch response[0] sacar 'GBSeq_feature-table', y de la parte de CDS 'GBFeature_location', después cortar la secuencia desde inicio - 1 al fin
    # # La secuencia está en 'GBSeq_sequence'
    return entez_response


def alignment_preparation(fasta_file):
    cd_hit_output_file = cd_hit(fasta_file)
    clustal_output = clustal(cd_hit_output_file)
    uniprot_to_entrez = map_uniprot_ids_to_entrez_ids(clustal_output)
    fetch_cds_from_entrez(uniprot_to_entrez)


if __name__ == '__main__':
    fasta_file = Blast().run('./query.fasta', '../swissprot/swissprot')

    # Checkpoint
    alignment_preparation(fasta_file)
