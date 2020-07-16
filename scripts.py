import datetime
import os
import subprocess
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
            'to': 'REFSEQ_NT_ID',
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


def uniprot_id(sequence_ids, ids):
    entrez_id = next((x for x in sequence_ids if x.startswith('ref|')), None).replace('ref|', '').replace('|', '')
    return next((ids_map['uniprot_id'] for ids_map in ids if entrez_id in ids_map['entrez_id']), None)


def cds_location(features):
    return next((feature['GBFeature_location'] for feature in features if feature['GBFeature_key'] == 'CDS'), None)


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
    return list(map(lambda element: {
        'uniprot_id': uniprot_id(element['GBSeq_other-seqids'], ids),
        'location': cds_location(element['GBSeq_feature-table']),
        'sequence': element['GBSeq_sequence']
    }, entez_response))


def protein_based_nucleotide_alignment(entrez_response, protein_alignment_path):
    alignments = {}
    current = None
    with open(protein_alignment_path, "r") as f:
        for line in f.readlines():
            if line.startswith(">"):
                gen_id = line.split(" ")[0].replace(">", "")
                alignments[gen_id] = ''
                current = gen_id
            else:
                alignments[current] += line.replace('\n', '')
    nucleotide_alignments = {}
    for entrez_row in entrez_response:
        alignment_id = next((x for x in alignments if x.startswith(entrez_row['uniprot_id'])), None)
        alignment = alignments[alignment_id]



def alignment_preparation(fasta_file):
    cd_hit_output_file = cd_hit(fasta_file)
    uniprot_to_entrez = map_uniprot_ids_to_entrez_ids(cd_hit_output_file)
    entrez_response = fetch_cds_from_entrez(uniprot_to_entrez)
    clustal_output = clustal(cd_hit_output_file)
    # clustal_output = 'results/clustalo_alignedclust100_P00784_blast_result.fasta'
    # entrez_response = [{'uniprot_id': 'O65493', 'location': '1..312', 'sequence': 'tttttttttttttttttggnatggtcagtgcattttattgaatcagcacagtacaaaaataaataaaaataagggaagggganttaaattacagccaaactgagcttcatgactttgtcagattataaaccacacatactcacaaacacactacacatacatacaaaatgagacaaatntaaattaatnttaacantacccacagtttgggtcaaaggattagnctacaggaaggaanttgcactaaaaanccaacatacatcacacgngtgtanttgggcgtttcaaatttacaggctntggattanttct'}, {'uniprot_id': 'Q9LM66', 'location': '1..470', 'sequence': 'cttctgtgtaggtgaccggagcactgagaggcagctctgatgcactattgtgtgtcagcagctcaaaggccctaaaacactgaaggttctgcatctgaagtattagattgttagcagcaaaatatgaaagatgaggtggacagtcctctaagccctatttagggaagcttttccaagccacaatcttaactacctacccaaaggatttgcattacccccagattctgtgccaacaaccttttaaggaaatacagtccttgggaaatgagttttgatggtgaattggggtgttaaggaagggaaagattgtcatagatggtagggctttgaaaaatgcagggtattcagcttgccactcctgggctttcaacacatttgagttcacttgcctaggacgggttctcttgggtctttatttccccatnctgggcccattgctttaaatactattttgtttgaaaattaatttt'}]
    protein_based_nucleotide_alignment(entrez_response, clustal_output)


if __name__ == '__main__':
    files = os.listdir('./fasta')
    print(f"Whole process init time: {datetime.datetime.now()}")
    failed_count = 0
    for file in [files[0]]:
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
