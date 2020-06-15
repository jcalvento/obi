import os
import subprocess

import requests
from Bio import Entrez
from Bio.Blast import NCBIXML


def most_similar_between(alignment, another_alignment):
    if identity_percentage(alignment) >= identity_percentage(another_alignment) or\
            evalue(alignment) <= evalue(another_alignment) or coverage(alignment) >= coverage(another_alignment):
        return alignment
    return another_alignment


def evalue(alignment):
    hsp = alignment.hsps[0]
    if hasattr(hsp, 'evalue'):
        return hsp.evalue
    else:
        return 0


def identity_percentage(alignment):
    hsp = alignment.hsps[0]
    return round(hsp.identities / hsp.align_length, 2)


def coverage(protein_chain, alignment):
    hsp = alignment.hsps[0]
    return ((hsp.query_end - hsp.query_start) / len(protein_chain)) * 100


def meets_expected_criteria(protein_chain, alignment):
    return identity_percentage(alignment) >= 0.4 and evalue(alignment) <= 0.005 and\
           coverage(protein_chain, alignment) >= 90
# parametrizable pero defaults:
# porcentaje de identidad > 40%
# coverage > 90% ((dfblst['q.end']-dfblst['q.start'])/query_length)*100.0
# evalue < 0.005


def fasta_info(alignment):
    accession = alignment.hit_id.split('|')[1]
    gen_name = alignment.hit_def.split(';')[0].split('RecName: Full=')[-1]
    raw_species = alignment.hit_def.split(';')[-1]
    species = raw_species[raw_species.find("[") + 1:raw_species.find("]")]
    header = f">{accession} {species}, gene for {gen_name}"
    sequence = alignment.hsps[0].sbjct

    return f"{header}\n{sequence}"


def make_fasta(alignments):
    content = "\n\n".join(list(map(lambda alignment: fasta_info(alignment), alignments)))
    fasta_file = open("blast_result.fasta", "w")
    fasta_file.write(content)
    fasta_file.close()
    return fasta_file.name


def blast(input_fasta, db_path):
    # Blast
    blast_result = os.popen(f'blastp -query {input_fasta} -db {db_path} -outfmt 5')
    blast_records = NCBIXML.read(blast_result)

    sequence = []
    with open(input_fasta, "r") as f:
        for line in f.readlines():
            if not line.startswith(">"):
                sequence.append(line)
    protein_chain = "".join(sequence).replace("\n", "")

    # Mejores hits segun criterio
    filtered_results = list(filter(
        lambda alignment: meets_expected_criteria(protein_chain, alignment), blast_records.alignments
    ))

    # Crear Fasta con hits
    return make_fasta(filtered_results)


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
    Entrez.email = 'juliancalvento@gmail.com'
    entrez_ids = list(map(lambda entrez_dict: entrez_dict['entrez_id'], ids))
    fetch_response = Entrez.efetch(db="nucleotide", id=entrez_ids, rettype="gb", retmode="xml")
    parsed_response = Entrez.read(fetch_response)
    # # del fetch response[0] sacar 'GBSeq_feature-table', y de la parte de CDS 'GBFeature_location', después cortar la secuencia desde inicio - 1 al fin
    # # La secuencia está en 'GBSeq_sequence'
    return parsed_response


if __name__ == '__main__':
    fasta_file = blast('./query.fasta', '../swissprot/swissprot')

    # Checkpoint
    cd_hit_output_file = cd_hit(fasta_file)
    clustal_output = clustal(cd_hit_output_file)
    uniprot_to_entrez = map_uniprot_ids_to_entrez_ids(clustal_output)
    fetch_cds_from_entrez(uniprot_to_entrez)
