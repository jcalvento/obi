import subprocess
from io import StringIO

import requests
from Bio.Blast import NCBIXML, NCBIWWW


def perform_request(path):
    return requests.get(
        url=f"https://rest.ensembl.org/{path}",
        headers={"Content-type": "text/x-fasta"}
    ).text


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


if __name__ == '__main__':
    # Blast
    protein_chain = 'MKWVTFISLLFLFSSAYSRGVFRRDAHKSEVAHRFKDLGEENFKALVLIAFAQYLQQCPFEDHVKLVNEVTEFAKTCVADESAENCDKSLHTLFGDKLCTVATLRETYGEMADCCAKQEPERNECFLQHKDDNPNLPRLVRPEVDVMCTAFHDNEETFLKKYLYEIARRHPYFYAPELLFFAKRYKAAFTECCQAADKAACLLPKLDELRDEGKASSAKQRLKCASLQKFGERAFKAWAVARLSQRFPKAEFAEVSKLVTDLTKVHTECCHGDLLECADDRADLAKYICENQDSISSKLKECCEKPLLEKSHCIAEVENDEMPADLPSLAADFVESKDVCKNYAEAKDVFLGMFLYEYARRHPDYSVVLLLRLAKTYETTLEKCCAAADPHECYAKVFDEFKPLVEEPQNLIKQNCELFKQLGEYKFQNALLVRYTKKVPQVSTPTLVEVSRNLGKVGSKCCKHPEAKRMPCAEDYLSVVLNQLCVLHEKTPVSDRVTKCCTESLVNRRPCFSALEVDETYVPKEFNAETFTFHADICTLSEKERQIKKQTALVELVKHKPKATKEQLKAVMDDFAAFVEKCCKADDKETCFAEEGKKLVAASQAALGL'
    # scan_result = NCBIWWW.qblast(
    #     "blastp", "swissprot", protein_chain, word_size=6,
    #     threshold=10, matrix_name="BLOSUM62", gapcosts="11 1"
    # )
    with open("blast_result.xml", "r") as f:
        file = f.read()
    scan_result = StringIO(file)
    blast_records = NCBIXML.read(scan_result)

    # Mejores hits segun criterio
    filtered_records = list(filter(
        lambda alignment: meets_expected_criteria(protein_chain, alignment), blast_records.alignments
    ))

    # Crear Fasta con hits
    fasta_file = make_fasta(filtered_records)

    # Clustering con cd-hit
    cd_hit_output_file = 'clust100_' + fasta_file
    subprocess.call(['./../cdhit/cd-hit', '-i', fasta_file, '-o', cd_hit_output_file, '-c', '1'])

    # Alineamiento con Clustal Omega
    clustal_output = "clustalo_aligned.fasta"
    subprocess.call(['clustalo', '-i', cd_hit_output_file, '-o', clustal_output, "--force"])

    # Mappeo de Uniprot a Ensembl
    uniprot_ids = []
    with open(clustal_output, "r") as clustal_file:
        for line in clustal_file.readlines():
            if line.startswith(">"):
                uniprot_ids.append(line.split(" ")[0].replace(">", ""))

    uniprot_response = requests.post(
        url='https://www.uniprot.org/uploadlists/',
        data={
            'from': 'ACC+ID',
            'to': 'ENSEMBL_ID',
            'format': 'tab',
            'query': " ".join(uniprot_ids)
        }
    )

    uniprot_to_ensembl = []
    for result_line in uniprot_response.text.split("\n")[1:-1]:
        uniprot_to_ensembl.append({"uniprot_id": result_line.split("\t")[0], "ensembl_id": result_line.split("\t")[1]})

    # Lookup en Ensembl
    lookup_responses = []
    for uniprot_mapping in uniprot_to_ensembl:
        lookup_response = requests.get(
            url=f"https://rest.ensembl.org/lookup/id/{uniprot_mapping['ensembl_id']}?expand=1",
            headers={"Content-type": "application/json"}
        )

        lookup_responses.append(lookup_response.json())
    lookup_responses
