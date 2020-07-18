import datetime
import os
import pathlib
import shutil
import subprocess

import requests

from src.blast import Blast, EmptyBlastResultError, BlastResultsError
from src.entrez import EntrezDB, InvalidEntrezIds


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


def cd_hit(input_fasta, results_dir):
    cd_hit_output_file = f'{results_dir}/cd_hit'
    subprocess.call(['./../cdhit/cd-hit', '-i', input_fasta, '-o', cd_hit_output_file, '-c', '1'])

    return cd_hit_output_file


def clustal(input_fasta, results_dir):
    clustal_output = f"{results_dir}/clustalo_aligned.fasta"
    subprocess.call(['clustalo', '-i', input_fasta, '-o', clustal_output, "--force"])

    return clustal_output


def protein_based_nucleotide_alignment(entrez_response, protein_alignment_path, results_dir):
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
        # TODO: Validar largo de la cadena con el resultado del blastcmd? Alignment podria agrega cosas al final
        init, end = entrez_row['location'].split('..')
        adn_sequence = entrez_row['sequence'][int(init) - 1:int(end) - 3]  # -3 removes stop codon
        adn_codons = [adn_sequence[index:index + 3] for index in range(0, len(adn_sequence), 3)]
        nucleotide_alignment = ''
        codon_index = 0
        for amino_acid in alignment:
            if amino_acid == '-':
                nucleotide_alignment += '---'
            else:
                nucleotide_alignment += adn_codons[codon_index]
                codon_index += 1
        nucleotide_alignments[alignment_id] = nucleotide_alignment

    nucleotide_fasta = open(f"{results_dir}/nucleotide_alignment.fasta", "w")
    nucleotide_fasta.write('\n'.join(list(map(lambda key: f'>{key}\n{nucleotide_alignments[key]}', nucleotide_alignments))))
    nucleotide_fasta.close()
    return nucleotide_fasta.name


def alignment_preparation(fasta_file, results_dir):
    cd_hit_output_file = cd_hit(fasta_file, results_dir)
    uniprot_to_entrez = map_uniprot_ids_to_entrez_ids(cd_hit_output_file)
    entrez_response = EntrezDB().fetch_cds(uniprot_to_entrez)
    clustal_output = clustal(cd_hit_output_file, results_dir)
    return protein_based_nucleotide_alignment(entrez_response, clustal_output, results_dir)


def generate_tree(nucleotide_alignment_path, boostrap):
    os.popen(f'iqtree -s {nucleotide_alignment_path} -bb {boostrap}').read()

    return f'{nucleotide_alignment_path}.treefile'


def hyphy(nucleotide_alignment_path, tree_path):
    subprocess.call(['hyphy', 'meme', '--alignment', nucleotide_alignment_path, '-bb', tree_path])
    # os.popen(f'hyphy meme --alignment {nucleotide_alignment_path} --tree {tree_path}').read()


def step_2(nucleotide_alignment_path, boostrap):
    tree_path = generate_tree(nucleotide_alignment_path, boostrap)
    hyphy(nucleotide_alignment_path, tree_path)


def create_results_dir(file):
    filename = file.split('/')[-1].split('.')[0]
    dir_path = f'{pathlib.Path(__file__).parent.parent.absolute()}/results/{filename}'
    if os.path.isdir(dir_path):
        shutil.rmtree(dir_path)
    os.mkdir(dir_path)
    return dir_path


if __name__ == '__main__':
    fastas_dir = f"{pathlib.Path(__file__).parent.parent.absolute()}/fasta"
    files = os.listdir(fastas_dir)
    print(f"Whole process init time: {datetime.datetime.now()}")
    failed_count = 0
    blast = Blast(f'{pathlib.Path(__file__).parent.parent.parent.absolute()}/swissprot/swissprot')

    for file in [files[0]]:
        print(f"Init time: {datetime.datetime.now()}")
        results_dir = create_results_dir(file)
        try:
            fasta_file = blast.run(f"{fastas_dir}/{file}", results_dir)
        except (EmptyBlastResultError, BlastResultsError) as e:
            print(e.message)
            failed_count += 1
            continue
        try:
            nucleotide_alignment = alignment_preparation(fasta_file, results_dir)
            # step_2(nucleotide_alignment, 1000)
        except InvalidEntrezIds:
            failed_count += 1
        print(f"Finish time: {datetime.datetime.now()}")
    print(f"Whole process finish time: {datetime.datetime.now()}, failed: {failed_count}")
