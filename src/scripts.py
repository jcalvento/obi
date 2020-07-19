import datetime
import os
import pathlib
import shutil
import subprocess

import requests

from src.blast import Blast, BlastResultsError
from src.entrez import EntrezDB, InvalidEntrezIds
from src.nucleotide_aligner import NucleotideAligner


def map_uniprot_ids_to_entrez_ids(fasta_file):
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


def alignment_preparation(fasta_file, results_dir):
    cd_hit_output_file = cd_hit(fasta_file, results_dir)
    uniprot_to_entrez = map_uniprot_ids_to_entrez_ids(cd_hit_output_file)
    entrez_response = EntrezDB().fetch_cds(uniprot_to_entrez)
    clustal_output = clustal(cd_hit_output_file, results_dir)
    return NucleotideAligner().protein_based_nucleotide_alignment(entrez_response, clustal_output, results_dir)


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
    blast = Blast(f'{pathlib.Path(__file__).parent.parent.absolute()}/swissprot/swissprot')

    for file in [files[0]]:
        print(f"Init time: {datetime.datetime.now()}")
        results_dir = create_results_dir(file)
        try:
            fasta_file = blast.run(f"{fastas_dir}/{file}", results_dir)
        except BlastResultsError as e:
            print(e.message)
            failed_count += 1
            continue
        try:
            nucleotide_alignment = alignment_preparation(fasta_file, results_dir)
            step_2(nucleotide_alignment, 1000)
        except InvalidEntrezIds:
            failed_count += 1
        print(f"Finish time: {datetime.datetime.now()}")
    print(f"Whole process finish time: {datetime.datetime.now()}, failed: {failed_count}")
