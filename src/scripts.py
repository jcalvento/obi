import datetime
import os
import pathlib
import shutil
import subprocess

from src.blast import BlastResultsError, Blast
from src.entrez import EntrezDB, InvalidEntrezIds
from src.hyphy import Hyphy
from src.nucleotide_aligner import NucleotideAligner
from src.uniprot_api import UniprotAPIClient


def cd_hit(input_fasta, results_dir):
    cd_hit_output_file = f'{results_dir}/cd_hit'
    subprocess.call(['cd-hit', '-i', input_fasta, '-o', cd_hit_output_file, '-c', '1'])

    return cd_hit_output_file


def clustal(input_fasta, results_dir):
    clustal_output = f"{results_dir}/clustalo_aligned.fasta"
    subprocess.call(['clustalo', '-i', input_fasta, '-o', clustal_output, "--force"])

    return clustal_output


def alignment_preparation(fasta_file, results_dir, email):
    cd_hit_output_file = cd_hit(fasta_file, results_dir)
    uniprot_to_entrez = UniprotAPIClient().refseq_ids(cd_hit_output_file)
    entrez_response = EntrezDB(email).fetch_cds(uniprot_to_entrez)
    clustal_output = clustal(cd_hit_output_file, results_dir)
    return NucleotideAligner().protein_based_nucleotide_alignment(entrez_response, clustal_output, results_dir)


def create_results_dir(file):
    filename = file.split('/')[-1].split('.')[0]
    dir_path = f'{root_path}/results/{filename}'
    if os.path.isdir(dir_path):
        shutil.rmtree(dir_path)
    os.mkdir(dir_path)
    return dir_path


if __name__ == '__main__':
    print(f"Whole process init time: {datetime.datetime.now()}")
    email = 'juliancalvento@gmail.com'
    root_path = pathlib.Path(__file__).parent.parent.absolute()
    fastas_dir = f"{root_path}/fasta"
    files = os.listdir(fastas_dir)
    failed_count = 0
    blast = Blast(f'{root_path}/swissprot/swissprot')

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
            nucleotide_alignment = alignment_preparation(fasta_file, results_dir, email)
            Hyphy(nucleotide_alignment).run(1000)
        except InvalidEntrezIds:
            failed_count += 1
        print(f"Finish time: {datetime.datetime.now()}")
    print(f"Whole process finish time: {datetime.datetime.now()}, failed: {failed_count}")
