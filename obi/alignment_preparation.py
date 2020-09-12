import csv
import subprocess
from shutil import copyfile

from obi.entrez import EntrezDB
from obi.nucleotide_aligner import NucleotideAligner
from obi.uniprot_api import UniprotAPIClient


def parse_uniprot_ids(fasta_file):
    uniprot_ids = []
    with open(fasta_file, "r") as f:
        for line in f.readlines():
            if line.startswith(">"):
                uniprot_ids.append(line.split(" ")[0].replace(">", ""))
    return uniprot_ids


class AlignmentPreparation:
    def __init__(self, fasta_file, results_dir, email, pdb_csv_file):
        self.__pdb_csv_file = pdb_csv_file
        self.__email = email
        self.__results_dir = results_dir
        self.__fasta_file = fasta_file

    def run(self):
        self.__map_uniprot_to_pdb()
        cd_hit_output_file = self.__cd_hit()
        uniprot_to_entrez = UniprotAPIClient().refseq_ids(self.__uniprot_ids)
        entrez_response = EntrezDB(self.__email).fetch_cds(uniprot_to_entrez)
        if len(self.__uniprot_ids) <= 1:
            print(f"Not enough records to align {self.__uniprot_ids}")
            clustal_output = f"{self.__results_dir}/clustalo_aligned.fasta"
            copyfile(self.__fasta_file, clustal_output)
        else:
            clustal_output = self.__clustal(cd_hit_output_file)

        return NucleotideAligner().protein_based_nucleotide_alignment(
            entrez_response, clustal_output, self.__results_dir
        )

    def __cd_hit(self):
        cd_hit_output_file = f'{self.__results_dir}/cd_hit'
        subprocess.call(['cd-hit', '-i', self.__fasta_file, '-o', cd_hit_output_file, '-c', '1'])

        return cd_hit_output_file

    def __clustal(self, cd_hit_results_file):
        clustal_output = f"{self.__results_dir}/clustalo_aligned.fasta"
        subprocess.call(['clustalo', '-i', cd_hit_results_file, '-o', clustal_output, "--force"])

        return clustal_output

    def __map_uniprot_to_pdb(self):
        self.__uniprot_ids = parse_uniprot_ids(self.__fasta_file)
        ids = list(map(lambda uniprot_id: uniprot_id.split('.')[0], self.__uniprot_ids))
        with open(self.__pdb_csv_file, newline='') as csv_file:
            reader = csv.DictReader(filter(lambda csv_row: csv_row[0] != '#', csv_file))
            pdbs = []
            for row in reader:
                if row.get('SP_PRIMARY') in ids:
                    pdbs.append({'uniprot': row['SP_PRIMARY'], 'pdb': row['PDB']})
        self.__pdbs = pdbs

    @property
    def pdb_mapping(self):
        return self.__pdbs
