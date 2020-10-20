import csv
import subprocess
from shutil import copyfile

from obi import cd_hit
from obi.cd_hit import replace_cluster_heads
from obi.entrez import EntrezDB
from obi.logger import info
from obi.nucleotide_aligner import NucleotideAligner
from obi.uniprot_api import UniprotAPIClient


def parse_uniprot_ids(fasta_file):
    uniprot_ids = []
    with open(fasta_file, "r") as f:
        for line in f.readlines():
            if line.startswith(">"):
                uniprot_ids.append(line.split(" ")[0].replace(">", ""))
    return uniprot_ids


def fasta_content(fasta_file):
    content = {}
    with open(fasta_file, "r") as f:
        current_id = None
        for line in f.readlines():
            if line.startswith(">"):
                current_id = line.split(" ")[0].replace(">", "")
                content[current_id] = {"header": line}
            else:
                content[current_id].setdefault("sequence", "")
                content[current_id]["sequence"] += line
    return content


class AlignmentPreparationResult:
    def __init__(self, nucleotide_alignment_result, uniprot_pdb_mapping, uniprot_entrez_mapping):
        self.uniprot_entrez_mapping = uniprot_entrez_mapping
        self.uniprot_pdb_mapping = uniprot_pdb_mapping
        self.__nucleotide_alignment_result = nucleotide_alignment_result

    @property
    def nucleotide_alignment(self):
        return self.__nucleotide_alignment_result.nucleotide_alignment

    @property
    def amino_acid_alignment(self):
        return self.__nucleotide_alignment_result.amino_acid_alignment

    @property
    def codons_and_translations(self):
        return self.__nucleotide_alignment_result.codons_and_translations

    @property
    def nucleotide_alignment_path(self):
        return self.__nucleotide_alignment_result.nucleotide_alignment_path


class AlignmentPreparation:
    def __init__(self, fasta_file, results_dir, email, pdb_csv_file):
        self.__pdb_csv_file = pdb_csv_file
        self.__email = email
        self.__results_dir = results_dir
        self.__fasta_file = fasta_file
        self.__fasta_content = fasta_content(self.__fasta_file)

    def run(self):
        uniprot_pds_mapping = self.__map_uniprot_to_pdb()
        cd_hit_output = self.__cd_hit(uniprot_pds_mapping)
        uniprot_entrez_mapping = UniprotAPIClient().refseq_ids(self.__uniprot_ids)
        entrez_response = EntrezDB(self.__email).fetch_cds(uniprot_entrez_mapping)
        clustal_output = self.__amino_acids_alignment(cd_hit_output)

        nucleotide_alignment_result = NucleotideAligner().protein_based_nucleotide_alignment(
            entrez_response, clustal_output, self.__results_dir
        )

        return AlignmentPreparationResult(nucleotide_alignment_result, uniprot_pds_mapping, uniprot_entrez_mapping)

    def __cd_hit(self, uniprot_pds_mapping):
        cd_hit_output_file = f'{self.__results_dir}/cd_hit'
        cd_hit_output = cd_hit.process(self.__fasta_file, cd_hit_output_file)

        return replace_cluster_heads(
            cd_hit_output, self.__fasta_content, uniprot_pds_mapping, f"{cd_hit_output.output_file}_replaced"
        )

    def __amino_acids_alignment(self, cd_hit_output_file):
        if len(self.__uniprot_ids) <= 1:
            info(f"Not enough records to align {self.__uniprot_ids}", self.__results_dir)
            clustal_output = f"{self.__results_dir}/clustalo_aligned.fasta"
            copyfile(self.__fasta_file, clustal_output)
        else:
            clustal_output = self.__clustal(cd_hit_output_file)
        return clustal_output

    def __clustal(self, cd_hit_results_file):
        clustal_output = f"{self.__results_dir}/clustalo_aligned.fasta"
        subprocess.call(['clustalo', '-i', cd_hit_results_file, '-o', clustal_output, "--force"])

        return clustal_output

    def __map_uniprot_to_pdb(self):
        self.__uniprot_ids = parse_uniprot_ids(self.__fasta_file)
        ids = list(map(lambda uniprot_id: uniprot_id.split('.')[0], self.__uniprot_ids))
        with open(self.__pdb_csv_file, newline='') as csv_file:
            reader = csv.DictReader(filter(lambda csv_row: csv_row[0] != '#', csv_file))
            mapping = {}
            for row in reader:
                if row.get('SP_PRIMARY') in ids:
                    if not mapping.get(row['SP_PRIMARY']):
                        mapping[row['SP_PRIMARY']] = [row['PDB']]
                    else:
                        mapping[row['SP_PRIMARY']].append(row['PDB'])
        return mapping
