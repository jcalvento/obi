from marshmallow import Schema, fields, post_load

from src.obi.logger import logger
from src.obi.utils import detect


class NucleotideAlignmentResultSchema(Schema):
    codons_and_translations = fields.Dict()
    nucleotide_alignment_path = fields.Str()
    nucleotide_alignment = fields.Dict()
    amino_acid_alignment = fields.Dict()

    @post_load
    def deserialize(self, data, **kwargs):
        return NucleotideAlignerResult(**data)


class NucleotideAlignerResult:
    def __init__(self, amino_acid_alignment, nucleotide_alignment, nucleotide_alignment_path, codons_and_translations):
        self.codons_and_translations = codons_and_translations
        self.nucleotide_alignment_path = nucleotide_alignment_path
        self.nucleotide_alignment = nucleotide_alignment
        self.amino_acid_alignment = amino_acid_alignment

    def codons_and_translations_of(self, uniprot_id):
        alignment_key = detect(
            lambda key: key.startswith(uniprot_id), self.codons_and_translations.keys()
        )
        return self.codons_and_translations[alignment_key]

    def nucleotide_alignment_of(self, uniprot_id):
        alignment_key = detect(
            lambda key: key.startswith(uniprot_id), self.nucleotide_alignment.keys()
        )
        return self.nucleotide_alignment[alignment_key]

    def amino_acid_alignment_of(self, uniprot_id):
        alignment_key = detect(
            lambda key: key.startswith(uniprot_id), self.amino_acid_alignment.keys()
        )
        return self.amino_acid_alignment[alignment_key]


class NucleotideAligner:
    def __init__(self):
        self.__logger = None

    def protein_based_nucleotide_alignment(self, entrez_response, protein_alignment_path, results_dir):
        self.__logger = logger(results_dir)
        alignments = self.__map_alignments_data(protein_alignment_path)
        nucleotide_alignments = self.__nucleotide_alignments(alignments, entrez_response)
        nucleotide_alignments_path = self.__write_results(nucleotide_alignments, results_dir)

        return NucleotideAlignerResult(
            alignments, nucleotide_alignments, nucleotide_alignments_path, self.__codons_and_translation
        )

    def __nucleotide_alignments(self, alignments, entrez_response):
        nucleotide_alignments = {}
        self.__codons_and_translation = {}
        for entrez_row in entrez_response:
            alignment_id = detect(lambda seq_id: seq_id.startswith(entrez_row.uniprot_id), alignments)
            alignment = alignments.get(alignment_id)
            if not alignment:
                continue  # FIXME: Mandar a entrez solo los cabeza de grupo de cd hit
            if not self.__is_the_same_sequence(alignment, entrez_row.translation):
                self.__logger.info(f"{alignment_id} differs on the sequence")
                self.__logger.info(f"Alignment: {alignment.replace('-', '')}")
                self.__logger.info(f"Translation: {entrez_row.translation}")
                self.__logger.info(f"NCBI ID: {entrez_row.locus_version}")
                continue
            adn_codons = self.__adn_codons(entrez_row)
            nucleotide_alignment = ''
            codon_index = 0
            self.__codons_and_translation[alignment_id] = {
                "codons": adn_codons, "translation": entrez_row.translation
            }
            for amino_acid in alignment:
                if amino_acid == '-':
                    nucleotide_alignment += '---'
                else:
                    nucleotide_alignment += adn_codons[codon_index]
                    codon_index += 1
            nucleotide_alignments[alignment_id] = nucleotide_alignment
        return nucleotide_alignments

    def __is_the_same_sequence(self, alignment, protein_sequence):
        alignment_sequence = alignment.replace("-", "")
        sequence_length = len(alignment_sequence)
        return sequence_length == len(protein_sequence) and\
            sum(1 for a, b in zip(alignment_sequence, protein_sequence) if a != b) * 100 / sequence_length <= 1

    def __adn_codons(self, entrez_row):
        init, end = entrez_row.location.split('..')
        adn_sequence = entrez_row.sequence[int(init) - 1:int(end) - 3]  # -3 removes stop codon
        adn_codons = [adn_sequence[index:index + 3] for index in range(0, len(adn_sequence), 3)]
        return adn_codons

    def __write_results(self, nucleotide_alignments, results_dir):
        nucleotide_fasta = open(f"{results_dir}/nucleotide_alignment.fasta", "w")
        nucleotide_fasta.write(
            '\n'.join(list(map(lambda key: f'>{key}\n{nucleotide_alignments[key]}', nucleotide_alignments))))
        nucleotide_fasta.close()
        return nucleotide_fasta.name

    def __map_alignments_data(self, protein_alignment_path):
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
        return alignments
