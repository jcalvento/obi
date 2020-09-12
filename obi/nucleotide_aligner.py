from obi.utils import detect


class NucleotideAligner:
    def __init__(self):
        self.alignments = []
        self.nucleotide_alignments = []

    def protein_based_nucleotide_alignment(self, entrez_response, protein_alignment_path, results_dir):
        self.alignments = self.__map_alignments_data(protein_alignment_path)

        self.nucleotide_alignments = self.__nucleotide_alignments(self.alignments, entrez_response)

        return self.__write_results(self.nucleotide_alignments, results_dir)

    def __nucleotide_alignments(self, alignments, entrez_response):
        nucleotide_alignments = {}
        for entrez_row in entrez_response:
            alignment_id = detect(lambda seq_id: seq_id.startswith(entrez_row.uniprot_id), alignments)
            alignment = alignments[alignment_id]
            adn_codons = self.__adn_codons(entrez_row)
            nucleotide_alignment = ''
            codon_index = 0
            for amino_acid in alignment:
                if amino_acid == '-':
                    nucleotide_alignment += '---'
                else:
                    nucleotide_alignment += adn_codons[codon_index]
                    codon_index += 1
            nucleotide_alignments[alignment_id] = nucleotide_alignment
        return nucleotide_alignments

    def __adn_codons(self, entrez_row):
        init, end = entrez_row.location.split('..')
        adn_sequence = entrez_row.sequence[int(init) - 1:int(end) - 3]  # -3 removes stop codon
        # TODO: Validar largo de la cadena con el resultado del blastcmd? Alignment podria agrega cosas al final
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
