import os

from Bio.Blast import NCBIXML


class EmptyBlastResultError(RuntimeError):
    def __init__(self, message):
        self.message = message


class BlastResultsError(RuntimeError):
    def __init__(self, message):
        self.message = message


class Blast:
    def run(self, input_fasta, db_path):
        # Blast
        try:
            blast_result = os.popen(f'blastp -query {input_fasta} -db {db_path} -outfmt 5')
            blast_records = NCBIXML.read(blast_result)
        except ValueError:
            raise EmptyBlastResultError(f"Couldn't find any blast result for file {input_fasta}")
        sequence = []
        with open(input_fasta, "r") as f:
            for line in f.readlines():
                if not line.startswith(">"):
                    sequence.append(line)
        protein_chain = "".join(sequence).replace("\n", "")

        # Mejores hits segun criterio
        filtered_results = list(filter(
            lambda alignment: self.__meets_expected_criteria(protein_chain, alignment), blast_records.alignments
        ))

        if not filtered_results:
            raise BlastResultsError("There is none blast result that meets filtering criteria")
        # Crear Fasta con hits
        return self.__make_fasta(input_fasta.split('/')[-1].split('.')[0], filtered_results)

    def __evalue(self, alignment):
        hsp = alignment.hsps[0]
        if hasattr(hsp, 'evalue'):
            return hsp.evalue
        else:
            return 0

    def __identity_percentage(self, alignment):
        hsp = alignment.hsps[0]
        return round(hsp.identities / hsp.align_length, 2)

    def __coverage(self, protein_chain, alignment):
        hsp = alignment.hsps[0]
        return ((hsp.query_end - hsp.query_start) / len(protein_chain)) * 100

    def __meets_expected_criteria(self, protein_chain, alignment):
        return self.__identity_percentage(alignment) >= 0.4 and self.__evalue(alignment) <= 0.005 and \
               self.__coverage(protein_chain, alignment) >= 90

    # parametrizable pero defaults:
    # porcentaje de identidad > 40%
    # coverage > 90% ((dfblst['q.end']-dfblst['q.start'])/query_length)*100.0
    # evalue < 0.005

    def __fasta_info(self, alignment):
        accession = alignment.hit_id.split('|')[1]
        gen_name = alignment.hit_def.split(';')[0].split('RecName: Full=')[-1]
        raw_species = alignment.hit_def.split(';')[-1]
        species = raw_species[raw_species.find("[") + 1:raw_species.find("]")]
        header = f">{accession} {species}, gene for {gen_name}"
        sequence = alignment.hsps[0].sbjct

        return f"{header}\n{sequence}"

    def __make_fasta(self, filename, alignments):
        content = "\n\n".join(list(map(lambda alignment: self.__fasta_info(alignment), alignments)))
        fasta_file = open(f"results/{filename}_blast_result.fasta", "w")
        fasta_file.write(content)
        fasta_file.close()
        return fasta_file.name