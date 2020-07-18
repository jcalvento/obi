import os

from Bio.Blast import NCBIXML


class EmptyBlastResultError(RuntimeError):
    def __init__(self, message):
        self.message = message


class BlastResultsError(RuntimeError):
    def __init__(self, message):
        self.message = message


class Blast:
    def __init__(self, db_path):
        self._db_path = db_path

    def run(self, input_fasta, results_dir):
        # Blast
        try:
            blast_result = os.popen(f'blastp -query {input_fasta} -db {self._db_path} -outfmt 5')
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
        return self.__make_fasta(results_dir, filtered_results)

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

    def __make_fasta(self, results_dir, alignments):
        ids = ",".join(list(map(lambda alignment: alignment.hit_id.split('|')[1], alignments)))
        blastdbcmd_result = os.popen(f'blastdbcmd -db {self._db_path} -entry {ids}').read()
        fasta_file = open(f"{results_dir}/blast_result.fasta", "w")
        fasta_file.write(blastdbcmd_result)
        fasta_file.close()
        return fasta_file.name
