import os

from Bio.Blast import NCBIXML

from obi.logger import logger


class BlastResultsError(RuntimeError):
    def __init__(self, message):
        self.message = message


class Blast:
    def __init__(self, db_path):
        self._db_path = db_path

    def run(self, input_fasta, results_dir):
        blast_records = self.__blast_records(input_fasta)

        protein_sequence = self.__query_sequence(input_fasta)

        filtered_results = self.__filter_results(blast_records, protein_sequence, logger(results_dir))

        return self.__make_fasta(results_dir, filtered_results)

    def __filter_results(self, blast_records, protein_sequence, obi_logger):
        filtered_results = list(filter(
            lambda alignment: self.__meets_expected_criteria(protein_sequence, alignment), blast_records.alignments
        ))
        if not filtered_results:
            error_message = "There is none blast result that meets filtering criteria"
            obi_logger.info(error_message)
            raise BlastResultsError(error_message)

        return filtered_results

    def __query_sequence(self, input_fasta):
        sequence = []
        with open(input_fasta, "r") as f:
            for line in f.readlines():
                if not line.startswith(">"):
                    sequence.append(line)

        return "".join(sequence).replace("\n", "")

    def __blast_records(self, input_fasta):
        blast_result = self.__blastp(input_fasta)
        blast_records = NCBIXML.read(blast_result)
        return blast_records

    def __blastp(self, input_fasta):
        return os.popen(f'blastp -query {input_fasta} -db {self._db_path} -outfmt 5')

    def __evalue(self, alignment):
        hsp = alignment.hsps[0]
        return hsp.expect

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
    # gaps < 6%

    def __make_fasta(self, results_dir, alignments):
        ids = ",".join(list(map(lambda alignment: alignment.hit_id.split('|')[1], alignments)))
        blastdbcmd_result = self.__blastdbcmd(ids)
        fasta_file = open(f"{results_dir}/blast_result.fasta", "w")
        fasta_file.write(blastdbcmd_result)
        fasta_file.close()
        return fasta_file.name

    def __blastdbcmd(self, ids):
        return os.popen(f'blastdbcmd -db {self._db_path} -entry {ids}').read()
