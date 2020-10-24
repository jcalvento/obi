import os

from Bio.Blast import NCBIXML

from obi.logger import logger


class BlastResultsError(RuntimeError):
    def __init__(self, message):
        self.message = message


class Blast:
    def __init__(self, db_path):
        self._db_path = db_path

    def run(self, input_fasta, results_dir, min_identity=0.4, max_evalue=0.005, min_coverage=90, max_gaps=6):
        blast_records = self.__blast_records(input_fasta)

        protein_sequence = self.__query_sequence(input_fasta)

        filtered_results = self.__filter_results(blast_records, protein_sequence, logger(results_dir), min_identity, max_evalue, min_coverage, max_gaps)

        return self.__make_fasta(results_dir, filtered_results)

    def __filter_results(self, blast_records, protein_sequence, obi_logger, min_identity, max_evalue, min_coverage, max_gaps):
        filtered_results = list(filter(
            lambda alignment: self.__meets_expected_criteria(protein_sequence, alignment, min_identity, max_evalue, min_coverage, max_gaps), blast_records.alignments
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

    def _evalue(self, alignment):
        hsp = alignment.hsps[0]
        return hsp.expect

    def _identity_percentage(self, alignment):
        hsp = alignment.hsps[0]
        return round(hsp.identities / hsp.align_length, 2)

    def _coverage(self, protein_chain, alignment):
        hsp = alignment.hsps[0]
        return ((hsp.query_end - hsp.query_start) / len(protein_chain)) * 100

    def _gaps_percentage(self, protein_chain, alignment):
        hsp = alignment.hsps[0]
        return (hsp.gaps / len(protein_chain)) * 100

    def __meets_expected_criteria(self, protein_chain, alignment, min_identity, max_evalue, min_coverage, max_gaps):
        return self._identity_percentage(alignment) >= min_identity and self._evalue(alignment) <= max_evalue and \
               self._coverage(protein_chain, alignment) >= min_coverage and\
               self._gaps_percentage(protein_chain, alignment) <= max_gaps

    def __make_fasta(self, results_dir, alignments):
        ids = ",".join(list(map(lambda alignment: alignment.hit_id.split('|')[1], alignments)))
        blastdbcmd_result = self.__blastdbcmd(ids)
        fasta_file = open(f"{results_dir}/blast_result.fasta", "w")
        fasta_file.write(blastdbcmd_result)
        fasta_file.close()
        return fasta_file.name

    def __blastdbcmd(self, ids):
        return os.popen(f'blastdbcmd -db {self._db_path} -entry {ids}').read()
