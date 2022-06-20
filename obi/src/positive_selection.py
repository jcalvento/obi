import json
from abc import abstractmethod
from functools import reduce

from obi.src.hyphy import Hyphy
from obi.src.logger import info
from obi.src.sifts import Sifts
from obi.src.utils import detect, get_element


class HyphyError(RuntimeError):
    def __init__(self, message):
        self.message = message


class PositiveSelectionAnalyzer:
    LOCAL_MODE = "local"
    REMOTE_MODE = "remote"

    @classmethod
    def remote(cls):
        return RemoteAnalyzer()

    @classmethod
    @abstractmethod
    def suits_mode(cls, mode):
        pass

    @classmethod
    def for_mode(cls, mode):
        return detect(lambda subclass: subclass.suits_mode(mode.lower()), cls.__subclasses__(), default=LocalAnalyzer)()

    def analyse(self, results_dir, alignment_preparation_result, email=None):
        self._validate_alignments(alignment_preparation_result, results_dir)
        return self._run_analysis(results_dir, alignment_preparation_result, email)

    @abstractmethod
    def _run_analysis(self, results_dir, alignment_preparation_result, email):
        pass

    def _generate_report(self, alignment_preparation_result, hyphy_result, pdb_mappings, results_dir):
        report = PositiveSelectionReport(alignment_preparation_result, hyphy_result, pdb_mappings).generate()

        with open(f"{results_dir}/positive_selection.json", "w") as f:
            f.write(json.dumps(report, indent=2))
        return report

    def _map_uniprot_to_pdb(self, alignment_preparation_result):
        sifts = Sifts()

        def get_pdb_data(result, item):
            uniprot_id, pdbs = item
            result[uniprot_id] = list(map(lambda pdb_id: sifts.map_to(pdb_id), pdbs))

            return result

        return reduce(get_pdb_data, alignment_preparation_result.uniprot_pdb_mapping.items(), {})

    def _validate_alignments(self, alignment_preparation_result, results_dir):
        if len(alignment_preparation_result.nucleotide_alignment) < 3:
            error_message = f"Not enough nucleotide alignments for {results_dir.split('/')[-1]}, alignments:" \
                            f" {len(alignment_preparation_result.nucleotide_alignment)}"
            info(error_message, results_dir)
            raise HyphyError(error_message)

    def _process_hyphy_result_and_generate_report(self, alignment_preparation_result, hyphy_result, results_dir):
        pdb_mappings = self._map_uniprot_to_pdb(alignment_preparation_result)
        report = self._generate_report(alignment_preparation_result, hyphy_result, pdb_mappings, results_dir)

        return report


class LocalAnalyzer(PositiveSelectionAnalyzer):
    @classmethod
    def suits_mode(cls, mode):
        return mode == cls.LOCAL_MODE

    def _run_analysis(self, results_dir, alignment_preparation_result, email):
        hyphy_result = Hyphy.local().run(alignment_preparation_result.nucleotide_alignment_path)

        return self._process_hyphy_result_and_generate_report(alignment_preparation_result, hyphy_result, results_dir)


class RemoteAnalyzer(PositiveSelectionAnalyzer):
    def __init__(self):
        self._hyphy = Hyphy.remote()

    @classmethod
    def suits_mode(cls, mode):
        return mode == cls.REMOTE_MODE

    def _run_analysis(self, results_dir, alignment_preparation_result, email):
        hyphy_result = self._hyphy.run(alignment_preparation_result.nucleotide_alignment_path, email)

        with open(f"{results_dir}/datamonkey_response.json", "w") as f:
            f.write(json.dumps(hyphy_result, indent=2))
        print(f"Job queued in Datamonkey: {hyphy_result}")

        return hyphy_result

    def resume(self, results_dir, alignment_preparation_result):
        with open(f"{results_dir}/datamonkey_response.json", "r") as f:
            job_data = json.load(f)
            hyphy_result = self._hyphy.job_result(results_dir, job_data['id'])

            return self._process_hyphy_result_and_generate_report(
                alignment_preparation_result, hyphy_result, results_dir
            )


class PositiveSelectionReport:
    def __init__(self, alignment_preparation_result, hyphy_result, pdb_mappings):
        self._pdb_mappings = pdb_mappings
        self._hyphy_result = hyphy_result
        self._alignment_preparation_result = alignment_preparation_result

    def generate(self):
        positive_selection_rows = self._select_positive_rows()
        result = {}
        for uniprot_id, alignment in self._alignment_preparation_result.nucleotide_alignment.items():
            id_prefix = uniprot_id.split('.')[0]
            if id_prefix not in self._pdb_mappings.keys():
                continue
            acc_number = self._accession_number(uniprot_id)
            rows = []
            codons = self._alignment_preparation_result.codons_and_translations[uniprot_id]['codons']
            alignment_codons = [alignment[index:index + 3] for index in range(0, len(alignment), 3)]
            translation = self._alignment_preparation_result.codons_and_translations[uniprot_id]['translation']
            amino_acid_alignment = self._alignment_preparation_result.amino_acid_alignment[uniprot_id]
            for selection_row in positive_selection_rows:
                row = self.map_data_to_row(
                    acc_number, alignment_codons, amino_acid_alignment, codons, selection_row, translation, id_prefix
                )
                rows.append(row)
            result[uniprot_id] = rows
        return result

    def map_data_to_row(self, acc_number, alignment_codons, amino_acid_alignment, codons, selection_row, translation,
                        id_prefix):
        index = selection_row['index']
        row = {
            'acc_number': acc_number,
            'index': index,
            'p-value': selection_row['p-value'],
            'codon': get_element(codons, index),
            'aa': get_element(translation, index),
            'al_codon': alignment_codons[index],
            'al_aa': amino_acid_alignment[index],
            'pdbs': []
        }
        for pdb_info in self._pdb_mappings[id_prefix]:
            residue = pdb_info[0]
            pdb_id = list(residue.keys())[0]
            chains = {}
            for chain, mapping in residue[pdb_id].items():
                chains[chain] = mapping.get(index + 1)
            row['pdbs'].append({'id': pdb_id, 'chains': chains})
        return row

    def _accession_number(self, uniprot_id):
        return detect(
            lambda mapping: uniprot_id.startswith(mapping.from_id),
            self._alignment_preparation_result.uniprot_entrez_mapping
        ).to_id.split('.')[0]

    def _select_positive_rows(self):
        positive_selection_rows = []
        for index, row in enumerate(self._hyphy_result):
            if row['p-value'] <= 0.1:
                row['index'] = index
                positive_selection_rows.append(row)
        return positive_selection_rows
