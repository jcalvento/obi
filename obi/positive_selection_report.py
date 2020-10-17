from obi.utils import detect, get_element


class PositiveSelectionReport:
    def __init__(self, alignment_preparation_result, hyphy_result, pdb_mappings):
        self._pdb_mappings = pdb_mappings
        self._hyphy_result = hyphy_result
        self._alignment_preparation_result = alignment_preparation_result

    def generate(self):
        positive_selection_rows = []
        for index, row in enumerate(self._hyphy_result):
            if row['p-value'] <= 0.1:
                row['index'] = index
                positive_selection_rows.append(row)
        result = {}
        for uniprot_id, alignment in self._alignment_preparation_result.nucleotide_alignment.items():
            uniprot_id_prefix = uniprot_id.split('.')[0]
            if uniprot_id_prefix not in self._pdb_mappings.keys():
                continue
            acc_number = detect(
                lambda mapping: uniprot_id.startswith(mapping.from_id),
                self._alignment_preparation_result.uniprot_entrez_mapping
            ).to_id.split('.')[0]
            rows = []
            codons = self._alignment_preparation_result.codons_and_translations[uniprot_id]['codons']
            alignment_codons = [alignment[index:index + 3] for index in range(0, len(alignment), 3)]
            translation = self._alignment_preparation_result.codons_and_translations[uniprot_id]['translation']
            amino_acid_alignment = self._alignment_preparation_result.amino_acid_alignment[uniprot_id]
            for selection_row in positive_selection_rows:
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
                for pdb_info in self._pdb_mappings[uniprot_id_prefix]:
                    residue = pdb_info[0]
                    pdb_id = list(residue.keys())[0]
                    chains = {}
                    for chain, mapping in residue[pdb_id].items():
                        chains[chain] = mapping.get(index + 1)
                    row['pdbs'].append({'id': pdb_id, 'chains': chains})
                rows.append(row)
            result[uniprot_id] = rows
        return result
