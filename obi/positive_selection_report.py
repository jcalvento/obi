from obi.utils import detect


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
        result = []
        for uniprot_id, alignment in self._alignment_preparation_result.nucleotide_alignment.items():
            acc_number = detect(
                lambda mapping: uniprot_id.startswith(mapping.from_id),
                self._alignment_preparation_result.uniprot_entrez_mapping
            ).to_id
            rows = []
            codons = [alignment[index:index + 3] for index in range(0, len(alignment), 3)]
            for selection_row in positive_selection_rows:
                row = {
                    'acc_number': acc_number,
                    'index': selection_row['index'],
                    'p-value': selection_row['p-value'],
                    'codon': codons[selection_row['index']],  # Creo que no se necesita del alineamiento si no de la cadena
                    'aa': self._alignment_preparation_result.amino_acid_alignment[uniprot_id][selection_row['index']]
                }
                rows.append(row)
            result.append({uniprot_id: rows})
        return result
