from pytest import fixture

from obi.src.entrez import EntrezDB
from obi.src.uniprot_api import UniprotIdMapping


class TestEntrezDB:
    @fixture
    def entrez_api(self):
        return EntrezDB('juliancalvento@gmail.com')

    def test_given_a_uniprot_refseq_mapping_when_it_fetches_the_cds_returns_the_location_transation_and_sequence(self, entrez_api):
        mapping = UniprotIdMapping.for_tab_separated_ids('O65493\tNM_119701.4')

        entrez_response = entrez_api.fetch_cds([mapping])

        entrez_mapping = entrez_response[0]
        sequence = 'accatctcttcttcttcttcctcttctctcagtgaacaaatttggctatggctttttctgcaccatcactttccaaattctctcttttggttgccat' \
                   'ttcagcatcagctctcctctgttgtgcttttgcccgtgatttctccattgttggatacacgccggagcatttgacaaacactgacaagcttctagag' \
                   'ctcttcgagtcatggatgtcagaacacagcaaggcttacaaaagcgtggaggagaaggtgcacaggtttgaggttttcagagagaatctgatgcata' \
                   'tagaccagaggaacaatgagatcaacagttactggctcggtttgaacgagtttgcggatttgacccatgaagagttcaaaggaagatatctaggact' \
                   'tgcaaagccacaattctctagaaagagacagccctcagctaacttcaggtacagagatatcacggacttgcctaaatccgtagactggagaaagaaa' \
                   'ggcgctgtggctcctgtcaaggaccagggtcaatgtggtagctgttgggcattttcaacagttgcagctgtcgaggggatcaaccagatcacaacag' \
                   'ggaatctgagttcgctttcagagcaagaactcatagactgtgacacaactttcaacagtggctgcaatggaggtctcatggactacgcattccagta' \
                   'cataatttcgaccggtggtctccacaaagaagatgattacccttatctcatggaggaaggaatttgtcaagagcagaaagaggatgtggaacgtgtg' \
                   'acaatcagcggctacgaagatgtccctgaaaatgatgacgaaagcctggtgaaggctttagctcatcagccagtcagtgtggctattgaggcttcag' \
                   'gaagagacttccagttctacaaagggggagtgtttaatggtaaatgtggaacagacctagaccacggtgtggcagcggttggatatggttcatcaaa' \
                   'gggatctgactatgttattgtcaagaactcatggggaccaagatggggagagaaagggtttattaggatgaagagaaacactggtaaaccagaggga' \
                   'ctctgtggaatcaacaagatggcctcatatcctaccaagaccaagtgatagatatacctttgttcctgcatcctacctttttttctttctttcaatg' \
                   'catttctattgcacattatattgtattgatgaatgtgagattttaaaataacatggtaagagttctcttaccatttttctgactgcacatcctccaa' \
                   'atgctcatccagcttcaagaaacctaaaaaaacatgaccatagagctacctctctc'
        translation = 'MAFSAPSLSKFSLLVAISASALLCCAFARDFSIVGYTPEHLTNTDKLLELFESWMSEHSKAYKSVEEKVHRFEVFRENLMHIDQRNNEINSYWL' \
                      'GLNEFADLTHEEFKGRYLGLAKPQFSRKRQPSANFRYRDITDLPKSVDWRKKGAVAPVKDQGQCGSCWAFSTVAAVEGINQITTGNLSSLSEQE' \
                      'LIDCDTTFNSGCNGGLMDYAFQYIISTGGLHKEDDYPYLMEEGICQEQKEDVERVTISGYEDVPENDDESLVKALAHQPVSVAIEASGRDFQFY' \
                      'KGGVFNGKCGTDLDHGVAAVGYGSSKGSDYVIVKNSWGPRWGEKGFIRMKRNTGKPEGLCGINKMASYPTKTK'
        assert 1 == len(entrez_response)
        assert 'O65493' == entrez_mapping.uniprot_id
        assert sequence == entrez_mapping.sequence
        assert translation == entrez_mapping.translation
        assert '48..1115' == entrez_mapping.location
        assert 'NM_119701.4' == entrez_mapping.locus_version

    def test_given_two_uniprot_refseq_mapping_when_it_fetches_the_cds_returns_results_for_both(self, entrez_api):
        mapping = UniprotIdMapping.for_tab_separated_ids('O65493\tNM_119701.4')
        another_mapping = UniprotIdMapping.for_tab_separated_ids('Q9LM66\tNM_101938.5')

        entrez_response = entrez_api.fetch_cds([mapping, another_mapping])

        entrez_mapping = entrez_response[0]
        another_entrez_mapping = entrez_response[1]
        assert 2 == len(entrez_response)
        assert another_entrez_mapping.uniprot_id != entrez_mapping.uniprot_id
        assert another_entrez_mapping.sequence != entrez_mapping.sequence
        assert another_entrez_mapping.translation != entrez_mapping.translation
        assert another_entrez_mapping.location != entrez_mapping.location
        assert another_entrez_mapping.locus_version != entrez_mapping.locus_version
