from obi.src.entrez import EntrezElement
from obi.src.nucleotide_aligner import NucleotideAligner
from obi.src.utils import detect
from tests.utils import get_resource, results_path


class TestNucleotideAligner:
    def test_given_nucleotide_data_from_entrez_and_protein_alignment_it_returns_a_nucleotide_alignment_combining_both(self):
        q9lm66_data = EntrezElement(
            uniprot_id="Q9LM66", location="98..1168", locus_version="NM_101938.5",
            sequence="aaaagagcgacctctatcgagtctccagcctatataaacaggctcttctcttccatttcctcttatcaagtgacaaaagaaaaacaaacgtactcaaatggctctttcttcaccttcaagaatcctctgttttgctcttgccttatccgctgcttctctctccctctctttcgcttcttcccacgattactccatcgttggatactcccccgaggatttggaatctcatgacaaactcatagaactcttcgaaaactggatctcaaattttgagaaagcttatgaaaccgttgaagagaagtttcttaggttcgaagttttcaaggataatctaaagcacatcgatgagactaacaagaaagggaaaagctactggctcgggctcaacgagtttgcggatttgagccatgaggagttcaagaaaatgtatttagggctcaagactgatatagtgagacgcgatgaagaaagatcttacgcagagttcgcttacagggacgtcgaagctgttcctaagtctgttgactggagaaagaaaggagctgtggcggaagttaagaaccagggctcttgtggaagttgttgggcgttttcgacagtagcagctgtcgaaggtataaacaagattgtgacaggaaacttgacaacattgtcagaacaagaactcatagactgtgacacgacctacaacaatggctgcaacggtggtctcatggactatgcctttgagtacattgttaagaacggaggtctacgcaaggaagaagattatccttactctatggaagaaggaacttgcgagatgcaaaaggatgaatctgaaacagtaaccattaatggacaccaagacgtacctactaatgatgagaagagtctcttgaaggcattggctcatcagcctctcagtgtcgccattgatgcatctggtagagagttccagttctatagcggcggcgtgtttgatgggcggtgcggggttgatcttgaccacggtgtggctgcggttgggtatggatcaagcaagggttcagattacatcattgtgaagaattcttggggaccaaaatggggagaaaaaggttacatcaggctgaagaggaacactgggaaaccagagggtctctgtggaatcaacaagatggcttctttccccaccaaaactaagtgatcatcagatcattctctcttatctgcatttgatcagcatcttacagtgtcttcgttttttctcttttgatttatcaacggaagtcttgtaacaactatgtgcaataatacctatatataagaagcagagaattgacaaatttctagtaccaaccatgaagaaattaaagcagagaatattaaaaataatgattattgaagaaaagaacag",
            translation="MALSSPSRILCFALALSAASLSLSFASSHDYSIVGYSPEDLESHDKLIELFENWISNFEKAYETVEEKFLRFEVFKDNLKHIDETNKKGKSYWLGLNEFADLSHEEFKKMYLGLKTDIVRRDEERSYAEFAYRDVEAVPKSVDWRKKGAVAEVKNQGSCGSCWAFSTVAAVEGINKIVTGNLTTLSEQELIDCDTTYNNGCNGGLMDYAFEYIVKNGGLRKEEDYPYSMEEGTCEMQKDESETVTINGHQDVPTNDEKSLLKALAHQPLSVAIDASGREFQFYSGGVFDGRCGVDLDHGVAAVGYGSSKGSDYIIVKNSWGPKWGEKGYIRLKRNTGKPEGLCGINKMASFPTKTK"
        )
        entrez_response = [q9lm66_data]

        result_path = results_path("")[:-1]
        nucleotide_alignment = NucleotideAligner().protein_based_nucleotide_alignment(
            entrez_response, get_resource("protein_alignment.fasta"), result_path
        )

        assert nucleotide_alignment.nucleotide_alignment_path == f"{result_path}/nucleotide_alignment.fasta"
        self.assert_there_is_an_alignment_for_the_given_nucleotide(nucleotide_alignment, q9lm66_data)

    def test_when_there_is_one_sequence_that_is_not_the_same_but_matches_more_than_99_percent_it_includes_it_on_the_alignment(self):
        o65493_data = EntrezElement(
            uniprot_id="O65493", location="48..1115", locus_version="NM_119701.4",
            sequence="accatctcttcttcttcttcctcttctctcagtgaacaaatttggctatggctttttctgcaccatcactttccaaattctctcttttggttgccatttcagcatcagctctcctctgttgtgcttttgcccgtgatttctccattgttggatacacgccggagcatttgacaaacactgacaagcttctagagctcttcgagtcatggatgtcagaacacagcaaggcttacaaaagcgtggaggagaaggtgcacaggtttgaggttttcagagagaatctgatgcatatagaccagaggaacaatgagatcaacagttactggctcggtttgaacgagtttgcggatttgacccatgaagagttcaaaggaagatatctaggacttgcaaagccacaattctctagaaagagacagccctcagctaacttcaggtacagagatatcacggacttgcctaaatccgtagactggagaaagaaaggcgctgtggctcctgtcaaggaccagggtcaatgtggtagctgttgggcattttcaacagttgcagctgtcgaggggatcaaccagatcacaacagggaatctgagttcgctttcagagcaagaactcatagactgtgacacaactttcaacagtggctgcaatggaggtctcatggactacgcattccagtacataatttcgaccggtggtctccacaaagaagatgattacccttatctcatggaggaaggaatttgtcaagagcagaaagaggatgtggaacgtgtgacaatcagcggctacgaagatgtccctgaaaatgatgacgaaagcctggtgaaggctttagctcatcagccagtcagtgtggctattgaggcttcaggaagagacttccagttctacaaagggggagtgtttaatggtaaatgtggaacagacctagaccacggtgtggcagcggttggatatggttcatcaaagggatctgactatgttattgtcaagaactcatggggaccaagatggggagagaaagggtttattaggatgaagagaaacactggtaaaccagagggactctgtggaatcaacaagatggcctcatatcctaccaagaccaagtgatagatatacctttgttcctgcatcctacctttttttctttctttcaatgcatttctattgcacattatattgtattgatgaatgtgagattttaaaataacatggtaagagttctcttaccatttttctgactgcacatcctccaaatgctcatccagcttcaagaaacctaaaaaaacatgaccatagagctacctctctc",
            translation="MAFSAPSLSKFSLLVAISASALLCCAFARDFSIVGYTPEHLTNTDKLLELFESWMSEHSKAYKSVEEKVHRFEVFRENLMHIDQRNNEINSYWLGLNEFADLTHEEFKGRYLGLAKPQFSRKRQPSANFRYRDITDLPKSVDWRKKGAVAPVKDQGQCGSCWAFSTVAAVEGINQITTGNLSSLSEQELIDCDTTFNSGCNGGLMDYAFQYIISTGGLHKEDDYPYLMEEGICQEQKEDVERVTISGYEDVPENDDESLVKALAHQPVSVAIEASGRDFQFYKGGVFNGKCGTDLDHGVAAVGYGSSKGSDYVIVKNSWGPRWGEKGFIRMKRNTGKPEGLCGINKMASYPTKTK"
        )
        q9lm66_data = EntrezElement(
            uniprot_id="Q9LM66", location="98..1168", locus_version="NM_101938.5",
            sequence="aaaagagcgacctctatcgagtctccagcctatataaacaggctcttctcttccatttcctcttatcaagtgacaaaagaaaaacaaacgtactcaaatggctctttcttcaccttcaagaatcctctgttttgctcttgccttatccgctgcttctctctccctctctttcgcttcttcccacgattactccatcgttggatactcccccgaggatttggaatctcatgacaaactcatagaactcttcgaaaactggatctcaaattttgagaaagcttatgaaaccgttgaagagaagtttcttaggttcgaagttttcaaggataatctaaagcacatcgatgagactaacaagaaagggaaaagctactggctcgggctcaacgagtttgcggatttgagccatgaggagttcaagaaaatgtatttagggctcaagactgatatagtgagacgcgatgaagaaagatcttacgcagagttcgcttacagggacgtcgaagctgttcctaagtctgttgactggagaaagaaaggagctgtggcggaagttaagaaccagggctcttgtggaagttgttgggcgttttcgacagtagcagctgtcgaaggtataaacaagattgtgacaggaaacttgacaacattgtcagaacaagaactcatagactgtgacacgacctacaacaatggctgcaacggtggtctcatggactatgcctttgagtacattgttaagaacggaggtctacgcaaggaagaagattatccttactctatggaagaaggaacttgcgagatgcaaaaggatgaatctgaaacagtaaccattaatggacaccaagacgtacctactaatgatgagaagagtctcttgaaggcattggctcatcagcctctcagtgtcgccattgatgcatctggtagagagttccagttctatagcggcggcgtgtttgatgggcggtgcggggttgatcttgaccacggtgtggctgcggttgggtatggatcaagcaagggttcagattacatcattgtgaagaattcttggggaccaaaatggggagaaaaaggttacatcaggctgaagaggaacactgggaaaccagagggtctctgtggaatcaacaagatggcttctttccccaccaaaactaagtgatcatcagatcattctctcttatctgcatttgatcagcatcttacagtgtcttcgttttttctcttttgatttatcaacggaagtcttgtaacaactatgtgcaataatacctatatataagaagcagagaattgacaaatttctagtaccaaccatgaagaaattaaagcagagaatattaaaaataatgattattgaagaaaagaacag",
            translation="MALSSPSRILCFALALSAASLSLSFASSHDYSIVGYSPEDLESHDKLIELFENWISNFEKAYETVEEKFLRFEVFKDNLKHIDETNKKGKSYWLGLNEFADLSHEEFKKMYLGLKTDIVRRDEERSYAEFAYRDVEAVPKSVDWRKKGAVAEVKNQGSCGSCWAFSTVAAVEGINKIVTGNLTTLSEQELIDCDTTYNNGCNGGLMDYAFEYIVKNGGLRKEEDYPYSMEEGTCEMQKDESETVTINGHQDVPTNDEKSLLKALAHQPLSVAIDASGREFQFYSGGVFDGRCGVDLDHGVAAVGYGSSKGSDYIIVKNSWGPKWGEKGYIRLKRNTGKPEGLCGINKMASFPTKTK"
        )
        entrez_response = [o65493_data, q9lm66_data]

        result_path = results_path("")[:-1]
        nucleotide_alignment = NucleotideAligner().protein_based_nucleotide_alignment(
            entrez_response, get_resource("protein_alignment.fasta"), result_path
        )

        assert nucleotide_alignment.nucleotide_alignment_path == f"{result_path}/nucleotide_alignment.fasta"
        self.assert_there_is_an_alignment_for_the_given_nucleotide(nucleotide_alignment, o65493_data)
        self.assert_there_is_an_alignment_for_the_given_nucleotide(nucleotide_alignment, q9lm66_data)

    def test_when_there_is_one_sequence_that_differs_with_the_alignments_one_in_more_than_1_percent_it_ignores_it_for_the_nucleotide_alignment(self):
        different_nucleotide_data = EntrezElement(
            uniprot_id="Q9LM66", location="98..1168", locus_version="NM_101938.5",
            sequence="isdoesnotmatter",
            translation="THISISWRONG"
        )

        result_path = results_path("")[:-1]
        nucleotide_alignment = NucleotideAligner().protein_based_nucleotide_alignment(
            [different_nucleotide_data], get_resource("protein_alignment.fasta"), result_path
        )

        assert nucleotide_alignment.nucleotide_alignment_path == f"{result_path}/nucleotide_alignment.fasta"
        assert not detect(lambda key: key.startswith("Q9LM66"), nucleotide_alignment.codons_and_translations.keys())

    def test_given_an_alignment_of_two_sequences_and_its_data_it_returns_the_nucleotide_alignment(self):
        o65493_data = EntrezElement(
            uniprot_id="O65493", location="1..21", locus_version="NM_119701.4",
            sequence="atggctttttctgcaccauag",
            translation="MAFSAP"
        )
        q9lm66_data = EntrezElement(
            uniprot_id="Q9LM66", location="1..21", locus_version="NM_101938.5",
            sequence="atggctctttcttcacctugadwadwa",
            translation="MALSSP"
        )

        result_path = results_path("")[:-1]
        nucleotide_alignment = NucleotideAligner().protein_based_nucleotide_alignment(
            [o65493_data, q9lm66_data], get_resource("short_protein_alignment.fasta"), result_path
        )

        assert nucleotide_alignment.nucleotide_alignment_path == f"{result_path}/nucleotide_alignment.fasta"
        o65493_nucleotide_alignment = nucleotide_alignment.nucleotide_alignment_of(o65493_data.uniprot_id)
        o65493_aa_alignment = nucleotide_alignment.amino_acid_alignment_of(o65493_data.uniprot_id)
        assert o65493_nucleotide_alignment == "atggct------ttttctgcacca"
        assert o65493_aa_alignment == "MA--FSAP"
        q9lm66_nucleotide_alignment = nucleotide_alignment.nucleotide_alignment_of(q9lm66_data.uniprot_id)
        q9lm66_aa_alignment = nucleotide_alignment.amino_acid_alignment_of(q9lm66_data.uniprot_id)
        assert q9lm66_nucleotide_alignment == "atggct------ctttcttcacct"
        assert q9lm66_aa_alignment == "MA--LSSP"

    def assert_there_is_an_alignment_for_the_given_nucleotide(self, nucleotide_alignment, nucleotide_data):
        init, end = nucleotide_data.location.split('..')
        alignment = nucleotide_alignment.codons_and_translations_of(nucleotide_data.uniprot_id)
        assert nucleotide_data.translation == alignment['translation']
        assert len(nucleotide_data.translation) == len(alignment['codons'])
        assert all(len(codon) == 3 for codon in alignment['codons'])
        assert nucleotide_data.sequence[int(init) - 1:int(end) - 3] == "".join(alignment['codons'])
