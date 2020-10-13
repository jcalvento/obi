import os
import pathlib

from obi import cd_hit
from obi.alignment_preparation import fasta_content
from obi.cd_hit import replace_cluster_heads, CdHitReport
from obi.utils import flat_map


def get_resource(file_path):
    return f"{str(pathlib.Path(__file__).parent.absolute())}/resources/{file_path}"


def results_path(file_path):
    return f"{str(pathlib.Path(__file__).parent.absolute())}/tests_result_files/{file_path}"


def uniprot_pdbs_mapping():
    return {
        'P02768': ['1ao6', '1ao6', '1bj5', '1bke'],
        'P02769': ['2l7u', '3v03', '3v03', '4f5s', '4f5s'],
        'P02771': ['3mrk'],
        'P35747': ['3v08', '4f5t', '4f5u', '4j2v', '4ot2', '4zbq'],
        'P49065': ['3v09', '6ock', '6ocl'],
        'P14639': ['4luf', '4luh', '5orf', '5orf', '5orf', '5orf', '6hn0'],
        'P49822': ['5ghk'],
        'P49064': ['5yxe']
    }


class TestCdHit:
    def test_when_processing_fasta_file_it_returns_clusters_report(self):
        fasta_file = get_resource("cd_hit_P02769_blast_result.fasta")
        fasta_file_content = fasta_content(fasta_file)
        cd_hit_output = results_path("cd_hit_P02769")

        report = cd_hit.process(fasta_file, cd_hit_output)

        report_sequences = flat_map(
            lambda x: list(map(lambda value: value['id'], x.values())), report.clusters.values()
        )
        assert all(sequence in report_sequences for sequence in fasta_file_content.keys())
        assert 19 == len(report.clusters)

    def test_when_there_is_a_cluster_which_head_has_no_pdb_mapping_but_other_in_the_cluster_does_it_replaces_it(self):
        cd_hit_report = CdHitReport(get_resource("cd_hit_to_replace"))
        fasta_file = get_resource("cd_hit_P02769_blast_result.fasta")
        fasta_file_content = fasta_content(fasta_file)
        output_file_path = results_path("cd_hit_replaced")

        replaced_output_file = replace_cluster_heads(
            cd_hit_report, fasta_file_content, uniprot_pdbs_mapping(), output_file_path
        )

        replaced_output = CdHitReport(replaced_output_file)
        assert output_file_path == replaced_output.output_file
        assert "P02768.2" == replaced_output.clusters['0']['0']['id']
        assert "A2V9Z4.1" == replaced_output.clusters['0']['1']['id']
        assert "P02769.4" == replaced_output.clusters['1']['0']['id']

    def test_when_there_is_a_cluster_with_nothing_to_change_it_returns_the_original_file_path(self):
        cd_hit_report = CdHitReport(get_resource("cd_hit_nothing_to_replace"))
        fasta_file = get_resource("cd_hit_P02769_blast_result.fasta")
        fasta_file_content = fasta_content(fasta_file)
        output_file_path = results_path("cd_hit_nothing_replaced")

        replaced_output_file = replace_cluster_heads(
            cd_hit_report, fasta_file_content, uniprot_pdbs_mapping(), output_file_path
        )

        replaced_output = CdHitReport(replaced_output_file)
        assert cd_hit_report.output_file == replaced_output.output_file
        assert not os.path.isfile(output_file_path)
