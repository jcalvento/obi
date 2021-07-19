import copy
import re
import subprocess

from src.obi.utils import detect


def process(fasta_file_path, output_file_path):
    subprocess.call(['cd-hit', '-i', fasta_file_path, '-o', output_file_path, '-c', '1'])

    return CdHitReport(output_file_path)


def replace_cluster_heads(cd_hit_report, fasta_content, uniprot_pds_mapping, output_file_path=None):
    to_replace = copy.deepcopy(cd_hit_report.clusters)
    replaced = False
    for index, cluster in cd_hit_report.clusters.items():
        uniprot_with_mapping = detect(
            lambda el_data: uniprot_pds_mapping.get(el_data['id'].split(".")[0]) is not None,
            list(cluster.values())
        )
        cluster_head = cluster['0']
        if uniprot_with_mapping and uniprot_with_mapping != cluster_head:
            index_to_replace = list(cluster.keys())[list(cluster.values()).index(uniprot_with_mapping)]
            to_replace[index][index_to_replace], to_replace[index]["0"] = cluster_head, uniprot_with_mapping
            replaced = True
    if replaced:
        result = ""
        for index, cluster in to_replace.items():
            head_info = fasta_content[cluster['0']['id']]
            result += f"{head_info['header']}{head_info['sequence']}"
        with open(output_file_path, 'w') as f:
            f.write(result)
        with open(f"{output_file_path}.clstr", 'w') as f:
            list(map(
                lambda index_cluster: f.writelines(
                    [
                        f">Cluster {index_cluster[0]}\n",
                        *list(map(lambda cindex_cluster: f"{cindex_cluster[0]}   {' '.join(cindex_cluster[1]['line'])}\n", index_cluster[1].items()))
                    ]
                ), to_replace.items()
            ))
        return output_file_path
    return cd_hit_report.output_file


class CdHitReport:
    def __init__(self, output_file_path):
        self.output_file = output_file_path
        self.clusters = self.__parse_clusters()

    def __parse_clusters(self):
        clusters = {}
        with open(f'{self.output_file}.clstr', "r") as f:
            for line in f.readlines():
                split_line = line.split()
                if line.startswith(">"):
                    current = split_line[-1]
                else:
                    clusters.setdefault(current, {})
                    clusters[current][split_line[0]] = {
                        "id": re.search(">(.*)...", split_line[2]).group(1), "line": split_line[1:]
                    }
        return clusters
