import csv


def parse_csv():
    with open('./proteina_protmiscuity.csv') as csv_file:
        csv_reader = csv.DictReader(csv_file, delimiter=',')
        for row in csv_reader:
            with open(f"./fasta/{row['Uniprot_id']}.fasta", 'w') as f:
                f.writelines([
                        f">{row['Uniprot_id']} {row['protein_name']} {row['organism']}\n",
                        f"{row['seq']}"
                    ]
                )
