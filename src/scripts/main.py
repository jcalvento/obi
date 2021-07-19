import argparse
import datetime

from src.obi.blast import Blast, BlastResultsError
from src.obi.entrez import InvalidEntrezIds
from src.obi.obi_1 import Obi
from src.obi.positive_selection import HyphyError, PositiveSelectionAnalyzer
from src.obi.utils import create_results_dir, root_path

if __name__ == "__main__":
    root_path = root_path()
    parser = argparse.ArgumentParser()
    parser.add_argument('--fasta', help='Path of the Fasta file containing the problem sequence')
    parser.add_argument(
        '--email', help="Email to register while using Entrez"
    )
    parser.add_argument(
        '--output-path', help='Directory where the program creates output files. Default: ./results/<fasta-filename>'
    )
    parser.add_argument(
        '--blast',
        help='Path to Blast Swissprot DB. Default: ./swissprot/swissprot',
        default=f'{root_path}/swissprot/swissprot'
    )
    parser.add_argument(
        '--min-identity',
        help='Minimum identity percentage for Blast query expressed as float. Default: 0.4',
        default=0.4
    )
    parser.add_argument(
        '--max-evalue',
        help='Max evalue permitted for Blast results. Default: 0.005',
        default=0.005
    )
    parser.add_argument(
        '--min-coverage',
        help='Min coverage permitted for Blast results. Default: 90',
        default=90
    )
    parser.add_argument(
        '--max-gaps',
        help='Maximum number of gaps permitted for Blast results. Default: 6',
        default=6
    )
    parser.add_argument(
        '--include-analysis', help='If present, runs positive selection analysis, which includes running Hyphy and'
                                   ' getting final results. If not present, result includes until nucleotide alignment',
        action='store_true', default=False
    )
    parser.add_argument(
        '--mode', help='Only used when --include-analysis is present. '
                       '"local" if you want to run Hyphy locally, "remote" to use Datamonkey instead. Default: local',
        default=PositiveSelectionAnalyzer.LOCAL_MODE
    )
    parser.add_argument(
        "--api-key", help='Required when running in remote mode. To get one go to http://datamonkey.org/apiKey'
    )
    args = parser.parse_args()

    results_dir = args.output_path or create_results_dir(args.fasta)
    blast = Blast(args.blast)
    pdb_uniprot_mapping = f"{root_path}/pdb_chain_uniprot.csv"
    obi_1 = Obi(
        blast=blast, email=args.email, uniprot_pdb_csv_path=pdb_uniprot_mapping, results_dir=results_dir
    )

    print(f"Init time: {datetime.datetime.now()}")
    try:
        obi_1.run(
            args.fasta, args.include_analysis, min_identity=float(args.min_identity),
            max_evalue=float(args.max_evalue), min_coverage=float(args.min_coverage),
            max_gaps=int(args.max_gaps), mode=args.mode, api_key=args.api_key
        )
    except (BlastResultsError, HyphyError, InvalidEntrezIds) as e:
        print(e.message)
    print(f"Finish time: {datetime.datetime.now()}")
    print(f"Results: {results_dir}")
