#!/usr/bin/env python3
import sys
from pathlib import Path
path_root = Path(__file__).parents[1]
sys.path.append(str(path_root))

import argparse

from monomer_aligner import run_multiple_sequence_alignment


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument('--fasta', required=True, help='path to polyketide fasta file')
    parser.add_argument('--scores', required=True, help='path to scoring matrix file')
    args = parser.parse_args()
    run_multiple_sequence_alignment(
        args.fasta,
        gap_cost=2,
        gap_end_cost=2,
        scoring_matrix=args.scores
    )


if __name__ == '__main__':
    main()
