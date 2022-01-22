#!/usr/bin/env python3
import sys
from pathlib import Path
path_root = Path(__file__).parents[1]
sys.path.append(str(path_root))

import argparse

from monomer_aligner import run_pairwise_alignment


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument('seq1', type=str, help='polyketide sequence 1')
    parser.add_argument('seq2', type=str, help='polyketide sequence 2')
    parser.add_argument('--gap', type=int, help='gap penalty (default=2)', required=False, default=2)
    parser.add_argument('--end', type=int, help='gap end penalty (default=2)', required=False, default=2)
    args = parser.parse_args()
    run_pairwise_alignment(args.seq1, args.seq2, gap_cost=args.gap, gap_end_cost=args.end)


if __name__ == '__main__':
    main()
