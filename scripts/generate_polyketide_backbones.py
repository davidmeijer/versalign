#!/usr/bin/env python3
import sys
from pathlib import Path
path_root = Path(__file__).parents[1]
sys.path.append(str(path_root))

import argparse

from monomer_aligner import generate_polyketide_backbones


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        'fasta',
        type=str,
        help='path to MSA polyketide fasta file'
    )
    parser.add_argument(
        'n',
        type=int,
        help='number of polyketide backbones to generate'
    )
    args = parser.parse_args()
    generate_polyketide_backbones(args.fasta, num=args.n)


if __name__ == '__main__':
    main()
