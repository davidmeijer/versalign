#!/usr/bin/env python3
import sys
from pathlib import Path
path_root = Path(__file__).parents[1]
sys.path.append(str(path_root))

import argparse

from monomer_aligner import make_logo


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument('fasta', type=str, help='path to MSA polyketide fasta file')
    parser.add_argument('--out', type=str, help='output path of logo', required=False, default=None)
    args = parser.parse_args()
    make_logo(args.fasta, args.out)


if __name__ == '__main__':
    main()
