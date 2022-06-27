#!/usr/bin/env python3
import sys
from pathlib import Path
path_root = Path(__file__).parents[1]
sys.path.append(str(path_root))

import argparse

from monomer_aligner import cluster_polyketide_backbones


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        'fasta',
        type=str,
        help='path to polyketide backbone MSA fasta file'
    )
    parser.add_argument(
        'out_dir',
        type=str,
        help='path to output dir'
    )
    args = parser.parse_args()
    cluster_polyketide_backbones(args.fasta, args.out_dir)


if __name__ == '__main__':
    main()
