#!/usr/bin/env python3
import sys
from pathlib import Path
path_root = Path(__file__).parents[1]
sys.path.append(str(path_root))

from monomer_aligner import run_pairwise_alignment


def main() -> None:
    seq1 = sys.argv[1]
    seq2 = sys.argv[2]
    run_pairwise_alignment(seq1, seq2, gap_cost=2, gap_end_cost=2)


if __name__ == '__main__':
    main()
