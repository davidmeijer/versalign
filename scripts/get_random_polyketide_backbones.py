#!/usr/bin/env python3
"""
NOTE: random sequences will not yield a logo
"""
from argparse import ArgumentParser
from random import randint, choices


def main() -> None:
    parser = ArgumentParser()
    parser.add_argument('--min', type=int, required=True, help='min length')
    parser.add_argument('--max', type=int, required=True, help='max length')
    parser.add_argument('--n', type=int, required=True, help='number of sequences')
    parser.add_argument('--scores', required=True, help='path to scoring matrix file')
    args = parser.parse_args()

    with open(args.scores, 'r') as handle:
        mods = handle.readline().strip().split()[:-1]  # Don't include gap

    for i in range(1, args.n + 1):
        print(f'>seq{i}')
        seq = ''.join(choices(mods, k=randint(args.min, args.max)))
        print(seq)


if __name__ == '__main__':
    main()
