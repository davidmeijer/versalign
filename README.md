[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

---

### Polyketide subunit monomer aligner

This repository contains a small program to perform multiple 
sequence alignment on polyketide backbone structures. 

The algorithm is a naive implementation of a progressive multiple sequence
alignment, using a guide tree constructed from pairwise alignments
(Needleman-Wunsch).

### Usage

An example input file is shown in `./data/pk_test.fasta`:
```text
>seq1
A1A2A1
>seq2
A1A2A1A1
>seq3
B2
>seq4
A2B2A1A1I1
>seq5
A1A2A1
```

A multiple sequence alignment of the sequences displayed above can be performed
by running `./scripts/run_multiple_sequence_alignment ./data/pk_test.fasta > ./out/pk_test_msa.fasta` 
from root.

The output file should look like this:
```text
>seq3
--/--/B2/--/--/--
>seq4
--/A2/B2/A1/A1/I1
>seq2
A1/A2/--/A1/A1/--
>seq1
A1/A2/--/--/A1/--
>seq5
A1/A2/--/--/A1/--

```

#### Optional settings

`--gap`: gap penalty (default: 2)

`--end`: end gap penalty (default: 2)

