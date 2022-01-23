from .parser import parse_fasta, Record
from .aligner import ModuleSequence


def make_logo(path: str) -> None:
    fasta = parse_fasta(path)

    seqs = []
    for record in fasta:
        seqs.append(ModuleSequence(record.name, record.seq))

    for seq in seqs:
        print(seq)

