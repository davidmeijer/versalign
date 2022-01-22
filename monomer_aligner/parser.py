from typing import List


class Record:
    def __init__(self, name: str, seq: str) -> None:
        self.name = name
        self.seq = seq


def parse_fasta(path: str) -> List[Record]:
    records = {}
    with open(path, 'r') as handle:
        header = None
        for line in handle:
            line = line.strip()
            if line.startswith('>'):
                header = line[1:]
                records[header] = []
            else:
                records[header].append(line)
    return [Record(header, ''.join(seqs)) for header, seqs in records.items()]
