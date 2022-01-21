from __future__ import annotations
from typing import Optional, List

from .scoring_matrix import PKModule, ScoringMatrix
from .alignment_matrix import AlignmentMatrix


class Alignment:
    def __init__(
        self,
        aligned_seq1: List[PKModule],
        aligned_seq2: List[PKModule],
        alignment_score: int,
        gap_cost: int,
        end_gap_cost: int
    ) -> None:
        self.seq1 = aligned_seq1
        self.seq2 = aligned_seq2
        self.score = alignment_score
        self.gap = gap_cost
        self.end_gap = end_gap_cost

    def percentage_identity(self, decimals: Optional[int] = 2):
        zipped_alignment = list(zip(self.seq1, self.seq2))
        same = sum([(ab == ba) for (ab, ba) in zipped_alignment])
        identity_score = (same / len(zipped_alignment)) * 100
        return round(identity_score, decimals)

    def display(self):
        identity = self.percentage_identity()
        seq1 = " ".join([m.display_in_alignment() for m in self.seq1])
        seq2 = " ".join([m.display_in_alignment() for m in self.seq2])
        msg = (
            f'>> Polyketide backbone alignment:'
            f'\n>> PK1: {seq1}'
            f'\n>> PK2: {seq2}'
            f'\n>> Results: score: {self.score}, identity: {identity} %'
            f'\n>> Options: gap cost: {self.gap}, end gap cost: {self.end_gap}'
        )
        print(msg)


class ModuleSequence:
    def __init__(self, module_sequence_string: str) -> None:
        self._seq = self.parse(module_sequence_string)

    def __str__(self) -> str:
        return ' '.join([m.display_name() for m in self._seq])

    def parse(self, module_sequence_string: str) -> List[PKModule]:
        module_list = []
        module = None
        for char in module_sequence_string:
            if char.isalpha() and not module:  # For starting module
                module = char
            elif char.isalpha():  # New module starts with letter
                module_list.append(PKModule[module])
                module = char
            else:  # Chars other than letters are assigned to existing module
                module += char
        module_list.append(PKModule[module])  # Make sure last module is added
        return module_list

    def alignment_matrix(
        self,
        other: ModuleSequence,
        gap_cost: int,
        end_gap_cost: int
    ) -> AlignmentMatrix:
        # Instantiate zero matrix
        nrows = len(self._seq) + 1
        ncols = len(other._seq) + 1
        mat = AlignmentMatrix(nrows, ncols)
        mat.build(0)

        # Fill in initial alignment for when there is zero alignment
        for row in range(1, nrows):
            mat.add(row, 0, (mat.get(row - 1, 0) - end_gap_cost))

        for col in range(1, ncols):
            mat.add(0, col, (mat.get(0, col - 1) - end_gap_cost))

        # Go through matrix consecutively and calculate the scores
        score = ScoringMatrix()
        for col in range(1, ncols):
            for row in range(1, nrows):

                # Calculate pairwise score
                cost = score(self._seq[row - 1], other._seq[col - 1])

                # Calculate penalty score
                if col == len(other._seq) or row == len(self._seq):
                    penalty = end_gap_cost
                else:
                    penalty = gap_cost

                # Calculate final score and fill in
                mat.add(row, col, max([
                    mat.get(row - 1, col) - penalty,
                    mat.get(row, col - 1) - penalty,
                    mat.get(row - 1, col - 1) + cost
                ]))
        return mat

    def traceback(
        self,
        other: ModuleSequence,
        row: int,
        col: int,
        mat: AlignmentMatrix,
        gap_cost: int,
        end_gap_cost: int
    ) -> List[PKModule]:
        # End stage for one seq when seq is depleted of modules
        if col == 0:
            return self._seq[:row]
        if row == 0:
            return [PKModule.GAP] * col

        # Gap penalty depends on location in matrix
        if row == len(self._seq) or col == len(other._seq):
            penalty = end_gap_cost
        else:
            penalty = gap_cost

        # Define current location and directions
        current = mat.get(row, col)
        horizontal = mat.get(row, col - 1)
        diagonal = mat.get(row - 1, col - 1)
        vertical = mat.get(row - 1, col)

        # Traceback defined by conditions for every direction
        if (
            self._seq[row - 1] == other._seq[col - 1] or
            diagonal > max([horizontal, vertical]) or
            (
                vertical - current != penalty and
                horizontal - current != penalty
            )
        ):
            return self.traceback(
                other,
                row - 1,
                col - 1,
                mat,
                gap_cost,
                end_gap_cost
            ) + [self._seq[row - 1]]

        if horizontal > vertical:
            return self.traceback(
                other,
                row,
                col - 1,
                mat,
                gap_cost,
                end_gap_cost
            ) + [PKModule.GAP]

        if vertical >= horizontal:
            return self.traceback(
                other,
                row - 1,
                col,
                mat,
                gap_cost,
                end_gap_cost
            ) + [self._seq[row - 1]]

    def optimal_alignment(
        self,
        other: ModuleSequence,
        gap_cost: int,
        end_gap_cost: int
    ) -> Alignment:
        mat = self.alignment_matrix(other, gap_cost, end_gap_cost)
        aligned_self = self.traceback(
            other,
            len(self._seq),
            len(other._seq),
            mat,
            gap_cost,
            end_gap_cost
        )
        mat.transpose()
        aligned_other = other.traceback(
            self,
            len(other._seq),
            len(self._seq),
            mat,
            gap_cost,
            end_gap_cost
        )
        alignment_score = mat.get(-1, -1)
        alignment = Alignment(
            aligned_self,
            aligned_other,
            alignment_score,
            gap_cost,
            end_gap_cost
        )
        return alignment


def run_pairwise_alignment(
    seq1: str,
    seq2: str,
    gap_cost: int,
    gap_end_cost: int
) -> None:
    seq1 = ModuleSequence(seq1)
    seq2 = ModuleSequence(seq2)
    alignment = seq1.optimal_alignment(seq2, gap_cost, gap_end_cost)
    alignment.display()
