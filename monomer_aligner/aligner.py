from __future__ import annotations
from typing import Optional, Union, List, Tuple

from scipy.spatial.distance import pdist
from scipy.cluster.hierarchy import linkage

from .scoring_matrix import PKModule, ScoringMatrix
from .alignment_matrix import AlignmentMatrix, PairwiseScoreMatrix
from .parser import parse_fasta, Record


Module = Tuple[PKModule, Union[int, None]]  # PKModule with tag


class PairwiseAlignment:
    def __init__(
        self,
        name_seq1: str,
        aligned_seq1: List[Module],
        name_seq2: str,
        aligned_seq2: List[Module],
        alignment_score: int,
        gap_cost: int,
        end_gap_cost: int
    ) -> None:
        self.name_seq1 = name_seq1
        self.seq1 = aligned_seq1
        self.name_seq2 = name_seq2
        self.seq2 = aligned_seq2
        self.score = alignment_score
        self.gap = gap_cost
        self.end_gap = end_gap_cost

    def aligned_sequences(self) -> Tuple[ModuleSequence, ModuleSequence]:
        return (
            ModuleSequence(self.name_seq1, self.seq1),
            ModuleSequence(self.name_seq2, self.seq2)
        )

    def percentage_identity(self, decimals: Optional[int] = 2):
        zipped_alignment = list(zip(self.seq1, self.seq2))
        same = sum([(ab == ba) for ((ab, _), (ba, _)) in zipped_alignment])
        identity_score = (same / len(zipped_alignment)) * 100
        return round(identity_score, decimals)

    def display(self):
        identity = self.percentage_identity()
        seq1 = "".join([m.display_in_alignment() for m, _ in self.seq1])
        seq2 = "".join([m.display_in_alignment() for m, _ in self.seq2])
        msg = (
            f'>> Polyketide backbone alignment:'
            f'\nPK1: {seq1}'
            f'\nPK2: {seq2}'
            f'\n>> Results: score: {self.score}, identity: {identity} %'
            f'\n>> Options: gap cost: {self.gap}, end gap cost: {self.end_gap}'
        )
        print(msg)


class ModuleSequence:
    def __init__(
        self,
        name: str,
        module_sequence: Union[str, List[Module]]
    ) -> None:
        self.name = name
        if isinstance(module_sequence, str):
            self._seq = self.parse(module_sequence)
        elif isinstance(module_sequence, list):
            if all([
                isinstance(m, PKModule) and
                (isinstance(m_idx, int) or m_idx is None)
                for m, m_idx in module_sequence
            ]):
                self._seq = module_sequence
            else:
                raise ValueError('modules are not of type PKModule')
        else:
            raise ValueError('module sequence is not of type str or '
                             'PKModule list')

    def __str__(self) -> str:
        return "".join([m.display_name() for m, _ in self._seq])

    def tag_idx(self) -> None:
        self._seq = [
            (module, module_idx)
            for module_idx, (module, _) in enumerate(self._seq)
        ]

    def clear_tags(self) -> None:
        self._seq = [
            (module, None)
            for module, _ in self._seq
        ]

    def insert(self, idx: int, module: Module) -> None:
        self._seq.insert(idx, module)

    def parse(
        self,
        module_sequence_string: str,
        gap: Optional[str] = '-'
    ) -> List[Module]:

        def get_polyketide_subunit(module: str, gap: Optional[str] = '-'):
            if module == gap:
                module = 'GAP'

            try:
                module = PKModule[module]
            except KeyError as err:
                print(f'{err}: unknown polyketide module in sequence')

            return module

        def starting(char: char) -> bool:
            return char.isalpha() or char == PKModule.GAP.display_name()

        module_list = []
        module = None
        tag = None
        for char in module_sequence_string:

            if starting(char) and module is None:  # For starting module
                module = char

            elif starting(char):  # New module starts with letter
                module_list.append(
                    (get_polyketide_subunit(module, gap=gap), None)
                )
                module = char

            else:
                # Chars other than starting chars are appended to existing
                # module
                module += char

        # Make sure last module is added:
        module_list.append(
            (get_polyketide_subunit(module, gap=gap), tag)
        )
        return module_list

    def alignment_matrix(
        self,
        other: ModuleSequence,
        gap_cost: int,
        end_gap_cost: int,
        scoring_matrix: str,
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
        score = ScoringMatrix(scoring_matrix)
        for col in range(1, ncols):
            for row in range(1, nrows):

                # Calculate pairwise score
                cost = score(self._seq[row - 1][0], other._seq[col - 1][0])

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
    ) -> List[Module]:
        # End stage for one seq when seq is depleted of modules
        if col == 0:
            return self._seq[:row]
        if row == 0:
            return [(PKModule.GAP, None)] * col

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
            self._seq[row - 1][0] == other._seq[col - 1][0] or
            diagonal > max([horizontal, vertical]) or
            (
                vertical - current != penalty and
                horizontal - current != penalty
            )
        ):
            return self.traceback(
                other, row - 1, col - 1, mat, gap_cost, end_gap_cost
            ) + [self._seq[row - 1]]

        if horizontal > vertical:
            return self.traceback(
                other, row, col - 1, mat, gap_cost, end_gap_cost
            ) + [(PKModule.GAP, None)]

        if vertical >= horizontal:
            return self.traceback(
                other, row - 1, col, mat, gap_cost, end_gap_cost
            ) + [self._seq[row - 1]]

    def optimal_alignment(
        self,
        other: ModuleSequence,
        gap_cost: int,
        end_gap_cost: int,
        scoring_matrix: str
    ) -> PairwiseAlignment:
        """Needleman-Wunsch pairwise sequence alignment
        """
        mat = self.alignment_matrix(other, gap_cost, end_gap_cost, scoring_matrix)
        aligned_self = self.traceback(
            other, len(self._seq), len(other._seq), mat, gap_cost, end_gap_cost
        )
        mat.transpose()
        aligned_other = other.traceback(
            self, len(other._seq), len(self._seq), mat, gap_cost, end_gap_cost
        )
        alignment_score = mat.get(-1, -1)
        alignment = PairwiseAlignment(
            self.name,
            aligned_self,
            other.name,
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
    seq1 = ModuleSequence('seq1', seq1)
    seq2 = ModuleSequence('seq1', seq2)
    alignment = seq1.optimal_alignment(seq2, gap_cost, gap_end_cost)
    alignment.display()


class MultipleSequenceAlignment:
    def __init__(
        self,
        records: List[Record],
        gap_cost: int,
        gap_end_cost: int,
        scoring_matrix: str
    ) -> None:
        self._records = [ModuleSequence(r.name, r.seq) for r in records]
        self.gap = gap_cost
        self.end = gap_end_cost
        self.sm = scoring_matrix
        self.msa = self._align()

    def _align(self) -> None:
        """Based on `Progressive Sequence Alignment as a Prerequisite to
        Correct Phylogenetic Trees` by Feng and Doolittle, 1987
        """
        def all_pairwise_scores(
            seqs1: List[ModuleSequence],
            seqs2: List[ModuleSequence]
        ) -> List[PKModule]:
            mat = PairwiseScoreMatrix(len(self._records), len(self._records))
            mat.build(0.0)
            for idx1, seq1 in enumerate(seqs1):
                for idx2, seq2 in enumerate(seqs2):
                    # Skip second calculation since matrix is mirrored around
                    # the diagonal:
                    if idx2 < idx1:
                        continue

                    if idx1 == idx2:
                        score = 0.0

                    else:
                        alignment = seq1.optimal_alignment(
                            seq2, self.gap, self.end, self.sm
                        )
                        score = 100.0 - alignment.percentage_identity()
                    mat.add(idx1, idx2, score)
                    mat.add(idx2, idx1, score)
            return mat

        if len(self._records) == 0:
            msa = {0: [ModuleSequence('', '')]}  # Return empty alignment
        elif len(self._records) == 1:
            msa = {0: self._records}  # Return seq if single
        else:
            # Identification of most closely related pair:
            mat = all_pairwise_scores(self._records, self._records)
            guide_tree = linkage(pdist(mat._matrix), method='ward')
            msa = {seq_idx: [seq] for seq_idx, seq in enumerate(self._records)}

            # Progressive insertion of neutral elements
            # (can create new gaps, but cannot remove existing gaps):
            for pair_idx, pair in enumerate(guide_tree):

                # Every pair in the guide tree can be a new pair that needs
                # seed alignment or it is a leaf connecting to existing
                # alignment
                j1, j2 = int(pair[0]), int(pair[1])
                new_idx = pair_idx + len(self._records)

                if len(msa[j1]) == 1 and len(msa[j2]) == 1:
                    seed1, seed2 = msa[j1][0], msa[j2][0]
                    alignment = seed1.optimal_alignment(
                        seed2, self.gap, self.end, self.sm
                    )
                    msa[new_idx] = list(alignment.aligned_sequences())

                elif len(msa[j1]) == 1 or len(msa[j2]) == 1:
                    if len(msa[j1]) == 1:
                        leaf, seqs = msa[j1][0], msa[j2]
                    else:
                        leaf, seqs = msa[j2][0], msa[j1]

                    # Tag already aligned sequences with original location
                    # for insertion possible gaps after annealing new seq
                    front_seq = seqs[0]
                    front_seq.tag_idx()
                    back_seq = seqs[-1]
                    back_seq.tag_idx()

                    front_alignment = front_seq.optimal_alignment(
                        leaf, self.gap, self.end, self.sm
                    )
                    back_alignment = back_seq.optimal_alignment(
                        leaf, self.gap, self.end, self.sm
                    )
                    front_score = front_alignment.percentage_identity()
                    back_score = back_alignment.percentage_identity()

                    if front_score >= back_score:
                        align_to = front_alignment
                        other_seqs = seqs[1:]
                    else:
                        align_to = back_alignment
                        other_seqs = seqs[:-1]

                    anchor, new = align_to.aligned_sequences()
                    for m_idx, (_, m_tag) in enumerate(anchor._seq):
                        if m_tag is None:  # New insertion
                            for seq in other_seqs:
                                seq.insert(m_idx, (PKModule.GAP, None))

                    for seq in other_seqs:
                        seq.clear_tags()
                    anchor.clear_tags()

                    if front_score >= back_score:
                        new_msa = [new, anchor] + other_seqs
                    else:
                        new_msa = other_seqs + [anchor, new]

                    msa[new_idx] = new_msa


                elif len(msa[j1]) != 1 and len(msa[j2]) != 1:
                    # First we need to decide which MSA comes on top. To
                    # determine this we need to score the top of msa1 with the
                    # bottom of msa2 and vica versa
                    msa1, msa2 = msa[j1], msa[j2]

                    # Tag already aligned sequences with original location
                    # for insertion possible gaps after annealing new seq
                    msa1_top, msa1_bottom = msa1[0], msa1[-1]
                    msa2_top, msa2_bottom = msa2[0], msa2[-1]
                    msa1_top.tag_idx()
                    msa1_bottom.tag_idx()
                    msa2_top.tag_idx()
                    msa2_bottom.tag_idx()

                    msa1_top_alignment = msa1_bottom.optimal_alignment(
                        msa2_top, self.gap, self.end, self.sm
                    )
                    msa2_top_alignment = msa2_bottom.optimal_alignment(
                        msa1_top, self.gap, self.end, self.sm
                    )
                    msa1_top_score = msa1_top_alignment.percentage_identity()
                    msa2_top_score = msa2_top_alignment.percentage_identity()

                    if msa1_top_score >= msa2_top_score:
                        msa1_align_to, msa2_align_to = \
                            msa1_top_alignment.aligned_sequences()
                        msa1_other_seqs = msa1[:-1]
                        msa2_other_seqs = msa2[1:]
                    else:
                        msa2_align_to, msa1_align_to = \
                            msa2_top_alignment.aligned_sequences()
                        msa1_other_seqs = msa1[1:]
                        msa2_other_seqs = msa2[:-1]

                    for m_idx, (_, m_tag) in enumerate(msa1_align_to._seq):
                        if m_tag is None:  # New insertion
                            for seq in msa1_other_seqs:
                                seq.insert(m_idx, (PKModule.GAP, None))

                    for m_idx, (_, m_tag) in enumerate(msa2_align_to._seq):
                        if m_tag is None:  # New insertion
                            for seq in msa2_other_seqs:
                                seq.insert(m_idx, (PKModule.GAP, None))

                    for seq in msa1_other_seqs:
                        seq.clear_tags()
                    for seq in msa2_other_seqs:
                        seq.clear_tags()
                    msa1_align_to.clear_tags()
                    msa2_align_to.clear_tags()

                    if msa1_top_score >= msa2_top_score:
                        new_msa = \
                            msa1_other_seqs + \
                            [msa1_align_to, msa2_align_to] + \
                            msa2_other_seqs
                    else:
                        new_msa = \
                            msa2_other_seqs + \
                            [msa2_align_to, msa1_align_to] + \
                            msa1_other_seqs

                    msa[new_idx] = new_msa

                else:
                    raise IndexError(
                        'aligned sequences from guide tree failed'
                    )

                del msa[j1]
                del msa[j2]
        return msa

    def display(self):
        # After merging all branches we should be left with a single cluster
        try:
            keys = list(self.msa.keys())
        except AttributeError as err:
            print(f'{err}: did you supply more than 1 sequence?')

        if len(keys) == 1:
            for seq in self.msa[keys[0]]:
                print(f'>{seq.name}\n{seq}')

        else:
            raise ValueError('MSA is unfinished')


def multiple_sequence_alignment(
    records: List[Record],
    gap_cost: int,
    gap_end_cost: int,
    scoring_matrix: str
) -> MultipleSequenceAlignment:
    msa = MultipleSequenceAlignment(records, gap_cost, gap_end_cost, scoring_matrix)
    return msa


def run_multiple_sequence_alignment(
    path: str,
    gap_cost: int,
    gap_end_cost: int,
    scoring_matrix: str
) -> None:
    records = parse_fasta(path)
    alignment = multiple_sequence_alignment(records, gap_cost, gap_end_cost, scoring_matrix)
    alignment.display()
