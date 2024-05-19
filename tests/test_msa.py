# -*- coding: utf-8 -*-

"""Contains unit tests for the versalign.msa module."""

import typing as ty
import unittest

from versalign.motif import Motif
from versalign.msa import multiple_sequence_alignment, pairwise_scoring_matrix
from versalign.sequence import Sequence


class A(Motif):
    """Dummy motif."""

    def __eq__(self, other: ty.Any) -> bool:
        """Compare equality of motifs.

        :param other: The other motif to compare.
        :type other: ty.Any
        :return: True if the motifs are equal, False otherwise.
        :rtype: bool
        """
        return isinstance(other, A)

    def __str__(self) -> str:
        """Convert motif to string representation."""
        return "A"


class B(Motif):
    """Dummy motif."""

    def __eq__(self, other: ty.Any) -> bool:
        """Compare equality of motifs.

        :param other: The other motif to compare.
        :type other: ty.Any
        :return: True if the motifs are equal, False otherwise.
        :rtype: bool
        """
        return isinstance(other, B)

    def __str__(self) -> str:
        """Convert motif to string representation."""
        return "B"


class C(Motif):
    """Dummy motif."""

    def __eq__(self, other: ty.Any) -> bool:
        """Compare equality of motifs.

        :param other: The other motif to compare.
        :type other: ty.Any
        :return: True if the motifs are equal, False otherwise.
        :rtype: bool
        """
        return isinstance(other, C)

    def __str__(self) -> str:
        """Convert motif to string representation."""
        return "C"


def score_func(a: Motif, b: Motif) -> int:
    """Score function for pairwise alignment."""
    if a == b:
        return 1

    return -1


class TestPairwiseScoringMatrix(unittest.TestCase):
    """Test the pairwise_scoring_matrix function."""

    def test_pairwise_scoring_matrix_1(self) -> None:
        """Test the pairwise_scoring_matrix function."""
        seq1 = Sequence("seq_a", [A(), A(), A(), A()])
        seq2 = Sequence("seq_b", [B(), B(), B(), B()])
        seq3 = Sequence("seq_c", [C(), C(), C(), C()])
        seqs = [seq1, seq2, seq3]
        matrix = pairwise_scoring_matrix(seqs, 2, 1, score_func)
        self.assertEqual(matrix.min_value, 0.0)
        self.assertEqual(matrix.max_value, 4.0)

    def test_pairwise_scoring_matrix_2(self) -> None:
        """Test the pairwise_scoring_matrix function."""
        seq1 = Sequence("seq_a", [A(), A(), A(), A()])
        seq2 = Sequence("seq_b", [A(), A(), A(), A()])
        seq3 = Sequence("seq_c", [B(), B(), B(), B()])
        seqs = [seq1, seq2, seq3]
        matrix = pairwise_scoring_matrix(seqs, 2, 1, score_func)
        self.assertEqual(matrix.min_value, 0.0)
        self.assertEqual(matrix.max_value, 4.0)

    def test_pairwise_scoring_matrix_3(self) -> None:
        """Test the pairwise_scoring_matrix function."""
        seq1 = Sequence("seq_a", [A(), A(), A(), A()])
        seq2 = Sequence("seq_b", [A(), A(), A(), A()])
        seq3 = Sequence("seq_c", [A(), A(), A(), A()])
        seqs = [seq1, seq2, seq3]
        matrix = pairwise_scoring_matrix(seqs, 2, 1, score_func)
        self.assertEqual(matrix.min_value, 4.0)
        self.assertEqual(matrix.max_value, 4.0)


class TestMultipleSequenceAlignment(unittest.TestCase):
    """Test the multiple_sequence_alignment function."""

    def test_multiple_sequence_alignment_1(self) -> None:
        """Test the multiple_sequence_alignment function."""
        seq1 = Sequence("seq_a", [A(), A(), A(), A()])
        seq2 = Sequence("seq_b", [B(), B(), B(), B()])
        seq3 = Sequence("seq_c", [C(), C(), C(), C()])
        seqs = [seq1, seq2, seq3]
        result = multiple_sequence_alignment(seqs, 2, 1, score_func)
        self.assertEqual(len(result), 3)
        self.assertTrue(all(isinstance(seq, Sequence) for seq in result))
        self.assertTrue(all([len(seq) == 12 for seq in result]))

    def test_multiple_sequence_alignment_2(self) -> None:
        """Test the multiple_sequence_alignment function."""
        seq1 = Sequence("seq_a", [A(), A(), A(), A()])
        seq2 = Sequence("seq_b", [B(), B(), B(), B()])
        seq3 = Sequence("seq_c", [B(), B(), B(), B()])
        seqs = [seq1, seq2, seq3]
        result = multiple_sequence_alignment(seqs, 2, 1, score_func)
        self.assertEqual(len(result), 3)
        self.assertTrue(all(isinstance(seq, Sequence) for seq in result))
        self.assertTrue(all([len(seq) == 8 for seq in result]))
