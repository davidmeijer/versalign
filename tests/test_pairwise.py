# -*- coding: utf-8 -*-

"""Contains unit tests for the versalign.pairwise module."""

import typing as ty
import unittest

from versalign.motif import Gap, Motif
from versalign.pairwise import PairwiseAlignment, align_pairwise
from versalign.sequence import Sequence

GLOBAL = PairwiseAlignment.NEEDLEMAN_WUNSCH


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


class TestPairwise(unittest.TestCase):
    """Test pairwise sequence alignment."""

    def test_align_sequence_1(self):
        """Test aligning two sequences."""
        seq_a = Sequence("seq_a", [A(), A(), A(), A()])
        seq_b = Sequence("seq_b", [B(), B(), B(), B()])
        result_a, result_b, _ = align_pairwise(seq_a, seq_b, 2, 1, score_func, GLOBAL)
        expected_a = [Gap(), Gap(), Gap(), Gap(), A(), A(), A(), A()]
        expected_b = [B(), B(), B(), B(), Gap(), Gap(), Gap(), Gap()]
        self.assertEqual(result_a._motifs, expected_a)
        self.assertEqual(result_b._motifs, expected_b)

    def test_align_sequence_2(self):
        """Test aligning two sequences."""
        seq_a = Sequence("seq_a", [A(), A(), A(), A()])
        seq_b = Sequence("seq_b", [A(), A(), A(), A()])
        result_a, result_b, _ = align_pairwise(seq_a, seq_b, 2, 1, score_func, GLOBAL)
        expected_a = [A(), A(), A(), A()]
        expected_b = [A(), A(), A(), A()]
        self.assertEqual(result_a._motifs, expected_a)
        self.assertEqual(result_b._motifs, expected_b)

    def test_align_sequence_3(self):
        """Test aligning two sequences."""
        seq_a = Sequence("seq_a", [A(), A(), B(), B()])
        seq_b = Sequence("seq_b", [B(), B(), C(), C()])
        result_a, result_b, _ = align_pairwise(seq_a, seq_b, 2, 1, score_func, GLOBAL)
        expected_a = [A(), A(), B(), B(), Gap(), Gap()]
        expected_b = [Gap(), Gap(), B(), B(), C(), C()]
        self.assertEqual(result_a._motifs, expected_a)
        self.assertEqual(result_b._motifs, expected_b)

    def test_align_sequence_4(self):
        """Test aligning two sequences."""
        seq_a = Sequence("seq_a", [A(), A(), A()])
        seq_b = Sequence("seq_b", [A(), A(), A(), A()])
        result_a, result_b, _ = align_pairwise(seq_a, seq_b, 1, 2, score_func, GLOBAL)
        expected_a = [Gap(), A(), A(), A()]
        expected_b = [A(), A(), A(), A()]
        self.assertEqual(result_a._motifs, expected_a)
        self.assertEqual(result_b._motifs, expected_b)
