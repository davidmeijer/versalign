# -*- coding: utf-8 -*-

"""Contains unit tests for the versalign.pairwise module."""

import typing as ty
import unittest

from versalign.motif import Gap, Motif
from versalign.pairwise import align_sequences
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


class TestPairwise(unittest.TestCase):
    """Test pairwise sequence alignment."""

    def test_align_sequence_1(self):
        """Test aligning two sequences."""
        seq_a = Sequence([A(), A(), A(), A()])
        seq_b = Sequence([B(), B(), B(), B()])
        alignment = align_sequences(seq_a, seq_b, 2, 1, score_func)
        expected = [
            (Gap(), B()),
            (Gap(), B()),
            (Gap(), B()),
            (Gap(), B()),
            (A(), Gap()),
            (A(), Gap()),
            (A(), Gap()),
            (A(), Gap()),
        ]
        self.assertEqual(alignment, expected)

    def test_align_sequence_2(self):
        """Test aligning two sequences."""
        seq_a = Sequence([A(), A(), A(), A()])
        seq_b = Sequence([A(), A(), A(), A()])
        alignment = align_sequences(seq_a, seq_b, 2, 1, score_func)
        expected = [
            (A(), A()),
            (A(), A()),
            (A(), A()),
            (A(), A()),
        ]
        self.assertEqual(alignment, expected)

    def test_align_sequence_3(self):
        """Test aligning two sequences."""
        seq_a = Sequence([A(), A(), B(), B()])
        seq_b = Sequence([B(), B(), C(), C()])
        alignment = align_sequences(seq_a, seq_b, 2, 1, score_func)
        expected = [
            (A(), Gap()),
            (A(), Gap()),
            (B(), B()),
            (B(), B()),
            (Gap(), C()),
            (Gap(), C()),
        ]
        self.assertEqual(alignment, expected)

    def test_align_sequence_4(self):
        """Test aligning two sequences."""
        seq_a = Sequence([A(), A(), A()])
        seq_b = Sequence([A(), A(), A(), A()])
        alignment = align_sequences(seq_a, seq_b, 1, 2, score_func)
        expected = [
            (Gap(), A()),
            (A(), A()),
            (A(), A()),
            (A(), A()),
        ]
        self.assertEqual(alignment, expected)
