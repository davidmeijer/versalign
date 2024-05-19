# -*- coding: utf-8 -*-

"""Implementation of sequence object comprising a list of motifs."""

import typing as ty

from .motif import Motif


class Sequence:
    """Class for representing a sequence object."""

    def __init__(self, motifs: ty.Optional[ty.List[Motif]] = None) -> None:
        """Initialize the sequence object.

        :param motifs: The list of motifs in the sequence.
        :type motifs: ty.Optional[ty.List[Motif]]
        :raises ValueError: If any element in the list is not a motif.
        """
        if motifs is None:
            motifs = []

        if not all(isinstance(motif, Motif) for motif in motifs):
            raise ValueError("All elements in the list must be motifs.")

        self._motifs = motifs

    def __len__(self) -> int:
        """Return the length of the sequence.

        :return: The length of the sequence.
        :rtype: int
        """
        return len(self._motifs)

    def __getitem__(self, index: int) -> Motif:
        """Return the motif at the specified index.

        :param index: The index of the motif to return.
        :type index: int
        :return: The motif at the specified index.
        :rtype: Motif
        """
        return self._motifs[index]

    def __setitem__(self, index: int, motif: Motif) -> None:
        """Set the motif at the specified index.

        :param index: The index of the motif to set.
        :type index: int
        :param motif: The motif to set.
        :type motif: Motif
        :raises ValueError: If the motif is not a motif object.
        """
        if not isinstance(motif, Motif):
            raise ValueError("The motif must be a motif object.")

        self._motifs[index] = motif

    def __iter__(self) -> ty.Iterator[Motif]:
        """Return an iterator over the motifs in the sequence.

        :return: An iterator over the motifs.
        :rtype: ty.Iterator[Motif]
        """
        return iter(self._motifs)

    def __str__(self) -> str:
        """Convert the sequence to a string representation.

        :return: The string representation of the sequence.
        :rtype: str
        """
        return "".join(str(motif) for motif in self._motifs)

    def insert(self, index: int, motif: Motif) -> None:
        """Insert the motif at the specified index.

        :param index: The index to insert the motif.
        :type index: int
        :param motif: The motif to insert.
        :type motif: Motif
        :raises ValueError: If the motif is not a motif object.
        """
        if not isinstance(motif, Motif):
            raise ValueError("The motif must be a motif object.")

        self._motifs.insert(index, motif)

    def tag(self) -> None:
        """Tag all motifs in the sequence."""
        for motif_idx, motif in enumerate(self._motifs):
            motif.set_tag(motif_idx)

    def clear_tags(self) -> None:
        """Clear the tags for all motifs in the sequence."""
        for motif in self._motifs:
            motif.clear_tag()
