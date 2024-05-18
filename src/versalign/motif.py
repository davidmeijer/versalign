# -*- coding: utf-8 -*-

"""Implementation of a sequence motif object."""

import typing as ty
from abc import ABC, abstractmethod


class Motif(ABC):
    """Base class for sequence motif objects."""

    def __init__(self, tag: ty.Optional[int] = None) -> None:
        """Initialize the motif object.

        :param tag: The tag for the motif.
        :type tag: ty.Optional[int]
        """
        self._tag = tag

    def get_tag(self) -> ty.Optional[int]:
        """Return the tag for the motif.

        :return: The tag for the motif.
        :rtype: ty.Optional[int]
        """
        return self._tag

    def set_tag(self, tag: int) -> None:
        """Set the tag for the motif.

        :param tag: The tag for the motif.
        :type tag: int
        """
        self._tag = tag

    def clear_tag(self) -> None:
        """Clear the tag for the motif."""
        self._tag = None

    @abstractmethod
    def __eq__(self, other: ty.Any) -> bool:
        """Compare the motif to another motif.

        :param other: The other motif to compare.
        :type other: ty.Any
        :return: True if the motifs are equal, False otherwise.
        :rtype: bool
        """
        pass

    @abstractmethod
    def __str__(self) -> str:
        """Return the string representation of the motif.

        :return: The string representation of the motif.
        :rtype: str
        """
        pass


class Gap(Motif):
    """Class for representing a gap in a sequence motif."""

    def __init__(self) -> None:
        """Initialize the gap object."""
        super().__init__()  # Gap is initialized with tag=None.

    def __eq__(self, other: ty.Any) -> bool:
        """Compare the gap to another motif.

        :param other: The other motif to compare.
        :type other: Motif
        :return: True if the motifs are equal, False otherwise.
        :rtype: bool
        """
        return isinstance(other, Gap)

    def __str__(self) -> str:
        """Return the string representation of the gap.

        :return: The string representation of the gap.
        :rtype: str
        """
        return "-"
