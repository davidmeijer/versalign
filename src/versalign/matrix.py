# -*- coding: utf-8 -*-

"""Implementation of the scoring matrix for sequence alignment."""

import logging
import typing as ty

import numpy as np


class Matrix:
    """Base class for matrix objects."""

    def __init__(self, nrows: int, ncols: int, fill: float) -> None:
        """Initialize the matrix object.

        :param nrows: The number of rows in the matrix.
        :type nrows: int
        :param ncols: The number of columns in the matrix.
        :type ncols: int
        :param fill: The value to fill the matrix with.
        :type fill: float
        """
        self._nrows = nrows
        self._ncols = ncols
        self._matrix = self._create_matrix(fill)

    @property
    def nrows(self) -> int:
        """Return the number of rows in the matrix.

        :return: The number of rows in the matrix.
        :rtype: int
        """
        return self._nrows

    @property
    def ncols(self) -> int:
        """Return the number of columns in the matrix.

        :return: The number of columns in the matrix.
        :rtype: int
        """
        return self._ncols

    @property
    def max_value(self) -> float:
        """Return the maximum value in the matrix.

        :return: The maximum value in the matrix.
        :rtype: float
        """
        return np.max(self._matrix)

    @property
    def min_value(self) -> float:
        """Return the minimum value in the matrix.

        :return: The minimum value in the matrix.
        :rtype: float
        """
        return np.min(self._matrix)

    def _create_matrix(self, fill: float) -> np.ndarray:
        """Create the matrix object with NumPy.

        :param fill: The value to fill the matrix with.
        :type fill: float
        :return: The matrix object.
        :rtype: np.ndarray
        :raises ValueError: If the fill value is not an integer or float.
        """
        logger = logging.getLogger(__name__)

        if isinstance(fill, float):
            return np.full((self._nrows, self._ncols), fill, dtype=np.float32)

        else:
            msg = "Fill value must be a float."
            logger.error(msg)
            raise ValueError(msg)

    def transpose(self) -> np.ndarray:
        """Transpose the matrix.

        :return: The transposed matrix.
        :rtype: np.ndarray
        """
        return self._matrix.T

    def set_value(self, row: int, col: int, value: ty.Union[int, float]) -> None:
        """Set the value of a cell in the matrix.

        :param row: The row index.
        :type row: int
        :param col: The column index.
        :type col: int
        :param value: The value to set.
        :type value: Union[int, float]
        :raises ValueError: If the row or column index is out of bounds.
        :raises ValueError: If the value is not an integer.
        """
        logger = logging.getLogger(__name__)

        if not (0 <= row < self._nrows):
            msg = "Row index out of bounds."
            logger.error(msg)
            raise ValueError(msg)

        if not (0 <= col < self._ncols):
            msg = "Column index out of bounds."
            logger.error(msg)
            raise ValueError(msg)

        if (
            not isinstance(value, int)
            and not isinstance(value, float)
            and not isinstance(value, np.float32)
        ):
            msg = "Value must be an integer or float."
            logger.error(msg)
            raise ValueError(msg)

        if isinstance(value, int):
            value = float(value)

        self._matrix[row, col] = value

    def get_value(self, row: int, col: int) -> float:
        """Get the value of a cell in the matrix.

        :param row: The row index.
        :type row: int
        :param col: The column index.
        :type col: int
        :return: The value of the cell.
        :rtype: float
        :raises ValueError: If the row or column index is out of bounds.
        """
        logger = logging.getLogger(__name__)

        if not (0 <= row < self._nrows):
            msg = "Row index out of bounds."
            logger.error(msg)
            raise ValueError(msg)

        if not (0 <= col < self._ncols):
            msg = "Column index out of bounds."
            logger.error(msg)
            raise ValueError(msg)

        return self._matrix[row, col]

    def normalize(self) -> None:
        """Normalize the matrix so that all values are between 0 and 1."""
        self._matrix = (self._matrix - self.min_value) / (self.max_value - self.min_value)

    def to_distances(self) -> np.ndarray:
        """Convert the matrix to a distance matrix.

        :return: The distance matrix.
        :rtype: np.ndarray
        """
        self.normalize()
        return 1 - self._matrix


class AlignmentMatrix(Matrix):
    """Class for alignment matrix objects."""

    def alignment_score(self) -> float:
        """Return the alignment score of the matrix.

        :return: The alignment score.
        :rtype: float
        """
        return self._matrix[-1, -1]
