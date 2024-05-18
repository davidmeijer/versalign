# -*- coding: utf-8 -*-

"""Implementation of the scoring matrix for sequence alignment."""

import numpy as np


class Matrix:
    """Base class for matrix objects."""

    def __init__(self, nrows: int, ncols: int) -> None:
        """Initialize the matrix object.

        :param nrows: The number of rows in the matrix.
        :type nrows: int
        :param ncols: The number of columns in the matrix.
        :type ncols: int
        """
        self._nrows = nrows
        self._ncols = ncols
        self._matrix = self._create_matrix()

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

    def _create_matrix(self) -> np.ndarray:
        """Create the matrix object with NumPy, initialized with zeros.

        :return: The matrix object.
        :rtype: np.ndarray
        """
        return np.zeros((self._nrows, self._ncols), dtype=np.int32)

    def transpose(self) -> np.ndarray:
        """Transpose the matrix.

        :return: The transposed matrix.
        :rtype: np.ndarray
        """
        return self._matrix.T

    def set_value(self, row: int, col: int, value: int) -> None:
        """Set the value of a cell in the matrix.

        :param row: The row index.
        :type row: int
        :param col: The column index.
        :type col: int
        :param value: The value to set.
        :type value: int
        :raises ValueError: If the row or column index is out of bounds.
        :raises ValueError: If the value is not an integer.
        """
        if not (0 <= row < self._nrows):
            raise ValueError("Row index out of bounds.")

        if not (0 <= col < self._ncols):
            raise ValueError("Column index out of bounds.")

        if not isinstance(value, int):
            raise ValueError("Value must be an integer.")

        self._matrix[row, col] = value

    def get_value(self, row: int, col: int) -> int:
        """Get the value of a cell in the matrix.

        :param row: The row index.
        :type row: int
        :param col: The column index.
        :type col: int
        :return: The value of the cell.
        :rtype: int
        :raises ValueError: If the row or column index is out of bounds.
        """
        if not (0 <= row < self._nrows):
            raise ValueError("Row index out of bounds.")

        if not (0 <= col < self._ncols):
            raise ValueError("Column index out of bounds.")

        return self._matrix[row, col]
