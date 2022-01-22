from typing import List, Union


class Matrix:
    def __init__(self, nrows: int, ncols: int) -> None:
        self._matrix = None
        self._nrows = nrows
        self._ncols = ncols

    def __repr__(self) -> str:
        flat_matrix = [str(column) for row in self._matrix for column in row]
        padding = len(max(flat_matrix, key=len))

        display = []
        for row in [list(map(str, row)) for row in self._matrix]:
            row = [(' ' * (padding - len(item)) + item) for item in row]
            display.append(' '.join([str(item) for item in row]))
        return '\n'.join(display)

    def build(self, fill: Union[int, float]) -> None:
        self._matrix = [
            [fill for _ in range(self._ncols)]
            for _ in range(self._nrows)
        ]

    def transpose(self) -> List[List[Union[int, float]]]:
        transposed_matrix = [[] for _ in range(self._ncols)]
        for row in self._matrix:
            for i, column in enumerate(row):
                transposed_matrix[i].append(column)
        self._nrows, self._ncols = self._ncols, self._nrows
        self._matrix = transposed_matrix
        return self._matrix

    def add(self, row: int, col: int, value: Union[int, float]) -> List[List[Union[int, float]]]:
        self._matrix[row][col] = value
        return self._matrix

    def get(self, row: int, col: int) -> Union[int, float]:
        return self._matrix[row][col]


class AlignmentMatrix(Matrix):
    def __init__(self, ncols: int, nrows: int):
        super().__init__(ncols, nrows)


class PairwiseScoreMatrix(Matrix):
    def __init__(self, ncols: int, nrows: int):
        super().__init__(ncols, nrows)
