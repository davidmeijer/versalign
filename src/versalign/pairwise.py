# -*- coding: utf-8 -*-

"""Implementation of pairwise alignment between two sequences."""

import logging
import typing as ty
from enum import Enum, IntEnum, auto

import numpy as np

from .matrix import AlignmentMatrix
from .motif import Gap, Motif
from .sequence import Sequence


class PairwiseAlignment(Enum):
    """Enum for the type of pairwise alignment."""

    NEEDLEMAN_WUNSCH = auto()  # Needleman-Wunsch for global alignment.
    SMITH_WATERMAN = auto()  # Smith-Waterman for local alignment.


def create_alignment_matrix_needelman_wunsch(
    seq_a: Sequence,
    seq_b: Sequence,
    gap_penalty: int,
    end_gap_penalty: int,
    score_func: ty.Callable[[Motif, Motif], int],
) -> AlignmentMatrix:
    """Create the alignment matrix for the two sequences.

    :param seq_a: The first sequence to align.
    :type seq_a: Sequence
    :param seq_b: The second sequence to align.
    :type seq_b: Sequence
    :param gap_penalty: The penalty for inserting a gap.
    :type gap_penalty: int
    :param end_gap_penalty: The penalty for inserting a gap at the end of a sequence.
    :type end_gap_penalty: int
    :param score_func: The function to score a pair of motifs.
    :type score_func: ty.Callable[[Motif, Motif], int]
    :return: The alignment matrix for the two sequences.
    :rtype: AlignmentMatrix
    :raises ValueError: If the sequences are not Sequence objects.
    :raises ValueError: If the gap penalties are not integers.
    :raises ValueError: If the score function is not callable.
    """
    logger = logging.getLogger(__name__)
    logger.debug("Creating alignment matrix...")

    # Check if seq_a is Sequence and seq_b is Sequence.
    if not isinstance(seq_a, Sequence) or not isinstance(seq_b, Sequence):
        msg = "seq_a and seq_b must be Sequence objects."
        logger.error(msg)
        raise ValueError(msg)

    # Check if gap_penalty is int and end_gap_penalty is int.
    if not isinstance(gap_penalty, int) or not isinstance(end_gap_penalty, int):
        msg = "gap_penalty and end_gap_penalty must be integers."
        logger.error(msg)
        raise ValueError(msg)

    # Check if score_func is callable.
    if not callable(score_func):
        msg = "score_func must be callable."
        logger.error(msg)
        raise ValueError(msg)

    # Instatiate the alignment matrix.
    nrows = len(seq_a) + 1
    ncols = len(seq_b) + 1
    matrix = AlignmentMatrix(nrows, ncols, fill=0.0)

    # Fill in initial alignment for when there is zero alignment.
    for ri in range(nrows):
        matrix.set_value(ri, 0, ri * end_gap_penalty)

    for ci in range(ncols):
        matrix.set_value(0, ci, ci * end_gap_penalty)

    # Fill in the rest of the alignment matrix.
    for ci in range(1, ncols):
        for ri in range(1, nrows):

            # Calculate the pairwise score.
            score = score_func(seq_a[ri - 1], seq_b[ci - 1])

            # Calculate the penalty for inserting a gap.
            if ci == len(seq_b) or ri == len(seq_a):
                penalty = end_gap_penalty
            else:
                penalty = gap_penalty

            # Calculate the alignment score.
            matrix.set_value(
                row=ri,
                col=ci,
                value=max(
                    matrix.get_value(ri - 1, ci) - penalty,  # Insert gap in seq_a.
                    matrix.get_value(ri, ci - 1) - penalty,  # Insert gap in seq_b.
                    matrix.get_value(ri - 1, ci - 1) + score,  # Match or mismatch.
                ),
            )

    logger.debug("Alignment matrix created.")
    return matrix


def traceback_alignment_needelman_wunsch(
    matrix: AlignmentMatrix,
    seq_a: Sequence,
    seq_b: Sequence,
    gap_penalty: int,
    end_gap_penalty: int,
) -> ty.List[ty.Tuple[Motif, Motif]]:
    """Traceback the alignment matrix to find the optimal alignment.

    :param matrix: The alignment matrix.
    :type matrix: AlignmentMatrix
    :param seq_a: The first sequence to align.
    :type seq_a: Sequence
    :param seq_b: The second sequence to align.
    :type seq_b: Sequence
    :param gap_penalty: The penalty for inserting a gap.
    :type gap_penalty: int
    :param end_gap_penalty: The penalty for inserting a gap at the end of a sequence.
    :type end_gap_penalty: int
    :return: The optimal alignment of the two sequences.
    :rtype: ty.List[ty.Tuple[Motif, Motif]]
    :raises ValueError: If the sequences and matrix are not compatible.
    :raises ValueError: If the sequences are not Sequence objects.
    :raises ValueError: If the gap penalties are not integers.
    :raises ValueError: If the matrix is not a Matrix object.
    :raises ValueError: If the traceback is invalid.
    """

    def traceback(ri: int, ci: int) -> ty.List[ty.Tuple[Motif, Motif]]:
        """Traceback the alignment matrix recursively.

        :param ri: The row index.
        :type ri: int
        :param ci: The column index.
        :type ci: int
        :return: The optimal alignment of the two sequences.
        :rtype: ty.List[ty.Tuple[Motif, Motif]]
        :raises ValueError: If the traceback is invalid.
        """
        # End of the traceback.
        if ri == 0 and ci == 0:
            logger.debug("End of traceback.")
            return []

        if ri == 0:
            return traceback(ri, ci - 1) + [(Gap(), seq_b[ci - 1])]

        if ci == 0:
            return traceback(ri - 1, ci) + [(seq_a[ri - 1], Gap())]

        # Gap penalty.
        if ri == len(seq_a) or ci == len(seq_b):
            penalty = end_gap_penalty
        else:
            penalty = gap_penalty

        # Calculate score of current cell and possible moves.
        current = matrix.get_value(ri, ci)
        horizontal = matrix.get_value(ri, ci - 1)
        diagonal = matrix.get_value(ri - 1, ci - 1)
        vertical = matrix.get_value(ri - 1, ci)

        # Traceback for match.
        if (
            seq_a[ri - 1] == seq_b[ci - 1]
            or diagonal > max(horizontal, vertical)
            or ((vertical - current) != penalty and (horizontal - current) != penalty)
        ):
            return traceback(ri - 1, ci - 1) + [(seq_a[ri - 1], seq_b[ci - 1])]

        # Traceback for gap.
        if horizontal > vertical:
            return traceback(ri, ci - 1) + [(Gap(), seq_b[ci - 1])]

        if vertical >= horizontal:
            return traceback(ri - 1, ci) + [(seq_a[ri - 1], Gap())]

        raise ValueError("Invalid traceback.")

    logger = logging.getLogger(__name__)
    logger.debug("Tracing back alignment...")

    # Check if matrix is a Matrix object.
    if not isinstance(matrix, AlignmentMatrix):
        msg = "matrix must be a Matrix object."
        logger.error(msg)
        raise ValueError(msg)

    # Check if seq_a is Sequence and seq_b is Sequence.
    if not isinstance(seq_a, Sequence) or not isinstance(seq_b, Sequence):
        msg = "seq_a and seq_b must be Sequence objects."
        logger.error(msg)
        raise ValueError(msg)

    # Check if gap_penalty is int and end_gap_penalty is int.
    if not isinstance(gap_penalty, int) or not isinstance(end_gap_penalty, int):
        msg = "gap_penalty and end_gap_penalty must be integers."
        logger.error(msg)
        raise ValueError(msg)

    # Check if seq_a and seq_b could have created the matrix.
    if len(seq_a) + 1 != matrix.nrows or len(seq_b) + 1 != matrix.ncols:
        msg = "The sequences and matrix are not compatible."
        logger.error(msg)
        raise ValueError(msg)

    return traceback(ri=matrix.nrows - 1, ci=matrix.ncols - 1)


def needleman_wunsch(
    seq_a: Sequence,
    seq_b: Sequence,
    score_func: ty.Callable[[Motif, Motif], int],
    gap_penalty: int = 2,
    end_gap_penalty: int = 1,
) -> ty.Tuple[Sequence, Sequence, float]:
    """Align two sequences using Needleman-Wunsch algorithm.

    :param seq_a: The first sequence to align.
    :type seq_a: Sequence
    :param seq_b: The second sequence to align.
    :type seq_b: Sequence
    :param score_func: The function to score a pair of motifs.
    :type score_func: ty.Callable[[Motif, Motif], int]
    :param gap_penalty: The penalty for inserting a gap. Default is 2.
    :type gap_penalty: int
    :param end_gap_penalty: The penalty for inserting a gap at the end of a sequence.
        Default is 1.
    :type end_gap_penalty: int
    :return: The optimal alignment of the two sequences and the alignment score.
    :rtype: ty.Tuple[Sequence, Sequence, float]
    :raises ValueError: If the sequences are not Sequence objects.
    :raises ValueError: If the gap penalties are not integers.
    :raises ValueError: If the score function is not callable.
    """
    logger = logging.getLogger(__name__)
    logger.debug("Aligning sequences using Needleman-Wunsch algorithm...")

    # Check if seq_a is Sequence and seq_b is Sequence.
    if not isinstance(seq_a, Sequence) or not isinstance(seq_b, Sequence):
        msg = "seq_a and seq_b must be Sequence objects."
        logger.error(msg)
        raise ValueError(msg)

    # Check if gap_penalty is int and end_gap_penalty is int.
    if not isinstance(gap_penalty, int) or not isinstance(end_gap_penalty, int):
        msg = "gap_penalty and end_gap_penalty must be integers."
        logger.error(msg)
        raise ValueError(msg)

    # Check if score_func is callable.
    if not callable(score_func):
        msg = "score_func must be callable."
        logger.error(msg)
        raise ValueError(msg)

    # Create the alignment matrix.
    matrix = create_alignment_matrix_needelman_wunsch(
        seq_a=seq_a,
        seq_b=seq_b,
        gap_penalty=gap_penalty,
        end_gap_penalty=end_gap_penalty,
        score_func=score_func,
    )

    alignment_score = matrix.alignment_score()

    # Traceback the alignment matrix.
    aligned = traceback_alignment_needelman_wunsch(
        matrix=matrix,
        seq_a=seq_a,
        seq_b=seq_b,
        gap_penalty=gap_penalty,
        end_gap_penalty=end_gap_penalty,
    )

    seq_a_aligned = Sequence(seq_a.identifier, [motif_a for motif_a, _ in aligned])
    seq_b_aligned = Sequence(seq_b.identifier, [motif_b for _, motif_b in aligned])

    logger.debug("Sequences aligned.")
    return seq_a_aligned, seq_b_aligned, alignment_score


class SmithWatermanDirection(IntEnum):
    """Enum for the direction of the Smith-Waterman algorithm."""

    STOP = 0
    DIAGONAL = 1
    VERTICAL = 2
    HORIZONTAL = 3


def smith_waterman(
    seq_a: Sequence,
    seq_b: Sequence,
    score_func: ty.Callable[[Motif, Motif], int],
    gap_penalty: int = 2,
) -> ty.Tuple[Sequence, Sequence, float]:
    """Create the alignment matrix for the two sequences.

    :param seq_a: The first sequence to align.
    :type seq_a: Sequence
    :param seq_b: The second sequence to align.
    :type seq_b: Sequence
    :param score_func: The function to score a pair of motifs.
    :type score_func: ty.Callable[[Motif, Motif], int]
    :param gap_penalty: The penalty for inserting a gap. Default is 2.
    :type gap_penalty: int
    :return: The aligned sequences and the alignment score.
    :rtype: ty.Tuple[Sequence, Sequence, float]
    :raises ValueError: If the sequences are not Sequence objects.
    :raises ValueError: If the gap penalty is not an integer.
    :raises ValueError: If the score function is not callable.
    """
    logger = logging.getLogger(__name__)
    logger.debug("Aligning sequences using Smith-Waterman algorithm...")

    # Check if seq_a is Sequence and seq_b is Sequence.
    if not isinstance(seq_a, Sequence) or not isinstance(seq_b, Sequence):
        msg = "seq_a and seq_b must be Sequence objects."
        logger.error(msg)
        raise ValueError(msg)

    # Check if gap_penalty is int.
    if not isinstance(gap_penalty, int):
        msg = "gap_penalty must be an integer."
        logger.error(msg)
        raise ValueError(msg)

    # Check if score_func is callable.
    if not callable(score_func):
        msg = "score_func must be callable."
        logger.error(msg)
        raise ValueError(msg)

    # Instatiate the alignment matrix.
    nrows = len(seq_a) + 1
    ncols = len(seq_b) + 1
    matrix = AlignmentMatrix(nrows, ncols, fill=0.0)

    # Initialize tracing matrix.
    tracing_matrix = np.zeros(shape=(nrows, ncols), dtype=int)

    # Initialize maximum score and its position.
    max_score: ty.Optional[float] = None
    max_index: ty.Optional[ty.Tuple[int, int]] = None

    logger.debug("Creating alignment matrix...")

    # Fill in the alignment matrix.
    for ri in range(1, nrows):
        for ci in range(1, ncols):

            # Calculate score of current cell and possible moves.
            diagonal = matrix.get_value(ri - 1, ci - 1) + score_func(seq_a[ri - 1], seq_b[ci - 1])
            vertical = matrix.get_value(ri - 1, ci) - gap_penalty
            horizontal = matrix.get_value(ri, ci - 1) - gap_penalty

            directions = [
                (SmithWatermanDirection.STOP, 0),
                (SmithWatermanDirection.DIAGONAL, diagonal),
                (SmithWatermanDirection.VERTICAL, vertical),
                (SmithWatermanDirection.HORIZONTAL, horizontal),
            ]

            # Calculate the alignment score.
            highest_score_index = np.argmax([score for _, score in directions])
            highest_score_direction = directions[highest_score_index][0]
            highest_score_value = directions[highest_score_index][1]
            matrix.set_value(row=ri, col=ci, value=highest_score_value)
            tracing_matrix[ri, ci] = highest_score_direction

            # Update maximum score and its position.
            if max_score is None or highest_score_value > max_score:
                max_score = highest_score_value
                max_index = (ri, ci)

    if max_score is None or max_index is None:
        msg = "Maximum score and its position not found."
        logger.error(msg)
        raise ValueError(msg)

    logger.debug("Alignment matrix created.")

    # Traceback the alignment matrix.
    aligned_a: ty.List[Motif] = []
    aligned_b: ty.List[Motif] = []
    current_aligned_a: ty.Optional[Motif] = None
    current_aligned_b: ty.Optional[Motif] = None

    logger.debug("Tracing back alignment...")
    ri, ci = max_index

    while tracing_matrix[ri, ci] != SmithWatermanDirection.STOP:
        if tracing_matrix[ri, ci] == SmithWatermanDirection.DIAGONAL:
            current_aligned_a = seq_a[ri - 1]
            current_aligned_b = seq_b[ci - 1]
            ri -= 1
            ci -= 1

        elif tracing_matrix[ri, ci] == SmithWatermanDirection.VERTICAL:
            current_aligned_a = seq_a[ri - 1]
            current_aligned_b = Gap()
            ri -= 1

        elif tracing_matrix[ri, ci] == SmithWatermanDirection.HORIZONTAL:
            current_aligned_a = Gap()
            current_aligned_b = seq_b[ci - 1]
            ci -= 1

        else:
            msg = "Invalid traceback. Should be diagonal, vertical, or horizontal."
            logger.error(msg)
            raise ValueError(msg)

        if current_aligned_a is not None and current_aligned_b is not None:
            aligned_a = aligned_a + [current_aligned_a]
            aligned_b = aligned_b + [current_aligned_b]
        else:
            msg = "Invalid traceback. Found motif as None."
            logger.error(msg)
            raise ValueError(msg)

    logger.debug("Alignment traced back.")

    # Reverse the aligned sequences.
    aligned_a = aligned_a[::-1]
    aligned_b = aligned_b[::-1]

    logger.debug("Sequences aligned.")

    return (
        Sequence(seq_a.identifier, aligned_a),
        Sequence(seq_b.identifier, aligned_b),
        max_score,
    )


def align_pairwise(
    seq_a: Sequence,
    seq_b: Sequence,
    score_func: ty.Callable[[Motif, Motif], int],
    algorithm: PairwiseAlignment = PairwiseAlignment.NEEDLEMAN_WUNSCH,
    options: ty.Optional[ty.Dict[str, ty.Any]] = None,
) -> ty.Tuple[Sequence, Sequence, float]:
    """Align two sequences using pairwise alignment.

    :param seq_a: The first sequence to align.
    :type seq_a: Sequence
    :param seq_b: The second sequence to align.
    :type seq_b: Sequence
    :param score_func: The function to score a pair of motifs.
    :type score_func: ty.Callable[[Motif, Motif], int]
    :param algorithm: The algorithm to use for pairwise alignment.
    :type algorithm: PairwiseAlignment
    :param options: Additional keyword arguments for the algorithm.
    :type options: ty.Optional[ty.Dict[str, ty.Any]]
    :return: The optimal alignment of the two sequences and the alignment score.
    :rtype: ty.Tuple[Sequence, Sequence, float]
    :raises ValueError: If the sequences are not Sequence objects.
    :raises ValueError: If the gap penalties are not integers.
    :raises ValueError: If the score function is not callable.
    :raises ValueError: If the algorithm is not PairwiseAlignment.
    :raises NotImplementedError: If the algorithm is not implemented.
    """
    logger = logging.getLogger(__name__)

    # Check if seq_a is Sequence and seq_b is Sequence.
    if not isinstance(seq_a, Sequence) or not isinstance(seq_b, Sequence):
        msg = "seq_a and seq_b must be Sequence objects."
        logger.error(msg)
        raise ValueError(msg)

    # Check if score_func is callable.
    if not callable(score_func):
        msg = "score_func must be callable."
        logger.error(msg)
        raise ValueError(msg)

    # Check if algorithm is PairwiseAlignment.
    if not isinstance(algorithm, PairwiseAlignment):
        msg = "algorithm must be PairwiseAlignment."
        logger.error(msg)
        raise ValueError(msg)

    if algorithm == PairwiseAlignment.NEEDLEMAN_WUNSCH:
        seq_a_aligned, seq_b_aligned, alignment_score = needleman_wunsch(
            seq_a=seq_a,
            seq_b=seq_b,
            score_func=score_func,
            **(options or {}),
        )
        return seq_a_aligned, seq_b_aligned, alignment_score

    elif algorithm == PairwiseAlignment.SMITH_WATERMAN:
        seq_a_aligned, seq_b_aligned, alignment_score = smith_waterman(
            seq_a=seq_a,
            seq_b=seq_b,
            score_func=score_func,
            **(options or {}),
        )
        return seq_a_aligned, seq_b_aligned, alignment_score

    else:
        msg = f"Pairwise alignment algorithm {algorithm} is not implemented."
        logger.error(msg)
        raise NotImplementedError(msg)
