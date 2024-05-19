# -*- coding: utf-8 -*-

"""Implementation of multiple sequence alignment."""

import logging
import typing as ty

from scipy.cluster.hierarchy import linkage
from scipy.spatial.distance import pdist

from .matrix import Matrix
from .motif import Gap, Motif
from .pairwise import PairwiseAlignment, align_pairwise, create_alignment_matrix_needelman_wunsch
from .sequence import Sequence


def pairwise_scoring_matrix(
    seqs: ty.List[Sequence],
    gap_penalty: int,
    end_gap_penalty: int,
    score_func: ty.Callable[[Motif, Motif], int],
) -> Matrix:
    """Compute the pairwise similarity matrix of the given sequences.

    :param seqs: List of sequences to compute the similarity matrix.
    :type seqs: List[Sequence]
    :param gap_penalty: Penalty for opening a gap.
    :type gap_penalty: int
    :param end_gap_penalty: Penalty for extending a gap.
    :type end_gap_penalty: int
    :param score_func: Function to score two motifs.
    :type score_func: Callable[[Motif, Motif], int]
    :return: Pairwise similarity matrix.
    :rtype: Matrix

    The similarity matrix is computed using the global Needleman-Wunsch
    algorithm. The similarity matrix is symmetric, so only the upper triangle is
    computed.
    """
    num_seqs = len(seqs)
    matrix = Matrix(num_seqs, num_seqs, fill=0.0)

    for i, seq_a in enumerate(seqs):
        for j, seq_b in enumerate(seqs):

            # Skip if the similarity has already been computed.
            if j < i:
                continue

            score = create_alignment_matrix_needelman_wunsch(
                seq_a=seq_a,
                seq_b=seq_b,
                gap_penalty=gap_penalty,
                end_gap_penalty=end_gap_penalty,
                score_func=score_func,
            ).alignment_score()

            matrix.set_value(i, j, score)
            matrix.set_value(j, i, score)

    return matrix


def merge_singles(
    single1: Sequence,
    single2: Sequence,
    gap_penalty: int,
    end_gap_penalty: int,
    score_func: ty.Callable[[Motif, Motif], int],
) -> ty.List[Sequence]:
    """Merge two single sequences into a single alignment.

    :param single1: First sequence to merge.
    :type single1: Sequence
    :param single2: Second sequence to merge.
    :type single2: Sequence
    :param gap_penalty: Penalty for opening a gap.
    :type gap_penalty: int
    :param end_gap_penalty: Penalty for extending a gap.
    :type end_gap_penalty: int
    :param score_func: Function to score two motifs.
    :type score_func: Callable[[Motif, Motif], int]
    :return: List with two aligned sequences.
    :rtype: ty.List[Sequence]
    """
    single1_aligned, single2_aligned, _ = align_pairwise(
        seq_a=single1,
        seq_b=single2,
        score_func=score_func,
        algorithm=PairwiseAlignment.NEEDLEMAN_WUNSCH,
        options={"gap_penalty": gap_penalty, "end_gap_penalty": end_gap_penalty},
    )

    # Clear all tags.
    single1_aligned.clear_tags()
    single2_aligned.clear_tags()

    return [single1_aligned, single2_aligned]


def merge_single_with_cluster(
    single: Sequence,
    cluster: ty.List[Sequence],
    gap_penalty: int,
    end_gap_penalty: int,
    score_func: ty.Callable[[Motif, Motif], int],
) -> ty.List[Sequence]:
    """Merge a single sequence with a cluster of sequences.

    :param single: Single sequence to merge.
    :type single: Sequence
    :param cluster: Cluster of sequences to merge with.
    :type cluster: List[Sequence]
    :param gap_penalty: Penalty for opening a gap.
    :type gap_penalty: int
    :param end_gap_penalty: Penalty for extending a gap.
    :type end_gap_penalty: int
    :param score_func: Function to score two motifs.
    :type score_func: Callable[[Motif, Motif], int]
    :return: Merged sequence.
    :rtype: ty.List[Sequence]
    """
    # Tag already aligned sequences with original location for insertion
    # of possible gaps after annealing the new sequence.
    first = cluster[0]
    first.tag()
    last = cluster[-1]
    last.tag()

    # Align the leaf with the first and the last sequence in the cluster.
    single_aligned_with_first, first_aligned, score_with_first = align_pairwise(
        seq_a=single,
        seq_b=first,
        score_func=score_func,
        algorithm=PairwiseAlignment.NEEDLEMAN_WUNSCH,
        options={"gap_penalty": gap_penalty, "end_gap_penalty": end_gap_penalty},
    )
    single_aligned_with_last, last_aligned, score_with_last = align_pairwise(
        seq_a=single,
        seq_b=last,
        score_func=score_func,
        algorithm=PairwiseAlignment.NEEDLEMAN_WUNSCH,
        options={"gap_penalty": gap_penalty, "end_gap_penalty": end_gap_penalty},
    )

    if score_with_first >= score_with_last:
        new, anchor = single_aligned_with_first, first_aligned
        others = cluster[1:]
    else:
        new, anchor = single_aligned_with_last, last_aligned
        others = cluster[:-1]

    for motif_idx, motif in enumerate(anchor):
        if motif.get_tag() is None:  # New insertion if tag is None.
            for seq in others:
                seq.insert(motif_idx, Gap())

    # Clear all tags.
    for seq in others:
        seq.clear_tags()
    anchor.clear_tags()
    new.clear_tags()

    # Insert the new sequence into the cluster.
    if score_with_first >= score_with_last:
        new_cluster = [new, anchor] + others
    else:
        new_cluster = others + [anchor, new]

    return new_cluster


def merge_clusters(
    cluster1: ty.List[Sequence],
    cluster2: ty.List[Sequence],
    gap_penalty: int,
    end_gap_penalty: int,
    score_func: ty.Callable[[Motif, Motif], int],
) -> ty.List[Sequence]:
    """Merge two clusters of sequences.

    :param cluster1: First cluster of sequences to merge.
    :type cluster1: List[Sequence]
    :param cluster2: Second cluster of sequences to merge.
    :type cluster2: List[Sequence]
    :param gap_penalty: Penalty for opening a gap.
    :type gap_penalty: int
    :param end_gap_penalty: Penalty for extending a gap.
    :type end_gap_penalty: int
    :param score_func: Function to score two motifs.
    :type score_func: Callable[[Motif, Motif], int]
    :return: Merged sequence.
    :rtype: ty.List[Sequence]
    """
    # Tag already aligned sequences with original location for insertion
    # of possible gaps after annealing the new sequence.
    msa1_top, msa1_bottom = cluster1[0], cluster1[-1]
    msa2_top, msa2_bottom = cluster2[0], cluster2[-1]
    msa1_top.tag()
    msa1_bottom.tag()
    msa2_top.tag()
    msa2_bottom.tag()

    # Align the top of msa1 with the bottom of msa2 and vice versa.
    msa1_bottom_aligned_with_msa2_top, msa2_top_aligned_with_msa1_bottom, msa1_top_score = (
        align_pairwise(
            seq_a=msa1_bottom,
            seq_b=msa2_top,
            score_func=score_func,
            algorithm=PairwiseAlignment.NEEDLEMAN_WUNSCH,
            options={"gap_penalty": gap_penalty, "end_gap_penalty": end_gap_penalty},
        )
    )
    msa1_top_aligned_to_msa2_bottom, msa2_bottom_aligned_to_msa1_top, msa2_top_score = (
        align_pairwise(
            seq_a=msa1_top,
            seq_b=msa2_bottom,
            score_func=score_func,
            algorithm=PairwiseAlignment.NEEDLEMAN_WUNSCH,
            options={"gap_penalty": gap_penalty, "end_gap_penalty": end_gap_penalty},
        )
    )

    if msa1_top_score >= msa2_top_score:
        msa1_align_to, msa2_align_to = (
            msa1_bottom_aligned_with_msa2_top,
            msa2_top_aligned_with_msa1_bottom,
        )
        msa1_others = cluster1[:-1]
        msa2_others = cluster2[1:]
    else:
        msa2_align_to, msa1_align_to = (
            msa1_top_aligned_to_msa2_bottom,
            msa2_bottom_aligned_to_msa1_top,
        )
        msa1_others = cluster1[1:]
        msa2_others = cluster2[:-1]

    for motif_idx, motif in enumerate(msa1_align_to):
        if motif.get_tag() is None:  # New insertion if tag is None.
            for seq in msa1_others:
                seq.insert(motif_idx, Gap())

    for motif_idx, motif in enumerate(msa2_align_to):
        if motif.get_tag() is None:  # New insertion if tag is None.
            for seq in msa2_others:
                seq.insert(motif_idx, Gap())

    # Clear all tags.
    for seq in msa1_others:
        seq.clear_tags()
    for seq in msa2_others:
        seq.clear_tags()
    msa1_align_to.clear_tags()
    msa2_align_to.clear_tags()

    # Insert the new sequence into the cluster.
    if msa1_top_score >= msa2_top_score:
        new_cluster = msa1_others + [msa1_align_to, msa2_align_to] + msa2_others
    else:
        new_cluster = msa2_others + [msa2_align_to, msa1_align_to] + msa1_others

    return new_cluster


def multiple_sequence_alignment(
    seqs: ty.List[Sequence],
    gap_penalty: int,
    end_gap_penalty: int,
    score_func: ty.Callable[[Motif, Motif], int],
) -> ty.List[Sequence]:
    """Perform multiple sequence alignment on the given sequences.

    :param seqs: List of sequences to align.
    :type seqs: List[Sequence]
    :param gap_penalty: Penalty for opening a gap.
    :type gap_penalty: int
    :param end_gap_penalty: Penalty for extending a gap.
    :type end_gap_penalty: int
    :param score_func: Function to score two motifs.
    :type score_func: Callable[[Motif, Motif], int]
    :return: List of aligned sequences.
    :rtype: List[Sequence]
    :raises ValueError: If the multiple sequence alignment is incomplete.

    Based on: 'Progressive Sequence Alignment as a Prerequisite to Correct
    Phylogenetic Trees' by Feng and Doolittle, 1987
    """
    logger = logging.getLogger(__name__)

    if not seqs:
        logger.error("No sequences to align.")
        return []
    elif len(seqs) == 1:
        logger.error("Only one sequence to align.")
        return seqs
    else:
        logger.debug("Computing pairwise similarity matrix...")
        scores = pairwise_scoring_matrix(seqs, gap_penalty, end_gap_penalty, score_func)
        scores.to_distances()
        logger.debug("Computed pairwise similarity matrix.")

    # Identification of most closely related pairs.
    guide_tree = linkage(pdist(scores._matrix), method="ward")
    msa = {seq_idx: [seq] for seq_idx, seq in enumerate(seqs)}

    # Progressive insertion of neutral elements (can create new gaps, but cannot
    # remove existing gaps).
    for pair_idx, pair in enumerate(guide_tree):

        # Every pair in the guide tree can be a new pair that needs seed alignment
        # or it is a leaf connecting to to an existing alignment.
        j1, j2 = int(pair[0]), int(pair[1])
        new_idx = pair_idx + len(seqs)

        if len(msa[j1]) == 1 and len(msa[j2]) == 1:
            seed1, seed2 = msa[j1][0], msa[j2][0]
            seq1, seq2 = merge_singles(seed1, seed2, gap_penalty, end_gap_penalty, score_func)
            msa[new_idx] = [seq1, seq2]

        elif len(msa[j1]) == 1 or len(msa[j2]) == 1:

            # One of the sequences is a leaf, the other is an aligned cluster.
            if len(msa[j1]) == 1:
                leaf, cluster = msa[j1][0], msa[j2]
            else:
                leaf, cluster = msa[j2][0], msa[j1]

            new_cluster = merge_single_with_cluster(
                single=leaf,
                cluster=cluster,
                gap_penalty=gap_penalty,
                end_gap_penalty=end_gap_penalty,
                score_func=score_func,
            )

            # Update the MSA.
            msa[new_idx] = new_cluster

        # Both items are already aligned cluste: len(msa[j1]) > 1 and len(msa[j2]) > 1
        else:
            # First we need to decide which MSA comes on top. To determine this
            # we need to score the top of msa1 with the bottom of msa2 and vice
            # versa.
            cluster1, cluster2 = msa[j1], msa[j2]

            new_cluster = merge_clusters(
                cluster1=cluster1,
                cluster2=cluster2,
                gap_penalty=gap_penalty,
                end_gap_penalty=end_gap_penalty,
                score_func=score_func,
            )

            # Update the MSA.
            msa[new_idx] = new_cluster

        del msa[j1]
        del msa[j2]

    # Retrieve final MSA.
    clusters = list(msa.keys())

    if len(clusters) == 1:
        return msa[clusters[0]]
    else:
        msg = f"MSA incomplete. Found {len(clusters)} unaligned clusters."
        logger.error(msg)
        raise ValueError(msg)
