import os
import typing as ty

import numpy as np

from scipy.spatial.distance import squareform
from scipy.cluster.hierarchy import linkage, to_tree


from .parser import parse_fasta
from .aligner import ModuleSequence, PairwiseAlignment


def to_newick(linkage_matrix: np.ndarray, labels: ty.List[str]) -> str:
    """
    Outputs linkage matrix as Newick file

    Parameters
    ----------
    linkage_matrix : np.ndarray
        condensed distance matrix
    labels : list of str, optional
        leaf labels 

    Returns
    -------
    newick : str
        linkage matrix in newick format tree

    Source:
    https://stackoverflow.com/questions/28222179/save-dendrogram-to-newick-format
    """
    tree = to_tree(linkage_matrix, rd=False)

    def get_newick(node, newick, parentdist, leaf_names):
        if node.is_leaf():
            return "%s:%.2f%s" % (
                leaf_names[node.id],
                parentdist - node.dist,
                newick
            )
        else:
            if len(newick) > 0:
                newick = "):%.2f%s" % (parentdist - node.dist, newick)
            else:
                newick = ");"
            newick = get_newick(
                node.get_left(),
                newick,
                node.dist,
                leaf_names
            )
            newick = get_newick(
                node.get_right(),
                ",%s" % (newick),
                node.dist,
                leaf_names
            )
            newick = "(%s" % (newick)
            return newick

    newick = get_newick(tree, "", tree.dist, labels)

    return newick



def unit_to_color(unit: str) -> str:
    colors = {
        "A1": "#E50000",
        "A2": "#DC143C",
        "A3": "#FF0000",
        "B1": "#0343DF",
        "B2": "#0000FF",
        "C": "#FF00FF",
        "D": "#FFD700",
        "E1": "#ED7014",
        "E2": "#8D4004",
        "F1": "#FDA172",
        "F2": "#ED7117",
        "G1": "#FC6A03",
        "G2": "#ED820E",
        "H1": "#008080",
        "H2": "#029386",
        "O": "#D3D3D3",
    }
    return colors[unit]



def cluster_polyketide_backbones(msa: str, out_dir: str) -> None:
    records = parse_fasta(msa)  # ty.List[Record]
    records = [ModuleSequence(r.name, r.seq) for r in records]

    dist_matrix = np.zeros((len(records), len(records)))
    labels = [r.name for r in records]

    for i, seq_a in enumerate(records):
        for j, seq_b in enumerate(records):
            if j > i:
                alignment = PairwiseAlignment(
                    name_seq1=seq_a.name,
                    aligned_seq1=seq_a._seq,
                    name_seq2=seq_b.name,
                    aligned_seq2=seq_b._seq,
                    alignment_score=float("inf"),
                    gap_cost=float("inf"),
                    end_gap_cost=float("inf")
                )
                dist = 1 - (alignment.percentage_identity() / 100)
                dist_matrix.itemset((i, j), dist)
                dist_matrix.itemset((j, i), dist)
            else:
                continue 
            
    condensed = squareform(dist_matrix)
    linkage_matrix = linkage(condensed, method="ward")
    newick = to_newick(linkage_matrix, labels)

    fo = open(os.path.join(out_dir, "_tree.nwk"), "w")
    fo.write(newick)
    fo.close()

    msa_length = len(records[0]._seq)

    fo = open(os.path.join(out_dir, "_tree_subunit_annotation.txt"), "w")
    fo.write("DATASET_DOMAINS\n")
    fo.write("SEPARATOR COMMA\n")
    fo.write("DATASET_LABEL,polyketide_backbones\n")
    fo.write("COLOR,#ff0000\n")
    # fo.write("WIDTH,400\n")
    fo.write("DATA\n")
    for r in records:
        label = r.name 
        seq = [item.name for item, _ in r._seq]
        subunits = ",".join([
            f"RE|{idx}|{idx + 1}|{unit_to_color(unit)}|{unit}" 
            for idx, unit in enumerate(seq)
            if unit != "GAP"
        ])
        fo.write(f"{label},{msa_length},{subunits}\n")
    fo.close()
