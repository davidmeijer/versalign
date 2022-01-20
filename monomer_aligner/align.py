"""
Author: David Meijer
Description: functionality for aligning polyketide monomers.
"""
import os
import subprocess

import numpy as np
import scipy.cluster.hierarchy as shc
import matplotlib.pyplot as plt


def run_mafft(
    fasta: str,
    out_dir: str = './out/',
    scoring_matrix: str = './data/monomer_scoring.txt'
):
    path_msa = os.path.join(out_dir, 'mafft_msa.fasta')
    return_code = subprocess.call(
        f'mafft --aamatrix {scoring_matrix} {fasta} > {path_msa}',
        shell=True
    )
    return return_code, path_msa


# Names are copied from GRAPE documentation.
# A = malonyl (Mal)
# R = methylmalonyl (MeM)
# N = methoxylmalonyl (OMeMal)
# D = benzoyl (Bz)
# C = ethylmalonyl (EtM)
# Q = isobutryl (IBu)
# E = methylbuteryl2 (MBu)
# G = Unknown (Unk)
def run_aligner():
    d = {
        'Mal': 'A',
        'MeM': 'R',
        'OMeMal': 'N',
        'Bz': 'D',
        'EtM': 'C',
        'IBu': 'Q',
        'MBu': 'E',
        'unknown': 'G'
    }
    encoded = open('./data/backbones_encoded.fasta', 'w')
    with open('./data/backbones.fasta', 'r') as handle:
        for line in handle:
            if line.startswith('>'):
                encoded.write(line)
            else:
                seq = line.strip().split('-')
                seq = [i.split('/')[0] for i in seq]
                seq = [d[name] for name in seq]
                seq = ''.join(seq)
                encoded.write(f'{seq}\n')
    return_code, path_msa = run_mafft('./data/backbones_encoded.fasta')

    e = {
        '-': 0,
        'A': 1,
        'R': 2,
        'N': 3,
        'D': 4,
        'C': 5,
        'Q': 6,
        'E': 7,
        'G': 8
    }
    d_r = {
        'A': 'Mal',
        'R': 'MeM',
        'N': 'OMeMal',
        'D': 'Bz',
        'C': 'EtM',
        'Q': 'IBu',
        'E': 'MBu',
        'G': 'Unk'
    }
    vecs = []
    labels = []
    with open(path_msa, 'r') as handle:
        rec = []
        for idx, line in enumerate(handle):
            if line.startswith('>') and idx > 0:
                rec = ''.join(rec)
                label = []
                vec = []
                for char in rec:
                    vec.append(e[char])
                    if char != '-':
                        label.append(d_r[char])
                labels.append('-'.join(label))
                vecs.append(np.array(vec))
                rec = []
            elif not line.startswith('>'):
                rec.append(line.strip())
            else:
                pass
    vecs = np.stack(vecs)

    linkage = shc.linkage(vecs, method='ward')
    dend = shc.dendrogram(
        linkage,
        color_threshold=15,
        orientation='right',
        labels=labels
    )
    ax = plt.gca()
    ax.tick_params(axis='y', which='major', labelsize=1)
    plt.savefig('./out/dend.png', dpi=3000)

    e_r = {
        0: '-',
        1: 'A',
        2: 'R',
        3: 'N',
        4: 'D',
        5: 'C',
        6: 'Q',
        7: 'E',
        8: 'G'
    }

    clusters = {}
    prev_c = None
    cs = []
    for idx, c in zip(dend['leaves'], dend['leaves_color_list']):
        if c == prev_c:
            pass
        elif c != prev_c and prev_c == None:
            cs.append(len(cs))
        else:
            cs.append(len(cs))
        prev_c = c

        cl = len(cs)
        if not cl in clusters:
            clusters[cl] = [idx]
        else:
            clusters[cl].append(idx)

    for cl, inds in clusters.items():

        lens = []

        fn = f'./out/seqs_cluster_{cl}.fasta'
        handle = open(fn, 'w')
        for idx, vec in enumerate(vecs[inds]):
            handle.write(f'>seq{idx}\n')
            seq = []
            for i in vec:
                if i != 0:
                    seq.append(e_r[i])
            handle.write(f"{''.join(seq)}\n")
            lens.append(len(seq))
        handle.close()

        avg_len = sum(lens)/len(lens)
        print(avg_len)

        path_msa = f'./out/mafft_msa_cluster_{cl}.fasta'
        command = f'mafft --aamatrix ./data/monomer_scoring.txt {fn} > {path_msa}'
        _ = subprocess.call(command, shell=True)
