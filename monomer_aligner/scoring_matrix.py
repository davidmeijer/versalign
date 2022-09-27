from __future__ import annotations
from enum import Enum, auto
from typing import Dict

import matplotlib.pyplot as plt
import numpy as np

from rdkit import Chem
from rdkit.Chem import rdmolops, AllChem, DataStructs


def mol_to_fingerprint(mol: Chem.Mol, n: int = 2048) -> np.array:
    fingerprint = np.zeros((0,), dtype=int)
    DataStructs.ConvertToNumpyArray(
        AllChem.GetMorganFingerprintAsBitVect(mol, radius=2, nBits=n), 
        fingerprint
    )
    return fingerprint

def tanimoto_similarity(fp1: np.array, fp2: np.array) -> float:
    return (np.logical_and(fp1, fp2).sum() / (float(np.logical_or(fp1, fp2).sum())))


def parse_scoring_matrix(path: str) -> Dict[PKModule, Dict[PKModule, int]]:
    pairwise_scores = {}
    with open(path, 'r') as handle:
        header = handle.readline().strip().split()
        for idx, line in enumerate(handle):
            k_a = PKModule[header[idx]]
            pairwise_scores[k_a] = {}
            k_b = [PKModule[pk] for pk in header]
            # Parse scores line -- first item is PK module name:
            scores = map(int, line.strip().split()[1:])
            for k, score in zip(k_b, scores):
                pairwise_scores[k_a][k] = score
    return pairwise_scores


class PKModule(Enum):
    A1, A2, A3, A4, A5 = auto(), auto(), auto(), auto(), auto()
    B1, B2, B3, B4, B5 = auto(), auto(), auto(), auto(), auto()
    C1 = auto()
    D1 = auto()
    E1, E1a, E1b, E1c = auto(), auto(), auto(), auto()
    E2, E2a, E2b, E2c = auto(), auto(), auto(), auto()
    E3, E3a, E3b, E3c = auto(), auto(), auto(), auto()
    F1, F2 = auto(), auto()
    R, Rr, Rs = auto(), auto(), auto()
    H, Hr, Hs = auto(), auto(), auto()
    K, Kr, Ks = auto(), auto(), auto()
    D, Dr, Ds = auto(), auto(), auto()
    E, Er, Es = auto(), auto(), auto()
    S, Sr, Ss = auto(), auto(), auto()
    T, Tr, Ts = auto(), auto(), auto()
    N, Nr, Ns = auto(), auto(), auto()
    Q, Qr, Qs = auto(), auto(), auto()
    C, Cr, Cs = auto(), auto(), auto()
    G = auto()
    P, Pr, Ps = auto(), auto(), auto()
    A, Ar, As = auto(), auto(), auto()
    V, Vr, Vs = auto(), auto(), auto()
    I, Ir, Is = auto(), auto(), auto()
    L, Lr, Ls = auto(), auto(), auto()
    M, Mr, Ms = auto(), auto(), auto()
    F, Fr, Fs = auto(), auto(), auto()
    Y, Yr, Ys = auto(), auto(), auto()
    W, Wr, Ws = auto(), auto(), auto()
    GAP = auto()  # Denotes gap in PK module sequence

    def display_name(self):
        if self == PKModule.GAP:
            return '-'
        else:
            return self.name

    def display_in_alignment(self):
        max_length = max([len(m.display_name()) for m in PKModule])
        padding = max_length - len(self.display_name())
        return (' ' * padding) + self.display_name()

    def logo_color(self):
        cm = plt.get_cmap('gist_rainbow')
        colors = {
            module.display_name(): cm(1. * module_idx/len(PKModule))
            for module_idx, module in enumerate(PKModule)
        }
        return colors[self.display_name()]


class ScoringMatrix:
    def __init__(self, path: str) -> None:
        self.pairwise_scores = parse_scoring_matrix(path)

    def __call__(self, module1: PKModule, module2: PKModule) -> float:
        try:  # Validate scoring matrix on first call
            if self.__class__.__call__._called:
                return self._score(module1, module2)  # Score modules
        except AttributeError:
            if not self._validate():  # Validate scoring matrix
                raise ValueError('ScoringMatrix is invalid')
            self.__class__.__call__._called = True
            # Call scoring matrix again after validation
            return self(module1, module2)

    def _validate(self) -> bool:
        # Validate if all modules are pairwise represented in scoring matrix
         return all([
            all([
                m2 in self.pairwise_scores[m1]
                for m2 in PKModule
            ])
            for m1 in PKModule
        ])

    def _score(self, module1: PKModule, module2: PKModule) -> float:
        try: 
            score = self.pairwise_scores[module1][module2]
        except KeyError: 
            if not isinstance(module1, PKModule) and not isinstance(module2, PKModule):
                fp1 = mol_to_fingerprint(Chem.MolFromSmiles(module1))
                fp2 = mol_to_fingerprint(Chem.MolFromSmiles(module2))
                return (tanimoto_similarity(fp1, fp2) - 0.5) * 20
            elif isinstance(module1, PKModule) or isinstance(module2, PKModule):
                score = -5.0
        return score
