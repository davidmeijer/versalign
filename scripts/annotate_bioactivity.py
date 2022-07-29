#!/usr/bin/env python3
from sys import argv 
import numpy as np
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem
from tqdm import tqdm

def mol_to_fp(mol, num_bits, radius) -> np.array:
    bit_fingerprint = np.zeros((0,), dtype=int)
    morgan_bit_vector = AllChem.GetMorganFingerprintAsBitVect(mol, radius, num_bits)
    DataStructs.ConvertToNumpyArray(morgan_bit_vector, bit_fingerprint)
    return bit_fingerprint

def sim(first_fingerprint, second_fingerprint):
    return (
        np.logical_and(first_fingerprint, second_fingerprint).sum() / (
            float(np.logical_or(first_fingerprint, second_fingerprint).sum())
        )
    )


fn_donphan = argv[1]

mols = {
    # "Epothilone_sequence_1": r"CC1CCCC2C(O2)CC(OC(=O)CC(C(C(=O)C(C1O)C)(C)C)O)C(=CC3=CSC(=N3)C)C,C[C@H]1CCC[C@@H]2[C@@H](O2)C[C@H](OC(=O)C[C@H](C(C(=O)[C@@H]([C@H]1O)C)(C)C)O)/C(=C/C3=CSC(=N3)C)/C",
    "Amphidinolide_P_sequence_1": r"C[C@@H]1C(=C)C[C@H]2[C@H]3[C@@H](O3)CC(=C)/C=C/[C@H](OC(=O)C[C@@]1(O2)O)[C@H](C)C(=C)C",
    "Lovastatin_sequence_1": r"CC[C@H](C)C(=O)O[C@H]1C[C@H](C=C2[C@H]1[C@H]([C@H](C=C2)C)CC[C@@H]3C[C@H](CC(=O)O3)O)C",
    "Lovastatin_sequence_2": r"CC[C@H](C)C(=O)O[C@H]1C[C@H](C=C2[C@H]1[C@H]([C@H](C=C2)C)CC[C@@H]3C[C@H](CC(=O)O3)O)C",
    "Dictyostatin_sequence_1": r"C[C@H]1CC[C@H]([C@@H]([C@@H](OC(=O)/C=C\C=C\[C@H]([C@H](C[C@@H](/C=C\[C@@H]([C@@H]([C@H](C1)C)O)C)O)O)C)[C@@H](C)/C=C\C=C)C)O",
    "Discodermolide_sequence_1": r"C[C@H]1[C@@H](OC(=O)[C@@H]([C@H]1O)C)C[C@@H](/C=C\[C@H](C)[C@@H]([C@@H](C)/C=C(/C)\C[C@H](C)[C@H]([C@H](C)[C@H]([C@@H](C)/C=C\C=C)OC(=O)N)O)O)O",
    "Amphidinolide_J_sequence_1": r"CCC/C=C/[C@@H](C)[C@H]1C(/C=C\C([C@H](C=CCCC(=C)[C@H](CC(=O)O1)C)O)C)O",
    "Macrolactin_A_sequence_1": r"CC1CCC/C=C/C=C/C(CC(C/C=C\C=C\C(C/C=C/C=C\C(=O)O1)O)O)O",
    "Thermolide_A_sequence_1": r"C[C@@H]1C[C@H]([C@@H](OC(=O)[C@H](NC(=O)C[C@H](C[C@@H]1O)O)C)[C@@H](C)C[C@H](C)[C@@H]([C@H](C)[C@@H](C[C@H](C)O)OC(=O)C)O)C",
    "Gephyornic_acid_1": r"C[C@@H]1[C@@H](O[C@@](C([C@H]1OC)(C)C)([C@@H](C)C[C@H](C)[C@@H]([C@@]2([C@H](O2)[C@@H](C)C=C(C)C)C)O)O)CC(=O)O"
}

mols = {k: mol_to_fp(Chem.MolFromSmiles(v), 2048, 3) for k, v in mols.items()}


fn_donphan = argv[1]
antibacterial_records = []
not_antibacterial_records = []
antifungal_records = []
not_antifungal_records = []
antiviral_records = []
not_antiviral_records = []
with open(fn_donphan, "r") as fo:
    fo.readline()
    for line in tqdm(fo):
        line = line.strip().split(",")
        mol_id,smiles,antibacterial,antifungal,antiviral,not_antibacterial,not_antifungal,not_antiviral = line
        try: 
            mol = Chem.MolFromSmiles(smiles)
            fp = mol_to_fp(mol, 2048, 3)
            antibacterial_records.append((mol_id, fp))
        except: continue
        if antibacterial != "":
            if int(antibacterial) > 0:
                antibacterial_records.append((mol_id, fp))
        if not_antibacterial != "":
            if int(not_antibacterial) > 0:
                not_antibacterial_records.append((mol_id, fp))
        if antifungal != "":
            if int(antifungal) > 0:
                antifungal_records.append((mol_id, fp))
        if not_antifungal != "":
            if int(not_antifungal) > 0:
                not_antifungal_records.append((mol_id, fp))
        if antiviral != "":
            if int(antiviral) > 0:
                antiviral_records.append((mol_id, fp))
        if not_antiviral != "":
            if int(not_antiviral) > 0:
                not_antiviral_records.append((mol_id, fp))


for k, v in mols.items():
    sims = [sim(v, fp) > 0.9 for mol_id, fp in antibacterial_records]
    print(f"antibacterial {k}:", sum(sims))
    sims = [sim(v, fp) > 0.9 for mol_id, fp in not_antibacterial_records]
    print(f"not_antibacterial {k}:", sum(sims))
    sims = [sim(v, fp) > 0.9 for mol_id, fp in antifungal_records]
    print(f"antifungal {k}:", sum(sims))
    sims = [sim(v, fp) > 0.9 for mol_id, fp in not_antifungal_records]
    print(f"not_antifungal {k}:", sum(sims))
    sims = [sim(v, fp) > 0.9 for mol_id, fp in antiviral_records]
    print(f"antiviral {k}:", sum(sims))
    sims = [sim(v, fp) > 0.9 for mol_id, fp in not_antiviral_records]
    print(f"not_antiviral {k}:", sum(sims))

