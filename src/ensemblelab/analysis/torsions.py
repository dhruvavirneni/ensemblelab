from collections.abc import Iterable
from dataclasses import dataclass

import numpy as np
from rdkit import Chem
from rdkit.Chem import rdMolTransforms

from ensemblelab.generators import Ensemble


@dataclass
class Torsion:
    """Class representing a torsion defined by four atom indices."""

    atom_indices: tuple[int, int, int, int]
    central_bond: tuple[int, int] | None = None

    def __post_init__(self):
        if len(self.atom_indices) != 4:
            raise ValueError("Torsion must be defined by exactly four atom indices.")


def find_heavy_atom_torsions(mol: Ensemble) -> dict[str, Torsion]:
    mol = mol.get_mol()
    torsions = {}

    for bond in mol.GetBonds():
        if bond.GetBondType() != Chem.BondType.SINGLE:
            continue

        if bond.IsInRing():
            continue

        atom_2 = bond.GetBeginAtom()
        atom_3 = bond.GetEndAtom()

        if atom_2.GetAtomicNum() == 1 or atom_3.GetAtomicNum() == 1:
            continue

        left_neighbors = [
            atom.GetIdx()
            for atom in atom_2.GetNeighbors()
            if atom.GetIdx() != atom_3.GetIdx() and atom.GetAtomicNum() > 1
        ]

        right_neighbors = [
            atom.GetIdx()
            for atom in atom_3.GetNeighbors()
            if atom.GetIdx() != atom_2.GetIdx() and atom.GetAtomicNum() > 1
        ]

        if not left_neighbors or not right_neighbors:
            continue

        atom_1 = left_neighbors[0]
        atom_4 = right_neighbors[0]

        name = f"torsion_{len(torsions) + 1}_{atom_1}_{atom_2.GetIdx()}_{atom_3.GetIdx()}_{atom_4}"
        torsions[name] = Torsion(atom_indices=(atom_1, atom_2.GetIdx(), atom_3.GetIdx(), atom_4))

    return torsions



def torsion_angles(
    ensemble: Ensemble, torsions: Iterable[Torsion] | None = None
) -> np.ndarray:
    """Compute torsion angle matrix for an Ensemble.
    Returns a 2d numpy array of shape (n_torsions, n_conformers) where each row corresponds to a torsion and each column corresponds to a conformer.
    """
    if torsions is None:
        torsions = find_heavy_atom_torsions(ensemble).values()
    torsions = tuple(torsions)
    result = np.zeros((len(torsions), len(ensemble.conformers)), dtype=float)

    for conformer_index, conformer in enumerate(ensemble.conformers):
        rdkit_conformer = ensemble.rdkit_conformer(conformer.id)
        for torsion_index, torsion in enumerate(torsions):
            result[torsion_index, conformer_index] = rdMolTransforms.GetDihedralDeg(
                rdkit_conformer, *torsion.atom_indices
            )

    return result
    

