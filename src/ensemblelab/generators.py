"""Conformer-ensemble generation backed by RDKit ETKDG."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np
from ase import Atoms
from rdkit import Chem, rdBase
from rdkit.Chem import AllChem


@dataclass(slots=True)
class Conformer:
    """One generated conformer and its conformer-level data.

    ``energy`` is ``None`` until an optimization backend assigns a value. It
    belongs to the conformer rather than a parallel ensemble-level array so a
    conformer's geometry and computed properties remain together.
    """

    id: int
    atoms: Atoms
    energy: float | None = None


@dataclass(slots=True)
class Ensemble:
    """A molecular conformational ensemble with aligned conformer-level data.

    Generation initializes structures and provenance only. Every conformer has
    ``energy=None`` until optimization assigns a physical energy.
    """

    smiles: str
    molecule: Chem.Mol
    conformers: tuple[Conformer, ...]
    metadata: dict[str, Any]

    def __post_init__(self) -> None:
        ids = [conformer.id for conformer in self.conformers]
        if len(ids) != len(set(ids)):
            raise ValueError("conformer IDs must be unique within an ensemble.")

    @classmethod
    def from_smiles(cls, smiles: str, n_confs: int = 25) -> "Ensemble":
        """Create an unoptimized ensemble from SMILES.

        This is the object oriented equivalent of :func:`generate` and keeps
        the public workflow ready for subsequent ``optimize`` and ``cluster``
        methods.
        """
        return generate(smiles, n_confs=n_confs)

    @property
    def conformer_ids(self) -> tuple[int, ...]:
        """RDKit conformer IDs in the alignment order used by ensemble data."""
        return tuple(conformer.id for conformer in self.conformers)

    def rdkit_conformer(self, conformer_id: int) -> Chem.Conformer:
        """Return an RDKit conformer by its stable ID."""
        if conformer_id not in self.conformer_ids:
            raise KeyError(f"Unknown conformer ID: {conformer_id}")
        return self.molecule.GetConformer(conformer_id)


def generate(smiles: str, n_confs: int = 25) -> Ensemble:
    """Embed an unoptimized conformational ensemble from a SMILES string.

    Uses the notebook's RDKit ETKDGv3 setup (explicit seed 42) on a
    hydrogenated molecule. Optimization and energies intentionally belong to
    the next module.
    """
    if not isinstance(smiles, str) or not smiles.strip():
        raise ValueError("smiles must be a non-empty SMILES string.")
    if isinstance(n_confs, bool) or not isinstance(n_confs, int) or n_confs < 1:
        raise ValueError("n_confs must be a positive integer.")

    molecule = Chem.MolFromSmiles(smiles)
    if molecule is None or molecule.GetNumAtoms() == 0:
        raise ValueError(f"Could not parse SMILES: {smiles!r}")
    molecule = Chem.AddHs(molecule)

    params = AllChem.ETKDGv3()
    params.randomSeed = 42
    conformer_ids = tuple(
        int(identifier)
        for identifier in AllChem.EmbedMultipleConfs(
            molecule, numConfs=n_confs, params=params
        )
    )
    if not conformer_ids:
        raise ValueError("RDKit could not embed any conformers for this molecule.")

    symbols = tuple(atom.GetSymbol() for atom in molecule.GetAtoms())
    conformers = tuple(
        Conformer(
            id=conformer_id,
            atoms=Atoms(
                symbols=symbols,
                positions=np.asarray(
                    molecule.GetConformer(conformer_id).GetPositions(), dtype=float
                ).copy(),
            ),
        )
        for conformer_id in conformer_ids
    )
    metadata: dict[str, Any] = {
        "generator": "rdkit.ETKDGv3",
        "requested_smiles": smiles,
        "canonical_smiles": Chem.MolToSmiles(Chem.RemoveHs(molecule), canonical=True),
        "n_conformers_requested": n_confs,
        "n_conformers_generated": len(conformers),
        "random_seed": 42,
        "optimization_status": "unoptimized",
        "energy_status": "uncomputed",
        "rdkit_version": rdBase.rdkitVersion,
    }
    return Ensemble(
        smiles=smiles,
        molecule=molecule,
        conformers=conformers,
        metadata=metadata,
    )
