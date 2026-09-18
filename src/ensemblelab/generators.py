"""Conformer and Ensemble generation backed by RDKit ETKDG."""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass, field
from typing import TYPE_CHECKING, Any

import numpy as np
from ase import Atoms
from rdkit import Chem, rdBase
from rdkit.Chem import AllChem

from .filters import BaseFilter, CompositeFilter, PopulationFilter

if TYPE_CHECKING:
    from .optimizers import BaseOptimizer


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
    energy_unit: str | None = None
    optimization_method: str | None = None
    optimization_converged: bool | None = None
    _molecule: Chem.Mol | None = field(default=None, repr=False, compare=False, hash=False)

    def show(self) -> None:
        """Display a concise, human-readable conformer summary."""
        from .display.summaries import conformer_summary

        print(conformer_summary(self))

# generation history generator helper
def _generation_history(
    smiles: str,
    canonical_smiles: str,
    n_requested: int,
    n_generated: int,
) -> dict[str, Any]:
    """Create the initial provenance record for an ensemble."""

    return {
        "process": "generation",
        "method": "rdkit.ETKDGv3",
        "requested_smiles": smiles,
        "canonical_smiles": canonical_smiles,
        "n_requested": n_requested,
        "n_generated": n_generated,
        "random_seed": 42,
        "rdkit_version": rdBase.rdkitVersion,
    }


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
    def from_smiles(cls, smiles: str, n_confs: int = 25) -> Ensemble:
        """Create an unoptimized ensemble from SMILES.

        This is the object oriented equivalent of "generate" function and keeps
        the public workflow ready for subsequent "optimize" and "cluster"
        methods.
        """
        return generate(smiles, n_confs=n_confs)

    @classmethod
    def get_mol(self) -> Chem.Mol:
        """Return the RDKit molecule associated with this ensemble."""
        return self.molecule

    @property
    def conformer_ids(self) -> tuple[int, ...]:
        """RDKit conformer IDs in the alignment order used by ensemble data."""
        return tuple(conformer.id for conformer in self.conformers)

    def rdkit_conformer(self, conformer_id: int) -> Chem.Conformer:
        """Return an RDKit conformer by its stable ID."""
        if conformer_id not in self.conformer_ids:
            raise KeyError(f"Unknown conformer ID: {conformer_id}")
        return self.molecule.GetConformer(conformer_id)

    def show(
        self,
        *,
        history: bool = False,
        metadata: bool = False,
        conformers: bool = True,
    ) -> None:
        """Display a human-readable view of this ensemble.

        The default view contains the ensemble summary and conformer table.
        Set ``history`` or ``metadata`` to include the corresponding stored
        provenance sections. Rendering is read-only and never triggers
        expensive analysis or optimization.
        """
        from .display.summaries import (
            ensemble_history,
            ensemble_metadata,
            ensemble_summary,
        )

        sections = [ensemble_summary(self, include_conformers=conformers)]
        if history:
            sections.append(ensemble_history(self))
        if metadata:
            sections.append(ensemble_metadata(self))
        print("\n\n".join(sections))


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

    canonical_smiles = Chem.MolToSmiles(Chem.RemoveHs(molecule), canonical=True)

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
            _molecule=molecule,
        )
        for conformer_id in conformer_ids
    )
    metadata = {
    "n_conformers": len(conformers),
    "optimization_status": "unoptimized",
    "energy_status": "uncomputed",
    "energy_unit": None,
    "rdkit_version": rdBase.rdkitVersion,
    "canonical_smiles": canonical_smiles,
    "history": [
        _generation_history(
            smiles,
            canonical_smiles,
            n_confs,
            len(conformers),
            )
        ],
    }
    return Ensemble(
        smiles=smiles,
        molecule=molecule,
        conformers=conformers,
        metadata=metadata,
    )

# later --- ASE to SMILES can be unpredictable
# if user supplies SMILES and Atoms, interpret Atoms in context of SMILES molecule, with no need for helper function
def setup(
    molecule: str | Atoms,
    n_confs: int,
    *,
    optimizers: Sequence[BaseOptimizer] | None = None,
    filters: Sequence[BaseFilter] | None = None,
    target_conformers: int | None = None,
) -> Ensemble:
    """Quick generation, optimization, and filtration for a new Ensemble."""

    from .optimizers import (
        BaseOptimizer,
        GFN2xTBOptimizer,
        HierarchicalOptimizer,
        MMFFOptimizer,
    )

    # Validate all inputs before beginning computational work.
    if not isinstance(molecule, (str, Atoms)):
        raise TypeError("molecule must be a SMILES string or ASE Atoms object.")

    if isinstance(molecule, str):
        if not molecule.strip():
            raise ValueError("SMILES string cannot be empty.")
        if Chem.MolFromSmiles(molecule) is None:
            raise ValueError(f"Invalid SMILES string: {molecule!r}")
    elif len(molecule) == 0:
        raise ValueError("ASE Atoms object cannot be empty.")

    if not isinstance(n_confs, int):
        raise TypeError("n_confs must be an integer.")
    if n_confs < 1:
        raise ValueError("n_confs must be at least 1.")

    if target_conformers is not None:
        if not isinstance(target_conformers, int):
            raise TypeError("target_conformers must be an integer.")
        if target_conformers < 1:
            raise ValueError("target_conformers must be at least 1.")
        if target_conformers > n_confs:
            raise ValueError("target_conformers cannot exceed n_confs.")

    if optimizers is not None:
        if not optimizers:
            raise ValueError("optimizers cannot be empty.")
        if not all(isinstance(optimizer, BaseOptimizer) for optimizer in optimizers):
            raise TypeError("optimizers must contain BaseOptimizer instances.")

    if filters is not None:
        if not filters:
            raise ValueError("filters cannot be empty.")
        if not all(isinstance(filter_, BaseFilter) for filter_ in filters):
            raise TypeError("filters must contain BaseFilter instances.")

    # Convert ASE input to SMILES so the existing generation pipeline can be used.
    if isinstance(molecule, Atoms):
        molecule = _atoms_to_smiles(molecule)

    # Use the standard pipeline unless the user supplies custom components.
    if optimizers is None:
        optimizers = (MMFFOptimizer(), GFN2xTBOptimizer())

    if filters is None:
        filters = (
            PopulationFilter(target_conformers=target_conformers),
        )

    ensemble = generate(molecule, n_confs=n_confs)

    # Run the configured optimizers as a hierarchical workflow.
    ensemble = HierarchicalOptimizer(
        optimizers=tuple(optimizers)
    ).optimize(ensemble)

    # Apply filters sequentially through the composite filter.
    ensemble = CompositeFilter(
        filters=tuple(filters)
    ).apply(ensemble)

    return ensemble

def _atoms_to_smiles(atoms: Atoms) -> str:
    """Helper for converting an ASE Atoms object to a canonical SMILES string."""

    if not isinstance(atoms, Atoms):
        raise TypeError("atoms must be an ASE Atoms object.")

    if len(atoms) == 0:
        raise ValueError("ASE Atoms object cannot be empty.")

    symbols = atoms.get_chemical_symbols()
    positions = atoms.get_positions()

    xyz = [
        str(len(atoms)),
        "Generated from ASE Atoms",
    ]

    xyz.extend(
        f"{symbol} {x:.8f} {y:.8f} {z:.8f}"
        for symbol, (x, y, z) in zip(symbols, positions)
    )

    molecule = Chem.MolFromXYZBlock("\n".join(xyz))

    if molecule is None:
        raise ValueError("Could not convert ASE Atoms to an RDKit molecule.")

    try:
        Chem.rdDetermineBonds.DetermineBonds(molecule)
        Chem.SanitizeMol(molecule)
    except Exception as exc:
        raise ValueError(
            "Could not determine valid molecular connectivity from the "
            "ASE Atoms object."
        ) from exc

    return Chem.MolToSmiles(molecule)