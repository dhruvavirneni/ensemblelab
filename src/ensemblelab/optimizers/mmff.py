"""RDKit MMFF and UFF conformer optimizers"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Literal

from ase import Atoms
from rdkit import Chem
from rdkit.Chem import AllChem

from ..generators import sync_atoms_from_rdkit
from .base import BaseOptimizer


@dataclass(frozen=True, slots=True)
class MMFFOptimizer(BaseOptimizer):
    """Optimize every conformer with RDKit's MMFF force field."""

    max_steps: int = 500

    def __post_init__(self) -> None:
        self._validate_max_steps(self.max_steps)

    @property
    def method(self) -> str:
        """Canonical optimization method name."""
        return "MMFF"

    def _optimize_conformer(
        self, molecule: Chem.Mol, conformer_id: int, atoms: Atoms
    ) -> tuple[float, bool]:
        properties = AllChem.MMFFGetMoleculeProperties(molecule)
        if properties is None:
            raise RuntimeError("MMFF parameters are unavailable for this molecule.")
        force_field = AllChem.MMFFGetMoleculeForceField(
            molecule, properties, confId=conformer_id
        )
        if force_field is None:
            raise RuntimeError(
                f"MMFF could not construct a force field for conformer {conformer_id}."
            )
        status = force_field.Minimize(maxIts=self.max_steps)
        if status not in {0, 1}:
            raise RuntimeError(
                f"MMFF optimization failed for conformer {conformer_id} (status {status})."
            )
        sync_atoms_from_rdkit(molecule, conformer_id, atoms)
        return float(force_field.CalcEnergy()), status == 0

    def _history_details(self) -> dict[str, Any]:
        return {
            "max_steps": self.max_steps,
        }


@dataclass(frozen=True, slots=True)
class UFFOptimizer(MMFFOptimizer):
    """Optimize every conformer with RDKit's UFF force field."""

    @property
    def method(self) -> Literal["UFF"]:
        """Canonical optimization method name."""
        return "UFF"

    def _optimize_conformer(self, molecule: Chem.Mol, conformer_id: int, atoms: Atoms) -> tuple[float, bool]:
        force_field = AllChem.UFFGetMoleculeForceField(molecule, confId=conformer_id)
        if force_field is None:
            raise RuntimeError(
                f"UFF could not construct a force field for conformer {conformer_id}."
            )
        status = force_field.Minimize(maxIts=self.max_steps)
        if status not in {0, 1}:
            raise RuntimeError(
                f"UFF optimization failed for conformer {conformer_id} (status {status})."
            )
        sync_atoms_from_rdkit(molecule, conformer_id, atoms)
        return float(force_field.CalcEnergy()), status == 0

    def _history_details(self) -> dict[str, Any]:
        return {
            "max_steps": self.max_steps,
        }