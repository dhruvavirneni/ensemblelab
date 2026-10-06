"""Shared interfaces and result handling for optimization backends."""

from __future__ import annotations

import warnings
from abc import ABC, abstractmethod
from copy import deepcopy
from typing import Any

import numpy as np
from ase import Atoms
from rdkit import Chem

from ..generators import Conformer, Ensemble


class BaseOptimizer(ABC):
    """Abstract optimizer that returns a new, synchronized ensemble"""

    @property
    @abstractmethod
    def method(self) -> str:
        """Canonical name recorded for this optimization backend"""


    @staticmethod
    def _append_history(
        metadata: dict[str, Any],
        process: str,
        method: str,
        **details: Any,
    ) -> None:
        """Append a workflow event to the ensemble history"""
        history = list(metadata.get("history", []))
        history.append(
            {
                "process": process,
                "method": method,
                **{key: value for key, value in details.items() if value is not None},
            })
        metadata["history"] = history

    def optimize(self, ensemble: Ensemble) -> Ensemble:
        """Optimize all conformers and return a new ensemble instance"""
        self._validate_ensemble(ensemble)

        optimized_molecule = Chem.Mol(ensemble.molecule)
        optimized_conformers: list[Conformer] = []
        convergence: dict[int, bool] = {}

        for conformer in ensemble.conformers:
            atoms = conformer.atoms.copy()
            self._set_rdkit_positions(optimized_molecule, conformer.id, atoms.positions)
            energy, converged = self._optimize_conformer(
                optimized_molecule, conformer.id, atoms
            )
            self._set_rdkit_positions(optimized_molecule, conformer.id, atoms.positions)

            optimized_conformers.append(
                Conformer(
                    id=conformer.id,
                    atoms=atoms,
                    _molecule=optimized_molecule,
                    energy=energy,
                    energy_unit="kcal/mol",
                    optimization_method=self.method,
                    optimization_converged=converged,
                )
            )
            convergence[conformer.id] = converged

        unconverged_ids = [
            identifier for identifier, converged in convergence.items() if not converged
        ]
        converged_ids = [
            identifier for identifier, converged in convergence.items() if converged
        ]
        if unconverged_ids:
            warnings.warn(
                "Optimization reached its maximum number of steps before convergence"
                f"for conformer IDs: {unconverged_ids}.",
                RuntimeWarning,
                stacklevel=2,
            )

        metadata = deepcopy(ensemble.metadata)
        # update history of ensemble with optimization details
        self._append_history(
            metadata,
            process="optimization",
            method=self.method,
            n_input_conformers=len(ensemble.conformers),
            n_output_conformers=len(optimized_conformers),
            n_converged=sum(convergence.values()),
            n_unconverged=len(unconverged_ids),
            converged_conformer_ids=converged_ids,
            unconverged_conformer_ids=unconverged_ids,
            energy_unit="kcal/mol",
            **self._history_details(),)
        
        # updated status metadata
        metadata.update({
            "optimization_status": ("optimized" if not unconverged_ids else "not_fully_converged"),
            "energy_status": "computed",
            "energy_unit": "kcal/mol",
            "n_conformers": len(optimized_conformers),
        })

        return Ensemble(
            smiles=ensemble.smiles,
            molecule=optimized_molecule,
            conformers=tuple(optimized_conformers),
            metadata=metadata,
        )

    @abstractmethod
    def _optimize_conformer(
        self, molecule: Chem.Mol, conformer_id: int, atoms: Atoms
    ) -> tuple[float, bool]:
        """Optimize one conformer, updating ``molecule`` or ``atoms`` in place"""

    @abstractmethod
    def _history_details(self) -> dict[str, Any]:
        """Return backend settings to store in optimization provenance"""
        ...

    @staticmethod
    def _validate_ensemble(ensemble: Ensemble) -> None:
        if not isinstance(ensemble, Ensemble):
            raise TypeError("ensemble must be an Ensemble instance")
        if not ensemble.conformers:
            raise ValueError("ensemble must contain at least one conformer")

    @staticmethod
    def _validate_max_steps(max_steps: int) -> None:
        if isinstance(max_steps, bool) or not isinstance(max_steps, int) or max_steps < 1:
            raise ValueError("max_steps must be a positive integer.")

    @staticmethod
    def _validate_fmax(fmax: float) -> None:
        if not isinstance(fmax, (float, int)) or isinstance(fmax, bool) or fmax <= 0:
            raise ValueError("fmax must be a positive number in eV/angstrom.")

    @staticmethod
    def _resolve_charge(ensemble: Ensemble, charge: int | None) -> int:
        if charge is None:
            return int(Chem.GetFormalCharge(ensemble.molecule))
        if isinstance(charge, bool) or not isinstance(charge, int):
            raise TypeError("charge must be an integer or None.")
        return charge

    @staticmethod
    def _resolve_multiplicity(ensemble: Ensemble, multiplicity: int | None) -> int:
        if multiplicity is None:
            return 1 + sum(
                atom.GetNumRadicalElectrons() for atom in ensemble.molecule.GetAtoms()
            )
        if (
            isinstance(multiplicity, bool)
            or not isinstance(multiplicity, int)
            or multiplicity < 1
        ):
            raise ValueError("multiplicity must be a positive integer or None.")
        return multiplicity

    @staticmethod
    def _set_rdkit_positions(
        molecule: Chem.Mol, conformer_id: int, positions: np.ndarray
    ) -> None:
        rdkit_conformer = molecule.GetConformer(conformer_id)
        for atom_index, (x, y, z) in enumerate(np.asarray(positions, dtype=float)):
            rdkit_conformer.SetAtomPosition(atom_index, (float(x), float(y), float(z)))
