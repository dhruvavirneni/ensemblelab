"""GFN2-xTB conformer optimization through TBLite and ASE"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

from ase import Atoms
from ase.optimize import BFGS
from rdkit import Chem

from ..generators import Ensemble
from .base import BaseOptimizer

_EV_TO_KCAL_PER_MOL = 23.0605478306


@dataclass(frozen=True, slots=True)
class GFN2xTBOptimizer(BaseOptimizer):
    """Optimize conformers with GFN2-xTB using TBLite's ASE calculator."""

    fmax: float = 0.05
    max_steps: int = 500
    charge: int | None = None
    multiplicity: int | None = None
    solvent: str | None = None
    _resolved_charge: int = field(init=False, repr=False, compare=False)
    _resolved_multiplicity: int = field(init=False, repr=False, compare=False)

    def __post_init__(self) -> None:
        self._validate_fmax(self.fmax)
        self._validate_max_steps(self.max_steps)
        if self.solvent is not None and (
            not isinstance(self.solvent, str) or not self.solvent.strip()
        ):
            raise ValueError("solvent must be a non-empty string or None.")

    @property
    def method(self) -> str:
        """Canonical optimization method name."""
        return "GFN2-xTB"

    def optimize(self, ensemble: Ensemble) -> Ensemble:
        """Optimize all conformers after resolving the ensemble's electronic state."""
        resolved_charge = self._resolve_charge(ensemble, self.charge)
        resolved_multiplicity = self._resolve_multiplicity(ensemble, self.multiplicity)
        object.__setattr__(self, "_resolved_charge", resolved_charge)
        object.__setattr__(self, "_resolved_multiplicity", resolved_multiplicity)
        try:
            return BaseOptimizer.optimize(self, ensemble)
        finally:
            object.__delattr__(self, "_resolved_charge")
            object.__delattr__(self, "_resolved_multiplicity")

    def _optimize_conformer(
        self, molecule: Chem.Mol, conformer_id: int, atoms: Atoms
    ) -> tuple[float, bool]:
        del molecule, conformer_id
        try:
            from tblite.ase import TBLite # type: ignore # noqa: I001
        except ImportError as error:
            raise ImportError(
                "GFN2-xTB optimization requires the optional 'tblite' package. "
                "Install ensemblelab with its xTB extra (upcoming)."
            ) from error

        calculator_options: dict[str, Any] = {
            "method": self.method,
            "charge": self._resolved_charge,
            "multiplicity": self._resolved_multiplicity,
        }
        if self.solvent is not None:
            calculator_options["solvation"] = ("alpb", self.solvent)
        atoms.calc = TBLite(**calculator_options)
        optimizer = BFGS(atoms, logfile=None)
        converged = bool(optimizer.run(fmax=self.fmax, steps=self.max_steps))
        return float(atoms.get_potential_energy()) * _EV_TO_KCAL_PER_MOL, converged

    def _history_details(self) -> dict[str, Any]:
        return {
            "fmax_eV_per_angstrom": self.fmax,
            "max_steps": self.max_steps,
            "solvent": self.solvent,
            "charge": self._resolved_charge,
            "multiplicity": self._resolved_multiplicity,
        }
