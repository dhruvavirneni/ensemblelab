from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING

from ensemblelab.filters.base import BaseFilter

if TYPE_CHECKING:
    from ensemblelab.generators import Ensemble


@dataclass(frozen=True, slots=True)
class EnergyFilter(BaseFilter):
    window: float | None = None
    top_n: int | None = None

    @property
    def method(self) -> str:
        return "energy"

    def __post_init__(self) -> None:
        if (self.window is None) == (self.top_n is None):
            raise ValueError("Specify exactly one of window or top_n.")

    def apply(self, ensemble: Ensemble) -> Ensemble:
        self._validate_ensemble(ensemble)
        if self.window is not None:
            ensemble_energies = [conformer.energy for conformer in ensemble.conformers if conformer.energy is not None]
            min_energy = min(ensemble_energies)
            filtered_conformers = [conformer for conformer in ensemble.conformers if conformer.energy is not None and conformer.energy - min_energy <= self.window]
            return self._build_filtered_ensemble(ensemble, filtered_conformers, {"window": self.window})
        elif self.top_n is not None:
            valid_conformers = [conformer for conformer in ensemble.conformers if conformer.energy is not None]
            sorted_conformers = sorted(valid_conformers, key=lambda c: c.energy)
            filtered_conformers = sorted_conformers[: self.top_n]
            return self._build_filtered_ensemble(ensemble, filtered_conformers, {"top_n": self.top_n})

    def _history_details(self) -> dict[str, float | int | None]:
        """Return the selected energy filter criterion for provenance."""
        return {"window": self.window, "top_n": self.top_n}
