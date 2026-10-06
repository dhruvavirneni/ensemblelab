from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING

from rdkit.Chem import rdMolAlign

from ensemblelab.filters.base import BaseFilter

if TYPE_CHECKING:
    from ensemblelab.generators import Ensemble


@dataclass(frozen=True, slots=True)
class RMSDFilter(BaseFilter):
    """
    RMSD threshold-based filtering:

    Specify RMSD cutoff in Angstroms (Å)
    Conformer RMSDs are determined via comparison to each previously approved conformer in ensemble, starting with lowest energy conformer.
    """
    cutoff: float = 0.5 # defaults to 50% RMSD similarity threshold

    def __post_init__(self) -> None:
        if self.cutoff is None:
            raise ValueError("RMSD filtering threshold (cutoff) must be specified")
        if self.cutoff <= 0:
            raise ValueError("RMSD cutoff must be greater than 0 Å.")

    def apply(self, ensemble: Ensemble):
        self._validate_ensemble(ensemble)
        if any(conformer.energy is None for conformer in ensemble.conformers):
            raise ValueError(
                "RMSD filtering requires all conformers to have computed energies.")
        energies = [float(conformer.energy) for conformer in ensemble.conformers]

        # list of Conformer objects sorted by energy
        sorted_conformers = [conformer for conformer, energy in sorted(zip(ensemble.conformers, energies), key=lambda item: item[1], reverse=False)]

        # RDKit molecule object
        approved_list = [sorted_conformers[0]]

        for conformer in sorted_conformers[1:]:
            is_duplicate = False

            for approved in approved_list:
                rmsd = rdMolAlign.GetBestRMS(
                    ensemble.molecule,
                    ensemble.molecule,
                    prbId=conformer.id,
                    refId=approved.id,
                )

                if rmsd <= self.cutoff:
                    is_duplicate = True
                    break
            if not is_duplicate:
                approved_list.append(conformer)
        return self._build_filtered_ensemble(
            ensemble,
            approved_list,)


    def _history_details(self) -> dict:
        return {
            "rmsd_cutoff_angstrom": self.cutoff,
        }
  