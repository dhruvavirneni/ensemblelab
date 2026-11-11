from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING, Any

import numpy as np
from rdkit import Chem
from rdkit.Chem import rdMolAlign

from ensemblelab.filters.base import BaseFilter

if TYPE_CHECKING:
    from ensemblelab.generators import Ensemble

# Max allowed difference (Å) between RDKit and ASE coordinates of a conformer.
_SYNC_TOLERANCE_ANGSTROM = 1e-4


@dataclass(frozen=True, slots=True)
class RMSDFilter(BaseFilter):
    """Remove geometric near-duplicates, keeping the lowest-energy representative.

    Conformers are visited in ascending energy order (ties keep input order).
    Each one is compared to every conformer already kept; if its RMSD to any of
    them is ``<= cutoff`` it is discarded as a duplicate, otherwise it is kept.
    The lowest-energy conformer is therefore always kept.

    RMSD conventions
    ----------------
    - Alignment: ``rdkit.Chem.rdMolAlign.GetBestRMS``, an optimal rigid-body
      superposition that also searches symmetry-equivalent atom mappings
      (e.g. permuted methyl hydrogens). All conformers share one RDKit
      molecule, so atom ordering is identical by construction.
    - Atom selection: all atoms by default. With ``heavy_atoms_only=True``,
      hydrogens are removed with ``Chem.RemoveHs`` first, so rotating an OH or
      other non-symmetric hydrogen does not inflate the RMSD and true heavy-atom
      duplicates merge. ``RemoveHs`` keeps a few special hydrogens (e.g.
      isotopic or stereo-defining ones).
    - Coordinates are read from the RDKit molecule, which must match each
      ``Conformer.atoms``; a mismatch raises ``ValueError``.

    The input ensemble is never modified. RMSDs are computed on a copy of the
    molecule because ``GetBestRMS`` leaves the probe conformer aligned in place.

    Parameters
    ----------
    cutoff
        RMSD threshold in Å. Pairs at or below it count as duplicates.
    heavy_atoms_only
        Compute RMSD over heavy atoms only. Defaults to ``False`` (all atoms).
    """

    cutoff: float = 0.5  # Å
    heavy_atoms_only: bool = False

    @property
    def method(self) -> str:
        return "rmsd"

    def __post_init__(self) -> None:
        if self.cutoff is None:
            raise ValueError("RMSD filtering threshold (cutoff) must be specified")
        if self.cutoff <= 0:
            raise ValueError("RMSD cutoff must be greater than 0 Å.")
        if not isinstance(self.heavy_atoms_only, bool):
            raise TypeError("heavy_atoms_only must be a bool.")

    def apply(self, ensemble: Ensemble) -> Ensemble:
        self._validate_ensemble(ensemble)
        if any(conformer.energy is None for conformer in ensemble.conformers):
            raise ValueError(
                "RMSD filtering requires all conformers to have computed energies.")
        self._check_coordinates_synced(ensemble)
        energies = [float(conformer.energy) for conformer in ensemble.conformers]

        # list of Conformer objects sorted by energy (stable: ties keep input order)
        sorted_conformers = [conformer for conformer, energy in sorted(zip(ensemble.conformers, energies), key=lambda item: item[1], reverse=False)]

        # Work on a copy: GetBestRMS moves the probe conformer, and the input
        # molecule is shared with the input ensemble. RemoveHs also returns a
        # new molecule and preserves conformer IDs.
        reference = (
            Chem.RemoveHs(ensemble.molecule)
            if self.heavy_atoms_only
            else Chem.Mol(ensemble.molecule)
        )

        approved_list = [sorted_conformers[0]]
        duplicate_of: dict[int, dict[str, float | int]] = {}

        for conformer in sorted_conformers[1:]:
            is_duplicate = False

            for approved in approved_list:
                rmsd = rdMolAlign.GetBestRMS(
                    reference,
                    reference,
                    prbId=conformer.id,
                    refId=approved.id,
                )

                if rmsd <= self.cutoff:
                    is_duplicate = True
                    duplicate_of[conformer.id] = {
                        "kept_conformer_id": approved.id,
                        "rmsd_angstrom": float(rmsd),
                    }
                    break
            if not is_duplicate:
                approved_list.append(conformer)

        return self._build_filtered_ensemble(
            ensemble,
            approved_list,
            {
                "kept_conformer_ids": [conformer.id for conformer in approved_list],
                "removed_conformer_ids": list(duplicate_of),
                "duplicate_of": duplicate_of,
            },
        )

    @staticmethod
    def _check_coordinates_synced(ensemble: Ensemble) -> None:
        """Ensure RDKit conformer coordinates match each Conformer's ASE atoms.

        RMSD is computed from the RDKit molecule, so stale RDKit coordinates
        (e.g. a backend that updated only ASE positions) would give silently
        wrong duplicate decisions.
        """
        for conformer in ensemble.conformers:
            rdkit_positions = ensemble.molecule.GetConformer(conformer.id).GetPositions()
            ase_positions = np.asarray(conformer.atoms.positions, dtype=float)
            if rdkit_positions.shape != ase_positions.shape or not np.allclose(
                rdkit_positions, ase_positions, rtol=0.0, atol=_SYNC_TOLERANCE_ANGSTROM
            ):
                raise ValueError(
                    f"Conformer {conformer.id}: RDKit coordinates do not match "
                    "Conformer.atoms. RMSD filtering reads RDKit coordinates, so "
                    "the ensemble must be synchronized first."
                )

    def _history_details(self) -> dict[str, Any]:
        return {
            "rmsd_cutoff_angstrom": self.cutoff,
            "heavy_atoms_only": self.heavy_atoms_only,
            "alignment": "rdkit.rdMolAlign.GetBestRMS (symmetry-aware)",
            "energy_ordering": "ascending",
        }
