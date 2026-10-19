"""Structural geometry metrics for molecular conformational ensembles."""

from __future__ import annotations

import math
from collections.abc import Sequence
from statistics import fmean, median, pstdev
from typing import Any

import numpy as np
from rdkit import Chem
from rdkit.Chem import Descriptors3D

from ensemblelab.generators import Conformer, Ensemble

try:
	from rdkit.Chem import rdFreeSASA
except ImportError:
	rdFreeSASA = None

RADIUS_OF_GYRATION_UNITS = "angstrom"
SASA_UNITS = "angstrom^2"
SASA_PROBE_RADIUS_A = 1.4
SASA_ALGORITHM = "LeeRichards"

_METRIC_ERRORS = (IndexError, OverflowError, RuntimeError, TypeError, ValueError)


def _calculation_molecule(molecule: Chem.Mol, conformer: Conformer) -> Chem.Mol:
	"""Build a private RDKit molecule using the conformer's ASE coordinates."""
	try:
		n_atoms = molecule.GetNumAtoms()
		atomic_numbers = np.asarray(conformer.atoms.get_atomic_numbers())
		positions = np.asarray(conformer.atoms.get_positions(), dtype=float)
	except (AttributeError, TypeError, ValueError) as exc:
		raise ValueError("conformer coordinates or molecular topology are unavailable.") from exc

	if n_atoms == 0:
		raise ValueError("molecule has no atoms.")
	if positions.shape != (n_atoms, 3):
		raise ValueError(
			f"expected coordinates with shape ({n_atoms}, 3), got {positions.shape}."
		)
	if atomic_numbers.shape != (n_atoms,):
		raise ValueError(
			f"expected {n_atoms} atomic numbers, got {atomic_numbers.shape}."
		)
	if not np.array_equal(
		atomic_numbers,
		np.asarray([atom.GetAtomicNum() for atom in molecule.GetAtoms()]),
	):
		raise ValueError("conformer atoms do not match the ensemble molecular topology.")
	if not np.isfinite(positions).all():
		raise ValueError("conformer coordinates contain non-finite values.")

	try:
		calculation_molecule = Chem.Mol(molecule)
		calculation_molecule.RemoveAllConformers()
		rdkit_conformer = Chem.Conformer(n_atoms)
		for atom_index, (x, y, z) in enumerate(positions):
			rdkit_conformer.SetAtomPosition(
				atom_index, (float(x), float(y), float(z))
			)
		calculation_molecule.AddConformer(rdkit_conformer, assignId=True)
	except (IndexError, OverflowError, RuntimeError, TypeError) as exc:
		raise ValueError(f"could not prepare coordinates for RDKit: {exc}") from exc
	return calculation_molecule


def _sasa_radii(molecule: Chem.Mol) -> list[float]:
	periodic_table = Chem.GetPeriodicTable()
	radii = [
		float(periodic_table.GetRvdw(atom.GetAtomicNum()))
		for atom in molecule.GetAtoms()
	]
	if any(not math.isfinite(radius) or radius <= 0.0 for radius in radii):
		raise ValueError("RDKit does not provide a positive van der Waals radius for every atom.")
	return radii


def _summary_statistics(values: list[float], prefix: str, units: str) -> dict[str, Any]:
	if not values:
		return {
			f"{prefix}_mean_{units}": None,
			f"{prefix}_median_{units}": None,
			f"{prefix}_std_{units}": None,
			f"{prefix}_min_{units}": None,
			f"{prefix}_max_{units}": None,
			f"{prefix}_range_{units}": None,
		}
	minimum = min(values)
	maximum = max(values)
	return {
		f"{prefix}_mean_{units}": fmean(values),
		f"{prefix}_median_{units}": median(values),
		f"{prefix}_std_{units}": pstdev(values),
		f"{prefix}_min_{units}": minimum,
		f"{prefix}_max_{units}": maximum,
		f"{prefix}_range_{units}": maximum - minimum,
	}


def analyze_geometry(
	ensemble: Ensemble,
	conformers: Sequence[Conformer] | None = None,
) -> dict[str, Any]:
	"""Calculate per-conformer radius of gyration and solvent-accessible area.

	Parameters
	----------
	ensemble
		Source molecular ensemble. By default, all conformers are analyzed in
		``ensemble.conformers`` order.
	conformers
		Optional ordered subset of the identical conformer objects in
		``ensemble``. Source indices always refer to their positions in the
		full ensemble.

	Returns
	-------
	dict
		A mapping with ordered ``conformers`` records and one ``summary``.
		Radius of gyration is mass-weighted over all atoms in the molecular
		graph, including explicit hydrogens, and is reported in angstroms.
		SASA is in square angstroms, using RDKit periodic-table van der Waals
		radii, a 1.4 angstrom probe, and the Lee-Richards FreeSASA algorithm.
		Metric values that cannot be calculated are ``None`` and carry an
		explicit validity flag, status, and error message. Summary standard
		deviation is the population standard deviation (denominator *n*), so
		it is zero for one valid value; statistics are ``None`` when there are
		no valid values.

	Notes
	-----
	The calculations use private RDKit molecule copies populated from each
	conformer's ASE coordinates, so neither representation in the input
	ensemble is changed. If the installed RDKit build does not provide
	``rdFreeSASA``, SASA is marked unavailable rather than approximated.
	Radius of gyration describes mass distribution about the center of mass;
	neither it nor SASA alone measures molecular shape or conformational
	diversity.

	Raises
	------
	ValueError
		If selected conformers do not have unique IDs, or are not the exact
		conformer objects owned by ``ensemble``.
	"""
	selected = list(ensemble.conformers if conformers is None else conformers)
	indices_by_id = {
		conformer.id: index for index, conformer in enumerate(ensemble.conformers)
	}
	selected_ids = [conformer.id for conformer in selected]
	if len(selected_ids) != len(set(selected_ids)):
		raise ValueError("conformers must have unique IDs.")
	if any(
		indices_by_id.get(conformer.id) is None
		or ensemble.conformers[indices_by_id[conformer.id]] is not conformer
		for conformer in selected
	):
		raise ValueError("conformers must belong to the supplied ensemble.")

	sasa_available = rdFreeSASA is not None
	rows: list[dict[str, Any]] = []
	valid_radius_values: list[float] = []
	valid_sasa_values: list[float] = []

	for conformer in selected:
		row: dict[str, Any] = {
			"conformer_id": conformer.id,
			"conformer_index": indices_by_id[conformer.id],
			"radius_of_gyration_A": None,
			"radius_of_gyration_valid": False,
			"radius_of_gyration_status": "invalid",
			"radius_of_gyration_error": None,
			"sasa_A2": None,
			"sasa_valid": False,
			"sasa_status": "invalid",
			"sasa_error": None,
		}
		try:
			calculation_molecule = _calculation_molecule(ensemble.molecule, conformer)
		except ValueError as exc:
			row["radius_of_gyration_error"] = str(exc)
			if sasa_available:
				row["sasa_error"] = str(exc)
			else:
				row["sasa_status"] = "unavailable"
				row["sasa_error"] = "RDKit FreeSASA is unavailable."
			rows.append(row)
			continue

		try:
			radius_of_gyration = float(
				Descriptors3D.RadiusOfGyration(
					calculation_molecule,
					confId=0,
					useAtomicMasses=True,
					force=True,
				)
			)
			if not math.isfinite(radius_of_gyration):
				raise ValueError("calculation returned a non-finite value.")
			row["radius_of_gyration_A"] = radius_of_gyration
			row["radius_of_gyration_valid"] = True
			row["radius_of_gyration_status"] = "valid"
			valid_radius_values.append(radius_of_gyration)
		except _METRIC_ERRORS as exc:
			row["radius_of_gyration_error"] = (
				f"{type(exc).__name__}: {exc}"
			)

		if not sasa_available:
			row["sasa_status"] = "unavailable"
			row["sasa_error"] = "RDKit FreeSASA is unavailable."
		else:
			try:
				options = rdFreeSASA.SASAOpts(
					rdFreeSASA.SASAAlgorithm.LeeRichards,
					rdFreeSASA.SASAClassifier.Protor,
					SASA_PROBE_RADIUS_A,
				)
				sasa = float(
					rdFreeSASA.CalcSASA(
						calculation_molecule,
						_sasa_radii(calculation_molecule),
						confIdx=0,
						opts=options,
					)
				)
				if not math.isfinite(sasa):
					raise ValueError("calculation returned a non-finite value.")
				row["sasa_A2"] = sasa
				row["sasa_valid"] = True
				row["sasa_status"] = "valid"
				valid_sasa_values.append(sasa)
			except _METRIC_ERRORS as exc:
				row["sasa_error"] = f"{type(exc).__name__}: {exc}"

		rows.append(row)

	summary: dict[str, Any] = {
		"n_conformers": len(selected),
		"n_valid_radius_of_gyration": len(valid_radius_values),
		"n_invalid_radius_of_gyration": len(selected) - len(valid_radius_values),
		"n_valid_sasa": len(valid_sasa_values),
		"n_invalid_sasa": len(selected) - len(valid_sasa_values),
		"radius_of_gyration_units": RADIUS_OF_GYRATION_UNITS,
		"radius_of_gyration_method": "RDKit Descriptors3D.RadiusOfGyration",
		"radius_of_gyration_mass_weighted": True,
		"radius_of_gyration_atom_scope": (
			"all atoms in the molecular graph, including explicit hydrogens"
		),
		"sasa_units": SASA_UNITS,
		"sasa_method": "RDKit FreeSASA",
		"sasa_available": sasa_available,
		"sasa_atomic_radius_model": "RDKit periodic-table van der Waals radii (GetRvdw)",
		"sasa_atom_scope": (
			"all atoms in the molecular graph, including explicit hydrogens"
		),
		"sasa_probe_radius_A": SASA_PROBE_RADIUS_A,
		"sasa_algorithm": SASA_ALGORITHM,
	}
	summary.update(
		_summary_statistics(
			valid_radius_values, "radius_of_gyration", "A"
		)
	)
	summary.update(_summary_statistics(valid_sasa_values, "sasa", "A2"))

	return {"conformers": rows, "summary": summary}


__all__ = ["analyze_geometry"]
