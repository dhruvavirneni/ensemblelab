"""Torsional angle, circular-statistic, and state-coverage metrics."""

from __future__ import annotations

import math
from bisect import bisect_right
from collections.abc import Sequence
from numbers import Real
from typing import Any

import numpy as np
from rdkit import Chem
from rdkit.Chem import rdMolTransforms

from ensemblelab.generators import Conformer, Ensemble

ANGLE_UNITS = "degrees"
DEFAULT_STATE_BIN_EDGES_DEG = tuple(float(edge) for edge in range(-180, 181, 30))
CIRCULAR_MEAN_RESULTANT_EPSILON = 1e-12

_CALCULATION_ERRORS = (
	IndexError,
	OverflowError,
	RuntimeError,
	TypeError,
	ValueError,
)


def _validate_conformers(
	ensemble: Ensemble,
	conformers: Sequence[Conformer] | None,
) -> tuple[list[Conformer], dict[int, int]]:
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
	return selected, indices_by_id


def _state_edges(edges: Sequence[Real] | None) -> tuple[float, ...]:
	if edges is None:
		return DEFAULT_STATE_BIN_EDGES_DEG
	if isinstance(edges, (str, bytes)):
		raise ValueError("state_bin_edges_deg must be a sequence of numeric edges.")
	try:
		values = tuple(edges)
	except TypeError as exc:
		raise ValueError("state_bin_edges_deg must be a sequence of numeric edges.") from exc
	if len(values) < 2:
		raise ValueError("state_bin_edges_deg must contain at least two edges.")
	if any(
		isinstance(edge, bool)
		or not isinstance(edge, Real)
		or not math.isfinite(float(edge))
		for edge in values
	):
		raise ValueError("state_bin_edges_deg must contain only finite real values.")
	normalized = tuple(float(edge) for edge in values)
	if normalized[0] != -180.0 or normalized[-1] != 180.0:
		raise ValueError("state_bin_edges_deg must start at -180 and end at 180.")
	if any(left >= right for left, right in zip(normalized, normalized[1:])):
		raise ValueError("state_bin_edges_deg must be strictly increasing.")
	return normalized


def _amide_or_ester_bond(bond: Chem.Bond) -> bool:
	"""Whether the bond joins a carbonyl carbon to N or O."""
	endpoints = (bond.GetBeginAtom(), bond.GetEndAtom())
	for carbonyl, heteroatom in (endpoints, endpoints[::-1]):
		if carbonyl.GetAtomicNum() != 6 or heteroatom.GetAtomicNum() not in (7, 8):
			continue
		if any(
			neighbor.GetAtomicNum() == 8
			and carbonyl.GetOwningMol()
			.GetBondBetweenAtoms(carbonyl.GetIdx(), neighbor.GetIdx())
			.GetBondType()
			== Chem.BondType.DOUBLE
			for neighbor in carbonyl.GetNeighbors()
			if neighbor.GetIdx() != heteroatom.GetIdx()
		):
			return True
	return False


def _eligible_torsions(molecule: Chem.Mol) -> list[dict[str, Any]]:
	torsions: list[dict[str, Any]] = []
	for bond in molecule.GetBonds():
		if bond.GetBondType() != Chem.BondType.SINGLE or bond.IsInRing():
			continue
		left_atom = bond.GetBeginAtom()
		right_atom = bond.GetEndAtom()
		if left_atom.GetAtomicNum() == 1 or right_atom.GetAtomicNum() == 1:
			continue
		if _amide_or_ester_bond(bond):
			continue

		left_index, right_index = sorted((left_atom.GetIdx(), right_atom.GetIdx()))
		left_neighbors = sorted(
			atom.GetIdx()
			for atom in molecule.GetAtomWithIdx(left_index).GetNeighbors()
			if atom.GetAtomicNum() > 1 and atom.GetIdx() != right_index
		)
		right_neighbors = sorted(
			atom.GetIdx()
			for atom in molecule.GetAtomWithIdx(right_index).GetNeighbors()
			if atom.GetAtomicNum() > 1 and atom.GetIdx() != left_index
		)
		if not left_neighbors or not right_neighbors:
			continue

		torsions.append(
			{
				"torsion_id": f"torsion_{left_index}_{right_index}",
				"central_bond": (left_index, right_index),
				"atom_indices": (
					left_neighbors[0],
					left_index,
					right_index,
					right_neighbors[0],
				),
			}
		)
	return torsions


def _calculation_molecule(
	molecule: Chem.Mol, conformer: Conformer
) -> Chem.Mol:
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


def _normalize_angle(angle: float) -> float:
	normalized = (angle + 180.0) % 360.0 - 180.0
	return 0.0 if normalized == 0.0 else normalized


def _state_index(angle: float, edges: tuple[float, ...]) -> int:
	# Bins are left-closed/right-open; an exact interior edge goes into the bin
	# beginning at that edge. Normalization maps +180 degrees to -180 degrees.
	return bisect_right(edges, angle) - 1


def _circular_statistics(angles: list[float]) -> dict[str, float | None]:
	if not angles:
		return {
			"circular_mean_deg": None,
			"circular_resultant_length": None,
			"circular_variance": None,
		}
	sine_mean = math.fsum(math.sin(math.radians(angle)) for angle in angles) / len(angles)
	cosine_mean = math.fsum(math.cos(math.radians(angle)) for angle in angles) / len(angles)
	resultant = min(1.0, math.hypot(sine_mean, cosine_mean))
	if resultant <= CIRCULAR_MEAN_RESULTANT_EPSILON:
		mean_angle = None
	else:
		mean_angle = _normalize_angle(math.degrees(math.atan2(sine_mean, cosine_mean)))
	return {
		"circular_mean_deg": mean_angle,
		"circular_resultant_length": resultant,
		"circular_variance": 1.0 - resultant,
	}


def _bin_metadata(edges: tuple[float, ...]) -> list[dict[str, Any]]:
	return [
		{
			"state_id": f"bin_{index:02d}",
			"lower_bound_deg": lower,
			"upper_bound_deg": upper,
			"lower_inclusive": True,
			"upper_inclusive": False,
		}
		for index, (lower, upper) in enumerate(zip(edges, edges[1:]))
	]


def analyze_torsions(
	ensemble: Ensemble,
	conformers: Sequence[Conformer] | None = None,
	*,
	state_bin_edges_deg: Sequence[Real] | None = None,
) -> dict[str, Any]:
	"""Analyze rotatable-bond angles, circular statistics, and state coverage.

	Eligible torsions use an acyclic single bond between two heavy atoms that
	each have at least one other heavy-atom neighbor. Bonds in rings and
	terminal heavy-atom bonds are excluded, as are amide/ester C-N and C-O
	bonds where the carbon is double-bonded to oxygen. Acyclic bonds linking
	ring systems are eligible. For each eligible bond, the lowest-index
	available heavy neighbor on each side is selected as the outer atom.
	Atom indices are the stable indices in ``ensemble.molecule`` and a
	torsion ID is ``torsion_<lower-central-index>_<higher-central-index>``.

	Angles are calculated by RDKit's ``GetDihedralDeg`` from copied ASE
	coordinates, normalized to [-180, 180). Circular means use the mean sine
	and cosine; the mean is ``None`` if the resultant length is at most
	``1e-12``. Circular variance is ``1 - resultant_length``.

	Parameters
	----------
	ensemble
		Source molecular ensemble.
	conformers
		Optional ordered subset of identical conformer objects owned by the
		ensemble. Output indices refer to the original ensemble.
	state_bin_edges_deg
		Strictly increasing bin boundaries spanning exactly -180 through 180.
		Bins are left-closed and right-open, so an angle exactly on an interior
		boundary belongs to the bin that starts at that boundary. Defaults to
		twelve equal-width 30-degree bins. These descriptive bins are not
		assumed to correspond to equal-energy states.

	Returns
	-------
	dict
		A mapping with ``conformers`` (one row per selected conformer and
		eligible torsion), ``torsions`` (one ensemble summary per torsion),
		and ``summary``. Missing angles remain ``None`` and have an explicit
		validity flag and error; valid zero-degree angles are retained as zero.
		Torsions with no valid observations have missing circular statistics
		and state coverage. Aggregate means use only torsions with at least one
		valid observation.

	Raises
	------
	ValueError
		If conformer selections are not uniquely identified by ensemble-owned
		objects, or if state-bin edges are invalid.
	"""
	selected, indices_by_id = _validate_conformers(ensemble, conformers)
	edges = _state_edges(state_bin_edges_deg)
	torsions = _eligible_torsions(ensemble.molecule)
	bins = _bin_metadata(edges)
	rows: list[dict[str, Any]] = []
	angles_by_torsion: dict[str, list[float]] = {
		torsion["torsion_id"]: [] for torsion in torsions
	}
	occupied_by_torsion: dict[str, set[int]] = {
		torsion["torsion_id"]: set() for torsion in torsions
	}
	n_invalid_observations = 0

	for conformer in selected:
		calculation_molecule: Chem.Mol | None = None
		coordinate_error: str | None = None
		try:
			calculation_molecule = _calculation_molecule(ensemble.molecule, conformer)
		except ValueError as exc:
			coordinate_error = str(exc)

		for torsion in torsions:
			row: dict[str, Any] = {
				"conformer_id": conformer.id,
				"conformer_index": indices_by_id[conformer.id],
				"torsion_id": torsion["torsion_id"],
				"central_bond": torsion["central_bond"],
				"atom_indices": torsion["atom_indices"],
				"angle_deg": None,
				"torsion_valid": False,
				"torsion_status": "invalid",
				"torsion_error": coordinate_error,
				"state_id": None,
			}
			if calculation_molecule is not None:
				try:
					angle = float(
						rdMolTransforms.GetDihedralDeg(
							calculation_molecule.GetConformer(0),
							*torsion["atom_indices"],
						)
					)
					if not math.isfinite(angle):
						raise ValueError("calculation returned a non-finite angle.")
					angle = _normalize_angle(angle)
					state_index = _state_index(angle, edges)
					state_id = bins[state_index]["state_id"]
					row.update(
						{
							"angle_deg": angle,
							"torsion_valid": True,
							"torsion_status": "valid",
							"torsion_error": None,
							"state_id": state_id,
						}
					)
					angles_by_torsion[torsion["torsion_id"]].append(angle)
					occupied_by_torsion[torsion["torsion_id"]].add(state_index)
				except _CALCULATION_ERRORS as exc:
					row["torsion_error"] = f"{type(exc).__name__}: {exc}"
			n_invalid_observations += not row["torsion_valid"]
			rows.append(row)

	torsion_rows: list[dict[str, Any]] = []
	for torsion in torsions:
		torsion_id = torsion["torsion_id"]
		angles = angles_by_torsion[torsion_id]
		occupied_states = sorted(occupied_by_torsion[torsion_id])
		circular = _circular_statistics(angles)
		torsion_rows.append(
			{
				"torsion_id": torsion_id,
				"central_bond": torsion["central_bond"],
				"atom_indices": torsion["atom_indices"],
				"n_valid_observations": len(angles),
				**circular,
				"n_occupied_states": len(occupied_states),
				"state_coverage": (
					len(occupied_states) / len(bins) if angles else None
				),
				"occupied_state_ids": [bins[index]["state_id"] for index in occupied_states],
				"state_bin_edges_deg": list(edges),
				"state_bins": bins,
			}
		)

	torsions_with_observations = [
		record for record in torsion_rows if record["n_valid_observations"] > 0
	]
	summary = {
		"n_conformers": len(selected),
		"n_eligible_torsions": len(torsions),
		"n_torsions_with_valid_observations": len(torsions_with_observations),
		"n_torsions_without_valid_observations": (
			len(torsions) - len(torsions_with_observations)
		),
		"n_torsions_with_undefined_circular_mean": sum(
			record["circular_mean_deg"] is None
			for record in torsions_with_observations
		),
		"n_valid_torsion_observations": len(rows) - n_invalid_observations,
		"n_invalid_torsion_observations": n_invalid_observations,
		"total_occupied_states": sum(
			record["n_occupied_states"] for record in torsions_with_observations
		),
		"mean_occupied_states": (
			math.fsum(
				record["n_occupied_states"] for record in torsions_with_observations
			)
			/ len(torsions_with_observations)
			if torsions_with_observations
			else None
		),
		"mean_state_coverage": (
			math.fsum(
				record["state_coverage"] for record in torsions_with_observations
			)
			/ len(torsions_with_observations)
			if torsions_with_observations
			else None
		),
		"angle_units": ANGLE_UNITS,
		"angle_interval": "[-180, 180)",
		"rotatable_bond_definition": (
			"acyclic single bond between nonterminal heavy atoms; exclude "
			"amide/ester C-N and C-O bonds"
		),
		"outer_atom_selection": "lowest-index heavy-atom neighbor on each side",
		"state_bin_edges_deg": list(edges),
		"state_bins": bins,
		"state_bin_boundary_convention": "left-closed, right-open",
		"circular_mean_undefined_resultant_threshold": (
			CIRCULAR_MEAN_RESULTANT_EPSILON
		),
	}
	return {"conformers": rows, "torsions": torsion_rows, "summary": summary}


__all__ = ["analyze_torsions"]
