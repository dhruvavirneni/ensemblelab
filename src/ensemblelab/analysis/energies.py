"""Energy metrics for conformational ensembles."""

from __future__ import annotations

import math
from collections.abc import Sequence
from numbers import Real
from statistics import fmean, median, pstdev
from typing import Any

from ensemblelab.generators import Conformer, Ensemble

ENERGY_UNITS = "kcal/mol"


def analyze_energies(
	ensemble: Ensemble,
	conformers: Sequence[Conformer] | None = None,
) -> dict[str, Any]:
	"""Return per-conformer energy records and summary statistics.

	Parameters
	----------
	ensemble
		Source molecular ensemble. Its conformer ordering is preserved by
		default.
	conformers
		Optional ordered subset of conformers from ``ensemble``. This supports
		prefix or subset analyses without rebuilding conformer data. Each item
		must be the identical object held by the ensemble, and IDs must be
		unique. ``conformer_index`` always refers to the position in the full
		ensemble.

	Returns
	-------
	dict
		A stable mapping with ``conformers`` (one record per selected conformer)
		and ``summary`` (one ensemble summary). Records can be passed to a
		DataFrame/CSV writer. Invalid records retain their source ``energy``;
		their ``relative_energy`` is ``None``. Summary statistics are ``None``
		when there are no valid energies. ``energy_std`` is the population
		standard deviation (denominator *n*); a single valid energy therefore
		has standard deviation zero.

	Raises
	------
	ValueError
		If a selected conformer is not the exact object in the ensemble, IDs
		are duplicated, or a finite numeric energy has a missing or unsupported
		unit. Energies must be in kcal/mol; no conversion is performed.
	"""
	selected = list(ensemble.conformers if conformers is None else conformers)
	indices_by_id = {
		conformer.id: index
		for index, conformer in enumerate(ensemble.conformers)
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

	valid_values: list[float] = []
	valid_positions: list[int] = []
	rows: list[dict[str, Any]] = []
	for conformer in selected:
		energy = conformer.energy
		energy_valid = (
			isinstance(energy, Real)
			and not isinstance(energy, bool)
			and math.isfinite(float(energy))
		)
		if energy_valid:
			unit = conformer.energy_unit
			if not isinstance(unit, str) or unit.lower().replace(" ", "") != ENERGY_UNITS:
				raise ValueError(
					f"conformer {conformer.id} has a finite energy with unsupported "
					f"unit {unit!r}; expected {ENERGY_UNITS!r}."
				)
			valid_positions.append(len(rows))
			valid_values.append(float(energy))
		rows.append(
			{
				"conformer_id": conformer.id,
				"conformer_index": indices_by_id[conformer.id],
				"energy": energy,
				"relative_energy": None,
				"energy_valid": energy_valid,
			}
		)

	if valid_values:
		minimum = min(valid_values)
		for position, energy in zip(valid_positions, valid_values):
			rows[position]["relative_energy"] = energy - minimum
		relative_values = [rows[position]["relative_energy"] for position in valid_positions]
		energy_min: float | None = minimum
		energy_max: float | None = max(valid_values)
		energy_range: float | None = energy_max - energy_min
		energy_mean: float | None = fmean(valid_values)
		energy_median: float | None = median(valid_values)
		energy_std: float | None = pstdev(valid_values)
		relative_energy_mean: float | None = fmean(relative_values)
		relative_energy_median: float | None = median(relative_values)
	else:
		energy_min = energy_max = energy_range = None
		energy_mean = energy_median = energy_std = None
		relative_energy_mean = relative_energy_median = None

	return {
		"conformers": rows,
		"summary": {
			"n_conformers": len(selected),
			"n_valid_energies": len(valid_values),
			"n_invalid_energies": len(selected) - len(valid_values),
			"energy_min": energy_min,
			"energy_max": energy_max,
			"energy_range": energy_range,
			"energy_mean": energy_mean,
			"energy_median": energy_median,
			"energy_std": energy_std,
			"relative_energy_mean": relative_energy_mean,
			"relative_energy_median": relative_energy_median,
			"energy_units": ENERGY_UNITS,
		},
	}