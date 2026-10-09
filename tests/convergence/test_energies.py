import csv
import io
import math
from dataclasses import replace

import pytest

from ensemblelab.convergence.energies import analyze_energies
from ensemblelab.generators import Ensemble, generate


def _ensemble_with_energies(energies, units=None):
	base = generate("C", n_confs=max(1, len(energies)))
	if units is None:
		units = ["kcal/mol" if energy is not None else None for energy in energies]
	conformers = tuple(
		replace(conformer, energy=energy, energy_unit=unit)
		for conformer, energy, unit in zip(base.conformers, energies, units)
	)
	return replace(base, conformers=conformers)


def test_analyze_different_energies_and_relative_energy_summary():
	ensemble = _ensemble_with_energies([1.0, 2.0, 3.0])
	result = analyze_energies(ensemble)

	assert set(result) == {"conformers", "summary"}
	assert set(result["conformers"][0]) == {
		"conformer_id",
		"conformer_index",
		"energy",
		"relative_energy",
		"energy_valid",
	}
	assert [row["relative_energy"] for row in result["conformers"]] == [0.0, 1.0, 2.0]
	assert [row["energy_valid"] for row in result["conformers"]] == [True, True, True]
	assert result["summary"] == {
		"n_conformers": 3,
		"n_valid_energies": 3,
		"n_invalid_energies": 0,
		"energy_min": 1.0,
		"energy_max": 3.0,
		"energy_range": 2.0,
		"energy_mean": 2.0,
		"energy_median": 2.0,
		"energy_std": math.sqrt(2 / 3),
		"relative_energy_mean": 1.0,
		"relative_energy_median": 1.0,
		"energy_units": "kcal/mol",
	}


def test_minimum_energy_conformer_has_zero_relative_energy():
	ensemble = _ensemble_with_energies([5.0, -2.0, 1.0])
	rows = analyze_energies(ensemble)["conformers"]

	assert rows[1]["relative_energy"] == 0.0
	assert rows[1]["energy"] == -2.0


def test_equal_energy_conformers_have_zero_relative_energy_and_range():
	result = analyze_energies(_ensemble_with_energies([4.5, 4.5]))

	assert [row["relative_energy"] for row in result["conformers"]] == [0.0, 0.0]
	assert result["summary"]["energy_range"] == 0.0
	assert result["summary"]["energy_std"] == 0.0


def test_single_conformer_has_zero_relative_energy_and_population_std():
	result = analyze_energies(_ensemble_with_energies([2.75]))

	assert result["conformers"][0]["relative_energy"] == 0.0
	assert result["summary"]["energy_std"] == 0.0
	assert result["summary"]["n_valid_energies"] == 1


def test_empty_ensemble_returns_empty_records_and_unavailable_statistics():
	ensemble = _ensemble_with_energies([])
	result = analyze_energies(ensemble)

	assert result["conformers"] == []
	assert result["summary"]["n_conformers"] == 0
	assert result["summary"]["n_valid_energies"] == 0
	assert result["summary"]["n_invalid_energies"] == 0
	assert result["summary"]["energy_min"] is None
	assert result["summary"]["energy_mean"] is None
	assert result["summary"]["energy_units"] == "kcal/mol"


def test_missing_and_non_finite_energies_are_invalid_and_counted():
	result = analyze_energies(_ensemble_with_energies([None, math.nan, math.inf]))

	assert [row["energy_valid"] for row in result["conformers"]] == [False] * 3
	assert all(row["relative_energy"] is None for row in result["conformers"])
	assert result["summary"]["n_valid_energies"] == 0
	assert result["summary"]["n_invalid_energies"] == 3
	assert result["summary"]["energy_min"] is None
	assert result["summary"]["relative_energy_median"] is None


def test_nonnumeric_energy_is_invalid_and_malformed_unit_raises_value_error():
	result = analyze_energies(_ensemble_with_energies(["unavailable", 1.0]))

	assert result["conformers"][0]["energy"] == "unavailable"
	assert result["conformers"][0]["energy_valid"] is False
	assert result["summary"]["n_invalid_energies"] == 1
	with pytest.raises(ValueError, match="expected 'kcal/mol'"):
		analyze_energies(_ensemble_with_energies([1.0], units=[2]))


def test_finite_energy_requires_explicit_supported_unit():
	with pytest.raises(ValueError, match="expected 'kcal/mol'"):
		analyze_energies(_ensemble_with_energies([1.0], units=[None]))
	with pytest.raises(ValueError, match="expected 'kcal/mol'"):
		analyze_energies(_ensemble_with_energies([1.0], units=["eV"]))


def test_analysis_preserves_order_identity_and_does_not_mutate_input():
	ensemble = _ensemble_with_energies([3.0, None, 1.0])
	original_energies = [conformer.energy for conformer in ensemble.conformers]
	original_positions = [conformer.atoms.positions.copy() for conformer in ensemble.conformers]
	result = analyze_energies(ensemble)

	assert [row["conformer_id"] for row in result["conformers"]] == list(ensemble.conformer_ids)
	assert [row["conformer_index"] for row in result["conformers"]] == [0, 1, 2]
	assert [row["energy"] for row in result["conformers"]] == original_energies
	assert [conformer.energy for conformer in ensemble.conformers] == original_energies
	for conformer, positions in zip(ensemble.conformers, original_positions):
		assert (conformer.atoms.positions == positions).all()


def test_subset_analysis_keeps_indices_from_source_ensemble():
	ensemble = _ensemble_with_energies([3.0, 1.0, 2.0])
	result = analyze_energies(ensemble, conformers=ensemble.conformers[:2])

	assert [row["conformer_index"] for row in result["conformers"]] == [0, 1]
	assert result["summary"]["n_conformers"] == 2
	assert result["summary"]["energy_min"] == 1.0


def test_output_records_are_csv_serializable_without_zero_filling_missing_values():
	result = analyze_energies(_ensemble_with_energies([0.0, None]))
	stream = io.StringIO()
	writer = csv.DictWriter(stream, fieldnames=result["conformers"][0].keys())
	writer.writeheader()
	writer.writerows(result["conformers"])
	rows = list(csv.DictReader(io.StringIO(stream.getvalue())))

	assert [row["conformer_id"] for row in rows] == ["0", "1"]
	assert rows[0]["relative_energy"] == "0.0"
	assert rows[0]["energy_valid"] == "True"
	assert rows[1]["energy"] == ""
	assert rows[1]["relative_energy"] == ""
	assert rows[1]["energy_valid"] == "False"