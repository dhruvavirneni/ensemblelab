import csv
import io
import math
from dataclasses import replace

import pytest

from ensemblelab.analysis.clusters import cluster_basins
from ensemblelab.analysis.populations import BOLTZMANN_CONSTANT_KCAL_MOL_K
from ensemblelab.convergence.populations import analyze_populations
from ensemblelab.generators import generate


def _ensemble_with_energies(energies):
	base = generate("C", n_confs=max(1, len(energies)))
	conformers = tuple(
		replace(
			conformer,
			energy=energy,
			energy_unit="kcal/mol" if energy is not None else None,
		)
		for conformer, energy in zip(base.conformers, energies)
	)
	return replace(base, conformers=conformers)


def test_known_energy_difference_produces_expected_normalized_populations():
	temperature = 300.0
	energy_difference = (
		BOLTZMANN_CONSTANT_KCAL_MOL_K * temperature * math.log(2.0)
	)
	result = analyze_populations(
		_ensemble_with_energies([5.0, 5.0 + energy_difference]), temperature
	)

	assert [row["boltzmann_weight"] for row in result["conformers"]] == pytest.approx([1.0, 0.5])
	assert [row["population"] for row in result["conformers"]] == pytest.approx([2 / 3, 1 / 3])
	assert result["summary"]["population_sum"] == pytest.approx(1.0)
	assert result["summary"]["max_population"] == pytest.approx(2 / 3)
	assert result["summary"]["population_entropy"] == pytest.approx(
		-(2 / 3) * math.log(2 / 3) - (1 / 3) * math.log(1 / 3)
	)
	assert result["summary"]["effective_ensemble_size"] == pytest.approx(
		math.exp(result["summary"]["population_entropy"])
	)
	assert result["summary"]["temperature_K"] == temperature
	assert result["summary"]["energy_units"] == "kcal/mol"


def test_equal_energy_conformers_have_equal_populations_and_maximum_effective_size():
	result = analyze_populations(_ensemble_with_energies([2.0, 2.0, 2.0]))

	assert [row["population"] for row in result["conformers"]] == pytest.approx([1 / 3] * 3)
	assert result["summary"]["effective_ensemble_size"] == pytest.approx(3.0)
	assert result["summary"]["population_entropy"] == pytest.approx(math.log(3.0))


def test_single_valid_conformer_has_population_one():
	result = analyze_populations(_ensemble_with_energies([3.25]))

	assert result["conformers"][0]["population"] == 1.0
	assert result["conformers"][0]["boltzmann_weight"] == 1.0
	assert result["summary"]["effective_ensemble_size"] == 1.0
	assert result["summary"]["population_entropy"] == 0.0


def test_empty_and_no_valid_energy_ensembles_have_explicit_missing_statistics():
	for energies in ([], [None, math.nan, math.inf]):
		result = analyze_populations(_ensemble_with_energies(energies))
		summary = result["summary"]
		assert summary["n_valid_energies"] == 0
		assert summary["n_invalid_energies"] == len(energies)
		assert summary["population_sum"] is None
		assert summary["max_population"] is None
		assert summary["effective_ensemble_size"] is None
		assert summary["population_entropy"] is None
		assert all(row["population"] is None for row in result["conformers"])


def test_invalid_energies_are_excluded_and_marked_missing():
	result = analyze_populations(_ensemble_with_energies([1.0, None, math.nan, 2.0]))

	assert [row["population_valid"] for row in result["conformers"]] == [True, False, False, True]
	assert result["conformers"][1]["boltzmann_weight"] is None
	assert result["conformers"][1]["population"] is None
	assert result["summary"]["n_invalid_energies"] == 2
	assert result["summary"]["population_sum"] == pytest.approx(1.0)


def test_extremely_high_relative_energy_underflows_without_breaking_normalization():
	result = analyze_populations(_ensemble_with_energies([0.0, 1_000_000.0]))

	assert result["conformers"][1]["population_valid"] is True
	assert result["conformers"][1]["boltzmann_weight"] == 0.0
	assert result["conformers"][1]["population"] == 0.0
	assert result["summary"]["population_sum"] == 1.0
	assert result["summary"]["effective_ensemble_size"] == 1.0


@pytest.mark.parametrize("temperature", [0.0, -1.0, math.nan, math.inf, -math.inf, True])
def test_invalid_temperature_is_rejected(temperature):
	with pytest.raises(ValueError, match="positive finite value in kelvin"):
		analyze_populations(_ensemble_with_energies([1.0]), temperature)


def test_default_temperature_is_recorded_and_input_is_not_mutated():
	ensemble = _ensemble_with_energies([2.0, 1.0])
	original_energies = [conformer.energy for conformer in ensemble.conformers]
	original_positions = [conformer.atoms.positions.copy() for conformer in ensemble.conformers]
	result = analyze_populations(ensemble)

	assert result["summary"]["temperature_K"] == 298.15
	assert [row["conformer_id"] for row in result["conformers"]] == list(ensemble.conformer_ids)
	assert [row["conformer_index"] for row in result["conformers"]] == [0, 1]
	assert [conformer.energy for conformer in ensemble.conformers] == original_energies
	for conformer, positions in zip(ensemble.conformers, original_positions):
		assert (conformer.atoms.positions == positions).all()


def test_population_records_are_csv_serializable_and_preserve_missing_values():
	result = analyze_populations(_ensemble_with_energies([0.0, None]))
	stream = io.StringIO()
	writer = csv.DictWriter(stream, fieldnames=result["conformers"][0].keys())
	writer.writeheader()
	writer.writerows(result["conformers"])
	rows = list(csv.DictReader(io.StringIO(stream.getvalue())))

	assert rows[0]["population_valid"] == "True"
	assert rows[0]["population"] == "1.0"
	assert rows[1]["boltzmann_weight"] == ""
	assert rows[1]["population"] == ""
	assert rows[1]["population_valid"] == "False"


def test_existing_basin_population_policy_is_preserved():
	complete_ensemble = _ensemble_with_energies([0.0, 1.0])
	complete_basins = cluster_basins(complete_ensemble)
	assert sum(basin.population for basin in complete_basins) == pytest.approx(1.0)

	partial_ensemble = _ensemble_with_energies([0.0, None])
	partial_basins = cluster_basins(partial_ensemble)
	assert all(basin.population is None for basin in partial_basins)