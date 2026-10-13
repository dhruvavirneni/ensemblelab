import csv
import io
import math
from dataclasses import replace

import pytest

from ensemblelab.analysis.clusters import Basin
from ensemblelab.analysis.entropy import entropy_statistics
from ensemblelab.analysis.populations import BOLTZMANN_CONSTANT_KCAL_MOL_K
from ensemblelab.convergence.entropy import analyze_entropy
from ensemblelab.generators import generate


def _ensemble_with_probabilities(probabilities, temperature=298.15):
	base = generate("C", n_confs=max(1, len(probabilities)))
	energies = [
		-BOLTZMANN_CONSTANT_KCAL_MOL_K * temperature * math.log(probability)
		for probability in probabilities
	]
	conformers = tuple(
		replace(conformer, energy=energy, energy_unit="kcal/mol")
		for conformer, energy in zip(base.conformers, energies)
	)
	return replace(base, conformers=conformers)


def _basin(identifier, conformers):
	return Basin(
		id=identifier,
		representative=conformers[0],
		members=tuple(conformers),
		population=None,
		min_energy=None,
	)


def test_entropy_statistics_known_distributions_and_validation():
	assert entropy_statistics([1.0])["entropy"] == 0.0
	assert entropy_statistics([1.0])["effective_size"] == 1.0
	assert entropy_statistics([0.5, 0.5])["entropy"] == pytest.approx(math.log(2.0))
	assert entropy_statistics([0.5, 0.5])["effective_size"] == pytest.approx(2.0)
	assert entropy_statistics([0.75, 0.25])["entropy"] == pytest.approx(
		-0.75 * math.log(0.75) - 0.25 * math.log(0.25)
	)
	assert entropy_statistics([1.0, 0.0])["entropy"] == 0.0
	assert entropy_statistics([])["entropy"] is None
	with pytest.raises(ValueError, match="normalized"):
		entropy_statistics([0.2, 0.2])
	with pytest.raises(ValueError, match="finite, non-negative"):
		entropy_statistics([math.nan])


def test_one_state_and_two_equal_conformers_have_expected_entropy():
	single = _ensemble_with_probabilities([1.0])
	single_result = analyze_entropy(
		single,
		basins=[_basin(4, single.conformers)],
	)
	assert single_result["summary"]["conformer_entropy"] == 0.0
	assert single_result["summary"]["conformer_effective_size"] == 1.0
	assert single_result["summary"]["cluster_entropy"] == 0.0
	assert single_result["summary"]["cluster_effective_size"] == 1.0

	equal = _ensemble_with_probabilities([0.5, 0.5])
	equal_result = analyze_entropy(
		equal,
		basins=[_basin(8, (equal.conformers[0],)), _basin(9, (equal.conformers[1],))],
	)
	assert equal_result["summary"]["conformer_entropy"] == pytest.approx(math.log(2))
	assert equal_result["summary"]["conformer_effective_size"] == pytest.approx(2.0)
	assert equal_result["summary"]["cluster_entropy"] == pytest.approx(math.log(2))


def test_conformers_in_one_cluster_aggregate_to_one_structural_state():
	ensemble = _ensemble_with_probabilities([0.5, 0.25, 0.25])
	result = analyze_entropy(ensemble, basins=[_basin(12, ensemble.conformers)])

	assert result["summary"]["conformer_effective_size"] == pytest.approx(2.8284271247461903)
	assert result["clusters"][0]["cluster_population"] == pytest.approx(1.0)
	assert result["clusters"][0]["cluster_entropy_contribution"] == 0.0
	assert result["summary"]["cluster_entropy"] == 0.0
	assert result["summary"]["cluster_effective_size"] == 1.0
	assert [record["cluster_id"] for record in result["conformers"]] == [12, 12, 12]


def test_unequal_populations_are_summed_over_cluster_members():
	ensemble = _ensemble_with_probabilities([0.5, 0.25, 0.25])
	basins = [
		_basin(20, ensemble.conformers[:2]),
		_basin(21, ensemble.conformers[2:]),
	]
	result = analyze_entropy(ensemble, basins=basins)
	cluster_entropy = -0.75 * math.log(0.75) - 0.25 * math.log(0.25)

	assert [row["cluster_population"] for row in result["clusters"]] == pytest.approx([0.75, 0.25])
	assert [row["cluster_entropy_contribution"] for row in result["clusters"]] == pytest.approx(
		[0.75 * -math.log(0.75), 0.25 * -math.log(0.25)]
	)
	assert result["summary"]["cluster_entropy"] == pytest.approx(cluster_entropy)
	assert result["summary"]["cluster_effective_size"] == pytest.approx(math.exp(cluster_entropy))
	assert result["summary"]["cluster_effective_size"] <= result["summary"]["conformer_effective_size"]


def test_invalid_populations_are_missing_not_zero_and_excluded():
	base = generate("C", n_confs=3)
	conformers = (
		replace(base.conformers[0], energy=0.0, energy_unit="kcal/mol"),
		replace(base.conformers[1], energy=None, energy_unit=None),
		replace(base.conformers[2], energy=1.0, energy_unit="kcal/mol"),
	)
	ensemble = replace(base, conformers=conformers)
	result = analyze_entropy(
		ensemble,
		basins=[_basin(30, (conformers[0], conformers[1])), _basin(31, (conformers[2],))],
	)

	assert result["summary"]["n_valid_conformer_populations"] == 2
	assert result["summary"]["n_invalid_populations"] == 1
	assert result["conformers"][1]["population"] is None
	assert result["conformers"][1]["cluster_id"] == 30
	assert result["clusters"][0]["cluster_population"] == pytest.approx(result["conformers"][0]["population"])
	assert result["clusters"][1]["cluster_population"] == pytest.approx(result["conformers"][2]["population"])
	assert result["summary"]["cluster_entropy"] is not None


def test_no_valid_populations_return_missing_conformer_and_cluster_entropy():
	base = generate("C", n_confs=2)
	conformers = tuple(
		replace(conformer, energy=None, energy_unit=None)
		for conformer in base.conformers
	)
	ensemble = replace(base, conformers=conformers)
	result = analyze_entropy(
		ensemble,
		basins=[_basin(35, (conformers[0],)), _basin(36, (conformers[1],))],
	)

	assert result["summary"]["n_valid_conformer_populations"] == 0
	assert result["summary"]["n_valid_cluster_populations"] == 0
	assert result["summary"]["conformer_entropy"] is None
	assert result["summary"]["cluster_entropy"] is None
	assert all(row["cluster_population"] is None for row in result["clusters"])


def test_underflowed_zero_population_cluster_has_zero_entropy_contribution():
	base = generate("C", n_confs=2)
	conformers = (
		replace(base.conformers[0], energy=0.0, energy_unit="kcal/mol"),
		replace(base.conformers[1], energy=1_000_000.0, energy_unit="kcal/mol"),
	)
	ensemble = replace(base, conformers=conformers)
	result = analyze_entropy(
		ensemble,
		basins=[_basin(40, (conformers[0],)), _basin(41, (conformers[1],))],
	)

	assert result["clusters"][1]["cluster_population"] == 0.0
	assert result["clusters"][1]["cluster_entropy_contribution"] == 0.0
	assert result["summary"]["cluster_entropy"] == 0.0


def test_incomplete_or_overlapping_assignments_make_cluster_entropy_unavailable():
	ensemble = _ensemble_with_probabilities([0.5, 0.3, 0.2])
	partial = analyze_entropy(ensemble, basins=[_basin(50, ensemble.conformers[:2])])
	assert partial["summary"]["cluster_assignment_complete"] is False
	assert partial["summary"]["cluster_entropy_available"] is False
	assert partial["summary"]["cluster_entropy"] is None
	assert partial["summary"]["n_unassigned_conformers"] == 1
	assert partial["conformers"][2]["cluster_id"] is None

	overlap = analyze_entropy(
		ensemble,
		basins=[
			_basin(51, ensemble.conformers[:2]),
			_basin(52, ensemble.conformers[1:]),
		],
	)
	assert overlap["summary"]["cluster_assignment_complete"] is False
	assert overlap["summary"]["cluster_entropy"] is None


def test_empty_ensemble_has_explicit_missing_entropy_and_no_clusters():
	base = generate("C", n_confs=1)
	empty = replace(base, conformers=())
	result = analyze_entropy(empty)

	assert result["conformers"] == []
	assert result["clusters"] == []
	assert result["summary"]["n_conformers"] == 0
	assert result["summary"]["n_clusters"] == 0
	assert result["summary"]["conformer_entropy"] is None
	assert result["summary"]["cluster_entropy"] is None
	assert result["summary"]["cluster_assignment_complete"] is True


def test_default_clustering_records_method_cutoff_and_complete_assignments():
	ensemble = _ensemble_with_probabilities([0.5, 0.5])
	result = analyze_entropy(ensemble, rmsd_threshold=0.9)

	assert result["summary"]["clustering_method"] == "greedy_best_fit_rmsd"
	assert result["summary"]["clustering_distance"] == "RDKit best-fit RMSD"
	assert result["summary"]["rmsd_threshold_A"] == 0.9
	assert result["summary"]["cluster_assignment_complete"] is True
	assert sum(row["cluster_population"] for row in result["clusters"]) == pytest.approx(1.0)


def test_output_is_csv_compatible_and_analysis_does_not_mutate_inputs():
	ensemble = _ensemble_with_probabilities([0.7, 0.3])
	basins = [_basin(60, ensemble.conformers)]
	energies_before = [conformer.energy for conformer in ensemble.conformers]
	positions_before = [conformer.atoms.positions.copy() for conformer in ensemble.conformers]
	basin_populations_before = [basin.population for basin in basins]
	result = analyze_entropy(ensemble, basins=basins)

	assert [row["conformer_index"] for row in result["conformers"]] == [0, 1]
	assert [row["conformer_id"] for row in result["conformers"]] == list(ensemble.conformer_ids)
	assert [conformer.energy for conformer in ensemble.conformers] == energies_before
	assert [basin.population for basin in basins] == basin_populations_before
	for conformer, positions in zip(ensemble.conformers, positions_before):
		assert (conformer.atoms.positions == positions).all()
	stream = io.StringIO()
	writer = csv.DictWriter(stream, fieldnames=result["clusters"][0].keys())
	writer.writeheader()
	writer.writerows(result["clusters"])
	assert "cluster_entropy_contribution" in stream.getvalue()