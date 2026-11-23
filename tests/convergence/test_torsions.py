import csv
import io
import math
from dataclasses import replace

import numpy as np
import pytest
from rdkit.Chem import rdMolTransforms

from ensemblelab.analysis.torsions import torsion_angles as legacy_torsion_angles
from ensemblelab.convergence.torsions import analyze_torsions
from ensemblelab.generators import generate


def _butane_with_angles(angles):
	ensemble = generate("CCCC", n_confs=len(angles))
	for conformer, angle in zip(ensemble.conformers, angles):
		rdkit_conformer = ensemble.molecule.GetConformer(conformer.id)
		rdMolTransforms.SetDihedralDeg(rdkit_conformer, 0, 1, 2, 3, angle)
		conformer.atoms.positions[:] = rdkit_conformer.GetPositions()
	return ensemble


def _only_torsion(result):
	assert len(result["torsions"]) == 1
	return result["torsions"][0]


def test_per_conformer_angles_match_rdkit_reference_and_have_stable_atom_mapping():
	ensemble = _butane_with_angles([-60.0, 90.0])
	result = analyze_torsions(ensemble)

	assert set(result) == {"conformers", "torsions", "summary"}
	assert [row["conformer_id"] for row in result["conformers"]] == [
		conformer.id for conformer in ensemble.conformers
	]
	assert [row["conformer_index"] for row in result["conformers"]] == [0, 1]
	assert [row["torsion_id"] for row in result["conformers"]] == [
		"torsion_1_2",
		"torsion_1_2",
	]
	assert [row["central_bond"] for row in result["conformers"]] == [(1, 2)] * 2
	assert [row["atom_indices"] for row in result["conformers"]] == [
		(0, 1, 2, 3),
		(0, 1, 2, 3),
	]
	legacy_angles = legacy_torsion_angles(ensemble)
	assert legacy_angles.shape == (1, 2)
	for row, conformer in zip(result["conformers"], ensemble.conformers):
		expected = rdMolTransforms.GetDihedralDeg(
			ensemble.molecule.GetConformer(conformer.id), 0, 1, 2, 3
		)
		assert row["angle_deg"] == pytest.approx(expected)
		assert row["torsion_valid"] is True
		assert row["torsion_status"] == "valid"
		assert row["angle_deg"] == pytest.approx(legacy_angles[0, row["conformer_index"]])
	assert result["summary"]["angle_units"] == "degrees"
	assert result["summary"]["angle_interval"] == "[-180, 180)"


def test_wraparound_angles_have_circular_mean_at_boundary_not_zero():
	result = analyze_torsions(_butane_with_angles([179.0, -179.0]))
	torsion = _only_torsion(result)

	assert abs(torsion["circular_mean_deg"]) == pytest.approx(180.0)
	assert torsion["circular_resultant_length"] == pytest.approx(
		math.cos(math.radians(1.0))
	)
	assert torsion["circular_variance"] == pytest.approx(
		1.0 - math.cos(math.radians(1.0))
	)
	assert torsion["n_valid_observations"] == 2


def test_identical_angles_have_unit_resultant_and_single_occupied_state():
	torsion = _only_torsion(analyze_torsions(_butane_with_angles([60.0] * 3)))

	assert torsion["circular_mean_deg"] == pytest.approx(60.0)
	assert torsion["circular_resultant_length"] == pytest.approx(1.0)
	assert torsion["circular_variance"] == pytest.approx(0.0)
	assert torsion["n_occupied_states"] == 1
	assert torsion["state_coverage"] == pytest.approx(1 / 12)
	assert torsion["occupied_state_ids"] == ["bin_08"]


def test_distributed_angles_have_undefined_circular_mean():
	torsion = _only_torsion(
		analyze_torsions(_butane_with_angles([0.0, 90.0, -90.0, 180.0]))
	)

	assert torsion["circular_mean_deg"] is None
	assert torsion["circular_resultant_length"] < 1e-12
	assert torsion["circular_variance"] == pytest.approx(1.0)


def test_exact_bin_boundaries_are_assigned_to_the_bin_starting_at_boundary():
	result = analyze_torsions(
		_butane_with_angles([-180.0, -90.0, 0.0, 90.0, 180.0]),
		state_bin_edges_deg=[-180, -90, 0, 90, 180],
	)

	assert [row["angle_deg"] for row in result["conformers"]] == pytest.approx(
		[-180.0, -90.0, 0.0, 90.0, -180.0]
	)
	assert [row["state_id"] for row in result["conformers"]] == [
		"bin_00",
		"bin_01",
		"bin_02",
		"bin_03",
		"bin_00",
	]
	torsion = _only_torsion(result)
	assert torsion["n_occupied_states"] == 4
	assert torsion["state_coverage"] == 1.0
	assert result["summary"]["state_bin_boundary_convention"] == (
		"left-closed, right-open"
	)


@pytest.mark.parametrize(
	"edges",
	[
		[-170, 180],
		[-180, 170],
		[-180, 0, 0, 180],
		[-180],
		[-180, math.inf, 180],
		[-180, False, 180],
	],
)
def test_invalid_bin_edges_are_rejected(edges):
	with pytest.raises(ValueError, match="state_bin_edges_deg"):
		analyze_torsions(generate("CCCC", n_confs=1), state_bin_edges_deg=edges)


def test_default_edges_are_equal_width_thirty_degree_bins():
	summary = analyze_torsions(generate("CCCC", n_confs=1))["summary"]

	assert summary["state_bin_edges_deg"] == list(range(-180, 181, 30))
	assert len(summary["state_bins"]) == 12
	assert summary["state_bins"][0] == {
		"state_id": "bin_00",
		"lower_bound_deg": -180.0,
		"upper_bound_deg": -150.0,
		"lower_inclusive": True,
		"upper_inclusive": False,
	}


def test_empty_ensemble_has_no_observations_and_missing_torsion_statistics():
	ensemble = replace(generate("CCCC", n_confs=1), conformers=())
	result = analyze_torsions(ensemble)
	torsion = _only_torsion(result)

	assert result["conformers"] == []
	assert result["summary"]["n_conformers"] == 0
	assert result["summary"]["n_eligible_torsions"] == 1
	assert result["summary"]["n_valid_torsion_observations"] == 0
	assert result["summary"]["n_torsions_without_valid_observations"] == 1
	assert result["summary"]["total_occupied_states"] == 0
	assert result["summary"]["mean_occupied_states"] is None
	assert result["summary"]["mean_state_coverage"] is None
	assert torsion["circular_mean_deg"] is None
	assert torsion["circular_resultant_length"] is None
	assert torsion["state_coverage"] is None


def test_single_conformer_produces_one_valid_observation():
	result = analyze_torsions(_butane_with_angles([0.0]))
	torsion = _only_torsion(result)

	assert result["conformers"][0]["angle_deg"] == pytest.approx(0.0, abs=1e-10)
	assert result["conformers"][0]["torsion_valid"] is True
	assert torsion["n_valid_observations"] == 1
	assert torsion["circular_mean_deg"] == pytest.approx(0.0, abs=1e-10)
	assert torsion["circular_resultant_length"] == pytest.approx(1.0)
	assert result["summary"]["mean_state_coverage"] == pytest.approx(1 / 12)


@pytest.mark.parametrize("smiles", ["C", "CC", "CC(=O)NC", "CC(=O)O", "C1CCCCC1"])
def test_molecules_without_eligible_torsions_return_empty_torsion_records(smiles):
	result = analyze_torsions(generate(smiles, n_confs=1))

	assert result["conformers"] == []
	assert result["torsions"] == []
	assert result["summary"]["n_conformers"] == 1
	assert result["summary"]["n_eligible_torsions"] == 0
	assert result["summary"]["n_valid_torsion_observations"] == 0
	assert result["summary"]["total_occupied_states"] == 0
	assert result["summary"]["mean_occupied_states"] is None
	assert result["summary"]["mean_state_coverage"] is None


def test_invalid_coordinates_and_partial_observations_are_reported_not_zero_filled():
	ensemble = _butane_with_angles([30.0, 60.0])
	ensemble.conformers[0].atoms.positions[0, 0] = math.nan
	result = analyze_torsions(ensemble)

	invalid, valid = result["conformers"]
	assert invalid["angle_deg"] is None
	assert invalid["torsion_valid"] is False
	assert "non-finite" in invalid["torsion_error"]
	assert invalid["state_id"] is None
	assert valid["torsion_valid"] is True
	torsion = _only_torsion(result)
	assert torsion["n_valid_observations"] == 1
	assert torsion["state_coverage"] == pytest.approx(1 / 12)
	assert result["summary"]["n_valid_torsion_observations"] == 1
	assert result["summary"]["n_invalid_torsion_observations"] == 1


def test_inconsistent_atom_mapping_invalidates_only_that_conformer():
	ensemble = _butane_with_angles([30.0, 60.0])
	atomic_numbers = ensemble.conformers[0].atoms.get_atomic_numbers()
	atomic_numbers[0] = 8
	ensemble.conformers[0].atoms.set_atomic_numbers(atomic_numbers)

	result = analyze_torsions(ensemble)

	assert result["conformers"][0]["torsion_valid"] is False
	assert "do not match" in result["conformers"][0]["torsion_error"]
	assert result["conformers"][1]["torsion_valid"] is True
	assert result["summary"]["n_valid_torsion_observations"] == 1


def test_one_dihedral_failure_does_not_drop_other_observations(monkeypatch):
	ensemble = _butane_with_angles([30.0, 60.0])
	original = rdMolTransforms.GetDihedralDeg
	calls = 0

	def fail_first(conformer, *atom_indices):
		nonlocal calls
		calls += 1
		if calls == 1:
			raise RuntimeError("degenerate dihedral")
		return original(conformer, *atom_indices)

	monkeypatch.setattr(rdMolTransforms, "GetDihedralDeg", fail_first)
	result = analyze_torsions(ensemble)

	assert result["conformers"][0]["torsion_valid"] is False
	assert "degenerate dihedral" in result["conformers"][0]["torsion_error"]
	assert result["conformers"][1]["torsion_valid"] is True
	assert result["summary"]["n_valid_torsion_observations"] == 1
	assert result["summary"]["n_invalid_torsion_observations"] == 1


def test_inter_ring_linkage_is_eligible_but_ring_bonds_are_not():
	ensemble = generate("c1ccccc1-c2ccccc2", n_confs=1)
	result = analyze_torsions(ensemble)

	assert result["summary"]["n_eligible_torsions"] == 1
	torsion = _only_torsion(result)
	bond_begin, bond_end = torsion["central_bond"]
	bond = ensemble.molecule.GetBondBetweenAtoms(bond_begin, bond_end)
	assert bond is not None
	assert bond.IsInRing() is False


def test_analysis_preserves_ase_and_rdkit_coordinates():
	ensemble = _butane_with_angles([-60.0, 90.0])
	ase_positions = [
		conformer.atoms.positions.copy() for conformer in ensemble.conformers
	]
	rdkit_positions = {
		conformer.id: np.asarray(
			ensemble.molecule.GetConformer(conformer.id).GetPositions()
		).copy()
		for conformer in ensemble.conformers
	}

	analyze_torsions(ensemble)

	for conformer, positions in zip(ensemble.conformers, ase_positions):
		assert np.array_equal(conformer.atoms.positions, positions)
	for conformer_id, positions in rdkit_positions.items():
		assert np.array_equal(
			ensemble.molecule.GetConformer(conformer_id).GetPositions(), positions
		)


def test_subset_preserves_order_and_full_ensemble_indices():
	ensemble = _butane_with_angles([-60.0, 0.0, 90.0])
	result = analyze_torsions(
		ensemble, conformers=(ensemble.conformers[2], ensemble.conformers[0])
	)

	assert [row["conformer_id"] for row in result["conformers"]] == [
		ensemble.conformers[2].id,
		ensemble.conformers[0].id,
	]
	assert [row["conformer_index"] for row in result["conformers"]] == [2, 0]
	assert [row["angle_deg"] for row in result["conformers"]] == pytest.approx(
		[90.0, -60.0]
	)


def test_observation_records_are_csv_serializable_with_missing_angles_blank():
	ensemble = _butane_with_angles([30.0, 60.0])
	ensemble.conformers[0].atoms.positions[0, 0] = math.nan
	rows = analyze_torsions(ensemble)["conformers"]
	stream = io.StringIO()
	writer = csv.DictWriter(stream, fieldnames=rows[0].keys())
	writer.writeheader()
	writer.writerows(rows)
	serialized_rows = list(csv.DictReader(io.StringIO(stream.getvalue())))

	assert serialized_rows[0]["angle_deg"] == ""
	assert serialized_rows[0]["torsion_valid"] == "False"
	assert serialized_rows[0]["state_id"] == ""
	assert serialized_rows[1]["angle_deg"] != ""
	assert serialized_rows[1]["torsion_valid"] == "True"
