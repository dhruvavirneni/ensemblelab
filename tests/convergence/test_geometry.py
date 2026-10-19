import csv
import io
import math
from dataclasses import replace

import numpy as np
import pytest
from rdkit import Chem
from rdkit.Chem import Descriptors3D

try:
	from rdkit.Chem import rdFreeSASA
except ImportError:
	rdFreeSASA = None

import ensemblelab.convergence.geometry as geometry
from ensemblelab.convergence.geometry import analyze_geometry
from ensemblelab.generators import generate


def _make_calculation_molecule(ensemble, conformer):
	molecule = Chem.Mol(ensemble.molecule)
	molecule.RemoveAllConformers()
	rdkit_conformer = Chem.Conformer(molecule.GetNumAtoms())
	for atom_index, position in enumerate(conformer.atoms.get_positions()):
		rdkit_conformer.SetAtomPosition(atom_index, tuple(position))
	molecule.AddConformer(rdkit_conformer, assignId=True)
	return molecule


def test_metrics_match_rdkit_reference_implementations_and_report_metadata():
	ensemble = generate("CCO", n_confs=1)
	conformer = ensemble.conformers[0]
	molecule = _make_calculation_molecule(ensemble, conformer)
	expected_radius = Descriptors3D.RadiusOfGyration(
		molecule, confId=0, useAtomicMasses=True, force=True
	)
	if rdFreeSASA is not None:
		radii = [
			Chem.GetPeriodicTable().GetRvdw(atom.GetAtomicNum())
			for atom in molecule.GetAtoms()
		]
		options = rdFreeSASA.SASAOpts(
			rdFreeSASA.SASAAlgorithm.LeeRichards,
			rdFreeSASA.SASAClassifier.Protor,
			1.4,
		)
		expected_sasa = rdFreeSASA.CalcSASA(
			molecule, radii, confIdx=0, opts=options
		)

	result = analyze_geometry(ensemble)
	record = result["conformers"][0]
	summary = result["summary"]

	assert set(result) == {"conformers", "summary"}
	assert record["conformer_id"] == conformer.id
	assert record["conformer_index"] == 0
	assert record["radius_of_gyration_A"] == pytest.approx(expected_radius)
	assert record["radius_of_gyration_valid"] is True
	assert record["radius_of_gyration_status"] == "valid"
	assert summary["n_conformers"] == 1
	assert summary["n_valid_radius_of_gyration"] == 1
	assert summary["radius_of_gyration_units"] == "angstrom"
	assert summary["radius_of_gyration_mass_weighted"] is True
	assert summary["sasa_units"] == "angstrom^2"
	if rdFreeSASA is None:
		assert record["sasa_A2"] is None
		assert record["sasa_status"] == "unavailable"
		assert summary["n_valid_sasa"] == 0
	else:
		assert record["sasa_A2"] == pytest.approx(expected_sasa)
		assert record["sasa_valid"] is True
		assert record["sasa_status"] == "valid"
		assert summary["n_valid_sasa"] == 1
		assert summary["sasa_probe_radius_A"] == 1.4
		assert summary["sasa_atomic_radius_model"] == (
			"RDKit periodic-table van der Waals radii (GetRvdw)"
		)


def test_multiple_geometries_produce_distinct_values_and_preserve_subset_order():
	ensemble = generate("CC", n_confs=3)
	ensemble.conformers[0].atoms.positions[:] *= 1.5
	result = analyze_geometry(
		ensemble, conformers=(ensemble.conformers[2], ensemble.conformers[0])
	)

	assert [row["conformer_id"] for row in result["conformers"]] == [
		ensemble.conformers[2].id,
		ensemble.conformers[0].id,
	]
	assert [row["conformer_index"] for row in result["conformers"]] == [2, 0]
	assert result["conformers"][0]["radius_of_gyration_A"] != pytest.approx(
		result["conformers"][1]["radius_of_gyration_A"]
	)
	assert result["summary"]["n_conformers"] == 2


def test_single_conformer_has_zero_population_standard_deviation():
	summary = analyze_geometry(generate("C", n_confs=1))["summary"]

	assert summary["radius_of_gyration_std_A"] == 0.0
	assert summary["radius_of_gyration_range_A"] == 0.0
	if rdFreeSASA is not None:
		assert summary["sasa_std_A2"] == 0.0
		assert summary["sasa_range_A2"] == 0.0
	else:
		assert summary["sasa_std_A2"] is None


def test_empty_ensemble_has_empty_records_and_missing_statistics():
	ensemble = replace(generate("C", n_confs=1), conformers=())
	result = analyze_geometry(ensemble)

	assert result["conformers"] == []
	assert result["summary"]["n_conformers"] == 0
	assert result["summary"]["n_valid_radius_of_gyration"] == 0
	assert result["summary"]["n_invalid_sasa"] == 0
	assert result["summary"]["radius_of_gyration_mean_A"] is None
	assert result["summary"]["sasa_mean_A2"] is None


def test_non_finite_coordinates_are_invalid_without_losing_other_conformers():
	ensemble = generate("CC", n_confs=2)
	ensemble.conformers[0].atoms.positions[0] = [math.nan, 0.0, 0.0]
	result = analyze_geometry(ensemble)

	invalid, valid = result["conformers"]
	assert invalid["radius_of_gyration_A"] is None
	assert invalid["radius_of_gyration_valid"] is False
	assert "non-finite" in invalid["radius_of_gyration_error"]
	assert invalid["sasa_A2"] is None
	assert invalid["sasa_valid"] is False
	assert "non-finite" in invalid["sasa_error"]
	assert valid["radius_of_gyration_valid"] is True
	assert result["summary"]["n_valid_radius_of_gyration"] == 1
	assert result["summary"]["n_invalid_radius_of_gyration"] == 1
	assert result["summary"]["n_invalid_sasa"] == 1
	if rdFreeSASA is not None:
		assert valid["sasa_valid"] is True
		assert result["summary"]["n_valid_sasa"] == 1
	else:
		assert valid["sasa_status"] == "unavailable"
		assert result["summary"]["n_valid_sasa"] == 0


def test_atom_order_mismatch_marks_available_metrics_invalid():
	ensemble = generate("CC", n_confs=1)
	atomic_numbers = ensemble.conformers[0].atoms.get_atomic_numbers()
	atomic_numbers[0] = 8
	ensemble.conformers[0].atoms.set_atomic_numbers(atomic_numbers)
	record = analyze_geometry(ensemble)["conformers"][0]

	assert record["radius_of_gyration_valid"] is False
	assert "do not match" in record["radius_of_gyration_error"]
	if rdFreeSASA is not None:
		assert record["sasa_valid"] is False
		assert "do not match" in record["sasa_error"]
	else:
		assert record["sasa_status"] == "unavailable"


def test_metric_failures_are_isolated_and_non_finite_results_are_rejected(monkeypatch):
	ensemble = generate("C", n_confs=1)

	def fail_radius(*args, **kwargs):
		raise RuntimeError("descriptor failed")

	monkeypatch.setattr(geometry.Descriptors3D, "RadiusOfGyration", fail_radius)
	result = analyze_geometry(ensemble)
	record = result["conformers"][0]
	assert record["radius_of_gyration_valid"] is False
	assert "descriptor failed" in record["radius_of_gyration_error"]
	if rdFreeSASA is not None:
		assert record["sasa_valid"] is True
	else:
		assert record["sasa_status"] == "unavailable"

	monkeypatch.undo()
	monkeypatch.setattr(
		geometry.Descriptors3D,
		"RadiusOfGyration",
		lambda *args, **kwargs: math.inf,
	)
	record = analyze_geometry(ensemble)["conformers"][0]
	assert record["radius_of_gyration_A"] is None
	assert record["radius_of_gyration_valid"] is False
	assert "non-finite" in record["radius_of_gyration_error"]


def test_sasa_failure_does_not_invalidate_radius_of_gyration(monkeypatch):
	if rdFreeSASA is None:
		pytest.skip("RDKit FreeSASA is not available.")
	ensemble = generate("C", n_confs=1)

	def fail_sasa(*args, **kwargs):
		raise RuntimeError("SASA failed")

	monkeypatch.setattr(geometry.rdFreeSASA, "CalcSASA", fail_sasa)
	record = analyze_geometry(ensemble)["conformers"][0]

	assert record["radius_of_gyration_valid"] is True
	assert record["sasa_A2"] is None
	assert record["sasa_valid"] is False
	assert "SASA failed" in record["sasa_error"]


def test_sasa_is_explicitly_unavailable_when_rdkit_build_lacks_freesasa(monkeypatch):
	monkeypatch.setattr(geometry, "rdFreeSASA", None)
	ensemble = generate("C", n_confs=1)
	ensemble.conformers[0].atoms.positions[0, 0] = math.nan
	result = analyze_geometry(ensemble)
	record = result["conformers"][0]

	assert result["summary"]["sasa_available"] is False
	assert result["summary"]["n_valid_sasa"] == 0
	assert result["summary"]["n_invalid_sasa"] == 1
	assert record["sasa_A2"] is None
	assert record["sasa_valid"] is False
	assert record["sasa_status"] == "unavailable"
	assert "unavailable" in record["sasa_error"]
	assert record["radius_of_gyration_valid"] is False


def test_analysis_does_not_mutate_ase_or_rdkit_coordinates():
	ensemble = generate("CCO", n_confs=2)
	ase_positions = [
		conformer.atoms.positions.copy() for conformer in ensemble.conformers
	]
	rdkit_positions = {
		conformer_id: np.asarray(
			ensemble.molecule.GetConformer(conformer_id).GetPositions()
		).copy()
		for conformer_id in ensemble.conformer_ids
	}

	analyze_geometry(ensemble)

	for conformer, positions in zip(ensemble.conformers, ase_positions):
		assert np.array_equal(conformer.atoms.positions, positions)
	for conformer_id, positions in rdkit_positions.items():
		assert np.array_equal(
			ensemble.molecule.GetConformer(conformer_id).GetPositions(), positions
		)


def test_csv_records_keep_missing_values_and_metric_validity_flags():
	ensemble = generate("CC", n_confs=2)
	ensemble.conformers[0].atoms.positions[0, 0] = math.nan

	rows = analyze_geometry(ensemble)["conformers"]
	stream = io.StringIO()
	writer = csv.DictWriter(stream, fieldnames=rows[0].keys())
	writer.writeheader()
	writer.writerows(rows)
	serialized_rows = list(csv.DictReader(io.StringIO(stream.getvalue())))

	assert rows[0]["radius_of_gyration_A"] is None
	assert rows[0]["radius_of_gyration_valid"] is False
	assert rows[1]["radius_of_gyration_A"] is not None
	assert rows[1]["radius_of_gyration_valid"] is True
	assert serialized_rows[0]["radius_of_gyration_A"] == ""
	assert serialized_rows[0]["radius_of_gyration_valid"] == "False"
	assert serialized_rows[1]["radius_of_gyration_A"] != ""
	assert serialized_rows[1]["radius_of_gyration_valid"] == "True"
