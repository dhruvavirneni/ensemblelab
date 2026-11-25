"""Charged molecules (ions) through generation, optimization, and filtering.

The charge is written into the SMILES; ``Ensemble.charge`` and
``Ensemble.multiplicity`` are computed from the RDKit molecule and must reach
generation metadata, the display, every optimizer's history entry, and the
GFN2-xTB calculator. Multi-fragment input (salts, metal adducts) is rejected.
"""

from __future__ import annotations

import math
import warnings

import numpy as np
import pytest

from ensemblelab import RMSDFilter, generate
from ensemblelab.optimizers import MMFFOptimizer

# Perfluorooctanoate (PFOA anion) and propylammonium (protonated amine).
PFOA_ANION = "[O-]C(=O)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F"
PROPYLAMMONIUM = "CCC[NH3+]"
IONS = [(PFOA_ANION, -1), (PROPYLAMMONIUM, +1)]
ION_IDS = ["PFOA-anion", "propylammonium"]

SYNC_TOLERANCE_ANGSTROM = 1e-6


def check_sync(ensemble) -> None:
    """Assert every conformer's ASE Atoms matches its RDKit conformer."""
    molecule = ensemble.molecule
    symbols = [atom.GetSymbol() for atom in molecule.GetAtoms()]
    for conformer in ensemble.conformers:
        assert conformer.get_mol() is molecule
        assert conformer.atoms.get_chemical_symbols() == symbols
        deviation = np.abs(
            molecule.GetConformer(conformer.id).GetPositions() - conformer.atoms.positions
        ).max()
        assert deviation < SYNC_TOLERANCE_ANGSTROM


def snapshot(ensemble) -> dict[int, tuple[np.ndarray, np.ndarray]]:
    return {
        conformer.id: (
            ensemble.molecule.GetConformer(conformer.id).GetPositions().copy(),
            conformer.atoms.positions.copy(),
        )
        for conformer in ensemble.conformers
    }


@pytest.fixture(scope="module", params=IONS, ids=ION_IDS)
def ion(request):
    smiles, charge = request.param
    return smiles, charge, generate(smiles, n_confs=6)


@pytest.fixture(scope="module")
def optimized_ion(ion):
    smiles, charge, generated = ion
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        return smiles, charge, MMFFOptimizer().optimize(generated)


def test_generate_reports_charge(ion, capsys):
    _, charge, ensemble = ion
    assert ensemble.charge == charge
    assert ensemble.multiplicity == 1
    assert ensemble.metadata["charge"] == charge
    assert ensemble.metadata["multiplicity"] == 1
    assert ensemble.metadata["history"][0]["charge"] == charge
    assert ensemble.metadata["history"][0]["multiplicity"] == 1

    ensemble.show(history=True, conformers=False)
    output = capsys.readouterr().out
    assert f"Charge: {charge:+d}" in output
    assert "Multiplicity: 1" in output


def test_neutral_molecule_has_zero_charge():
    ensemble = generate("CC(=O)O", n_confs=1)
    assert ensemble.charge == 0
    assert ensemble.metadata["charge"] == 0


def test_charge_is_computed_from_molecule_not_stored(ion):
    _, _, ensemble = ion
    with pytest.raises(AttributeError):
        ensemble.charge = 0


def test_mmff_runs_with_finite_energies_and_records_charge(optimized_ion):
    _, charge, ensemble = optimized_ion
    energies = [conformer.energy for conformer in ensemble.conformers]
    assert energies and all(math.isfinite(energy) for energy in energies)
    assert all(conformer.energy_unit == "kcal/mol" for conformer in ensemble.conformers)
    event = ensemble.metadata["history"][-1]
    assert event["method"] == "MMFF"
    assert event["charge"] == charge
    assert ensemble.charge == charge


@pytest.mark.parametrize(("smiles", "charge"), IONS, ids=ION_IDS)
def test_gfn2_receives_ion_charge(smiles, charge, monkeypatch):
    pytest.importorskip("tblite")
    from ensemblelab.optimizers import GFN2xTBOptimizer, xtb

    real_loader = xtb._load_tblite_calculator
    received: list[dict] = []

    def spying_loader():
        calculator_class = real_loader()

        def spy(**options):
            received.append(options)
            return calculator_class(**options)

        return spy

    monkeypatch.setattr(xtb, "_load_tblite_calculator", spying_loader)
    ensemble = generate(smiles, n_confs=1)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        result = GFN2xTBOptimizer(max_steps=3).optimize(ensemble)

    assert received and all(options["charge"] == charge for options in received)
    assert all(options["multiplicity"] == 1 for options in received)
    event = result.metadata["history"][-1]
    assert event["method"] == "GFN2-xTB"
    assert event["charge"] == charge
    assert math.isfinite(result.conformers[0].energy)


def test_generated_ion_is_synced(ion):
    check_sync(ion[2])


def test_mmff_ion_is_synced_and_independent(ion, optimized_ion):
    _, _, generated = ion
    _, _, optimized = optimized_ion
    check_sync(optimized)
    assert optimized.molecule is not generated.molecule
    assert optimized.conformer_ids == generated.conformer_ids


@pytest.mark.parametrize("heavy_atoms_only", [False, True], ids=["all-atoms", "heavy-atoms"])
def test_rmsd_filter_on_ion(optimized_ion, heavy_atoms_only):
    _, charge, optimized = optimized_ion
    before = snapshot(optimized)
    result = RMSDFilter(cutoff=0.5, heavy_atoms_only=heavy_atoms_only).apply(optimized)

    check_sync(result)
    assert result.charge == charge
    assert result.molecule is not optimized.molecule
    assert set(result.conformer_ids) <= set(optimized.conformer_ids)
    event = result.metadata["history"][-1]
    assert sorted(event["kept_conformer_ids"] + event["removed_conformer_ids"]) == sorted(
        optimized.conformer_ids
    )
    # The input ensemble is untouched.
    for conformer_id, (rdkit_before, ase_before) in before.items():
        np.testing.assert_array_equal(
            optimized.molecule.GetConformer(conformer_id).GetPositions(), rdkit_before
        )
        conformer = next(c for c in optimized.conformers if c.id == conformer_id)
        np.testing.assert_array_equal(conformer.atoms.positions, ase_before)


@pytest.mark.parametrize(
    "smiles",
    ["CC(=O)O.[Na+]", "CC(=O)[O-].[Na+]", "[Na+].[Cl-]", "CCO.O"],
    ids=["acid+Na", "salt", "NaCl", "solvate"],
)
def test_multi_fragment_smiles_raises_clear_error(smiles):
    with pytest.raises(ValueError, match="disconnected fragments"):
        generate(smiles, n_confs=1)
