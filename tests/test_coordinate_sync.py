"""RDKit <-> ASE coordinate synchronization at every geometry boundary.

Each ``Ensemble`` holds an RDKit molecule (with its own conformer coordinates)
and per-conformer ASE ``Atoms``. Every operation that returns an ensemble must
return one where both copies agree, owns its own RDKit molecule, and leaves the
input ensemble untouched.
"""

from __future__ import annotations

import warnings

import numpy as np
import pytest

from ensemblelab import (
    CompositeFilter,
    EnergyFilter,
    PopulationFilter,
    RMSDFilter,
    generate,
)
from ensemblelab.generators import sync_atoms_from_rdkit, sync_rdkit_from_atoms
from ensemblelab.optimizers import HierarchicalOptimizer, MMFFOptimizer
from ensemblelab.optimizers.mmff import UFFOptimizer

SYNC_TOLERANCE_ANGSTROM = 1e-6


def check_sync(ensemble) -> None:
    """Assert every conformer's ASE Atoms matches its RDKit conformer."""
    molecule = ensemble.molecule
    rdkit_ids = {conformer.GetId() for conformer in molecule.GetConformers()}
    rdkit_symbols = [atom.GetSymbol() for atom in molecule.GetAtoms()]
    for conformer in ensemble.conformers:
        assert conformer.id in rdkit_ids, f"conformer {conformer.id} missing from RDKit Mol"
        assert conformer.get_mol() is molecule, f"conformer {conformer.id} bound to another Mol"
        assert conformer.atoms.get_chemical_symbols() == rdkit_symbols
        rdkit_positions = molecule.GetConformer(conformer.id).GetPositions()
        deviation = np.abs(rdkit_positions - conformer.atoms.positions).max()
        assert deviation < SYNC_TOLERANCE_ANGSTROM, (
            f"conformer {conformer.id}: max |RDKit - ASE| = {deviation:.2e} angstrom"
        )


def snapshot(ensemble) -> dict[int, tuple[np.ndarray, np.ndarray]]:
    """Copy both coordinate sets of every conformer."""
    return {
        conformer.id: (
            ensemble.molecule.GetConformer(conformer.id).GetPositions().copy(),
            conformer.atoms.positions.copy(),
        )
        for conformer in ensemble.conformers
    }


def assert_unchanged(ensemble, before) -> None:
    assert snapshot(ensemble).keys() == before.keys()
    for conformer_id, (rdkit_after, ase_after) in snapshot(ensemble).items():
        rdkit_before, ase_before = before[conformer_id]
        np.testing.assert_array_equal(rdkit_after, rdkit_before)
        np.testing.assert_array_equal(ase_after, ase_before)


def assert_independent(result, source) -> None:
    """The result must not share the RDKit Mol or Atoms objects with its source."""
    assert result.molecule is not source.molecule
    source_atoms = {id(conformer.atoms) for conformer in source.conformers}
    assert not any(id(conformer.atoms) in source_atoms for conformer in result.conformers)


@pytest.fixture(scope="module")
def generated():
    # Butanol: several rotatable bonds, an OH, and enough conformers to dedupe.
    return generate("CCCCO", n_confs=8)


@pytest.fixture(scope="module")
def optimized(generated):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        return MMFFOptimizer().optimize(generated)


def test_generate_is_synced(generated):
    check_sync(generated)


@pytest.mark.parametrize(
    "optimizer",
    [
        MMFFOptimizer(),
        UFFOptimizer(),
        HierarchicalOptimizer([MMFFOptimizer(), UFFOptimizer()]),
    ],
    ids=["MMFF", "UFF", "MMFF->UFF"],
)
def test_optimizer_returns_synced_independent_copy(generated, optimizer):
    before = snapshot(generated)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        result = optimizer.optimize(generated)

    check_sync(result)
    assert_independent(result, generated)
    assert_unchanged(generated, before)
    assert result.conformer_ids == generated.conformer_ids
    # Optimization really moved atoms, so the sync check is not vacuous.
    assert any(
        not np.allclose(result.conformers[i].atoms.positions, before[c.id][1])
        for i, c in enumerate(generated.conformers)
    )


def test_hierarchical_is_synced_after_every_stage(generated):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        current = generated
        for stage in (MMFFOptimizer(), UFFOptimizer(), MMFFOptimizer()):
            current = stage.optimize(current)
            check_sync(current)


@pytest.mark.parametrize(
    "filter_",
    [
        RMSDFilter(cutoff=0.5),
        RMSDFilter(cutoff=0.5, heavy_atoms_only=True),
        EnergyFilter(top_n=3),
        EnergyFilter(window=2.0),
        PopulationFilter(cumulative_cutoff=0.9),
        CompositeFilter([EnergyFilter(window=5.0), RMSDFilter(cutoff=0.5)]),
    ],
    ids=["rmsd", "rmsd-heavy", "energy-top_n", "energy-window", "population", "composite"],
)
def test_filter_returns_synced_independent_copy(optimized, filter_):
    before = snapshot(optimized)
    result = filter_.apply(optimized)

    check_sync(result)
    assert_independent(result, optimized)
    assert_unchanged(optimized, before)
    assert set(result.conformer_ids) <= set(optimized.conformer_ids)
    for conformer in result.conformers:
        source = next(c for c in optimized.conformers if c.id == conformer.id)
        assert conformer.energy == source.energy
        assert conformer.energy_unit == source.energy_unit


def test_filtered_ensemble_edits_do_not_leak(optimized):
    result = EnergyFilter(top_n=2).apply(optimized)
    before = snapshot(optimized)
    target = result.conformers[0]

    target.atoms.positions += 1.0
    sync_rdkit_from_atoms(result.molecule, target.id, target.atoms)

    check_sync(result)
    assert_unchanged(optimized, before)


def test_rmsd_filter_does_not_align_real_coordinates(optimized):
    before = snapshot(optimized)
    result = RMSDFilter(cutoff=0.5).apply(optimized)
    assert_unchanged(optimized, before)
    for conformer in result.conformers:
        np.testing.assert_array_equal(conformer.atoms.positions, before[conformer.id][1])


def test_sync_helpers_round_trip(generated):
    from rdkit import Chem

    molecule = Chem.Mol(generated.molecule)
    conformer = generated.conformers[0]
    atoms = conformer.atoms.copy()

    atoms.positions += np.array([0.5, -0.25, 1.0])
    sync_rdkit_from_atoms(molecule, conformer.id, atoms)
    np.testing.assert_allclose(molecule.GetConformer(conformer.id).GetPositions(), atoms.positions)

    fresh = conformer.atoms.copy()
    sync_atoms_from_rdkit(molecule, conformer.id, fresh)
    np.testing.assert_allclose(fresh.positions, atoms.positions)
    # The source ensemble is untouched.
    assert not np.allclose(conformer.atoms.positions, atoms.positions)


def test_sync_helpers_reject_mismatched_pairs(generated):
    from rdkit import Chem

    molecule = Chem.Mol(generated.molecule)
    atoms = generated.conformers[0].atoms.copy()

    with pytest.raises(ValueError, match="no conformer with ID"):
        sync_rdkit_from_atoms(molecule, 10_000, atoms)
    with pytest.raises(ValueError, match="atoms but the RDKit molecule has"):
        sync_atoms_from_rdkit(molecule, generated.conformers[0].id, atoms[:-1])
    swapped = atoms.copy()
    symbols = swapped.get_chemical_symbols()
    symbols[0], symbols[-1] = symbols[-1], symbols[0]
    swapped.set_chemical_symbols(symbols)
    with pytest.raises(ValueError, match="element order"):
        sync_rdkit_from_atoms(molecule, generated.conformers[0].id, swapped)


def test_gfn2_returns_synced_independent_copy():
    pytest.importorskip("tblite")
    from ensemblelab.optimizers import GFN2xTBOptimizer

    ensemble = generate("CCO", n_confs=2)
    before = snapshot(ensemble)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        result = HierarchicalOptimizer([MMFFOptimizer(), GFN2xTBOptimizer()]).optimize(ensemble)

    check_sync(result)
    assert_independent(result, ensemble)
    assert_unchanged(ensemble, before)
