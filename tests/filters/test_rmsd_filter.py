"""RMSDFilter: energy ordering, atom selection, coordinate sync, immutability."""

from __future__ import annotations

from dataclasses import replace

import numpy as np
import pytest
from rdkit import Chem
from rdkit.Chem import rdMolAlign, rdMolTransforms

from ensemblelab import EnergyFilter, RMSDFilter
from ensemblelab.generators import Conformer, Ensemble, generate
from ensemblelab.optimizers import MMFFOptimizer


def _ensemble_from_molecule(molecule: Chem.Mol, energies: list[float]) -> Ensemble:
    """Build a synced ensemble whose ASE atoms mirror the RDKit conformers."""
    template = generate("CCO", n_confs=1).conformers[0].atoms
    conformers = []
    for rd_conformer, energy in zip(molecule.GetConformers(), energies):
        atoms = template.copy()
        atoms.positions[:] = rd_conformer.GetPositions()
        conformers.append(
            Conformer(
                id=rd_conformer.GetId(),
                atoms=atoms,
                energy=energy,
                energy_unit="kcal/mol",
                optimization_method="MMFF",
                optimization_converged=True,
                _molecule=molecule,
            )
        )
    return Ensemble(smiles="CCO", molecule=molecule, conformers=tuple(conformers), metadata={})


@pytest.fixture
def ethanol_oh_rotamers() -> Ensemble:
    """Two ethanol conformers identical except for a 120° rotation of the OH hydrogen."""
    base = MMFFOptimizer().optimize(generate("CCO", n_confs=1))
    molecule = Chem.Mol(base.molecule)
    first = molecule.GetConformer(base.conformers[0].id)

    second = Chem.Conformer(first)
    second.SetId(first.GetId() + 1)
    molecule.AddConformer(second, assignId=False)

    o_index = next(a.GetIdx() for a in molecule.GetAtoms() if a.GetSymbol() == "O")
    o_atom = molecule.GetAtomWithIdx(o_index)
    h_index = next(n.GetIdx() for n in o_atom.GetNeighbors() if n.GetSymbol() == "H")
    c_index = next(n.GetIdx() for n in o_atom.GetNeighbors() if n.GetSymbol() == "C")
    cc_index = next(
        n.GetIdx()
        for n in molecule.GetAtomWithIdx(c_index).GetNeighbors()
        if n.GetSymbol() == "C"
    )
    rotated = molecule.GetConformer(second.GetId())
    angle = rdMolTransforms.GetDihedralDeg(rotated, cc_index, c_index, o_index, h_index)
    rdMolTransforms.SetDihedralDeg(rotated, cc_index, c_index, o_index, h_index, angle + 120.0)

    return _ensemble_from_molecule(molecule, energies=[0.0, 0.5])


@pytest.fixture(scope="module")
def butanol() -> Ensemble:
    return MMFFOptimizer().optimize(generate("CCCCO", n_confs=20))


def test_filter_constructs_and_names_itself():
    assert RMSDFilter().method == "rmsd"
    assert RMSDFilter().heavy_atoms_only is False


@pytest.mark.parametrize("kwargs", [{"cutoff": 0}, {"cutoff": -1.0}])
def test_invalid_cutoff_rejected(kwargs):
    with pytest.raises(ValueError):
        RMSDFilter(**kwargs)


def test_heavy_atoms_only_must_be_bool():
    with pytest.raises(TypeError):
        RMSDFilter(heavy_atoms_only="yes")


def test_keeps_lowest_energy_and_sorts_ascending(butanol):
    filtered = RMSDFilter(cutoff=0.5).apply(butanol)

    lowest = min(butanol.conformers, key=lambda c: c.energy)
    assert filtered.conformers[0].id == lowest.id
    energies = [c.energy for c in filtered.conformers]
    assert energies == sorted(energies)
    assert set(filtered.conformer_ids) <= set(butanol.conformer_ids)


def test_kept_conformers_are_pairwise_beyond_cutoff(butanol):
    cutoff = 0.5
    filtered = RMSDFilter(cutoff=cutoff).apply(butanol)
    molecule = Chem.Mol(butanol.molecule)
    ids = filtered.conformer_ids
    for i, ref in enumerate(ids):
        for prb in ids[i + 1 :]:
            assert rdMolAlign.GetBestRMS(molecule, molecule, prbId=prb, refId=ref) > cutoff


def test_all_atom_rmsd_misses_oh_rotamer_duplicate(ethanol_oh_rotamers):
    filtered = RMSDFilter(cutoff=0.1).apply(ethanol_oh_rotamers)
    assert len(filtered.conformers) == 2


def test_heavy_atoms_only_merges_oh_rotamer_duplicate(ethanol_oh_rotamers):
    filtered = RMSDFilter(cutoff=0.1, heavy_atoms_only=True).apply(ethanol_oh_rotamers)

    kept_id, removed_id = ethanol_oh_rotamers.conformer_ids
    assert filtered.conformer_ids == (kept_id,)

    event = filtered.metadata["history"][-1]
    assert event["process"] == "filter"
    assert event["method"] == "rmsd"
    assert event["heavy_atoms_only"] is True
    assert event["kept_conformer_ids"] == [kept_id]
    assert event["removed_conformer_ids"] == [removed_id]
    assert event["duplicate_of"][removed_id]["kept_conformer_id"] == kept_id
    assert event["duplicate_of"][removed_id]["rmsd_angstrom"] < 1e-6


def test_equal_energies_keep_input_order(ethanol_oh_rotamers):
    tied = _ensemble_from_molecule(Chem.Mol(ethanol_oh_rotamers.molecule), energies=[1.0, 1.0])
    filtered = RMSDFilter(cutoff=0.1, heavy_atoms_only=True).apply(tied)
    assert filtered.conformer_ids == (tied.conformer_ids[0],)


def test_input_ensemble_is_not_modified(butanol):
    positions_before = {
        cid: butanol.molecule.GetConformer(cid).GetPositions().copy()
        for cid in butanol.conformer_ids
    }
    ids_before = butanol.conformer_ids
    history_before = list(butanol.metadata.get("history", []))

    for heavy in (False, True):
        RMSDFilter(cutoff=0.5, heavy_atoms_only=heavy).apply(butanol)

    assert butanol.conformer_ids == ids_before
    assert butanol.metadata.get("history", []) == history_before
    for conformer in butanol.conformers:
        rdkit_positions = butanol.molecule.GetConformer(conformer.id).GetPositions()
        # GetBestRMS aligns its probe in place; the filter must use a copy.
        np.testing.assert_array_equal(rdkit_positions, positions_before[conformer.id])
        np.testing.assert_allclose(rdkit_positions, conformer.atoms.positions, atol=1e-8)


def test_same_pool_gives_same_result_when_reused(butanol):
    first = RMSDFilter(cutoff=0.5).apply(butanol).conformer_ids
    second = RMSDFilter(cutoff=0.5).apply(butanol).conformer_ids
    assert first == second


def test_single_conformer_is_kept():
    single = MMFFOptimizer().optimize(generate("CCO", n_confs=1))
    filtered = RMSDFilter().apply(single)
    assert filtered.conformer_ids == single.conformer_ids
    assert filtered.metadata["history"][-1]["removed_conformer_ids"] == []


def test_missing_energies_rejected():
    with pytest.raises(ValueError, match="computed energies"):
        RMSDFilter().apply(generate("CCO", n_confs=2))


def test_desynced_coordinates_rejected(butanol):
    conformers = list(butanol.conformers)
    shifted = conformers[1].atoms.copy()
    shifted.positions += 0.5
    conformers[1] = replace(conformers[1], atoms=shifted)
    desynced = replace(butanol, conformers=tuple(conformers))

    with pytest.raises(ValueError, match="RDKit coordinates do not match"):
        RMSDFilter().apply(desynced)


def test_coordinates_synced_after_gfn2():
    """GFN2 updates ASE positions; RMSDFilter needs the RDKit copy synced too."""
    pytest.importorskip("tblite.ase")
    from ensemblelab.optimizers import GFN2xTBOptimizer

    optimized = GFN2xTBOptimizer().optimize(generate("CCO", n_confs=2))
    for conformer in optimized.conformers:
        np.testing.assert_allclose(
            optimized.molecule.GetConformer(conformer.id).GetPositions(),
            conformer.atoms.positions,
            atol=1e-8,
        )
    RMSDFilter(cutoff=0.1, heavy_atoms_only=True).apply(optimized)


def test_shared_filter_builder_works_for_other_filters(butanol):
    """Regression: _build_filtered_ensemble raised NameError for every filter."""
    filtered = EnergyFilter(window=100.0).apply(butanol)
    assert filtered.conformer_ids
    assert filtered.metadata["history"][-1]["method"] == "energy"
