"""analysis.rmsd functions must not move coordinates in the shared ensemble molecule."""

from __future__ import annotations

import numpy as np
import pytest

from ensemblelab import RMSDFilter
from ensemblelab.analysis.rmsd import rmsd, rmsd_matrix
from ensemblelab.generators import generate
from ensemblelab.optimizers import MMFFOptimizer


@pytest.fixture
def pool():
    return MMFFOptimizer().optimize(generate("CCCCO", n_confs=8))


def _rdkit_positions(ensemble):
    return {c.id: ensemble.molecule.GetConformer(c.id).GetPositions().copy() for c in ensemble.conformers}


def _assert_unchanged_and_synced(ensemble, before):
    for conformer in ensemble.conformers:
        positions = ensemble.molecule.GetConformer(conformer.id).GetPositions()
        np.testing.assert_array_equal(positions, before[conformer.id])
        np.testing.assert_allclose(positions, conformer.atoms.positions, atol=1e-8)


def test_rmsd_matrix_is_read_only(pool):
    before = _rdkit_positions(pool)
    matrix, ids = rmsd_matrix(pool)
    assert matrix.shape == (len(ids), len(ids))
    np.testing.assert_allclose(matrix, matrix.T)
    _assert_unchanged_and_synced(pool, before)


def test_pairwise_rmsd_is_read_only(pool):
    before = _rdkit_positions(pool)
    value = rmsd(pool.conformers[1], pool.conformers[0], pool.molecule)
    assert value >= 0.0
    _assert_unchanged_and_synced(pool, before)


def test_pool_can_be_analyzed_then_filtered(pool):
    """Analyzing a pool must not desync it for a later RMSDFilter pass."""
    rmsd_matrix(pool)
    RMSDFilter(cutoff=0.5).apply(pool)
