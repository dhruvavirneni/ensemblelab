"""Packaging tests: the core works without optional backends (tblite, ORCA).

tblite absence is simulated with ``sys.modules["tblite"] = None``, which makes
any ``import tblite`` raise ImportError even when tblite is installed.
"""

from __future__ import annotations

import importlib
import math
import os
import sys
from pathlib import Path

import pytest

import ensemblelab
from ensemblelab.generators import generate

PACKAGE_DIR = Path(ensemblelab.__file__).parent


@pytest.fixture
def fresh_ensemblelab_without_tblite(monkeypatch):
    """Hide tblite and re-import ensemblelab from scratch.

    All ``ensemblelab*`` modules are restored afterwards so later tests keep
    using the original class objects.
    """
    saved = {
        name: module
        for name, module in sys.modules.items()
        if name == "ensemblelab" or name.startswith("ensemblelab.")
    }
    for name in saved:
        del sys.modules[name]
    monkeypatch.setitem(sys.modules, "tblite", None)
    monkeypatch.setitem(sys.modules, "tblite.ase", None)
    try:
        yield importlib.import_module("ensemblelab")
    finally:
        for name in [n for n in sys.modules if n == "ensemblelab" or n.startswith("ensemblelab.")]:
            del sys.modules[name]
        sys.modules.update(saved)


@pytest.fixture(scope="module")
def ethanol():
    return generate("CCO", n_confs=3)


def _package_module_names() -> list[str]:
    names = []
    for path in sorted(PACKAGE_DIR.rglob("*.py")):
        parts = path.relative_to(PACKAGE_DIR).with_suffix("").parts
        if parts[-1] == "__init__":
            parts = parts[:-1]
        names.append(".".join(("ensemblelab", *parts)))
    return names


# --- tblite / GFN2-xTB -----------------------------------------------------


def test_import_ensemblelab_without_tblite(fresh_ensemblelab_without_tblite):
    ens = fresh_ensemblelab_without_tblite
    assert ens.Ensemble is not None
    assert ens.generate is not None
    # Lazy optimizer re-exports resolve without touching tblite.
    assert ens.MMFFOptimizer is not None
    assert ens.GFN2xTBOptimizer is not None
    assert sys.modules["tblite"] is None


def test_no_module_imports_tblite_at_top_level(fresh_ensemblelab_without_tblite):
    """Import every module in the package; none may require tblite to import.

    Modules that fail for an unrelated missing optional dependency (e.g.
    matplotlib in plotting helpers) are tolerated here; only tblite is checked.
    """
    for name in _package_module_names():
        try:
            importlib.import_module(name)
        except ImportError as error:
            assert "tblite" not in str(error), f"{name} imports tblite at import time"
            assert not (error.name or "").startswith("tblite"), name


@pytest.mark.parametrize("via", ["package", "submodule"])
def test_gfn2_optimizer_without_tblite_raises_helpful_error(
    fresh_ensemblelab_without_tblite, via
):
    if via == "package":
        optimizer_cls = fresh_ensemblelab_without_tblite.GFN2xTBOptimizer
    else:
        optimizer_cls = importlib.import_module("ensemblelab.optimizers.xtb").GFN2xTBOptimizer

    with pytest.raises(ImportError) as excinfo:
        optimizer_cls()

    message = str(excinfo.value)
    assert "GFN2xTBOptimizer requires tblite" in message
    assert "pip install ensemblelab[xtb]" in message
    assert "conda install -c conda-forge tblite-python" in message


def test_legacy_optimize_gfn2_without_tblite_raises_helpful_error(
    fresh_ensemblelab_without_tblite,
):
    legacy = importlib.import_module("ensemblelab.optimize")
    ensemble = fresh_ensemblelab_without_tblite.generate("C", n_confs=1)
    with pytest.raises(ImportError, match="conda install -c conda-forge tblite-python"):
        legacy.optimize(ensemble, method="GFN2-xTB")


def test_gfn2_optimization_with_tblite(ethanol):
    pytest.importorskip("tblite")
    pytest.importorskip("tblite.ase")
    from ensemblelab.optimizers import GFN2xTBOptimizer

    small = generate("CCO", n_confs=1)
    optimized = GFN2xTBOptimizer().optimize(small)

    (conformer,) = optimized.conformers
    assert conformer.optimization_method == "GFN2-xTB"
    assert conformer.energy_unit == "kcal/mol"
    assert math.isfinite(conformer.energy)
    # Ethanol total GFN2 energy is roughly -11.4 Eh, i.e. about -7150 kcal/mol.
    assert -7300 < conformer.energy < -7000


# --- MMFF (core only) ------------------------------------------------------


def test_mmff_runs_with_core_dependencies_only(fresh_ensemblelab_without_tblite):
    ens = fresh_ensemblelab_without_tblite
    ensemble = ens.generate("CCO", n_confs=3)
    original_positions = [c.atoms.positions.copy() for c in ensemble.conformers]

    optimized = ens.MMFFOptimizer().optimize(ensemble)

    assert [c.id for c in optimized.conformers] == [c.id for c in ensemble.conformers]
    for conformer in optimized.conformers:
        assert conformer.optimization_method == "MMFF"
        assert conformer.energy_unit == "kcal/mol"
        assert math.isfinite(conformer.energy)
        assert conformer.optimization_converged
    # The input ensemble is not mutated.
    for before, conformer in zip(original_positions, ensemble.conformers):
        assert (before == conformer.atoms.positions).all()
        assert conformer.energy is None


# --- ORCA (external executable) --------------------------------------------


@pytest.fixture
def no_ase_orca_config(monkeypatch):
    orca_module = importlib.import_module("ensemblelab.optimizers.orca")
    monkeypatch.setattr(orca_module, "_orca_command_from_ase_config", lambda: None)
    return orca_module


def test_orca_optimizer_constructs_without_orca(no_ase_orca_config):
    # The executable is checked on use, not at construction or import.
    no_ase_orca_config.ORCAOptimizer()


def test_orca_without_command_raises_setup_error(no_ase_orca_config, ethanol):
    with pytest.raises(FileNotFoundError) as excinfo:
        no_ase_orca_config.ORCAOptimizer().optimize(ethanol)
    message = str(excinfo.value)
    assert "No ORCA executable is configured" in message
    assert "not installed with ensemblelab" in message
    assert "orca_command" in message
    assert "[orca]" in message


def test_orca_with_missing_executable_raises_setup_error(no_ase_orca_config, ethanol, tmp_path):
    missing = str(tmp_path / "no-such-orca")
    with pytest.raises(FileNotFoundError, match="not found or is not executable"):
        no_ase_orca_config.ORCAOptimizer(orca_command=missing).optimize(ethanol)


def test_orca_command_falls_back_to_ase_config(monkeypatch):
    orca_module = importlib.import_module("ensemblelab.optimizers.orca")
    # Any real executable stands in for ORCA; it is resolved, not run.
    monkeypatch.setattr(orca_module, "_orca_command_from_ase_config", lambda: sys.executable)
    assert orca_module._resolve_orca_command(None) == os.path.abspath(sys.executable)
