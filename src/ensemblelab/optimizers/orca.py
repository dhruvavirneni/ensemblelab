"""ORCA conformer optimization through ASE."""

from __future__ import annotations

import os
import shutil
from dataclasses import dataclass, field
from pathlib import Path
from tempfile import TemporaryDirectory
from typing import Any

from ase import Atoms
from ase.optimize import BFGS
from rdkit import Chem

from ..generators import Ensemble
from .base import BaseOptimizer

_EV_TO_KCAL_PER_MOL = 23.0605478306

_ORCA_SETUP_HINT = (
    "ORCA is an external program and is not installed with ensemblelab. "
    "Download it separately (free for academic use) from the ORCA forum, "
    "https://orcaforum.kofo.mpg.de, then either pass the full path to the "
    "executable, e.g. ORCAOptimizer(orca_command='/opt/orca/orca'), or set it "
    "once in ASE's config file (~/.config/ase/config.ini, or the file named by "
    "$ASE_CONFIG_FILE):\n\n    [orca]\n    command = /opt/orca/orca\n\n"
    "Use the full path: ORCA requires it for parallel runs, and a bare 'orca' "
    "on PATH may be an unrelated program (e.g. the GNOME screen reader)."
)


def _orca_command_from_ase_config() -> str | None:
    """Return the ORCA command from ASE's config file, if one is set."""
    from ase.config import cfg

    try:
        command = cfg.parser["orca"]["command"]
    except KeyError:
        return None
    return command.strip() or None


def _resolve_orca_command(command: str | None) -> str:
    """Locate the ORCA executable, raising an actionable error if it is missing.

    ``command`` takes precedence; otherwise ASE's ``[orca] command`` config
    entry is used, which is how ASE's own ORCA calculator locates ORCA.
    """
    source = "orca_command"
    if not command:
        command = _orca_command_from_ase_config()
        source = "ASE config [orca] command"
    if not command:
        raise FileNotFoundError(
            "No ORCA executable is configured. " + _ORCA_SETUP_HINT
        )
    resolved = shutil.which(command)
    if resolved is None:
        raise FileNotFoundError(
            f"ORCA executable {command!r} (from {source}) was not found or is "
            "not executable. " + _ORCA_SETUP_HINT
        )
    return os.path.abspath(resolved)


@dataclass(frozen=True, slots=True)
class ORCAOptimizer(BaseOptimizer):
    """Optimize conformers using ASE's ORCA calculator."""

    orca_command: str | None = None
    orca_simple_input: str = "r2SCAN-3c"
    orca_blocks: str | None = None
    workdir: str | Path | None = None
    fmax: float = 0.05
    max_steps: int = 500
    charge: int | None = None
    multiplicity: int | None = None
    _resolved_charge: int = field(init=False, repr=False, compare=False)
    _resolved_multiplicity: int = field(init=False, repr=False, compare=False)
    _resolved_command: str = field(init=False, repr=False, compare=False)

    def __post_init__(self) -> None:
        self._validate_fmax(self.fmax)
        self._validate_max_steps(self.max_steps)

    @property
    def method(self) -> str:
        """Canonical optimization method name."""
        return "ORCA"

    def optimize(self, ensemble: Ensemble) -> Ensemble:
        """Optimize all conformers after resolving the ensemble's electronic state."""
        # Locate ORCA before any work starts; checked on use, not at import.
        resolved_command = _resolve_orca_command(self.orca_command)
        resolved_charge = self._resolve_charge(ensemble, self.charge)
        resolved_multiplicity = self._resolve_multiplicity(ensemble, self.multiplicity)
        # Frozen dataclass: per-run state is set via object.__setattr__.
        object.__setattr__(self, "_resolved_command", resolved_command)
        object.__setattr__(self, "_resolved_charge", resolved_charge)
        object.__setattr__(self, "_resolved_multiplicity", resolved_multiplicity)
        try:
            return BaseOptimizer.optimize(self, ensemble)
        finally:
            object.__delattr__(self, "_resolved_command")
            object.__delattr__(self, "_resolved_charge")
            object.__delattr__(self, "_resolved_multiplicity")

    def _optimize_conformer(
        self, molecule: Chem.Mol, conformer_id: int, atoms: Atoms
    ) -> tuple[float, bool]:
        del molecule
        # Part of core ASE; imported lazily only to keep package import light.
        from ase.calculators.orca import ORCA, OrcaProfile

        def run(directory: str) -> tuple[float, bool]:
            atoms.calc = ORCA(
                profile=OrcaProfile(command=self._resolved_command),
                directory=directory,
                charge=self._resolved_charge,
                mult=self._resolved_multiplicity,
                orcasimpleinput=self.orca_simple_input,
                orcablocks=self.orca_blocks or "%pal nprocs 1 end",
            )
            optimizer = BFGS(atoms, logfile=None)
            converged = bool(optimizer.run(fmax=self.fmax, steps=self.max_steps))
            return float(atoms.get_potential_energy()) * _EV_TO_KCAL_PER_MOL, converged

        if self.workdir is None:
            with TemporaryDirectory(prefix=f"ensemblelab-orca-{conformer_id}-") as directory:
                return run(directory)

        directory = Path(self.workdir) / f"conformer-{conformer_id}"
        directory.mkdir(parents=True, exist_ok=True)
        return run(str(directory))

    def _history_details(self) -> dict[str, Any]:
        return {
        "max_steps": self.max_steps,
        "fmax": self.fmax,
        "charge": self._resolved_charge,
        "multiplicity": self._resolved_multiplicity,
        "orca_simple_input": self.orca_simple_input,
        "orca_blocks": self.orca_blocks,
        "orca_command": self._resolved_command,
        }
    # add n_cores and memory through orca? affect reproducibility and performance on other machines
