# Optimize

## Overview
`optimize()` accepts a method name and returns a new `Ensemble`; it does not change the input ensemble. The returned RDKit and ASE coordinates remain synchronized.

All optimizer classes return a new `Ensemble`. Optimized energies are stored as kcal/mol on each conformer, and the optimization method, convergence status, and settings are recorded in `metadata["history"]`. If a calculation reaches its maximum number of steps, the returned conformer is still retained with `optimization_converged=False`, and the conformer ID is listed in the optimization history event. Optional backends raise actionable errors when their dependency or executable is unavailable.

**Convergence status.** Each optimized conformer stores `optimization_converged`, which `Ensemble.show()` displays as ✓/✗ (N/A before optimization).
- MMFF/UFF: `True` when RDKit's `Minimize()` returns 0, i.e. RDKit's default force and energy tolerances were met within `max_steps`. `fmax` is not used.
- GFN2-xTB/ORCA: `True` when ASE's BFGS reaches a maximum per-atom force below `fmax` (eV/Å) within `max_steps`. A geometry that is already converged returns `True` after zero steps.
- `False` only means the step limit was reached. A failed calculation raises an exception and stops the run instead.
- The criteria differ between backends, so ✓ from MMFF and ✓ from GFN2-xTB are not equivalent statements.

## Functions
- **MMFF / UFF:** RDKit force field optimization. 
- **GFN2-xTB:**  uses [TBLite through ASE](https://tblite.readthedocs.io/en/latest/api/python.html).
- **Hierarchical:** runs several optimizers in sequence, e.g. MMFF followed by GFN2-xTB.
- **ORCA:** uses [ASE's ORCA calculator](https://docs.ase-lib.org/ase/calculators/orca.html). Supply `orca_command` or configure ASE (see below); optional `orca_simple_input`, `orca_blocks`, and `workdir` configure the calculation.


### `MMFFOptimizer` / `UFFOptimizer`
**Overview:** 
`MMFFOptimizer` optimizes all conformers with RDKit's MMFF94 force field. `UFFOptimizer` (in the same module, not re-exported from `ensemblelab.optimizers`) uses RDKit's Universal Force Field (UFF) instead. Each conformer receives an optimized geometry, an energy in kcal/mol, and a convergence status. MMFF raises `RuntimeError` if MMFF parameters are unavailable for the molecule.

Metadata recorded: `max_steps`


**Usage**
```python
from ensemblelab.optimizers.mmff import MMFFOptimizer, UFFOptimizer

optimized = MMFFOptimizer(max_steps=500).optimize(ensemble)
optimized_uff = UFFOptimizer(max_steps=500).optimize(ensemble)
```



### `GFN2xTBOptimizer`
**Overview:** 
Optimizes conformers using the GFN2-xTB semiempirical quantum chemical method through TBLite and ASE. Each conformer receives an optimized geometry, an energy in kcal/mol, and a convergence status. **GFN2-xTB optimization is an optional backend powered by TBLite.** Install the core EnsembleLab package normally. TBLite is required only when
using `GFN2xTBOptimizer`.

Because TBLite contains compiled native code, installation availability
depends on the user's Python version and operating system. On Linux, use
`pip install ensemblelab[xtb]`. On macOS/Windows, use
`conda install -c conda-forge tblite-python`. Constructing a
`GFN2xTBOptimizer` without TBLite raises an `ImportError` with these
instructions, before any conformers are processed.

Metadata recorded: `fmax`, `max_steps`, `solvent`, `charge`, `multiplicity`


**Usage**
```python
from ensemblelab.optimizers.xtb import GFN2xTBOptimizer

optimizer = GFN2xTBOptimizer(
    fmax=0.05,
    max_steps=500,
    solvent=None,
)
optimized = optimizer.optimize(ensemble)
```



### `ORCAOptimizer`

**Overview:** 
Optimizes conformers using ORCA through ASE's ORCA calculator. ORCA performs the underlying quantum chemical calculation while ASE manages the optimization workflow. Each conformer receives an optimized geometry, an energy in kcal/mol, and a convergence status.

Metadata recorded: `fmax`, `max_steps`, `charge`, `multiplicity`, `orca_simple_input`, `orca_blocks`, `orca_command` (resolved absolute path)

ORCA is an external program and is not installed with EnsembleLab. Pass the
full path to the executable as `orca_command`, or leave it unset and configure
ASE's config file (`~/.config/ase/config.ini` or `$ASE_CONFIG_FILE`):

```ini
[orca]
command = /opt/orca/orca
```

The executable is located when `optimize()` is called, not at import or
construction. If it can't be found, `optimize()` raises `FileNotFoundError`
with setup instructions. A bare `orca` on `PATH` is never assumed, because on
Linux it is often the GNOME screen reader.


**Usage**
```python
from ensemblelab.optimizers.orca import ORCAOptimizer

optimizer = ORCAOptimizer(
    orca_command="/path/to/orca",
    orca_simple_input="r2SCAN-3c",
    max_steps=500,
    fmax=0.05,
)

optimized = optimizer.optimize(ensemble)
```



### `HierarchicalOptimizer`
**Overview:** 
Applies a sequence of optimizers in order, passing each stage's output ensemble to the next, e.g. a fast MMFF pre-optimization followed by GFN2-xTB refinement. Each stage returns a new `Ensemble`, and the input is never modified. Every stage appends its own optimization event to `metadata["history"]`, so per-stage convergence remains traceable.

- Each conformer's `energy`, `optimization_method` and `optimization_converged` reflect the **last** stage only. A conformer that hit MMFF's step limit but then converged under GFN2-xTB shows ✓. The MMFF result is still in that stage's history entry.
- No pruning or deduplication happens between stages. To prune, run the stages and filters yourself, e.g. MMFF → `EnergyFilter` / `RMSDFilter` → GFN2-xTB.
- Construction fails immediately if a stage's optional backend is missing (e.g. `GFN2xTBOptimizer` without tblite), before any stage runs.


**Usage**
```python
from ensemblelab.optimizers import GFN2xTBOptimizer, HierarchicalOptimizer, MMFFOptimizer

optimizer = HierarchicalOptimizer([MMFFOptimizer(), GFN2xTBOptimizer()])
optimized = optimizer.optimize(ensemble)
```
