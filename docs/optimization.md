# Optimize

## Overview
`optimize()` accepts a method name and returns a new `Ensemble`; it does not change the input ensemble. The returned RDKit and ASE coordinates remain synchronized.

All optimizer classes return a new `Ensemble`. Optimized energies are stored as kcal/mol on each conformer, and the optimization method, convergence status, and settings are recorded in `metadata["history"]`. If a calculation reaches its maximum number of steps, the returned conformer is still retained with `optimization_converged=False`, and the conformer ID is listed in the optimization history event. Optional backends raise actionable errors when their dependency or executable is unavailable.

## Functions
- **MMFF / UFF:** RDKit force field optimization. 
- **GFN2-xTB:**  uses [TBLite through ASE](https://tblite.readthedocs.io/en/latest/api/python.html).
- **ORCA:** uses [ASE's ORCA calculator](https://docs.ase-lib.org/ase/calculators/orca.html). Supply `orca_command`; optional `orca_simple_input`, `orca_blocks`, and `workdir` configure the calculation.


### `MMFFOptimizer`
**Overview:** 
Optimizes all conformers using RDKit's Universal Force Field (UFF). Each conformer receives an optimized geometry, an energy in kcal/mol, and a convergence status.

Metadata recorded: `max_steps`


**Usage**
```python
from ensemblelab.optimizers.mmff import UFFOptimizer

optimizer = UFFOptimizer(max_steps=500)
optimized = optimizer.optimize(ensemble)
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




