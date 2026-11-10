# EnsembleLab

`ensemblelab` is an open source Python library for quantitative analysis of molecular conformational ensembles. It provides functionality for ensemble analysis (basin visualization, RMSD heatmapping, etc.) with simple object-based workflows.

## The Focus
**How many conformers are enough?**

Computational conformer generation often requires choosing a sampling size in advance. EnsembleLab is being developed to investigate how sampling size affects ensemble level structural and energetic properties, with the eventual goal of developing adaptive methods for choosing an appropriate sampling size for a given molecule.

**Current work**

EnsembleLab currently provides tools for conformer generation, optimization, filtering, RMSD analysis, torsional analysis, clustering, and molecular descriptors. These components form the foundation for studying convergence of conformational ensembles as sampling size increases.

## Installation

EnsembleLab requires Python 3.10 or newer. It isn't published on PyPI yet, so install it from a clone of this repository.

### pip

```bash
pip install -e .              # core: generation, MMFF/UFF optimization, filtering, analysis
pip install -e ".[dev]"       # core + pytest, for development
pip install -e ".[xtb]"       # core + GFN2-xTB backend (Linux only, see below)
pytest
```

`tblite`, which provides GFN2-xTB, publishes PyPI wheels for Linux only. On macOS and Windows, install it from conda-forge (`conda install -c conda-forge tblite-python`) or use the conda environment below.

### conda

From the repository root:

```bash
conda env create -f environment.yml
conda activate ensemblelab
pytest
```

This environment installs every backend available from conda-forge, including tblite, and installs EnsembleLab in editable mode.

### Optimizer requirements

| Optimizer | Needs | Install |
| --- | --- | --- |
| `MMFFOptimizer`, MMFF/UFF in `optimize()` | RDKit (core) | included in `pip install -e .` |
| `GFN2xTBOptimizer` | `tblite` (optional) | `pip install -e ".[xtb]"` on Linux; `conda install -c conda-forge tblite-python` on macOS/Windows |
| `ORCAOptimizer` | ORCA executable (external) | download from the [ORCA forum](https://orcaforum.kofo.mpg.de), then pass `orca_command="/full/path/to/orca"` or set `command` under `[orca]` in ASE's config file (`~/.config/ase/config.ini`) |

Optional backends are only imported when you use them. `import ensemblelab` works without tblite or ORCA. Constructing a `GFN2xTBOptimizer` without tblite raises an `ImportError` that tells you how to install it. Running an `ORCAOptimizer` that can't find ORCA raises a `FileNotFoundError`.

> Note: `generators.setup()` uses MMFF → GFN2-xTB by default, so it needs tblite unless you pass your own `optimizers=`.

## Datatypes

- `Conformer`: EnsembleLab native object type storing ASE `Atoms`-based molecules; records conformer energy in kcal/mol and optimization type.
- `Ensemble`: EnsembleLab native object storing `Conformer` objects.
- [`Mol`](https://www.rdkit.org/docs/source/rdkit.Chem.rdchem.html#rdkit.Chem.rdchem.Mol): RDKit-native molecule object.
- [`Atoms`](https://ase-lib.org/ase/atoms.html): ASE-native atomic structure object.

## Current API

The package currently provides conformer generation, single-stage optimization, filtering, and RMSD-based analysis. Import the public functions directly from `ensemblelab`:

```python
from ensemblelab import Ensemble, generate
from ensemblelab.optimizers import BaseOptimizer, MMFFOptimizer, GFN2xTBOptimizer, HierarchicalOptimizer
from ensemblelab.analysis.rmsd import rmsd_matrix, rmsd_heatmap
```

## Documentation

[Visit the wiki!](https://github.com/dhruvavirneni/ensemblelab/wiki)

You can also find updated drafts of documentation in /docs, where docs will be pushed during development. The wiki will receive periodic batch updates.

## Modules

- **Generate:** implemented and documented.
- **Optimize:** implemented and documented.
- **Filter:** implemented and documented.
- **Analysis:** implemented for RMSD-based geometry analysis and documented.
- **Display:** implemented and documented for inspection summaries.
- **Descriptors:** in progress.
- **Clustering:** initial functionality complete.
- **Utils:** in progress.
