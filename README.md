# EnsembleLab

`ensemblelab` is an open source Python library for quantitative analysis of molecular conformational ensembles. It provides functionality for ensemble analysis (basin visualization, RMSD heatmapping, etc.) with simple object-based workflows.

## The Focus
**How many conformers are enough?**

Computational conformer generation often requires choosing a sampling size in advance. EnsembleLab is being developed to investigate how sampling size affects ensemble level structural and energetic properties, with the eventual goal of developing adaptive methods for choosing an appropriate sampling size for a given molecule.

**Current work**

EnsembleLab currently provides tools for conformer generation, optimization, filtering, RMSD analysis, and clustering, plus a `convergence` module of ensemble metrics (energies, Boltzmann populations, conformer- and cluster-level entropy, radius of gyration and SASA, and torsional statistics) that can be computed on any subset of an ensemble's conformers. These components form the foundation for studying convergence of conformational ensembles as sampling size increases. Convergence criteria and sampling-size recommendations are not implemented yet.

## Installation

EnsembleLab requires Python 3.10 or newer. It isn't published on PyPI yet, so install it from a clone of this repository.

### pip

```bash
pip install -e .                  # core: generation, MMFF/UFF optimization, filtering, analysis, convergence metrics
pip install -e ".[plot]"          # core + matplotlib, for rmsd_heatmap and plot_basin_analysis
pip install -e ".[notebooks]"     # core + matplotlib, pandas, ipykernel, for the repository notebooks
pip install -e ".[dev]"           # core + pytest, for development
pip install -e ".[xtb]"           # core + GFN2-xTB backend (Linux only, see below)
pip install -e ".[notebooks,dev]" # extras combine
pytest
```

The core install is RDKit, ASE, and NumPy. SASA in `convergence.geometry` uses RDKit's FreeSASA module, which ships with the RDKit pip and conda-forge builds; if it is missing, SASA is reported as unavailable rather than approximated.

Plotting is optional: matplotlib is imported only when you call a plotting helper, which raises an `ImportError` with install instructions if it is missing.

To run the notebooks (e.g. `src/experiments/sampling_convergence/convergence_study.ipynb`), install the `notebooks` extra and select that environment as the Jupyter kernel in VS Code or Jupyter.

`tblite`, which provides GFN2-xTB, publishes PyPI wheels for Linux only. On macOS and Windows, install it from conda-forge (`conda install -c conda-forge tblite-python`) or use the conda environment below.

### conda

From the repository root:

```bash
conda env create -f environment.yml
conda activate ensemblelab
pytest
```

This environment installs every backend available from conda-forge (including tblite), plotting and notebook support (matplotlib, pandas, ipykernel), and pytest, and installs EnsembleLab in editable mode. It is the easiest way to get everything, including GFN2-xTB on macOS/Windows.

### Optimizer requirements

| Optimizer | Needs | Install |
| --- | --- | --- |
| `MMFFOptimizer`, `UFFOptimizer`, MMFF/UFF in `optimize()` | RDKit (core) | included in `pip install -e .` |
| `GFN2xTBOptimizer` | `tblite` (optional) | `pip install -e ".[xtb]"` on Linux; `conda install -c conda-forge tblite-python` on macOS/Windows |
| `ORCAOptimizer` | ORCA executable (external) | download from the [ORCA forum](https://orcaforum.kofo.mpg.de), then pass `orca_command="/full/path/to/orca"` or set `command` under `[orca]` in ASE's config file (`~/.config/ase/config.ini`) |

Optional backends are only imported when you use them. `import ensemblelab` works without tblite or ORCA. Constructing a `GFN2xTBOptimizer` without tblite raises an `ImportError` that tells you how to install it. Running an `ORCAOptimizer` that can't find ORCA raises a `FileNotFoundError`.

> Note: `generators.setup()` uses MMFF → GFN2-xTB by default, so it needs tblite unless you pass your own `optimizers=`.

## Datatypes

- `Conformer`: EnsembleLab native object type storing ASE `Atoms`-based molecules; records conformer energy in kcal/mol and optimization type.
- `Ensemble`: EnsembleLab native object storing `Conformer` objects and the RDKit `Mol` they belong to.
- [`Mol`](https://www.rdkit.org/docs/source/rdkit.Chem.rdchem.html#rdkit.Chem.rdchem.Mol): RDKit-native molecule object.
- [`Atoms`](https://ase-lib.org/ase/atoms.html): ASE-native atomic structure object.

Each conformer's coordinates live in two places: its RDKit conformer and its ASE `Atoms`. Every optimizer and filter returns a new `Ensemble` in which both copies match and which owns its own `Mol` and `Atoms`, so the input ensemble is never modified. RMSD calculations (`RMSDFilter`, `analysis.rmsd`) align a copy of the molecule and never move stored coordinates. If you edit coordinates yourself, use `ensemblelab.generators.sync_rdkit_from_atoms` / `sync_atoms_from_rdkit` to keep the copies matched.

## Current API

The package currently provides conformer generation, single- and multi-stage (hierarchical) optimization, filtering, RMSD-based analysis, and convergence-facing ensemble metrics:

```python
from ensemblelab import Ensemble, generate
from ensemblelab import EnergyFilter, PopulationFilter, RMSDFilter, CompositeFilter
from ensemblelab.optimizers import BaseOptimizer, MMFFOptimizer, GFN2xTBOptimizer, HierarchicalOptimizer
from ensemblelab.analysis.rmsd import rmsd_matrix, rmsd_heatmap
from ensemblelab.convergence.energies import analyze_energies
from ensemblelab.convergence.populations import analyze_populations
from ensemblelab.convergence.entropy import analyze_entropy
from ensemblelab.convergence.geometry import analyze_geometry
from ensemblelab.convergence.torsions import analyze_torsions
```

A minimal MMFF workflow:

```python
ens = generate("CCCCO", n_confs=50)
ens = MMFFOptimizer().optimize(ens)
ens = RMSDFilter(cutoff=0.5, heavy_atoms_only=True).apply(ens)  # keeps lowest-energy representative
analyze_populations(ens)["summary"]["effective_ensemble_size"]
```

## Documentation

[Visit the wiki!](https://github.com/dhruvavirneni/ensemblelab/wiki)

You can also find updated drafts of documentation in /docs, where docs will be pushed during development. The wiki will receive periodic batch updates.

## Modules

- **Generate:** implemented and documented.
- **Optimize:** implemented and documented.
- **Filter:** implemented and documented. Energy, Boltzmann population, RMSD (with optional heavy-atom selection and recorded duplicate mapping), and composite filters.
- **Analysis:** implemented for RMSD-based geometry analysis and documented.
- **Display:** implemented and documented for inspection summaries.
- **Descriptors:** in progress.
- **Clustering:** initial functionality complete.
- **Convergence:** ensemble metrics implemented and tested (energies, populations, entropy, geometry, torsions). Conformer properties (dipoles, hydrogen bonds), convergence curves, and convergence criteria are not implemented yet.
- **Utils:** in progress.
