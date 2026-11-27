# EnsembleLab

[![tests](https://github.com/dhruvavirneni/ensemblelab/actions/workflows/tests.yml/badge.svg)](https://github.com/dhruvavirneni/ensemblelab/actions/workflows/tests.yml)
![version](https://img.shields.io/badge/version-0.3.0-blue)
![python](https://img.shields.io/badge/python-3.10%2B-blue)
[![license](https://img.shields.io/badge/license-MIT-green)](LICENSE)

**EnsembleLab** is a Python library for generating, optimizing, filtering, and analyzing molecular conformational ensembles. It provides one `Ensemble` object that keeps conformer geometries, energies, and a full provenance history together, built on [RDKit](https://www.rdkit.org) and [ASE](https://ase-lib.org).

It is developed alongside a research question: **how many conformers are enough?** Conformer generation usually asks you to fix a sample size in advance. EnsembleLab is being used to study how ensemble-level properties converge as sampling increases, with the long-term goal of predicting a suitable sample size for a given molecule.


## Features

- **Generation:** conformers from SMILES with RDKit ETKDGv3, using a fixed seed for reproducibility.
- **Optimization:** MMFF and UFF (built in), GFN2-xTB (optional, via tblite), and ORCA (external), plus multi-stage workflows such as MMFF → GFN2-xTB.
- **Filtering:** energy window, Boltzmann population, and symmetry-aware RMSD deduplication, which records the IDs of removed conformers.
- **Analysis:** pairwise RMSD matrices and heatmaps, and RMSD-based clustering.
- **Convergence metrics:** energies, Boltzmann populations, conformer- and cluster-level entropy, radius of gyration, SASA, and torsional statistics. Each can be computed on any subset of an ensemble.
- **Ions:** charged molecules written as SMILES (e.g. `CC(=O)[O-]`), with the charge tracked through every step.
- **Provenance:** every step returns a new `Ensemble` and appends settings, units, seeds, and convergence status to `metadata["history"]`. The input ensemble is never modified.

## Installation

EnsembleLab requires Python 3.10+. It is not on PyPI yet, so install it from a clone:

```bash
git clone https://github.com/dhruvavirneni/ensemblelab.git
cd ensemblelab
pip install -e .
```

The core install needs only RDKit, ASE, and NumPy. Optional extras:

| Extra | Adds |
| --- | --- |
| `pip install -e ".[plot]"` | matplotlib, for RMSD heatmaps and basin plots |
| `pip install -e ".[notebooks]"` | matplotlib, pandas, ipykernel, for the example notebooks |
| `pip install -e ".[xtb]"` | tblite, for GFN2-xTB (**Linux only**; see below) |
| `pip install -e ".[dev]"` | pytest |

**Everything at once (recommended on macOS/Windows):** tblite publishes pip wheels for Linux only, so the conda environment is the easiest way to get GFN2-xTB elsewhere:

```bash
conda env create -f environment.yml
conda activate ensemblelab
```

**ORCA** is a separate program. Download it from the [ORCA forum](https://orcaforum.kofo.mpg.de), then pass `orca_command="/full/path/to/orca"` to `ORCAOptimizer`. See [docs/optimization.md](docs/optimization.md) for setup.

Optional backends are imported only when used, so `import ensemblelab` always works, and a missing backend raises an error that explains how to install it.

## Quickstart

```python
from ensemblelab import RMSDFilter, generate
from ensemblelab.optimizers import MMFFOptimizer
from ensemblelab.convergence.populations import analyze_populations

ens = generate("CCCCO", n_confs=50)                              # 1-butanol
ens = MMFFOptimizer().optimize(ens)                              # energies in kcal/mol
ens = RMSDFilter(cutoff=0.5, heavy_atoms_only=True).apply(ens)   # remove near-duplicates
ens.show()                                                       # summary + conformer table

summary = analyze_populations(ens)["summary"]                    # Boltzmann weights at 298.15 K
summary["effective_ensemble_size"], summary["population_entropy"]
```

For a two-stage workflow (needs tblite):

```python
from ensemblelab.optimizers import GFN2xTBOptimizer, HierarchicalOptimizer

ens = HierarchicalOptimizer([MMFFOptimizer(), GFN2xTBOptimizer()]).optimize(generate("CCCCO", n_confs=50))
```

Ions work the same way. Put the charge in the SMILES, and it is passed to every backend:

```python
ion = generate("CC(=O)[O-]", n_confs=25)
ion.charge   # -1
```

## Documentation

| Topic | |
| --- | --- |
| Generation, `Ensemble`/`Conformer`, ions | [docs/generation.md](docs/generation.md) |
| Optimizers and backend setup | [docs/optimization.md](docs/optimization.md) |
| Filters | [docs/filtering.md](docs/filtering.md) |
| RMSD and clustering analysis | [docs/analysis.md](docs/analysis.md) |
| Convergence metrics | [docs/convergence.md](docs/convergence.md) |
| Display (`show()`) | [docs/display.md](docs/display.md) |

The [`docs/`](docs/) folder is updated during development and is the most current reference. The [wiki](https://github.com/dhruvavirneni/ensemblelab/wiki) receives periodic batch updates and may lag behind.

Worked examples are in the notebooks under [`tests/`](tests/) (for example, [`tests/ions/test_ions.ipynb`](tests/ions/test_ions.ipynb)), and the pilot convergence study is in [`experiments/sampling_convergence/`](experiments/sampling_convergence/).

Find the writeup [here](https://github.com/dhruvavirneni/ensemblelab/blob/main/experiments/sampling_convergence/docs/convergence_results.md)

## Roadmap

Implemented: generation, optimization, filtering, RMSD analysis and clustering, convergence metrics, and ion support.

Not yet implemented:
- conformer-level properties (dipole moments, hydrogen bonds) and collision cross sections;
- convergence curves, convergence criteria, and sample-size recommendations;
- a versioned conformational fingerprint (ensemble descriptor set);
- ML models that predict sampling requirements.

See [ROADMAP.md](ROADMAP.md) and the [CHANGELOG](CHANGELOG.md).

## Development

```bash
pip install -e ".[dev]"
pytest
```

GFN2-xTB tests are skipped automatically when tblite is not installed. Contributor guidelines and scientific conventions are in [AGENTS.md](AGENTS.md).

## Citation

If you use EnsembleLab in your research, please cite it using [CITATION.cff](CITATION.cff) (GitHub's "Cite this repository" button generates BibTeX and APA from it).

## License

[MIT](LICENSE)
