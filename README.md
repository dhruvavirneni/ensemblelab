# EnsembleLab

`ensemblelab` is an open source Python library for quantitative analysis of molecular conformational ensembles. It provides functionality for ensemble analysis (basin visualization, RMSD heatmapping, etc.) with simple object-based workflows.

## The Focus
**How many conformers are enough?**

Computational conformer generation often requires choosing a sampling size in advance. EnsembleLab is being developed to investigate how sampling size affects ensemble level structural and energetic properties, with the eventual goal of developing adaptive methods for choosing an appropriate sampling size for a given molecule.

**Current work**

EnsembleLab currently provides tools for conformer generation, optimization, filtering, RMSD analysis, torsional analysis, clustering, and molecular descriptors. These components form the foundation for studying convergence of conformational ensembles as sampling size increases.

## Install for development

```bash
pip install -e ".[dev]"
pytest
```

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
