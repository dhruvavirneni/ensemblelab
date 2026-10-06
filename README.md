# ensemblelab

`ensemblelab` is an open source Python library for quantitative analysis of molecular conformational ensembles. It provides functionality for ensemble analysis (basin visualization, RMSD heatmapping, etc.) with simple object-based workflows.

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

- [Generation](docs/generation.md): conformer generation, `Conformer`, `Ensemble`, and `generate()`
- [Optimization](docs/optimization.md): optimizer backends and optimization provenance
- [Filtering](docs/filtering.md): energy, population, and RMSD-based filtering workflows
- [Analysis](docs/analysis.md): geometric analysis, RMSD matrix, and heatmap utilities
- [Display](docs/display.md): ensemble/conformer display and inspection helpers

## Modules

- **Generate:** implemented and documented.
- **Optimize:** implemented and documented.
- **Filter:** implemented and documented.
- **Analysis:** implemented for RMSD-based geometry analysis and documented.
- **Display:** implemented and documented for inspection summaries.
- **Descriptors:** in progress.
- **Clustering:** in progress.
- **Utils:** in progress.