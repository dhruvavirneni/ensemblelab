# ensemblelab

`ensemblelab` is an open-source Python library for quantitative analysis of molecular conformational ensembles. It fills the gap between creating conformers and describing flexibility at the ensemble level.

## Current capability

Generate a reproducible, unoptimized conformer ensemble with aligned RDKit and ASE representations:

```python
from ensemblelab import generate

ensemble = generate("CCO", n_confs=20)
print(ensemble.metadata["n_conformers_generated"])
assert ensemble.conformers[0].energy is None  # assigned by optimize() later
```

The returned `Ensemble` stores the input SMILES, a hydrogenated RDKit molecule containing all conformers, ASE `Atoms` snapshots, and generation provenance. Each conformer owns its own placeholder `energy` value (`ensemble.conformers[0].energy`), which remains `None` until optimization.

## Install for development

```bash
pip install -e ".[dev]"
pytest
```

## Roadmap

1. Generate, optimize, filter, and cluster conformer ensembles.
2. Define and validate an interpretable conformational fingerprint.
3. Benchmark across molecule classes and solvent models.
4. Investigate experimental-property correlations and molecular-ML embeddings.
