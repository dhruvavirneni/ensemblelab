# Changelog

All notable changes to EnsembleLab are recorded here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and the project uses
[semantic versioning](https://semver.org/).

## [0.2.0] - 2026-11-24

### Added

- **Ions:** `Ensemble.charge` and `Ensemble.multiplicity`, computed from the
  RDKit molecule. Charged SMILES (e.g. `CC(=O)[O-]`, `CCC[NH3+]`) are recorded
  in generation metadata and history, shown by `Ensemble.show()`, and used as
  the default charge/multiplicity by GFN2-xTB and ORCA. Every optimizer history
  entry, including MMFF/UFF, now records `charge`.

### Changed

- `generate()` raises `ValueError` for multi-fragment SMILES (salts, metal
  adducts such as `CC(=O)O.[Na+]`), whose ETKDG conformers place fragments at
  arbitrary distances.

[0.2.0]: https://github.com/dhruvavirneni/ensemblelab/releases/tag/v0.2.0

## [0.1.0] - 2026-11-17

First tagged release. The API is still pre-1.0 and may change.

### Added

- **Data model:** `Conformer` and `Ensemble` objects that hold an RDKit `Mol`
  alongside per-conformer ASE `Atoms`, energies in kcal/mol, optimization
  method and convergence status, and a provenance history. Optimizers and
  filters return a new `Ensemble` with RDKit and ASE coordinates kept in sync
  (`sync_rdkit_from_atoms`, `sync_atoms_from_rdkit` for manual edits).
- **Generation:** `generate()` embeds conformers from SMILES with RDKit ETKDGv3
  (fixed seed); `setup()` runs generation, optimization, and filtering in one call.
- **Optimization:** `MMFFOptimizer`, `UFFOptimizer`, `GFN2xTBOptimizer` (optional
  `tblite` backend, with implicit solvent), `ORCAOptimizer` (external ORCA
  executable), and `HierarchicalOptimizer` for multi-stage workflows, plus the
  functional `optimize()` / `hierarchical_optimize()` helpers. Optional backends
  are imported only when used.
- **Filtering:** `EnergyFilter`, `PopulationFilter` (Boltzmann population),
  `RMSDFilter` (optional heavy-atom selection, keeps the lowest-energy
  representative, records duplicate mapping), and `CompositeFilter`.
- **Analysis:** pairwise `rmsd`, `rmsd_matrix`, and `rmsd_heatmap`; RMSD basin
  clustering with basin entropy and `plot_basin_analysis`; heavy-atom torsion
  detection and torsion angles.
- **Convergence metrics** (computable on any subset of an ensemble's conformers):
  energies and relative energies, Boltzmann populations and effective ensemble
  size, conformer- and cluster-level entropy, radius of gyration and SASA, and
  torsional statistics (circular statistics and discrete state bins).
- **Display:** `Ensemble.show()` and `Conformer.show()` text summaries, conformer
  and metadata tables, and an initial (untested) 3D `view()` helper.
- **Packaging:** optional extras `xtb`, `plot`, `notebooks`, `dev`, and `orca`;
  a conda `environment.yml`; pytest suite and GitHub Actions CI.
- **Convergence study:** subsampling experiment (n-octanol, MMFF) with figures and writeup in experiments/sampling_convergence/docs/convergence_results.md

### Not yet included

- Conformer property metrics (dipole moments, intramolecular hydrogen bonds).
- Descriptors and utility modules (in progress).

[0.1.0]: https://github.com/dhruvavirneni/ensemblelab/releases/tag/v0.1.0
