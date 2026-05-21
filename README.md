# Molecular Flexibility in Solvent Environments

This project studies how organic molecules' conformational landscapes change under environmental perturbation.

The pipeline will compare baseline vacuum conformer ensembles against solvent-perturbed ensembles using structural, thermodynamic, basin, torsional, and visualization analyses.

## Planned Pipeline

1. Generate conformers from SMILES.
2. Optimize geometries and compute conformer energies.
3. Filter high-energy or invalid conformers.
4. Compute pairwise RMSD and remove duplicates.
5. Cluster conformers into structural basins.
6. Compute Boltzmann populations and entropy-like ensemble metrics.
7. Analyze torsional variance and conformational accessibility.
8. Repeat the same workflow under implicit solvent conditions.
9. Compare vacuum and solvent ensembles.

## Repository Layout

```text
src/            Pipeline modules
notebooks/      Exploratory notebooks
data/raw/       Input molecule sets
data/processed/ Intermediate data products
data/results/   Computed metrics and tables
data/figures/   Generated plots and visualizations
docs/           Project notes and research log
```

