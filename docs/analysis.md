# Analysis

## Overview
Analysis helpers inspect geometric similarity and conformational basins in an `Ensemble`. RMSD routines are defined in `ensemblelab.analysis.rmsd`; basin routines are defined in `ensemblelab.analysis.clusters`. They use RDKit conformer geometry stored on the ensemble's molecule and accept `Conformer` objects rather than notebook-specific dictionaries. Basin calculations leave the source ensemble unchanged.

## Functions
- **rmsd:** compute the best-fit RMSD between two conformers.
- **rmsd_matrix:** compute the symmetric pairwise RMSD matrix for an ensemble or a selected conformer subset.
- **rmsd_heatmap:** plot and return a pairwise RMSD heatmap for visualization.
- **cluster_basins:** group conformers by RMSD to each basin's representative.
- **basin_rmsd_matrix:** return pairwise RMSD data ordered by basin membership.
- **basin_entropy:** calculate the Shannon entropy of basin populations.
- **plot_basin_analysis:** plot basin occupancy and the basin-ordered RMSD matrix.

### Basin Analysis
`cluster_basins` uses a greedy representative-based RMSD threshold (0.75 Å by default) and returns `Basin` values containing the representative and member `Conformer` objects. When every selected conformer has an energy in kcal/mol, basin populations are calculated from Boltzmann weights and normalized over the selected conformers. For unoptimized ensembles, clustering still works and population and entropy results are `None`.

```python
from ensemblelab.analysis.clusters import (
    basin_entropy,
    basin_rmsd_matrix,
    cluster_basins,
    plot_basin_analysis,
)

basins = cluster_basins(ensemble, rmsd_threshold=0.75, temperature=298.15)
matrix, conformer_ids, basin_ids = basin_rmsd_matrix(ensemble, basins)
entropy = basin_entropy(basins)
figure, axes = plot_basin_analysis(ensemble, basins, show=False)
```

Plotting requires Matplotlib; it is imported only when `plot_basin_analysis` is called.


### `rmsd`
**Overview**
Computes the best-fit root-mean-square deviation between two conformers using RDKit's RMSD routine. Values are reported in Angstroms and are appropriate for comparing conformer geometry after alignment.

**Usage**
```python
from ensemblelab.analysis.rmsd import rmsd

value = rmsd(conformer_a, conformer_b, ensemble.molecule)
```


### `rmsd_matrix`
**Overview**
Constructs the full pairwise RMSD matrix for an ensemble, or for an explicitly supplied subset of conformers. The output is a symmetric NumPy array with zero values on the diagonal. The function validates that the requested conformers are unique and belong to the same ensemble before computing distances.

The returned tuple is:

- `matrix`: NumPy array of pairwise RMSD values in Angstroms.
- `conformer_ids`: ordered conformer IDs matching the matrix axes.

**Usage**
```python
from ensemblelab.analysis.rmsd import rmsd_matrix

matrix, conformer_ids = rmsd_matrix(ensemble)
```


### `rmsd_heatmap`
**Overview**
Plots a heatmap of the pairwise RMSD matrix and returns the matrix plus the Matplotlib axes object used for rendering. The axes are labeled by conformer ID, and a colorbar indicates RMSD in Angstroms. When `reference_conformer` is provided, the reference conformer is moved to the first row and column.

**Usage**
```python
from ensemblelab.analysis.rmsd import rmsd_heatmap

matrix, axes = rmsd_heatmap(
    ensemble,
    conformers=ensemble.conformers,
    reference_conformer=ensemble.conformers[0],
    show=True,
)
```
