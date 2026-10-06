# Analysis

## Overview
Analysis helpers inspect geometric similarity between conformers in an `Ensemble`. The current implementation is focused on conformer-to-conformer RMSD and is defined in `ensemblelab.analysis.rmsd`. These routines read structure data without mutating the input ensemble, and they use the RDKit conformer geometry already stored on the ensemble's molecule.

## Functions
- **rmsd:** compute the best-fit RMSD between two conformers.
- **rmsd_matrix:** compute the symmetric pairwise RMSD matrix for an ensemble or a selected conformer subset.
- **rmsd_heatmap:** plot and return a pairwise RMSD heatmap for visualization.


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
