from collections.abc import Sequence

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.axes import Axes
from rdkit.Chem import rdMolAlign

from ensemblelab.generators import Conformer, Ensemble


def rmsd(conformer1: Conformer, conformer2: Conformer, molecule) -> float:
    molecule1 = conformer1.get_mol()
    molecule2 = conformer2.get_mol()

    return rdMolAlign.GetBestRMS(
        molecule1,
        molecule2,
        prbId=conformer1.id,
        refId=conformer2.id,
    )

def rmsd_matrix(ensemble: Ensemble, conformers: Sequence[Conformer] | None = None) -> tuple[np.ndarray, list[int]]:
    """Compute the pairwise RMSD matrix for an ensemble.

    Parameters
    ----------
    ensemble
        Ensemble whose RDKit molecule contains the conformer coordinates.
    conformers
        Optional ordered subset of conformers to include. The default uses
        ``ensemble.conformers``.

    Returns
    -------
    numpy.ndarray
        The RMSD matrix. Matrix labels use conformer IDs.
    """
    selected_conformers = list(
        ensemble.conformers if conformers is None else conformers
    )
    if not selected_conformers:
        raise ValueError("ensemble must contain at least one conformer.")

    ensemble_ids = set(ensemble.conformer_ids)
    selected_ids = [conformer.id for conformer in selected_conformers]
    if len(selected_ids) != len(set(selected_ids)):
        raise ValueError("conformers must have unique IDs.")
    if not set(selected_ids).issubset(ensemble_ids):
        raise ValueError("conformers must belong to the supplied ensemble.")

    molecule = ensemble.molecule
    conformer_ids = [conformer.id for conformer in selected_conformers]
    n_conformers = len(conformer_ids)
    rmsd = np.zeros((n_conformers, n_conformers), dtype=float)

    for row_index, conformer_id in enumerate(conformer_ids):
        for column_index in range(row_index + 1, n_conformers):
            rmsd[row_index, column_index] = rdMolAlign.GetBestRMS(
                molecule,
                molecule,
                prbId=conformer_id,
                refId=conformer_ids[column_index],
            )
            rmsd[column_index, row_index] = rmsd[row_index, column_index]
    return rmsd, conformer_ids

def rmsd_heatmap(
    ensemble: Ensemble,
    conformers: Sequence[Conformer] | None = None,
    *,
    reference_conformer: Conformer | None = None,
    show: bool = True,
) -> tuple[np.ndarray, Axes]:
    """Plot and return the pairwise RMSD matrix for an ensemble.

    Parameters
    ----------
    ensemble
        Ensemble whose RDKit molecule contains the conformer coordinates.
    conformers
        Optional ordered subset of conformers to include. The default uses
        ``ensemble.conformers``.
    reference_conformer
        Optional conformer used as the first matrix row and column. It must be
        included in ``conformers``; the pairwise values remain symmetric.
    show
        Display the figure immediately when ``True``.

    Returns
    -------
    tuple[numpy.ndarray, matplotlib.axes.Axes]
        The RMSD matrix and the heatmap axes. Matrix labels use conformer IDs.
    """
    selected_conformers = list(
        ensemble.conformers if conformers is None else conformers
    )
    if not selected_conformers:
        raise ValueError("ensemble must contain at least one conformer.")

    ensemble_ids = set(ensemble.conformer_ids)
    selected_ids = [conformer.id for conformer in selected_conformers]
    if len(selected_ids) != len(set(selected_ids)):
        raise ValueError("conformers must have unique IDs.")
    if not set(selected_ids).issubset(ensemble_ids):
        raise ValueError("conformers must belong to the supplied ensemble.")

    if reference_conformer is not None:
        if reference_conformer.id not in selected_ids:
            raise ValueError("reference_conformer must be included in conformers.")
        selected_conformers.insert(
            0,
            selected_conformers.pop(selected_ids.index(reference_conformer.id)),
        )
        selected_ids.insert(0, selected_ids.pop(selected_ids.index(reference_conformer.id)))

    molecule = ensemble.molecule
    conformer_ids = [conformer.id for conformer in selected_conformers]
    n_conformers = len(conformer_ids)
    rmsd_matrix = np.zeros((n_conformers, n_conformers), dtype=float)

    for row_index, conformer_id in enumerate(conformer_ids):
        for column_index in range(row_index + 1, n_conformers):
            rmsd_matrix[row_index, column_index] = rdMolAlign.GetBestRMS(
                molecule,
                molecule,
                prbId=conformer_id,
                refId=conformer_ids[column_index],
            )
            rmsd_matrix[column_index, row_index] = rmsd_matrix[row_index, column_index]

    figure, axes = plt.subplots(figsize=(10, 8))
    image = axes.imshow(rmsd_matrix, cmap="magma", vmin=0)
    figure.colorbar(image, ax=axes, label="RMSD (Å)")
    axes.set_xticks(range(n_conformers), labels=conformer_ids)
    axes.set_yticks(range(n_conformers), labels=conformer_ids)
    axes.set_xlabel("Conformer ID")
    axes.set_ylabel("Conformer ID")
    axes.set_title("Pairwise RMSD Heatmap")
    figure.tight_layout()

    if show:
        plt.show()

    return rmsd_matrix, axes


rmsd_headmap = rmsd_heatmap