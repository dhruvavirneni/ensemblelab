from __future__ import annotations

import math
from collections.abc import Sequence
from dataclasses import dataclass
from typing import Any

import numpy as np
from rdkit import Chem
from rdkit.Chem import rdMolAlign

from ensemblelab.analysis.entropy import entropy_statistics
from ensemblelab.analysis.populations import (
	BOLTZMANN_CONSTANT_KCAL_MOL_K,
	_boltzmann_weights,
)
from ensemblelab.generators import Conformer, Ensemble

R_KCAL = BOLTZMANN_CONSTANT_KCAL_MOL_K


@dataclass(frozen=True, slots=True)
class Basin:
	"""An RMSD-defined basin and its conformer-level summary."""

	id: int
	representative: Conformer
	members: tuple[Conformer, ...]
	population: float | None
	min_energy: float | None

	@property
	def representative_id(self) -> int:
		return self.representative.id

	@property
	def member_ids(self) -> tuple[int, ...]:
		return tuple(member.id for member in self.members)

	@property
	def size(self) -> int:
		return len(self.members)


def _select_conformers(
	ensemble: Ensemble, conformers: Sequence[Conformer] | None
) -> list[Conformer]:
	selected = list(ensemble.conformers if conformers is None else conformers)
	if not selected:
		raise ValueError("At least one conformer is required for basin analysis.")

	ensemble_conformers = {conformer.id: conformer for conformer in ensemble.conformers}
	selected_ids = [conformer.id for conformer in selected]
	if len(selected_ids) != len(set(selected_ids)):
		raise ValueError("conformers must have unique IDs.")
	if any(ensemble_conformers.get(conformer.id) is not conformer for conformer in selected):
		raise ValueError("conformers must belong to the supplied ensemble.")
	return selected


def _boltzmann_populations(
	conformers: Sequence[Conformer], temperature: float
) -> dict[int, float] | None:
	if any(conformer.energy is None for conformer in conformers):
		return None

	energy_units = {
		conformer.energy_unit.lower().replace(" ", "")
		for conformer in conformers
		if conformer.energy_unit is not None
	}
	if len(energy_units) > 1 or (energy_units and energy_units != {"kcal/mol"}):
		raise ValueError("Boltzmann populations require energies in kcal/mol.")

	energies = [float(conformer.energy) for conformer in conformers]
	if not all(math.isfinite(energy) for energy in energies):
		raise ValueError("conformer energies must be finite.")

	minimum_energy = min(energies)
	weights = _boltzmann_weights(
		[energy - minimum_energy for energy in energies], temperature
	)
	partition_function = sum(weights)
	return {
		conformer.id: weight / partition_function
		for conformer, weight in zip(conformers, weights)
	}


def cluster_basins(
	ensemble: Ensemble,
	conformers: Sequence[Conformer] | None = None,
	*,
	rmsd_threshold: float = 0.75,
	temperature: float = 298.15,
) -> list[Basin]:
	"""Group conformers by best-fit RMSD to each basin's first member.

	Basin representatives follow the input conformer order, matching the
	greedy grouping used in the analysis notebook. Populations are Boltzmann
	weights normalized over the selected conformers when all have energies;
	otherwise populations are ``None``.
	"""
	if not math.isfinite(rmsd_threshold) or rmsd_threshold <= 0:
		raise ValueError("rmsd_threshold must be a finite value greater than 0.")
	if not math.isfinite(temperature) or temperature <= 0:
		raise ValueError("temperature must be a finite value greater than 0 K.")

	selected = _select_conformers(ensemble, conformers)
	populations = _boltzmann_populations(selected, temperature)
	molecule = ensemble.molecule
	grouped: list[tuple[Conformer, list[Conformer]]] = []

	for conformer in selected:
		probe_molecule = Chem.Mol(molecule)
		for representative, members in grouped:
			distance = rdMolAlign.GetBestRMS(
				probe_molecule,
				molecule,
				prbId=conformer.id,
				refId=representative.id,
			)
			if distance < rmsd_threshold:
				members.append(conformer)
				break
		else:
			grouped.append((conformer, [conformer]))

	basins = []
	for basin_id, (representative, members) in enumerate(grouped, start=1):
		member_energies = [
			float(member.energy) for member in members if member.energy is not None
		]
		basins.append(
			Basin(
				id=basin_id,
				representative=representative,
				members=tuple(members),
				population=(
					sum(populations[member.id] for member in members)
					if populations is not None
					else None
				),
				min_energy=min(member_energies) if member_energies else None,
			)
		)
	return basins


def basin_rmsd_matrix(
	ensemble: Ensemble, basins: Sequence[Basin]
) -> tuple[np.ndarray, tuple[int, ...], tuple[int, ...]]:
	"""Return an RMSD matrix ordered by basin, plus conformer and basin IDs."""
	if not basins:
		raise ValueError("At least one basin is required.")
	basin_ids = [basin.id for basin in basins]
	if len(basin_ids) != len(set(basin_ids)):
		raise ValueError("basin IDs must be unique.")

	ordered_members = [member for basin in basins for member in basin.members]
	selected = _select_conformers(ensemble, ordered_members)
	conformer_ids = tuple(conformer.id for conformer in selected)
	ordered_basin_ids = tuple(
		basin.id for basin in basins for _ in basin.members
	)
	matrix = np.zeros((len(selected), len(selected)), dtype=float)

	for row_index, conformer in enumerate(selected):
		probe_molecule = Chem.Mol(ensemble.molecule)
		for column_index in range(row_index + 1, len(selected)):
			matrix[row_index, column_index] = rdMolAlign.GetBestRMS(
				probe_molecule,
				ensemble.molecule,
				prbId=conformer.id,
				refId=selected[column_index].id,
			)
			matrix[column_index, row_index] = matrix[row_index, column_index]

	return matrix, conformer_ids, ordered_basin_ids


def basin_entropy(basins: Sequence[Basin]) -> float | None:
	"""Calculate Shannon entropy (nats) from basin populations."""
	if not basins:
		raise ValueError("At least one basin is required.")
	if any(basin.population is None for basin in basins):
		return None

	populations = [float(basin.population) for basin in basins]
	if any(not math.isfinite(population) or population < 0 for population in populations):
		raise ValueError("basin populations must be finite and non-negative.")
	total_population = sum(populations)
	if total_population <= 0:
		raise ValueError("basin populations must sum to a positive value.")

	normalized_populations = [population / total_population for population in populations]
	return entropy_statistics(normalized_populations)["entropy"]


def plot_basin_analysis(
	ensemble: Ensemble, basins: Sequence[Basin], *, show: bool = True
) -> tuple[Any, tuple[Any, Any]]:
	"""Plot basin occupancy and the RMSD matrix ordered by basin.

	Matplotlib is imported only when plotting is requested. Occupancy uses
	population percentages when available, or basin sizes for unoptimized data.
	"""
	import matplotlib.pyplot as plt

	matrix, conformer_ids, _ = basin_rmsd_matrix(ensemble, basins)
	colors = plt.get_cmap("tab20")
	basin_colors = [colors(index % colors.N) for index in range(len(basins))]
	has_populations = all(basin.population is not None for basin in basins)
	occupancy = [
		basin.population * 100 if has_populations else basin.size
		for basin in basins
	]
	occupancy_label = "Population (%)" if has_populations else "Conformers"

	figure, (occupancy_axes, heatmap_axes) = plt.subplots(
		2,
		1,
		figsize=(max(9, len(basins) * 0.55), 9),
		gridspec_kw={"height_ratios": [1, 5]},
		constrained_layout=True,
	)
	positions = np.arange(len(basins))
	occupancy_axes.bar(positions, occupancy, color=basin_colors, edgecolor="black")
	occupancy_axes.set_xticks(positions, labels=[str(basin.id) for basin in basins])
	occupancy_axes.set_xlabel("Basin ID")
	occupancy_axes.set_ylabel(occupancy_label)
	occupancy_axes.set_title("Basin Occupancy")

	vmax = max(float(np.max(matrix)), 1e-12)
	image = heatmap_axes.imshow(matrix, cmap="magma", vmin=0, vmax=vmax)
	boundaries = np.cumsum([basin.size for basin in basins])[:-1]
	for boundary in boundaries:
		heatmap_axes.axhline(boundary - 0.5, color="white", linewidth=1.2)
		heatmap_axes.axvline(boundary - 0.5, color="white", linewidth=1.2)

	if len(conformer_ids) <= 40:
		tick_step = 1 if len(conformer_ids) <= 24 else 2
		ticks = np.arange(0, len(conformer_ids), tick_step)
		heatmap_axes.set_xticks(ticks, labels=[str(conformer_ids[i]) for i in ticks])
		heatmap_axes.set_yticks(ticks, labels=[str(conformer_ids[i]) for i in ticks])
	else:
		heatmap_axes.set_xticks([])
		heatmap_axes.set_yticks([])

	heatmap_axes.set_xlabel("Conformer ID")
	heatmap_axes.set_ylabel("Conformer ID")
	heatmap_axes.set_title("Pairwise RMSD Reordered by Basin")
	figure.colorbar(image, ax=heatmap_axes, label="RMSD (Å)")
	if show:
		plt.show()
	return figure, (occupancy_axes, heatmap_axes)
