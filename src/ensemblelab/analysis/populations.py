"""Energy-derived conformer populations for molecular ensembles."""

from __future__ import annotations

import math
from collections.abc import Sequence
from numbers import Real
from typing import Any

from ensemblelab.generators import Conformer, Ensemble

from .energies import ENERGY_UNITS, analyze_energies

BOLTZMANN_CONSTANT_KCAL_MOL_K = 0.00198720425864083
DEFAULT_TEMPERATURE_K = 298.15


def _boltzmann_weights(
	relative_energies: Sequence[float], temperature: float
) -> list[float]:
	thermal_energy = BOLTZMANN_CONSTANT_KCAL_MOL_K * temperature
	if thermal_energy == 0.0:
		return [1.0 if energy == 0.0 else 0.0 for energy in relative_energies]
	return [math.exp(-energy / thermal_energy) for energy in relative_energies]


def analyze_populations(
	ensemble: Ensemble,
	temperature: float = DEFAULT_TEMPERATURE_K,
	*,
	conformers: Sequence[Conformer] | None = None,
) -> dict[str, Any]:
	"""Calculate normalized Boltzmann weights and populations for conformers.

	Parameters
	----------
	ensemble
		Source molecular ensemble.
	temperature
		Positive, finite temperature in kelvin. Defaults to the existing
		EnsembleLab convention of 298.15 K; the temperature used is always
		recorded in the returned summary.
	conformers
		Optional ordered subset of conformers owned by ``ensemble``. Selection
		and source indices follow :func:`analyze_energies`.

	Returns
	-------
	dict
		A mapping with ordered ``conformers`` records and a ``summary``. Invalid
		energies retain their energy fields, have ``population_valid=False`` and
		``None`` weights/populations, and are excluded from normalization. When
		there are no valid energies, population statistics are ``None``.

	Raises
	------
	ValueError
		If temperature is not a positive finite real number, or if a finite
		energy has missing or unsupported units. Energy validation and relative
		energies are provided by :func:`analyze_energies`.

	Notes
	-----
	Weights are calculated from minimum-shifted energies in kcal/mol. The
	calculation assumes equal statistical degeneracy. These are energy-based
	 estimates, not necessarily solution populations: optimized energies alone
	 do not account for conformational entropy, vibrational contributions,
	 solvent effects, sampling completeness, or duplicate representations of a
	 physical state. Entropy uses natural logarithms and is reported in nats;
	 effective ensemble size is its exponential and is conformer-level, not a
	 count of distinct structural states.
	"""
	if (
		isinstance(temperature, bool)
		or not isinstance(temperature, Real)
		or not math.isfinite(float(temperature))
		or temperature <= 0
	):
		raise ValueError("temperature must be a positive finite value in kelvin.")
	temperature = float(temperature)

	energy_result = analyze_energies(ensemble, conformers=conformers)
	rows = [
		{
			**record,
			"boltzmann_weight": None,
			"population": None,
			"population_valid": False,
		}
		for record in energy_result["conformers"]
	]
	valid_positions = [
		index for index, record in enumerate(energy_result["conformers"])
		if record["energy_valid"]
	]
	if valid_positions:
		relative_energies = [
			energy_result["conformers"][index]["relative_energy"]
			for index in valid_positions
		]
		weights = _boltzmann_weights(relative_energies, temperature)
		weight_sum = math.fsum(weights)
		populations = [weight / weight_sum for weight in weights]
		for index, weight, population in zip(valid_positions, weights, populations):
			rows[index]["boltzmann_weight"] = weight
			rows[index]["population"] = population
			rows[index]["population_valid"] = True
		population_sum: float | None = math.fsum(populations)
		max_population: float | None = max(populations)
		population_entropy: float | None = -math.fsum(
			population * math.log(population)
			for population in populations
			if population > 0.0
		)
		effective_ensemble_size: float | None = math.exp(population_entropy)
	else:
		population_sum = max_population = None
		population_entropy = effective_ensemble_size = None

	energy_summary = energy_result["summary"]
	return {
		"conformers": rows,
		"summary": {
			"temperature_K": temperature,
			"n_conformers": energy_summary["n_conformers"],
			"n_valid_energies": energy_summary["n_valid_energies"],
			"n_invalid_energies": energy_summary["n_invalid_energies"],
			"population_sum": population_sum,
			"max_population": max_population,
			"effective_ensemble_size": effective_ensemble_size,
			"population_entropy": population_entropy,
			"population_entropy_units": "nats",
			"energy_units": ENERGY_UNITS,
			"boltzmann_constant_kcal_mol_K": BOLTZMANN_CONSTANT_KCAL_MOL_K,
		},
	}