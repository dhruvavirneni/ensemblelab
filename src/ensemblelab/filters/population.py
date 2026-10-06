from __future__ import annotations

# population filtration via boltzmann distribution
import math
from dataclasses import dataclass
from typing import TYPE_CHECKING

from ensemblelab.filters.base import BaseFilter

if TYPE_CHECKING:
    from ensemblelab.generators import Ensemble

"""
Boltzmann population filter


Two filter classes in Boltzmann population filter:

Individual population cutoff - return Conformers in Ensemble with minimum Boltzmann population % threshold
Cumulative population cutoff - sort Conformers in Ensemble by energy and maintain conformers (high to low population) until collective % population reaches inputted threshold

"""

R_KCAL = 0.00198720425864083


@dataclass(frozen=True, slots=True)
class PopulationFilter(BaseFilter):
    """Filter conformers using Boltzmann populations."""

    population_cutoff: float | None = None
    cumulative_cutoff: float | None = None
    temperature: float = 298.15

    @property
    def method(self) -> str:
        return "boltzmann_population"

    def __post_init__(self) -> None:
        # error handling for population filter type selection
        if (self.population_cutoff is None) == (self.cumulative_cutoff is None):
            raise ValueError("Specify exactly one of population_cutoff or cumulative_cutoff.")

        if self.temperature <= 0:
            raise ValueError("temperature must be greater than 0 K.")

        if self.population_cutoff is not None and not (0 < self.population_cutoff <= 1):
            raise ValueError("population_cutoff must be between 0 and 1.")

        if self.cumulative_cutoff is not None and not (0 < self.cumulative_cutoff <= 1):
            raise ValueError("cumulative_cutoff must be between 0 and 1.")

    def apply(self, ensemble: Ensemble) -> Ensemble:
        self._validate_ensemble(ensemble)
        if any(
            conformer.energy is None
            for conformer in ensemble.conformers):
            raise ValueError("Boltzmann filtering requires every conformer to have an energy.")

        energies = [float(conformer.energy) for conformer in ensemble.conformers]

        min_energy = min(energies)

        populations = self._calculate_populations(
            energies,
            min_energy,)

        # population cutoff filtering
        if self.population_cutoff is not None:
            filtered_conformers = [conformer for conformer, population in zip(ensemble.conformers, populations,) if population >= self.population_cutoff]

        # cumulative population cutoff
        else:
            # conformers sorted by energy
            ranked = sorted(zip(ensemble.conformers, populations), key=lambda item: item[1], reverse=True)
            filtered_conformers = []
            cumulative_population = 0

            for conformer, population in ranked:
                filtered_conformers.append(conformer)
                cumulative_population += population

                if cumulative_population >= self.cumulative_cutoff:
                    break
        return self._build_filtered_ensemble(
            ensemble,
            filtered_conformers,)
            

    def _calculate_populations(self, energies: list[float], min_energy: float) -> list[float]:
        # calculate normalized Boltzmann populations with equation: % population = (e^(-ΔEi/(RT)))/(∑j e^(-ΔEj/(RT)))
        weights = [
            math.exp(-(energy - min_energy) / (R_KCAL * self.temperature))
            for energy in energies
        ]

        # partition function (∑j e^(-ΔEj/(RT)))
        partition_function = sum(weights)

        return [weight / partition_function for weight in weights]

    def _history_details(self) -> dict:
        return {
            "temperature": self.temperature,
            "population_cutoff": self.population_cutoff,
            "cumulative_cutoff": self.cumulative_cutoff,
        }
