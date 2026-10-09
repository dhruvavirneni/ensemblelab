"""Shannon entropy utilities for normalized probability distributions."""

from __future__ import annotations

import math
from collections.abc import Sequence
from numbers import Real
from typing import Any


def entropy_statistics(probabilities: Sequence[float]) -> dict[str, Any]:
	"""Calculate Shannon entropy and effective size from normalized probabilities.

	Zero probabilities are valid and contribute zero under the convention
	``0 * ln(0) = 0``. Non-finite, negative, boolean, or non-normalized values
	raise ``ValueError``. An empty sequence returns missing statistics.

	Returns ``n_valid_populations``, ``population_sum``, ``entropy`` (nats), and
	``effective_size``. The caller decides which source values are valid and
	must normalize them before calling this function.
	"""
	if len(probabilities) == 0:
		return {
			"n_valid_populations": 0,
			"population_sum": 0.0,
			"entropy": None,
			"effective_size": None,
		}

	values: list[float] = []
	for probability in probabilities:
		if (
			isinstance(probability, bool)
			or not isinstance(probability, Real)
			or not math.isfinite(float(probability))
			or probability < 0
		):
			raise ValueError("probabilities must be finite, non-negative real values.")
		values.append(float(probability))

	population_sum = math.fsum(values)
	if not math.isclose(population_sum, 1.0, rel_tol=1e-9, abs_tol=1e-12):
		raise ValueError("probabilities must be normalized to sum to 1.")

	entropy = -math.fsum(
		probability * math.log(probability)
		for probability in values
		if probability > 0.0
	)
	return {
		"n_valid_populations": len(values),
		"population_sum": population_sum,
		"entropy": entropy,
		"effective_size": math.exp(entropy),
	}