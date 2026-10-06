"""Composite optimizer workflow orchestration."""

from __future__ import annotations

from collections.abc import Sequence
from typing import Any

from ..generators import Ensemble
from .base import BaseOptimizer


class HierarchicalOptimizer(BaseOptimizer):
    """Apply multiple optimizers sequentially to a conformational ensemble."""

    def __init__(self, optimizers: Sequence[BaseOptimizer]):
        if not optimizers:
            raise ValueError("optimizers cannot be empty.")
        if not all(isinstance(optimizer, BaseOptimizer) for optimizer in optimizers):
            raise TypeError("optimizers must contain BaseOptimizer instances.")

        self.optimizers = tuple(optimizers)

    @property
    def method(self) -> str:
        return "hierarchical"

    def _optimize_conformer(self, molecule, conformer_id, atoms):
        raise NotImplementedError(
            "HierarchicalOptimizer does not optimize conformers directly."
        )

    def _history_details(self) -> dict[str, Any]:
        return {}

    def optimize(self, ensemble: Ensemble) -> Ensemble:
        current = ensemble

        for optimizer in self.optimizers:
            current = optimizer.optimize(current)

        return current
