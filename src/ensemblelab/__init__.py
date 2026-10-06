"""Tools for reproducible molecular conformational-ensemble analysis."""

from .filters import (
    BaseFilter,
    CompositeFilter,
    EnergyFilter,
    PopulationFilter,
    RMSDFilter,
)
from .generators import Conformer, Ensemble, generate

__all__ = [
    "BaseFilter",
    "BaseOptimizer",
    "CompositeFilter",
    "Conformer",
    "EnergyFilter",
    "Ensemble",
    "GFN2xTBOptimizer",
    "HierarchicalOptimizer",
    "MMFFOptimizer",
    "ORCAOptimizer",
    "PopulationFilter",
    "RMSDFilter",
    "generate",
]


def __getattr__(name: str):
    """Lazily import optional optimizer classes only when explicitly requested."""
    if name in {"BaseOptimizer", "GFN2xTBOptimizer", "HierarchicalOptimizer", "MMFFOptimizer", "ORCAOptimizer"}:
        from . import optimizers

        return getattr(optimizers, name)
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")


def __dir__() -> list[str]:
    return sorted(set(globals()) | set(__all__))
