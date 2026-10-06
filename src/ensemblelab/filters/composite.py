from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass
from typing import TYPE_CHECKING

from ensemblelab.filters.base import BaseFilter

if TYPE_CHECKING:
    from ensemblelab.generators import Ensemble


@dataclass(frozen=True, slots=True)
class CompositeFilter(BaseFilter):
    """Apply multiple filters sequentially to an ensemble."""

    filters: Sequence[BaseFilter]

    def __post_init__(self) -> None:
        if not self.filters:
            raise ValueError("CompositeFilter requires at least one filter.")

        for filter_ in self.filters:
            if not isinstance(filter_, BaseFilter):
                raise TypeError(
                    "All items in filters must be BaseFilter instances."
                )

    @property
    def method(self) -> str:
        """Canonical name recorded for the composite filter."""
        return "CompositeFilter"

    def apply(self, ensemble: Ensemble) -> Ensemble:
        """Apply each filter sequentially and return the final ensemble."""
        self._validate_ensemble(ensemble)

        result = ensemble

        for filter_ in self.filters:
            result = filter_.apply(result)

        return result

    def _history_details(self) -> dict:
        """Return the filters included in the composite operation."""
        return {
            "filters": [filter_.method for filter_ in self.filters],
        }