"""Section-based, read-only text renderers for ensemble objects."""

from __future__ import annotations

from collections.abc import Mapping
from typing import Any

from .formatters import format_energy, format_value
from .tables import conformer_table, metadata_table


def conformer_summary(conformer: Any) -> str:
    """Return the stored summary for one conformer."""
    lines = [
        f"Conformer {conformer.id}",
        "=" * 32,
        f"Energy          {format_energy(conformer.energy, conformer.energy_unit)}",
        f"Optimization    {format_value(conformer.optimization_method)}",
        f"Converged       {format_value(conformer.optimization_converged) if conformer.optimization_converged is None else str(conformer.optimization_converged)}",
        f"Atoms           {len(conformer.atoms)}",
    ]
    return "\n".join(lines)


def ensemble_summary(ensemble: Any, *, include_conformers: bool = True) -> str:
    """Return the default ensemble section and, optionally, its conformer table."""
    conformers = ensemble.conformers
    energies = [conformer.energy for conformer in conformers]
    if all(energy is not None for energy in energies):
        energy_status = "computed"
    elif any(energy is not None for energy in energies):
        energy_status = "partially computed"
    else:
        energy_status = "uncomputed"
    methods = {conformer.optimization_method for conformer in conformers if conformer.optimization_method}
    optimization = next(iter(methods)) if len(methods) == 1 else ("mixed" if methods else "not run")

    section = "\n".join(
        [
            "Ensemble",
            "-" * 28,
            f"SMILES: {ensemble.smiles}",
            f"Atoms: {ensemble.molecule.GetNumAtoms()}",
            f"Charge: {ensemble.charge:+d}" if ensemble.charge else "Charge: 0",
            f"Multiplicity: {ensemble.multiplicity}",
            f"Conformers: {len(conformers)}",
            f"Energy: {energy_status}",
            f"Optimization: {optimization}",
        ]
    )
    if not include_conformers:
        return section
    return f"{section}\n\nConformers\n{'-' * 36}\n{conformer_table(conformers)}"


def ensemble_history(ensemble: Any) -> str:
    """Render provenance records from all currently supported history keys."""
    entries = _history_entries(ensemble.metadata)
    lines = ["Workflow History", "=" * 32]
    if not entries:
        return "\n".join(lines + ["No workflow history recorded."])
    for index, entry in enumerate(entries, start=1):
        lines.extend(["" if index > 1 else "", _history_entry(index, entry)])
    return "\n".join(lines).rstrip()


def ensemble_metadata(ensemble: Any) -> str:
    """Render raw ensemble metadata in a separate inspection section."""
    return f"Metadata\n{'=' * 32}\n{metadata_table(ensemble.metadata)}"


def _history_entries(metadata: Mapping[str, Any]) -> list[Mapping[str, Any]]:
    """Return canonical workflow history without mutating metadata."""
    value = metadata.get("history", [])
    return [entry for entry in value if isinstance(entry, Mapping)] if isinstance(value, list) else []


def _history_entry(index: int, entry: Mapping[str, Any]) -> str:
    """Render one heterogeneous provenance event using relevant populated fields."""
    process = format_value(entry.get("process", "workflow"))
    lines = [f"{index}. {process}"]
    fields = (
        ("method", "Method"),
        ("n_requested", "Requested"),
        ("n_generated", "Generated"),
        ("n_input_conformers", "Input"),
        ("n_output_conformers", "Output"),
        ("n_converged", "Converged"),
        ("n_unconverged", "Unconverged"),
        ("window", "Window"),
        ("top_n", "Top n"),
        ("rmsd_cutoff_angstrom", "RMSD cutoff (angstrom)"),
        ("temperature", "Temperature (K)"),
        ("population_cutoff", "Population cutoff"),
        ("cumulative_cutoff", "Cumulative cutoff"),
        ("max_steps", "Max steps"),
        ("fmax_eV_per_angstrom", "fmax (eV/angstrom)"),
        ("charge", "Charge"),
        ("multiplicity", "Multiplicity"),
        ("solvent", "Solvent"),
        ("energy_unit", "Energy unit"),
        ("random_seed", "Random seed"),
    )
    for key, label in fields:
        if entry.get(key) is not None:
            lines.append(f"   {label}: {format_value(entry[key])}")
    return "\n".join(lines)
