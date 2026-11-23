# Filter

## Overview
Filters select or remove conformers from an `Ensemble` and return a new `Ensemble`; the input is not modified. Each filter preserves the molecule and SMILES, recording the input/output conformer counts and filter settings in `metadata["history"]`. Filters that use energies require energies to be present on input Conformers (via optimization).

## Functions
- **EnergyFilter:** keeps conformers within an energy window of the lowest-energy conformer, or keeps the lowest-energy `top_n` conformers.
- **PopulationFilter:** keeps conformers by individual or cumulative Boltzmann population.
- **RMSDFilter:** removes geometrically redundant conformers using an RMSD cutoff, prioritizing lower-energy conformers.


### `EnergyFilter`
**Overview**
Filters conformers by computed energy. Specify exactly one of `window` or `top_n`. With `window`, all conformers whose energy is within the specified kcal/mol window of the minimum energy are retained. With `top_n`, conformers are sorted by energy and the lowest-energy entries are retained. Conformers without an energy are excluded.

Metadata recorded: `window` and `top_n` (the unused criterion is omitted)


**Usage**
```python
from ensemblelab.filters.energy import EnergyFilter

energy_filter = EnergyFilter(window=5.0)
filtered = energy_filter.apply(ensemble)
```


### `PopulationFilter`
**Overview**
Filters conformers using normalized Boltzmann populations. Specify exactly one of `population_cutoff` or `cumulative_cutoff`; both values are fractions between 0 and 1. `population_cutoff` retains every conformer whose individual population meets the cutoff. `cumulative_cutoff` ranks conformers by population and retains entries until the cumulative population reaches the cutoff. Every conformer must have an energy.

The default temperature is 298.15 K. Populations are calculated using the kcal/mol gas constant and the supplied temperature.

Metadata recorded: `temperature`, `population_cutoff`, `cumulative_cutoff`


**Usage**
```python
from ensemblelab.filters.population import PopulationFilter

population_filter = PopulationFilter(
	cumulative_cutoff=0.95,
	temperature=298.15,
)
filtered = population_filter.apply(ensemble)
```


### `RMSDFilter`
**Overview**
Removes conformers that are within the RMSD cutoff of an already approved conformer. Conformers are sorted from lowest to highest energy, so the lowest energy representative is considered first. Conformers with equal energies keep their input order. A conformer is a duplicate if its RMSD to any approved conformer is `<= cutoff`. The default cutoff is 0.5 Angstroms. Every conformer must have a computed energy.

- **Alignment:** `rdkit.Chem.rdMolAlign.GetBestRMS`, an optimal rigid-body superposition that also searches symmetry-equivalent atom mappings (e.g. permuted methyl hydrogens). All conformers share one RDKit molecule, so atom ordering is identical by construction.
- **Atom selection:** all atoms by default. With `heavy_atoms_only=True`, hydrogens are removed with `Chem.RemoveHs` before computing RMSD. Rotating an OH (or another non-symmetric hydrogen) then no longer inflates the RMSD, so true heavy-atom duplicates merge. `RemoveHs` keeps a few special hydrogens, such as isotopic or stereo-defining ones.
- **Coordinates:** read from the RDKit molecule, which must match each `Conformer.atoms`. A mismatch raises `ValueError` instead of producing silently wrong duplicate decisions. EnsembleLab's optimizers keep both copies synchronized.
- **Immutability:** the input ensemble, including its RDKit coordinates, is not modified. RMSDs are computed on a copy of the molecule, because `GetBestRMS` leaves the probe conformer aligned in place. The same pool can be filtered repeatedly with identical results.
- **Edge cases:** a single conformer is always kept. Missing energies raise `ValueError`.

Metadata recorded: `rmsd_cutoff_angstrom`, `heavy_atoms_only`, `alignment`, `energy_ordering`, `kept_conformer_ids`, `removed_conformer_ids`, and `duplicate_of`. `duplicate_of` maps each removed conformer ID to the kept conformer it matched (`kept_conformer_id`) and their `rmsd_angstrom`.


**Usage**
```python
from ensemblelab.filters.rmsd import RMSDFilter

rmsd_filter = RMSDFilter(cutoff=0.5)
filtered = rmsd_filter.apply(ensemble)

# Ignore hydrogen positions when deciding duplicates
heavy_filter = RMSDFilter(cutoff=0.5, heavy_atoms_only=True)
filtered = heavy_filter.apply(ensemble)
filtered.metadata["history"][-1]["duplicate_of"]
```
