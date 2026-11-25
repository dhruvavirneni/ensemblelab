# Generate

## Overview
`generate()` creates a new unoptimized `Ensemble` from a SMILES string using RDKit ETKDGv3 with a fixed random seed of 42. Each generated conformer is stored as a `Conformer` with an ASE `Atoms` geometry and optional energy metadata that is populated later by optimization steps. The returned ensemble retains the original SMILES, the RDKit molecule object, conformer IDs, and chronological workflow provenance in `metadata["history"]`.

## Functions
- **Conformer:** immutable per-conformer object with geometry, energy, and optimization metadata.
- **Ensemble:** aligned conformer collection with the original molecule and provenance.
- **generate():** generates a new ensemble from a SMILES string.
- **Ensemble.from_smiles():** convenience constructor equivalent to `generate()`.


### `Conformer`
**Overview:** 
A single generated conformer and its conformer-level data. The `energy` field is `None` until an optimization backend assigns a value. The geometry is stored as an ASE `Atoms` object, while `energy_unit`, `optimization_method`, and `optimization_converged` are attached to the same conformer rather than to a separate ensemble-level array.

Metadata recorded: `id`, `atoms`, `energy`, `energy_unit`, `optimization_method`, `optimization_converged`


**Usage**
```python
from ase import Atoms
from ensemblelab import Conformer

conformer = Conformer(
    id=0,
    atoms=Atoms("H2O", positions=[[0.0, 0.0, 0.0], [0.0, 0.0, 0.9], [0.0, 0.9, 0.0]]),
    energy=None,
    energy_unit=None,
    optimization_method=None,
    optimization_converged=None,
)
```


### `Ensemble`
**Overview:** 
An aligned conformational ensemble containing the input SMILES, the RDKit molecule, and a tuple of `Conformer` objects. The ensemble enforces unique conformer IDs and preserves one explicit conformer index for all downstream analyses. The molecular geometry and energies remain synchronized with the conformer-level metadata.

`Ensemble.charge` (net formal charge) and `Ensemble.multiplicity` (1 + radical electrons, i.e. radicals are assumed high-spin) are read-only properties computed from `molecule` on every access, so they cannot disagree with it. Optimizers use them as their default charge and multiplicity.

Metadata recorded: `n_conformers`, `optimization_status`, `energy_status`, `energy_unit`, `rdkit_version`, `canonical_smiles`, `charge`, `multiplicity`, `history`


**Usage**
```python
from ensemblelab import Ensemble

ensemble = Ensemble(
    smiles="CCO",
    molecule=molecule,
    conformers=(conformer_1, conformer_2),
    metadata={
        "n_conformers": 2,
        "history": [],
    },
)
```


### `generate`
**Overview:** 
Creates an unoptimized ensemble by parsing a SMILES string, adding explicit hydrogens, and embedding multiple conformers with RDKit ETKDGv3. The default random seed is fixed at 42 for reproducibility, and the ensemble stores the canonicalized molecule and generation provenance in `metadata["history"]`. Validation rejects empty or invalid SMILES strings, rejects SMILES with more than one disconnected fragment (see [Ions](#ions)), and requires a positive integer value for `n_confs`. The net charge and multiplicity are recorded in the metadata and in the generation history entry.

Metadata recorded: `n_conformers`, `optimization_status`, `energy_status`, `energy_unit`, `rdkit_version`, `canonical_smiles`, `charge`, `multiplicity`, `history`


**Usage**
```python
from ensemblelab import generate

ensemble = generate("CCO", n_confs=25)
```


### `Ensemble.from_smiles`
**Overview:** 
Class method equivalent to calling `generate()`. It is the object-oriented entry point for creating an unoptimized `Ensemble` from a SMILES string and is intended to fit naturally into the broader `Ensemble` workflow for optimization, filtering, and analysis.

Metadata recorded: `n_conformers`, `optimization_status`, `energy_status`, `energy_unit`, `rdkit_version`, `canonical_smiles`, `charge`, `multiplicity`, `history`


**Usage**
```python
from ensemblelab import Ensemble

ensemble = Ensemble.from_smiles("CCO", n_confs=25)
```


## Ions

Charged molecules are supported. Write the charge directly into the SMILES with bracket atoms:

| Species | SMILES | `ensemble.charge` |
|---|---|---|
| Acetate (deprotonated acid, [M-H]-) | `CC(=O)[O-]` | -1 |
| Perfluorooctanoate (PFOA anion) | `[O-]C(=O)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F` | -1 |
| Propylammonium (protonated amine, [M+H]+) | `CCC[NH3+]` | +1 |

```python
from ensemblelab import generate
from ensemblelab.optimizers import MMFFOptimizer

ensemble = generate("CC(=O)[O-]", n_confs=25)
ensemble.charge                     # -1, computed from the RDKit molecule
optimized = MMFFOptimizer().optimize(ensemble)
optimized.show()                    # summary includes "Charge: -1"
```

**What is supported**
- `generate()` / `Ensemble.from_smiles()` record `charge` and `multiplicity` in `metadata` and in the first history entry.
- `MMFFOptimizer` / `UFFOptimizer` use the molecule's formal charges and record the charge in their history entry. `GFN2xTBOptimizer` and `ORCAOptimizer` pass `ensemble.charge` and `ensemble.multiplicity` to the calculator unless you supply `charge=` / `multiplicity=` explicitly.
- Filters, RMSD analysis, and geometry metrics work unchanged. Formal charges do not change atomic masses, so mass-weighted quantities (center of mass, radius of gyration) are unaffected. `RMSDFilter` treats the two carboxylate oxygens as equivalent: RDKit's `GetBestRMS` symmetrizes conjugated terminal groups, even though the SMILES writes one oxygen as `=O` and the other as `[O-]`.

**Caveat: compare energies only within one charge state.** An anion and its neutral parent are different chemical species, so the difference between their energies is not a conformational energy. Relative energies, Boltzmann populations, and energy-based filters are meaningful only among conformers of the same charged species. Do not compare anion energies to neutral energies.

**Out of scope for now**
- *Multi-fragment input, including metal adducts* such as [M+Na]+. In SMILES these are separate fragments (`CC(=O)O.[Na+]`). ETKDG places the fragments at arbitrary distances from each other, so the conformers would be physically meaningless. `generate()` raises `ValueError` for any SMILES with `.`-separated fragments (salts, adducts, solvates). Handling metal adducts properly is a much harder problem.
- *Automatic ionization.* You supply the charged SMILES yourself. There is no helper that turns a neutral molecule into [M-H]- or [M+H]+, because that requires choosing which proton to remove or add (the protomer question).
