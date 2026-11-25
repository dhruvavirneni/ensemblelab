# ensemblelab contributor guide

This is the single agent/contributor guide for the repository. Part 1 covers the library as a whole; Part 2 covers the current research focus.

## Current focus

**Work exclusively on the conformational sampling convergence project (Part 2) until instructed otherwise.** Library-wide rules in Part 1 still apply to that work.

## Agent workflow

- Before each task, read this file and inspect relevant code, tests, and interfaces. Inspect the actual repository before assuming filenames, schemas, or APIs.
- Implement one module or well-defined pipeline component at a time. Avoid unrelated refactors, unnecessary dependencies, large rewrites, and unnecessary file moves; preserve existing functionality.
- Use focused tests, verify scientific calculations, and run relevant tests after changes.
- Ask before making major scientific or architectural changes.
- At the end of each task, report changes made, tests and outcomes, unresolved issues, and the next recommended step.
- **Update this file whenever a significant research decision, metric definition, convergence criterion, data schema, interface, completed phase, or limitation changes.** Keep it concise and factual. Reflect actual code and verified results; mark unfinished work clearly. Remove superseded plans rather than accumulating contradictory instructions.

---

# Part 1 — The library

## Mission

`ensemblelab` is an open-source Python library for computing **ensemble descriptors** that quantitatively describe molecular flexibility from conformational ensembles.

The project thesis is: *cheminformatics tooling makes conformers easy to generate but ensembles hard to analyze.* RDKit, xTB, and quantum-chemistry tools provide structures and energies; `ensemblelab` provides the consistent, reproducible analysis layer after generation.

Plots, reports, solvent comparisons, clustering, and benchmarks are supporting features. The central product is a scientifically sound `Ensemble` object and a robust, documented conformational fingerprint (CFP).

## Product principles

- Design for a short, discoverable user workflow. A useful first interaction should look like:

  ```python
  from ensemblelab import Ensemble

  ens = Ensemble.from_smiles(smiles, n_confs=300)
  ens.optimize(method="gfn2-xtb")
  ens.cluster()
  ens.descriptors
  ```

- Keep an actual object as the source of truth. Reports and figures are derived views, never the only output.
- Prefer safe defaults and clear escape hatches. Infer conformer counts from rotatable bonds when the user does not specify one, but record every chosen setting.
- Make computational provenance explicit: software versions, random seeds, force field or optimizer, solvent model, units, thresholds, and failed conformers.
- Do not silently turn an approximate result into a definitive scientific claim. Surface assumptions, warnings, missing optional dependencies, and non-convergence.
- Avoid web-app concerns. This repository is a reusable Python package first; notebook examples are secondary.

## Project structure

Actual layout (inspect before assuming anything finer-grained):

```text
src/ensemblelab/
  __init__.py        # public API: Ensemble, Conformer, generate, filters; optimizers lazily
  generators.py      # Conformer/Ensemble data model, ETKDG generation, setup()
  optimizers/        # BaseOptimizer, MMFF/UFF, GFN2-xTB, ORCA, HierarchicalOptimizer
  optimize.py        # older function-based optimize(), kept for compatibility
  filters/           # Energy, Population, RMSD, Composite filters
  analysis/          # general-purpose metrics: energies, populations, entropy, rmsd, clusters, torsions, geometry
  convergence/       # convergence-facing APIs (Part 2)
  display/           # read-only text summaries and viewers used by show()
  descriptors.py, io.py, report.py, utils.py   # placeholders (empty)
src/experiments/     # experiment orchestration, separate from library code (Part 2)
tests/  docs/  notebooks/
```

Planned additions not yet present: `descriptors` with a versioned CFP, solvent workflows, plotting/report modules, serialization (`io`), and `exceptions.py`. Only expose stable, documented names from `ensemblelab.__init__`. Keep backend-specific and provisional helpers private until their behavior is settled.

## Data model and invariants

`Ensemble` owns aligned conformer-level data. All public arrays/tables must preserve one explicit conformer index and cannot be reordered independently.

Expected core properties include:

- `conformers`: records containing molecular geometries, conformer identifiers, and their own energy values with explicit unit, method, and provenance;
- `populations`: normalized Boltzmann weights derived from conformer energies, temperature, and unit-aware constants;
- `rmsd`: a symmetric pairwise matrix with documented alignment and atom-selection rules;
- `torsions`: atom quartets, angle values, and fingerprints with circular-angle handling;
- `clusters`: assignments, method, parameters, and representative conformers;
- `descriptors`: named scalar ensemble descriptors plus metadata.

When filtering, deduplicating, or clustering changes membership, create a traceable mapping to original conformer IDs. Never discard failed calculations without recording why.

## Scientific conventions

- Use explicit units at every API boundary. Store and report energy units; convert before Boltzmann calculations. Temperature defaults must be visible and documented.
- RMSD must state whether hydrogens are included, the alignment method, and the atom mapping assumption. Do not compare molecules with incompatible atom ordering without an explicit mapping strategy.
- `rdMolAlign.GetBestRMS` aligns its probe in place. Every RMSD computation must run on a copy of the molecule (`Chem.Mol(...)`) so the shared ensemble molecule stays synchronized with `Conformer.atoms`.
- Treat torsions as circular variables. Linear mean/variance is incorrect near the -180/180 degree boundary. Preserve raw angles and support discrete states for comparison.
- Duplicate removal is a composite decision: energy proximity, geometry/RMSD similarity, and torsion-fingerprint similarity. Each criterion must be configurable and reported.
- Descriptors must define their formula, population weighting, unit, applicable ensemble state, and behavior for one-conformer or failed/empty ensembles.
- Keep deterministic behavior practical: accept a seed for stochastic generation and record it in metadata.
- Validation targets are relative/reproducible behavior first. Do not present force-field/xTB/implicit-solvent output as experimental truth.

## Conformational Fingerprint (CFP)

The initial CFP is an interpretable vector, not a black box. Candidate fields are:

1. Boltzmann entropy
2. effective conformer count
3. largest basin occupancy
4. mean pairwise RMSD
5. pairwise RMSD standard deviation
6. mean torsional variance
7. rotatable-bond count
8. energy spread
9. population Gini coefficient
10. radius of conformational space

Implement each descriptor independently with unit tests before bundling it into `cfp()`. Preserve a versioned descriptor schema so datasets remain comparable as the fingerprint evolves.

## Dependency and backend policy

- RDKit is the primary generation and molecular-graph dependency.
- MMFF and UFF should work without external executables where RDKit supports them.
- xTB and ORCA are optional integrations. Import lazily, detect availability early, and raise actionable installation/configuration errors.
- Keep optional plotting, report, and ML dependencies out of the minimal analysis install when feasible.
- ASE may be used as an interoperability layer where helpful, but the public API should not require users to understand ASE internals.

## Implementation workflow

1. Start with the smallest stable public API and an end-to-end RDKit/MMFF path.
2. Add one scientific capability at a time with a test molecule fixture and a documented definition.
3. Keep methods composable: users may supply coordinates/energies instead of invoking generation or optimization.
4. Put exploratory work in `notebooks/`; promote validated logic into `src/ensemblelab/` with tests.
5. Update docs and examples whenever public behavior changes.

Prefer typed, side-effect-light functions under the object methods. Make expensive work explicit; cached derived values must have clear invalidation when conformers, conformer energies, temperature, or settings change.

## Tests and quality bar

- Add tests for every descriptor formula, circular torsion behavior, Boltzmann normalization, duplicate-removal decisions, and conformer-index alignment.
- Include edge cases: a single conformer, degenerate energies, missing energies, failed optimization, no rotatable bonds, and empty post-filtered output.
- Use small deterministic molecules in unit tests; keep xTB/ORCA integration tests opt-in and skip cleanly when executables are unavailable.
- Test public examples as smoke tests once packaging exists.
- Do not add functionality based solely on notebook output. Capture expected numeric behavior in tests.

## Documentation expectations

Document methods in user terms: what scientific question they answer, inputs/units/defaults, outputs, assumptions, failure modes, and a minimal example. Distinguish implemented behavior from the roadmap.

For comparisons across solvents, record the physical model and computational settings; describe results as model-dependent ensemble changes. For ML features, focus on embeddings and molecular comparison, not unsupported property prediction claims.

## Scope roadmap

- **v1:** conformer generation, optimization, filtering, clustering, and core analysis.
- **v2:** validated ensemble-descriptor set and versioned CFP.
- **v3:** benchmark across molecule classes and solvents.
- **v4:** investigate correlations with experimental properties and usefulness in ML models.

Prioritize v1 foundations and descriptor correctness over an extensive visualization or report surface.

## Library implementation state

- **`Conformer` (`generators.py`)**: Mutable `slots` dataclass (not frozen) storing `id`, ASE `Atoms`, optional `energy` (float), `energy_unit` (e.g., `"kcal/mol"`), `optimization_method`, `optimization_converged`, and a private reference to its owning RDKit molecule. Library code treats it as immutable by convention: operations build new `Conformer`s rather than editing existing ones.
- **`Ensemble` (`generators.py`)**: Mutable `slots` dataclass holding `smiles`, RDKit `Mol` (`molecule`), tuple of `Conformer`, and `metadata`. `from_smiles()` delegates to `generate()`. `show()` renders read-only summaries via `display/`. `setup()` is a quick generate → hierarchical optimize (default MMFF → GFN2-xTB, so it needs tblite) → filter helper.
- **Coordinate sync (`generators.py`)**: RDKit conformer coordinates and `Conformer.atoms` are separate arrays. `sync_rdkit_from_atoms` / `sync_atoms_from_rdkit` are the only sync paths; they require an existing conformer ID, equal atom counts, and identical element order, and raise `ValueError` otherwise. Every operation returning an `Ensemble` returns one whose copies match and which owns its RDKit `Mol` and `Atoms` (never shared with the input). Tests: `tests/test_coordinate_sync.py`. Legacy `optimize.py` keeps its own sync code and is not covered.
- **Conformer Generation (`generators.py`)**: `generate()` uses RDKit ETKDGv3 with seed 42 to produce explicit-hydrogen conformers aligned between RDKit and ASE `Atoms`. It raises `ValueError` for multi-fragment SMILES (salts, metal adducts such as `CC(=O)O.[Na+]`), because ETKDG places fragments arbitrarily.
- **Ions / electronic state (`generators.py`)**: `Ensemble.charge` (net formal charge) and `Ensemble.multiplicity` (1 + radical electrons, high-spin assumption) are read-only properties computed from `molecule`, never stored. They are the single source of truth: `generate()` records both in `metadata` and the first history entry, `show()` displays them, every optimizer history entry records `charge`, and GFN2-xTB/ORCA (and legacy `optimize()`) default to them unless `charge=`/`multiplicity=` is passed. The user supplies charged SMILES; automatic (de)protonation (protomer choice) and metal adducts are out of scope. Energies are comparable only within one charge state. Tests: `tests/test_ions.py` (PFOA anion, propylammonium; GFN2 tests need tblite).
- **Optimizers (`optimizers/`)**: the primary API. `MMFFOptimizer`, `UFFOptimizer`, `GFN2xTBOptimizer` (tblite, imported lazily and checked at construction) and `ORCAOptimizer` (external executable, located on `optimize()`) share `BaseOptimizer.optimize()`. It returns a new `Ensemble` with RDKit and ASE coordinates written back in sync, per-conformer energies in kcal/mol, per-conformer `optimization_converged`, and an appended history event. `optimize.py` is the older function-based equivalent and is kept for compatibility.
- **Hierarchical Optimization (`optimizers/hierarchical.py`)**: `HierarchicalOptimizer` chains stages sequentially, and each stage appends its own history event. Conformer-level fields reflect the last stage only.
- **Filters (`filters/`)**: `EnergyFilter`, `PopulationFilter`, `RMSDFilter` and `CompositeFilter` return new ensembles holding a subset of the input conformers (on a copied `Mol` that still contains the removed RDKit conformers; conformers rebound with copied `Atoms`) and append a filter history event. `RMSDFilter` keeps the lowest-energy representative (ties keep input order), optionally uses heavy atoms only (`heavy_atoms_only=True`, via `Chem.RemoveHs`), refuses ensembles whose RDKit and ASE coordinates disagree, never mutates its input, and records `kept_conformer_ids`, `removed_conformer_ids`, and a `duplicate_of` mapping (removed ID → kept ID, RMSD). Tests: `tests/filters/test_rmsd_filter.py`.
- **RMSD analysis (`analysis/rmsd.py`)**: `rmsd`, `rmsd_matrix`, and `rmsd_heatmap` align a copy of the molecule, so they are read-only. Tests: `tests/analysis/test_rmsd_readonly.py`.

## Hierarchical optimization guidelines

`HierarchicalOptimizer` already provides multi-stage cascading refinement (e.g. fast force-field screening followed by quantum-chemical refinement). Remaining work should follow these guidelines:

1. **Multi-Stage Workflows** (implemented): ordered optimization tiers (e.g., Stage 1: `MMFF` -> Stage 2: `GFN2-xTB`).
2. **Intermediate Pruning & Filtering** (not implemented): optional relative energy windows ($\Delta E$ cutoffs in kcal/mol) or RMSD deduplication between stages, to prune high-energy/duplicate conformers before expensive quantum backend runs. Reuse `EnergyFilter`/`RMSDFilter` rather than reimplementing them.
3. **Traceability & History**: Accumulate full provenance across all stages in the chronological `metadata["history"]` list and preserve original conformer IDs. Discarded conformers should be recorded in provenance metadata (filters such as `RMSDFilter` already record removed IDs).
4. **Immutability & Safety**: Each tier returns a fresh `Ensemble` instance; input ensembles are never mutated in place.

---

# Part 2 — Research focus: conformational sampling convergence

## Research question

**How many conformers are needed to adequately characterize a molecule's conformational ensemble, and can that sampling requirement be predicted from molecular properties?**

The goal is to develop a reproducible experimental pipeline, identify meaningful convergence criteria, build a dataset relating molecular descriptors to sampling requirements, and eventually train a model that can guide conformer generation for new molecules.

This is an exploratory investigation. Do not assume the hypothesis is correct or claim novelty before establishing evidence.

## Research plan

### Phase A — Implement Ensemble Metrics

Build convergence-facing APIs under `src/ensemblelab/convergence/`, reusing general-purpose metric calculations and existing RMSD/clustering functionality where appropriate.

1. **Energies:** Relative energies and energy distributions.
2. **Populations:** Boltzmann weights, conformer and cluster populations, Shannon entropy, effective ensemble size.
3. **Geometry:** Radius of gyration and solvent-accessible surface area (SASA).
4. **Torsions:** Raw dihedral angles and discrete torsional states (e.g., gauche+, gauche−, anti).
5. **Properties:** Additional conformer-level properties, including dipole moments and intramolecular hydrogen bonds when suitable methods and data are available.
6. **Convergence:** Metric-versus-sample-size curves, cluster discovery rate, and metric-specific convergence criteria.

Each module must have documented assumptions, consistent interfaces, and focused tests. Implement incrementally; inspect existing code before duplicating functionality.

### Phase B — Build a Pilot Dataset

Experiment orchestration belongs in `src/experiments/sampling_convergence/`, separate from reusable library functionality.

Start with a small, chemically diverse test set, targeting 25–50 molecules if feasible.

For each molecule:

1. Record its identity, canonical representation, and relevant molecular descriptors.
2. Generate conformers using documented, reproducible settings.
3. Optimize structures consistently and record energies, failures, and runtimes.
4. Analyze ensembles at increasing sample sizes.
5. Calculate ensemble metrics at each sampling size and generate convergence curves.
6. Repeat independent sampling runs where feasible to assess reproducibility.
7. Save raw structures, per-conformer metrics, cluster-level metrics, convergence results, and molecule-level summaries.

Use checkpoints so interrupted runs can resume without losing completed work. Record random seeds, software versions, optimization methods, energy units, and configuration settings.

### Phase C — Define Sampling Sufficiency

Determine how to translate convergence curves into defensible sampling requirements.

* Define explicit, metric-specific convergence criteria.
* Expect energy, structural diversity, populations, torsions, and ensemble properties to converge at different rates.
* Investigate how thresholds and independent runs affect estimated sampling requirements.
* Define an overall recommended sample size only after establishing how individual criteria should be combined.
* Distinguish apparent convergence from evidence of adequate sampling; a flat curve does not guarantee that important conformational states were discovered.

The resulting labels should record the estimated sample requirement for each metric and, if justified, an overall requirement. Document the definition and uncertainty of every label.

### Phase D — Prepare for Machine Learning

Build a structured dataset linking molecular descriptors to the sampling requirements established in Phase C.

Potential input features include atom count, molecular weight, formal charge, rotatable-bond count, and other justified structural descriptors. Avoid selecting features solely because they are easy to calculate.

Potential prediction targets include metric-specific convergence sizes and an overall recommended sampling size.

Once sufficient reliable labels exist:

1. Establish simple baseline models.
2. Split data by molecule, preventing conformers or repeated runs of the same molecule from leaking across training and test sets.
3. Evaluate prediction error and whether predicted sampling sizes meet the intended convergence criteria.
4. Compare predictions against more extensive reference sampling.
5. Assess generalization, uncertainty, and sensitivity to the convergence definitions.

Do not prioritize ML training before the data-generation pipeline and labels are trustworthy.

### Phase E — Integrate into EnsembleLab

Eventually, accept a molecule such as a SMILES string, calculate its descriptors, predict an appropriate sampling size, and use that estimate to guide conformer generation.

Evaluate this adaptive workflow against a fixed sampling strategy. Treat predictions as estimates with limitations, not guarantees of complete conformational-space coverage.

## Experiment files

```text
src/experiments/sampling_convergence/
├── config.yaml
├── convergence_study.ipynb
├── data/
│   ├── checkpoints.csv
│   └── molecules.csv
└── results/
    ├── cluster_metrics.csv
    ├── conformer_metrics.csv
    ├── convergence_metrics.csv
    └── molecule_summary.csv
```

## Research standards

These add to the library-wide scientific conventions in Part 1.

* Distinguish observed conformer frequencies from energy-derived Boltzmann populations.
* Validate energy units, temperature assumptions, weighting, and cluster-population calculations.
* Calculate entropy at both conformer and cluster levels; near-duplicate conformers can inflate conformer-level effective ensemble size.
* Keep raw structures and detailed results separate from aggregated metrics.
* Ensure reproducibility, checkpointing, and recovery from interrupted jobs.
* Never fabricate results, runtime estimates, or scientific conclusions.
* Clearly distinguish implemented features, experimental findings, hypotheses, and future plans.

## Convergence API and implementation placement

* Keep public convergence entry points simple and ensemble-facing: accept the existing `Ensemble` API and return structured, dataset-friendly results rather than printed output or plots.
* Put general-purpose metric calculations in `src/ensemblelab/analysis/` when they are useful independently of convergence. A small function in `convergence/` may re-export that implementation to provide a clear, stable convergence-facing API; `convergence/energies.py` is the current example.
* Put logic in the convergence module when it is inherently about sampling size, selecting or comparing prefixes/subsets, calculating successive-sample curves, or applying convergence-specific criteria. That logic may compose functions from `analysis/`; do not duplicate their metric calculations.
* Keep the boundary pragmatic. Avoid wrappers that add no useful public API, but do not force convergence-specific orchestration into general analysis modules. Document parameter semantics, result fields, units/assumptions, and missing-data behavior, and add focused tests at the public entry point.

## Verified convergence implementation status

* **Energy metrics:** Implemented in `ensemblelab.analysis.energies.analyze_energies` and re-exported from `ensemblelab.convergence.energies`. Accepts an `Ensemble` and optional ordered subset of its owned conformers; returns ordered conformer records and one summary mapping. Source IDs and full-ensemble indices are retained.
* **Scientific contract:** Uses kcal/mol, matching the optimizer interfaces. Finite numeric energies must explicitly declare kcal/mol; unsupported or missing units raise `ValueError`. Missing, nonnumeric, and non-finite energies are retained, marked invalid, and excluded from statistics. Standard deviation is population standard deviation (`ddof=0`). No values are rounded or mutated.
* **Boltzmann conformer populations:** Implemented in `ensemblelab.analysis.populations.analyze_populations` and re-exported from `ensemblelab.convergence.populations`. Reuses `analyze_energies`, defaults to 298.15 K, records the temperature and kcal/mol Boltzmann constant, and reports normalized weights/populations, entropy (nats), and conformer-level effective size. Invalid energies are excluded and marked missing; no-valid-energy statistics are `None`. Existing basin weights reuse the same stable weight calculation.
* **Validation:** `tests/convergence/test_energies.py` covers hand-verifiable values, invalid/missing data, empty/equal/single cases, stable ordering/subsets, non-mutation, units, schema, and CSV serialization (11 passed). `tests/convergence/test_populations.py` covers known distributions, entropy/effective size, invalid inputs, underflow, temperatures, immutability, CSV serialization, and existing basin-population compatibility (15 passed).
* **Conformer/cluster entropy:** Implemented in `ensemblelab.convergence.entropy.analyze_entropy`, using `analysis.entropy.entropy_statistics`, Boltzmann populations, and `analysis.clusters.cluster_basins`/`Basin`. Cluster entropy is reported only for complete, non-overlapping membership of the analyzed conformers. Computed RMSD method, metric, and cutoff are recorded; supplied basins do not carry provenance, so their method/cutoff are reported as supplied/unavailable. Existing `basin_entropy` delegates to the shared entropy helper.
* **Entropy validation:** `tests/convergence/test_entropy.py` covers known distributions, basin aggregation, invalid and zero populations, empty/all-invalid data, incomplete/overlapping assignments, generated clustering metadata, CSV compatibility, and non-mutation. Entropy-only run: 11 passed; full pytest suite after structural geometry implementation: 49 passed.
* **Structural geometry:** Implemented in `ensemblelab.convergence.geometry.analyze_geometry`. It calculates mass-weighted radius of gyration (RDKit `Descriptors3D.RadiusOfGyration`, angstroms) and SASA (RDKit FreeSASA Lee-Richards, RDKit periodic-table van der Waals radii, 1.4 angstrom probe, square angstroms). It accepts ordered subsets of ensemble-owned conformers, reads ASE coordinates, validates atom order, calculates on private RDKit copies, and reports metric-specific status/error and distribution summaries. FreeSASA absence is explicitly reported as unavailable; no approximation is substituted.
* **Geometry validation:** `tests/convergence/test_geometry.py` covers reference calculations, metadata and units, differing geometries, empty/single ensembles, invalid coordinates, independent metric failures, non-finite outputs, unavailable FreeSASA, ordering, CSV serialization, and input immutability. Full pytest suite after structural geometry implementation: 49 passed.
* **Torsional metrics:** Implemented in `ensemblelab.convergence.torsions.analyze_torsions`. Selects acyclic single bonds between nonterminal heavy atoms, excludes ring bonds and amide/ester C-N/O bonds, and defines one reproducible dihedral per bond using lowest-index outer neighbors. Reports per-conformer angles normalized to [-180, 180), per-torsion circular mean/resultant/variance, and coverage of configurable angular bins (default 12 equal-width 30-degree bins). Failed observations remain missing and explicitly marked; aggregate means exclude torsions without valid observations.
* **Torsion validation:** `tests/convergence/test_torsions.py` covers RDKit and existing analysis API reference angles, wraparound and undefined circular means, state occupancy/boundaries, invalid edges and coordinates, partial failures, eligibility, ordering, CSV serialization, and input immutability. Full pytest suite after torsion implementation: 75 passed.
* **Deduplication (not yet integrated into convergence):** `RMSDFilter` (see Part 1, Library implementation state) is ready to act as a deduplicator for a shared reference pool. Open decisions: whether deduplication runs once on the pool or per sampling-size subset, and which atom selection to use. `cluster_basins` uses all atoms with a strict `<` threshold, while the filter uses `<=` with optional heavy-atom selection. Full suite after the RMSDFilter work: 103 passed, 2 skipped (tblite).
* **Limitations / next step:** Unit conversion and degeneracy factors are unsupported. Boltzmann estimates do not account automatically for conformational entropy, vibrational effects, solvent, sampling incompleteness, or duplicate physical states. Entropy is in nats and effective size is conformer- or cluster-state count conditional on population and clustering assumptions. Radius of gyration and SASA alone do not measure conformational diversity. RDKit asphericity is available but deferred pending a specific research justification. Equal-width torsion bins are not chemically equivalent or equal-energy states. Additional ensemble properties and overall convergence modules remain unimplemented; implement properties next.
