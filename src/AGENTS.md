# AGENT.md — EnsembleLab Conformational Sampling Research

## 1. Project Mission

Work exclusively on the conformational sampling convergence project until instructed otherwise.

EnsembleLab is a Python library for generating and analyzing molecular conformational ensembles. This research investigates:

**How many conformers are needed to adequately characterize a molecule's conformational ensemble, and can that sampling requirement be predicted from molecular properties?**

The goal is to develop a reproducible experimental pipeline, identify meaningful convergence criteria, build a dataset relating molecular descriptors to sampling requirements, and eventually train a model that can guide conformer generation for new molecules.

This is an exploratory investigation. Do not assume the hypothesis is correct or claim novelty before establishing evidence.

## 2. Research Plan

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

## 3. Current Project Structure

Existing package code is under `src/ensemblelab/`, which includes analysis, display, filters, optimizers, and other existing functionality.

Convergence modules are under `src/ensemblelab/convergence/`.

Experiment files currently live under:

```text
src/experiments/sampling_convergence/
├── config.yaml
├── data/
│   ├── checkpoints.csv
│   └── molecules.csv
└── results/
    ├── cluster_metrics.csv
    ├── conformer_metrics.csv
    ├── convergence_metrics.csv
    └── molecule_summary.csv
```

Inspect the actual repository before assuming filenames, schemas, or APIs. Preserve existing functionality and avoid unnecessary file moves.

## 4. Scientific and Engineering Standards

* Distinguish observed conformer frequencies from energy-derived Boltzmann populations.
* Validate energy units, temperature assumptions, weighting, and cluster-population calculations.
* Calculate entropy at both conformer and cluster levels; near-duplicate conformers can inflate conformer-level effective ensemble size.
* Treat torsion angles as circular variables. Preserve raw angles and support discrete states for comparison.
* Record failures explicitly rather than silently discarding them.
* Keep raw structures and detailed results separate from aggregated metrics.
* Ensure reproducibility, checkpointing, and recovery from interrupted jobs.
* Never fabricate results, runtime estimates, or scientific conclusions.
* Clearly distinguish implemented features, experimental findings, hypotheses, and future plans.

## 5. Agent Workflow and Persistent Context

### Convergence API and implementation placement

* Keep public convergence entry points simple and ensemble-facing: accept the existing `Ensemble` API and return structured, dataset-friendly results rather than printed output or plots.
* Put general-purpose metric calculations in `src/ensemblelab/analysis/` when they are useful independently of convergence. A small function in `convergence/` may re-export that implementation to provide a clear, stable convergence-facing API; `convergence/energies.py` is the current example.
* Put logic in the convergence module when it is inherently about sampling size, selecting or comparing prefixes/subsets, calculating successive-sample curves, or applying convergence-specific criteria. That logic may compose functions from `analysis/`; do not duplicate their metric calculations.
* Keep the boundary pragmatic. Avoid wrappers that add no useful public API, but do not force convergence-specific orchestration into general analysis modules. Document parameter semantics, result fields, units/assumptions, and missing-data behavior, and add focused tests at the public entry point.

Before each task, read this file and inspect relevant code, tests, and interfaces. Implement one module or well-defined pipeline component at a time.

Use focused tests, verify scientific calculations, and run relevant tests after changes. Avoid unrelated refactors, unnecessary dependencies, and large rewrites.

**Update this file whenever a significant research decision, metric definition, convergence criterion, data schema, interface, completed phase, or limitation changes.** Keep it concise and factual. Reflect actual code and verified results; mark unfinished work clearly. Remove superseded plans rather than accumulating contradictory instructions.

At the end of each task, report changes made, tests and outcomes, unresolved issues, and the next recommended step. Ask before making major scientific or architectural changes.

## Verified Implementation Status

* **Energy metrics:** Implemented in `ensemblelab.analysis.energies.analyze_energies` and re-exported from `ensemblelab.convergence.energies`. Accepts an `Ensemble` and optional ordered subset of its owned conformers; returns ordered conformer records and one summary mapping. Source IDs and full-ensemble indices are retained.
* **Scientific contract:** Uses kcal/mol, matching the optimizer interfaces. Finite numeric energies must explicitly declare kcal/mol; unsupported or missing units raise `ValueError`. Missing, nonnumeric, and non-finite energies are retained, marked invalid, and excluded from statistics. Standard deviation is population standard deviation (`ddof=0`). No values are rounded or mutated.
* **Validation:** `tests/convergence/test_energies.py` covers hand-verifiable values, invalid/missing data, empty/equal/single cases, stable ordering/subsets, non-mutation, units, schema, and CSV serialization. Focused result: 11 passed.
* **Limitations / next step:** Unit conversion is intentionally unsupported. The energy module does not implement populations, geometry, torsions, properties, clustering, or the convergence pipeline; implement and validate the population metrics as the next independent component.
