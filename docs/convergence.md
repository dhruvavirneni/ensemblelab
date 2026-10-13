# Convergence
* incomplete documentation section...just contains rough notes

## Central question: How does increasing conformer sampling affect the structural and energetic characterization of a molecule, and how can those changes help determine optimized conformer sampling sizes in ensemble analysis?
 * other stuff i'm wondering: in a given Ensemble with n conformer sample size, will all metrics converge relatively at the same n?

## goals and stuff
Investigate how many conformers are needed to adequately characterize a molecule's conformational ensemble, and whether that number can be predicted from the molecule's properties. The first stage is to study how ensemble-level metrics change as the number of generated conformers increases. Starting with a set of molecules spanning different sizes and levels of structural complexity, I will generate conformational ensembles at progressively larger sampling sizes, optimize their geometries, and standardize the resulting structures and energies for consistent comparison.

Ensemble wide metrics worth calculating/investing convergence on:
 * relative energy distribution: (for each conformer, calculate energy relative to lowest energy conformer) --> good for comparing conformational stability
 * boltzmann populations: assign each conformer a statistical weight based on its energy and a chosen temperature, then normalize the weights into estimated populations --> tells us which conformers contribute most to population weighted ensemble
 * conformational entropy and effective ensemble size: calculate Shannon entropy to measure how 'spread out' the ensemble is --> whether an ensemble is dominated by a few conformers or distribute across many...remember to calculate populations at cluster level too
    * according to google: calculating populations at the conformer cluster level is necessary because Shannon entropy is highly sensitive to how states are divided. If your ensemble contains many near-duplicate structures, treating them as distinct conformers artificially inflates the entropy $H$ and distorts the effective ensemble size ($N_{eff}$)
 * torsional state diversity: track distribution of torsions around rotatable bonds --> which internal rotations account for conformational variations
  * can track by simply tracking dihedral angles around given bonds, or assigning torsions to their 'discrete states'(gauche+, gauche-, anti) **investigate if one will be easier for later ML use cases**
 * radius of gyration
  * how spread out atoms are around conformer's center of mass --> whether ensemble is more compact or extended --> output mean, std, range, etc.
 * ensemble molecular properties: solvent accessible surface area, dipole moment, etc...intramolecular hbonds --> whether conformational changes affect properties beyond solely geometry

 * another cool thing (suggested by AI): cluster discovery rate: num of newly discovered clusters per additional batch of conformers --> could complement minimum energy and structural convergence curves


For each sampling size, I will calculate metrics such as minimum energy, pairwise RMSD, structural clustering, and other measures of conformational diversity. Plotting these metrics against the number of conformers should reveal how quickly each stabilizes and whether additional sampling continues to uncover meaningful structural or energetic differences. Since different metrics may converge at different rates, I will investigate how to define both metric-specific convergence points and an overall sampling size that satisfies the chosen criteria.


## pipeline + plan
 1. define experiement 
    * pick molecule of choice, sampling size tests, random seeds, optimization settings, etc.
 2. generate conformers
    * generate large reference pool in Ensemble...use same random seed to generate ensembles with different nconformers
 3. optimize and standardize conformers
    * consistent optimization, deduplication, etc.
 4. analyze each sampling size
 5. quantify convergence per sample size
    * measure metric changes between sample sizes (metrics to look at listed above)
 6. train ML model with generated dataset of multiple molecules, having run this analysis
 7. test ML model, using molecule input and other quick molecular metrics like logP, etc. etc. as parameters and outputting optimal n (sample size) for effectively using this molecule for ensemble analysis
 8. maybe build a web app interface or something that maintains records of ML analyses conducted and runs new ones to output n if an inputted molecule hasn't been analyzed by anyone yet (much later)

# Ensemble-wide analysis modules (to be tracked by convergenced analysis)
## Conformer Energies

### Metric identity and scope

* **Metric:** Per-conformer relative energy and ensemble energy distribution/statistics.
* **Implementation:** `ensemblelab.convergence.energies.analyze_energies`.
* **Scope:** One supplied `Ensemble`, or an ordered subset of its conformers. This function analyzes energies only; it does not calculate Boltzmann populations, entropy, cluster statistics, plots, or convergence criteria.
* **Purpose in convergence runs:** Re-run the analysis on successive ensemble prefixes or subsets to obtain comparable minimum-energy and distribution summaries. A stable minimum is only one convergence signal and does not establish adequate conformational-space coverage.

### Inputs and reproducibility

* `ensemble` (required): Existing `Ensemble` object. By default, all conformers are analyzed in `ensemble.conformers` order.
* `conformers` (optional): Ordered sequence of the identical `Conformer` objects in that ensemble. It enables subset or prefix calculations without copying energy arrays. Duplicate IDs or objects not owned by the ensemble raise `ValueError`.
* Conformer IDs come from `Conformer.id`. `conformer_index` is the zero-based position in the full source ensemble, including for subset analysis.
* No parallel explicit energy sequence is accepted, so there is no separate energy-array length to reconcile.
* Analysis is read-only and performs no rounding. Run, molecule, sampling-size, optimizer, and configuration identifiers remain experiment-layer metadata.

### Definitions and statistical convention

For each valid conformer $i$, relative energy is calculated against the minimum valid energy in the analyzed records:

$$\Delta E_i = E_i - \min_j(E_j)$$

The result includes minimum, maximum, range, arithmetic mean, median, population standard deviation, mean relative energy, and median relative energy. The population standard deviation uses denominator $n$ (equivalent to `statistics.pstdev`); it is zero for one valid energy. Equal energies therefore have zero relative energy and zero range.

### Units and energy validity

* Energies and energy-derived statistics are in **kcal/mol**, the unit assigned by EnsembleLab's optimization backends.
* Finite numeric energies must explicitly declare `energy_unit` as `kcal/mol` (case and whitespace normalized). Unit conversion is not performed. A finite energy with a missing or unsupported unit raises `ValueError`, preventing unitless or mixed-unit calculations.
* `None`, non-numeric, and non-finite values are retained as the original per-conformer `energy`, marked `energy_valid=False`, assigned `relative_energy=None`, counted in `n_invalid_energies`, and excluded from all statistics.
* If no energies are valid, statistics are `None`. For an empty ensemble, records are empty, both energy counts are zero, and statistics are `None`.

### Output schema

The function returns a dictionary with two keys:

* `conformers`: ordered list of records with `conformer_id`, `conformer_index`, original `energy`, `relative_energy`, and `energy_valid`.
* `summary`: one record with `n_conformers`, `n_valid_energies`, `n_invalid_energies`, `energy_min`, `energy_max`, `energy_range`, `energy_mean`, `energy_median`, `energy_std`, `relative_energy_mean`, `relative_energy_median`, and `energy_units`.

The lists and mappings are DataFrame/CSV-compatible. Missing values remain `None` in Python and serialize as blank CSV fields; `energy_valid` and the invalid count distinguish these from numeric zero. Preserve the conformer records as the detailed table and the summary as a separate ensemble-level row/table when aggregating datasets.

### Known limits

Only explicitly declared kcal/mol energies are analyzed. Other units require a deliberate, documented conversion layer before calling this function. This module makes no claim that the minimum energy, or any other single metric, demonstrates sufficient sampling.

## Boltzmann Conformer Populations

### Metric identity and scope

* **Metric:** Energy-derived conformer weights, normalized populations, population entropy, and conformer-level effective ensemble size.
* **Implementation:** `ensemblelab.analysis.populations.analyze_populations`, re-exported from `ensemblelab.convergence.populations`.
* **Scope:** One `Ensemble`, or an optional ordered subset of its conformers. This does not calculate cluster-level populations, geometry, torsions, or convergence criteria.
* **Interpretation:** These are Boltzmann estimates from the supplied conformer energies, not necessarily experimentally accurate solution populations. Effective size counts weighted conformer records, not distinct structural states.

### Inputs and reproducibility

* `ensemble` (required): Existing `Ensemble`; source conformer order and IDs are preserved.
* `temperature` (optional): Positive finite temperature in kelvin; defaults to **298.15 K**, matching the existing clustering convention. The actual value is always returned as `temperature_K`.
* `conformers` (optional keyword): Ordered subset of the identical conformer objects held by the ensemble. Source indices refer to positions in the full ensemble. Energy validity, ownership, relative energies, and units follow `analyze_energies`.
* The implementation records the Boltzmann constant as `boltzmann_constant_kcal_mol_K = 0.00198720425864083` and reports energy units as `kcal/mol`.

### Definitions and statistical convention

Weights use the minimum-shifted energy to improve numerical stability:

$$w_i = \exp\left(-\frac{E_i-E_{\min}}{k_B T}\right), \qquad p_i = \frac{w_i}{\sum_{j \in V} w_j}$$

Here $V$ is the set of valid-energy conformers in the analyzed selection. The denominator excludes invalid energies. Entropy uses natural logarithms and is reported in nats; effective ensemble size is its exponential:

$$H = -\sum_{i \in V,\ p_i>0} p_i \ln(p_i), \qquad N_{\mathrm{eff}} = \exp(H)$$

Zero weights from floating-point underflow are retained as valid zero populations and omitted from the $p_i\ln(p_i)$ sum. The minimum-energy conformer retains weight 1, so normalization remains defined whenever at least one valid energy exists. Equal statistical degeneracy is assumed; no degeneracy factor is represented in the current model.

### Output schema and invalid-energy policy

The result is a mapping with `conformers` and `summary`:

* Each ordered conformer record carries the energy-analysis fields `conformer_id`, `conformer_index`, `energy`, `relative_energy`, and `energy_valid`, plus `boltzmann_weight`, `population`, and `population_valid`.
* The summary contains `temperature_K`, `n_conformers`, `n_valid_energies`, `n_invalid_energies`, `population_sum`, `max_population`, `effective_ensemble_size`, `population_entropy`, `population_entropy_units`, `energy_units`, and `boltzmann_constant_kcal_mol_K`.
* `None`, nonnumeric, and non-finite energies follow the energy-analysis invalid-value policy: they remain traceable in the records, have `population_valid=False` and missing (`None`) weights/populations, are excluded from normalization, and contribute to `n_invalid_energies`. Finite energies with missing or unsupported units raise `ValueError`; units are not converted.
* For empty selections or selections with no valid energies, counts and temperature are returned, while population statistics are `None`. For one valid conformer, population and weight are 1, entropy is 0, and effective size is 1.
* Conformer records are CSV-compatible; invalid weights/populations serialize as blank fields and are disambiguated from valid zero populations by `population_valid`.

### Assumptions and limitations

Ordinary Boltzmann weighting of optimized conformer energies assumes equal degeneracy and does not automatically account for conformational entropy, vibrational contributions, solvent effects, sampling incompleteness, or duplicate representations of a physical state. Cluster-level populations and entropy are a separate analysis; do not interpret conformer-level effective size as a count of distinct states. Population stability alone does not establish adequate sampling.

## Conformer and Cluster Entropy

### Metric identity and scope

* **Metric:** Shannon entropy and effective ensemble size from conformer populations, compared with entropy and effective size after aggregating populations by structural cluster.
* **Implementation:** `ensemblelab.convergence.entropy.analyze_entropy`, using the reusable `ensemblelab.analysis.entropy.entropy_statistics` function.
* **Population source:** Boltzmann populations from `ensemblelab.analysis.populations.analyze_populations`; no separate entropy-specific weighting implementation is used.
* **Scope:** One `Ensemble` or an ordered conformer subset. This analysis does not generate a new clustering algorithm, population model, or convergence curve.

### Inputs and clustering provenance

* `ensemble` (required): Existing EnsembleLab ensemble.
* `temperature` (optional): Positive finite kelvin temperature passed to Boltzmann population analysis; defaults to 298.15 K and is returned in the summary.
* `conformers` (optional keyword): Ordered subset of ensemble-owned conformers. Existing energy/population selection and full-ensemble conformer indices are preserved.
* `basins` (optional): Existing `Basin` assignments from `analysis.clusters`. If omitted, `cluster_basins` is called using the configured `rmsd_threshold` (default 0.75 Å) and temperature. Supplied assignments are not changed or regenerated.
* Automatically generated assignments report method `greedy_best_fit_rmsd`, distance `RDKit best-fit RMSD`, and `rmsd_threshold_A`. Supplied basins report `provided_basin_assignments`; the `Basin` API does not retain the generating method or cutoff, so these are reported as unavailable rather than inferred.

### Definitions and normalization

For normalized probabilities $p_i$, entropy uses natural logarithms and is measured in nats. Effective size is the exponential of entropy:

$$H = -\sum_{i:p_i>0} p_i\ln(p_i), \qquad N_{\mathrm{eff}}=\exp(H)$$

Terms at exactly zero probability contribute zero (`0 ln 0 = 0`). `analysis.entropy.entropy_statistics` accepts only finite, nonnegative, already-normalized probabilities (sum within numerical tolerance); it rejects invalid or unnormalized input. Empty input produces missing entropy and effective size. Boltzmann conformer populations are already normalized over valid-energy conformers by the population API.

Cluster probability is the sum of the populations of valid member conformers. Cluster entropy is calculated from those aggregated values only when basin membership is a complete, non-overlapping partition of all selected conformers. Invalid conformer populations remain missing and are not treated as zero; clusters with no valid population have `cluster_population=None`. A valid cluster with an underflowed population of zero is retained and has zero entropy contribution. If assignments omit, overlap, or include a conformer outside the analyzed selection, cluster entropy and contributions are unavailable; conformer entropy remains calculable from valid populations.

### Output schema

The result has `conformers`, `clusters`, and `summary` records:

* Per-conformer records preserve energy-analysis and population-analysis fields, add `cluster_id` when uniquely assigned, and retain `None` for unassigned or ambiguous IDs.
* Per-cluster records contain `cluster_id`, `n_conformers`, `n_valid_populations`, `cluster_population`, and `cluster_entropy_contribution`. Join members through `cluster_id` in the per-conformer records. Per-cluster effective-size contribution is intentionally omitted: effective size is nonlinear and has no additive cluster contribution.
* The summary includes conformer/cluster counts and valid-population counts; conformer and cluster entropy/effective size; `population_source`; `temperature_K`; energy/entropy units; clustering method, distance, and cutoff metadata; assignment completeness/availability; and invalid/unassigned conformer counts.
* `conformer_entropy` can be calculated even if cluster assignments are incomplete. `cluster_entropy` and `cluster_effective_size` are `None` if assignment is incomplete or no valid populations exist. For an empty ensemble, records are empty and entropy values are `None`.

### Interpretation and limitations

Conformer-level entropy can be inflated by redundant near-duplicate structures. Cluster entropy is conditional on the clustering algorithm, distance metric, cutoff, and population model; it is not an absolute measure independent of those choices. Neither quantity alone establishes sampling completeness or should be interpreted as physical thermodynamic conformational entropy without the necessary statistical-mechanical assumptions. Boltzmann populations inherit the limitations documented above.
