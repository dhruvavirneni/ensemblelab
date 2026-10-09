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
