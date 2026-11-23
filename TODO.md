# TODO

## 1. Convergence definitions

Decide these before generating new data. Record them in `docs/convergence_definitions.md`.

- [ ] Physical observables that define N: Boltzmann-weighted radius of gyration, asphericity, lowest relative energy
- [ ] Sampling diagnostics reported but not used for N: `n_effective`, `n_within_3kcal`, `n_unique`
- [ ] Tolerance per observable (e.g. 1% for shape metrics, 0.1 kcal/mol for energy)
- [ ] Confidence level: N90 as the main label, with N80 and N95 reported
- [ ] Number of random draws per N (raise from 10 to ~30)
- [ ] Fixed pipeline settings: ETKDGv3 seed, force field, RMSD dedup cutoff (0.5 Å, heavy atoms), temperature
- [ ] Censoring rule for observables that never converge within the largest pool

## 2. Packaging and release

- [ ] Create PyPI and TestPyPI accounts with 2FA
- [ ] Add a GitHub Actions publish workflow using Trusted Publishing
- [ ] Upload to TestPyPI, install in a fresh environment, run a quick MMFF example
- [ ] Release v0.1.1 on GitHub to publish to PyPI
- [ ] Update README install instructions to `pip install ensemblelab`

## 3. Convergence analysis code

- [ ] Add N90 / NXX calculation from per-draw results
- [ ] Add per-observable convergence summary (one N per observable, plus max over chosen observables)
- [ ] Handle censored results (`> pool size`) explicitly in outputs
- [ ] Tests for convergence calculations on synthetic data with known answers

## 4. Larger pools

- [ ] Run 2,000 to 5,000 conformer pools for n-octanol and decanol
- [ ] Check whether energy and population metrics flatten
- [ ] Add results to `docs/convergence_results.md` and a research log entry
- [ ] Tag v0.2.0

## 5. Convergence dataset

- [ ] Build a diverse set of ~30 molecules (spread across rotatable bonds, rings, polarity); record SMILES and source
- [ ] Generate N90 labels per observable using the fixed pipeline
- [ ] Compute RDKit features: rotatable bonds, heavy atoms, rings, fraction sp3, HBD, HBA, logP, TPSA
- [ ] Save dataset as a versioned CSV in `experiments/`

## 6. N prediction (baseline first)

- [ ] Baseline: log(N90) vs. rotatable bonds
- [ ] Simple models: linear regression, then random forest / gradient boosting
- [ ] Evaluate with scaffold split (MurckoScaffold) and cross-validation
- [ ] Write up whether any model beats the baseline

## 7. Later (only if step 6 shows signal)

- [ ] Scale dataset to 100 to 500 molecules
- [ ] Fit each molecule's error curve, E(N) ~ a * N^(-b), and predict a and b
- [ ] Cost-vs-error selection of N: minimize C(N) + lambda * E(N)
- [ ] Ion support (charge, multiplicity) and CCS as an observable
