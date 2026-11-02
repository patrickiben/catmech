# Scoring the recount under corrected variances, and related checks

These are three checks from 2026-11-02, after a literature search. None changes package code or an earlier verdict. No outside reviewer has read them. A [later test](VARIANCE_MI_SYNTHESIS.md) carried these variances through imputation and changed the package's default.

## The recount under corrected variances

Under the [reading fixed in advance](VARIANCE_CHECK_PLAN.md), **`cr2`, `cv3` and `wcr` all qualify as candidates to try through imputation, and the current variance does not**. The current variance is what `catmech_pool_lm(cluster = )` uses and what `clubSandwich` calls CR1S. There were 24,500 requests on the [replacement study's](CLUSTERED_FIX_PLAN.md) cohorts, with no imputation. [Scored output](variance-check-decisions.txt) comes from a separate script.

| Mean coverage | Current | `cr2` with Satterthwaite df | `cv3` | `wcr` |
|---|---|---|---|---|
| 24 clustered simulated cells | 92.1 | 95.1 | 94.8 | 94.9 |
| 6 rare-count simulated cells | 83.2 | 87.5 | 87.5 | 89.5 |
| 8 row-only simulated cells | 93.0 | 95.5 | 95.1 | 95.2 |

| 11 real items | Current | `cr2` | `cv3` | `wcr` |
|---|---|---|---|---|
| Mean distance from 95%, points | 2.49 | 1.53 | 1.45 | 1.42 |
| Lowest coverage | 88.4 | 91.6 | 93.2 | 92.8 |
| Items at 92% or more | 7 | 10 | 11 | 11 |
| Simulated cells within 92% to 98%, of 38 | 22 | 32 | 31 | 32 |

What this shows:

- **The analysis itself was part of the shortfall.** With every row recounted, the current variance undercovers, and all three corrections largely repair that outside rare counts. Earlier verdicts read coverage against this low benchmark, so a method covering above the oracle may have been offsetting a variance that was too small. On `lbw`, the oracle covers 88.4% under the current variance and 91.6% to 94.0% under the corrections.
- **Not all of the gap is the analysis.** On `anc1_early` and `syphilis`, the corrected recount covers about 96% and 95% to 97%, while the replacement option covered 91.0% and 89.0% under the current variance. Those gaps still point to the imputation, though the option has not been scored under the corrections.
- **The corrections cost width.** On real items, `cr2` intervals are about 10% to 15% wider than the current ones and `cv3` about 14% to 23% wider, by the rounded widths in the scored output. `wcr` has no width here, because it tests the target directly.
- **69 clusters behave more like 20.** The mean Satterthwaite degrees of freedom on the real items are 18.2 to 23.0, against the 68 the current variance assumes. The covariate is skewed, and the largest cluster's partial leverage averages 0.11 to 0.14.
- **Rare counts are not repaired.** Every correction stays below 90% on average there, with the oracle's own data and no imputation. That problem is unsolved here.

Predictions: four of six held. `lbw` under `cr2` reached 91.6%, short of the predicted 92%. The current variance's mean standard error fell below the empirical standard deviation on only five of 11 real items, not nine, so on several items it undercovers despite a ratio above 1. Why is not established.

Limits: 500 cohorts give a Monte Carlo standard error of about 1 point, and three methods sit close together, so this does not rank them. The real items' results under the current variance had been seen. Complete-data performance need not carry through imputation, and no source found tests that. `cr2` was checked against explicit matrices and the HC2 form, not against `clubSandwich`, which is not installed here. The plan was written before the run, and its hash was recorded then in [`variance-check/plan-hashes.txt`](variance-check/plan-hashes.txt), but it was not committed first, so the order rests on that record and file times. A pilot of three cohorts per cell was read for errors only.

## Independent numerical checks

[All 50 scored checks passed](cluster-algebra-check.txt) on fixed datasets, seven of them added when `cv3` entered the package. Variance components match a one-way analysis of variance and `lme4` REML in a balanced layout. At REML components, generalized least squares matches `lme4` fixed effects and covariance and an explicit inverse. Conditional intercept means and variances match `lme4`. The CR1S variance matches a hand-written double sum, and the package's `cv3` matches delete-one-cluster refits. Rubin's rules and Barnard-Rubin degrees of freedom match a hand calculation, and the estimate and standard error match `mitools`, whose degrees of freedom differ slightly (21.5928 against 21.2917 in one case), which is reported and not scored. The sampler's draws match a separate implementation with explicit matrices, within Monte Carlo limits.

This catches coding errors. It does not validate the method. The separate implementation uses the same moment formula, so a wrong choice of formula would pass, and `pan`, `jomo` and `brms` have not been compared. Whether the draw of the intercept variance is well calibrated remains untested.

## The Zambia file

The data are from a performance-based financing study, [Grover, Bauhoff and Friedman, PLOS ONE 2019](https://doi.org/10.1371/journal.pone.0211262), and the record is at [Harvard Dataverse](https://doi.org/10.7910/DVN/LHUIBO) under CC0. Two files were read, `Zambia_data.xlsx` and its dictionary, and neither is in this repository. No model was fitted, so the data remain available for a planned test.

- **Structure.** There are 140 clinics with an identifier, in 18 districts (1 to 27 clinics each) and 9 provinces, with four intervention arms.
- **Paired values.** Nine items each have a tally-sheet value (`_ts`, verified) and a reported form value (`_hia`), for the first quarter only. Over-reporting flags cover all four quarters.
- **Seen so far.** Complete pairs per item number 106 to 125, and exact agreement 42 to 77. Mean verified values run from 5.3 (`exposedbabies`) to about 3,372 (`opd`).
- **Not settled.** Zeros may partly mean missing: up to 43 pairs per item are zero on both sides, and up to 15 are reported as positive with a verified zero. The two README files in the record were not read.

So this is a usable second health system, but not a repeat of Rwanda. Each clinic has one row per item, so a clinic-level audit has independent rows. Clustering can be tested only by district, with 18 unequal clusters, which is a hard case for every method here.
