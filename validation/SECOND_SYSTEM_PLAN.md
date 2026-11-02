# Prospective test: a second health system

Written 2026-11-02 and committed before the runner exists, so version history records the order. Amendments, if any, go in the section at the end.

## Why

Every clustered test so far uses one recount, from Rwanda, and every item in it has been seen. The evidence for `catmech_reference(cluster = )`, and for the shortfall of linear imputation on real clustered items, could be specific to that file. This test uses a health system that played no part in building the option.

## Data

Grover, Bauhoff and Friedman, *PLOS ONE* 2019 ([doi:10.1371/journal.pone.0211262](https://doi.org/10.1371/journal.pone.0211262)), data at [Harvard Dataverse](https://doi.org/10.7910/DVN/LHUIBO) under CC0: 140 Zambian clinics in 18 districts, with reported values from Health Information Aggregation forms (`_hia`) and verified values from tally sheets, activity sheets and registers (`_ts`), for nine incentivized items in the first quarter of 2013. Facilities were selected by district with probability proportional to catchment size.

**What has already been seen.** Column names, row counts, missing counts, the covariate's range, counts of pairs, of exact agreement and of zeros, and each item's mean verified value. Nothing has related an item to a covariate, and no method has been fitted.

**The pool for each item.** Clinics with both values and a catchment size. Rows where both values are zero are dropped, since the paper does not say whether a zero means no service or no record. That leaves 78 to 118 clinics per item. A positive reported value with a verified zero stays, as a possible real discrepancy. **A sensitivity run keeps the both-zero rows.**

**The analysis.** `verified value ~ covariate`, with the covariate the base-10 logarithm of catchment population, because catchment runs from 700 to 98,000. The target is that coefficient on the item's whole pool. The audited value is the outcome.

## Cohorts and audits

Each cohort draws the 18 districts with replacement, with each draw a new cluster. Four audits are drawn from it:

- **`clinic_25` and `clinic_50`.** A random 25% or 50% of rows, rounded, so most districts are partly audited. This exercises the option's conditional intercept, which the Rwanda tests barely used.
- **`district_33` and `district_50`.** A random third or half of the 18 draws, 6 or 9, with all their rows audited. A 25% audit would give 4 or 5 clusters, at or under the option's refusal limit, so it is not used.

There are 1,000 cohorts per item and audit, 36,000 in all. Seeds: cohort 2,010,000,000 + item x 1,000,000 + replicate, audit 2,020,000,000 + item x 1,000,000 + audit x 100,000 + replicate, imputation 2,030,000,000 plus the same, and `mice` 2,040,000,000 plus the same.

## Methods

All intervals use `cv3` with t(G - 1) as complete-data degrees of freedom and, for imputation, Barnard-Rubin:

- `oracle`: the verified value on every row.
- `recorded`: the reported value on every row.
- `audited_only`: the audited rows alone.
- `linear_mi`: `catmech_reference(use_recorded = TRUE)`, 20 imputations, pooled by `catmech_pool_lm(cluster = )`.
- `random_intercept`: the same with `cluster`.
- `mice_2l_lmer`: `mice` with `2l.lmer`, the recorded value kept by the wrapper from the [earlier rerun](OUTSIDE_METHODS_SENSITIVITY_PLAN.md), same pooling.

Raking and the difference estimator are left to the separate test of independently written setups. Requests run in a child process with a 60-second limit, and one over the limit counts as no interval.

**Sensitivity:** the both-zero rows kept, 500 cohorts, `oracle`, `linear_mi` and `random_intercept` only.

## Decision rule

Let D be the mean over the nine items of the distance of coverage from 95%, and R the geometric mean over the items of the RMSE of `random_intercept` over that of `linear_mi`, paired on cohorts where both gave an interval. The rule uses two cells: `clinic_25` and `district_50`. A cell **supports** the option if D for the option is at least 0.5 points below D for linear MI, R is at most 1.02, and the option gives an interval on at least 95% of requests. It **goes against** the option if D for the option is not below D for linear MI, or R is above 1.05, or that availability fails.

- **REPLICATES** if both cells support the option. The documents say the second system supports the Rwanda result.
- **CONTRADICTS** if both cells go against it. `catmech_reference(cluster = )` then gets a deprecation warning in the next version, since its case rests on real items, and the documents change to match. The code stays, so earlier results can be reproduced.
- **MIXED** otherwise. The documents say so.

The option stays experimental whatever happens. No earlier verdict changes, and thresholds do not change after viewing.

**Description only:** the other two audits; linear MI's and the oracle's coverage by item; whether a partial audit improved on the recorded value; `mice_2l_lmer` against the option; the error's intraclass correlation by district, the first such figure outside Rwanda; the sensitivity run.

## Predictions fixed in advance

1. The outcome is MIXED.
2. The oracle reaches 92% on at most six of nine items at `district_50`.
3. The option gives an interval on at least 95% of requests in every audit.
4. `recorded` has a lower RMSE than `linear_mi` on at least five of nine items at `clinic_25`.
5. Linear MI reaches 92% on at least seven of nine items at `clinic_50`.

## Limits known in advance

Eighteen unequal clusters, with one holding 27 clinics and three holding one, is a hard case for any cluster-robust interval, so the oracle may fall short and the benchmark is low. There is one quarter, so a clinic has one row per item, and clustering exists only by district. Facilities were chosen with probability proportional to size, which these unweighted analyses ignore, so the target is the pool's slope and not a national one. Items within a clinic are correlated, so nine items are not nine independent tests. The same assistant wrote the plan, runner and scoring, and no reviewer has read it.

## Process

The data file and dictionary enter the repository with their hashes. An engineering pilot, on separate seeds, is read for statuses and timing only. Any change it forces is listed below before the study. Hashes of the package, runner, scorer, plan and pool are checked before and after.

## Amendments
