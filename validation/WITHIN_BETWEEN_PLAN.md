# Prospective test: separate within-cluster and between-cluster slopes in the option

Written 2026-11-02 and committed before the runner, the scorer, the new simulation generator or any analysis of the new data exists, so version history records the order. Amendments, if any, go in the section at the end.

## Why

[A fresh agent's setup](INDEPENDENT_SETUPS_SYNTHESIS.md), `ind_mice`, reached 92% on all 11 real items at 25% and nearly removed the bias with the audited value as a predictor. Separate within-cluster and between-cluster slopes are its visible difference from `catmech_reference(cluster = )`, but it also differs in other ways, including a heuristic shrinkage of the imputation variance for an imputed predictor. Whether the slopes themselves help has not been tested. Rwanda, Zambia and the existing simulated cells have now all been seen, so they cannot be held out.

## The candidate

`within_between` is the option, unchanged, with cluster means added to `predictors`: the mean within each cluster, over all its rows, of the recorded value and of every predictor that varies within clusters. This is the same model as separate slopes, written another way. There is no shrinkage, no strata pooling and no change to the draws, the pooling or the analysis. It can be run with the package as it stands, so no package code changes for this test.

Means aliased on the audited rows are dropped. If the audited clusters would leave fewer than 3 degrees of freedom between clusters, all means are dropped for that request, so it is the option itself. Such requests are counted and reported.

## Comparators

`random_intercept` (the option), `linear_mi`, `ind_mice` as committed, `oracle` and `recorded`, each with the cluster-robust CV3 variance through the package, or `ind_mice`'s own analysis.

## Data that decides

**SNAP.** USDA Supplemental Nutrition Assistance Program Quality Control Database, fiscal year 2020, period 1, October 2019 through February 2020, doi:10.15482/USDA.ADC/1528542, U.S. Public Domain. The file has 18,319 household reviews from 53 state-level units with 104 to 477 reviews each. Period 2 was collected under pandemic waivers and is not used. No analysis of this file has been run, and only its columns, row counts, review statuses and cluster sizes have been read.

- Recorded value: `RAWBEN`, the benefit the unit was certified to receive. Reference: what the reviewer determined after re-interviewing the household, which is `RAWBEN` when `STATUS` is 1, `RAWBEN - AMTERR` when `STATUS` is 2 and `RAWBEN + AMTERR` when `STATUS` is 3. Rows are dropped when `CASE` is not 1, or when `STATUS` and `AMTERR` disagree. Dollar amounts are in hundreds.
- Analyses: `benefit_on_size`, benefit on `CERTHHSZ`; `benefit_on_gross`, benefit on `RAWGROSS`; and `gross_on_benefit`, `RAWGROSS` on benefit and `CERTHHSZ`, with the audited value as a predictor. The imputation model's predictors are the analysis's other variables.
- Cohorts: draw 53 states with replacement and then k reviews without replacement within each, k of 10 or 40. Audit whole states, 13 or 26 of the 53. The target is the slope on the whole pool using the reference value, with each state weighted equally, to match this sampling.
- That is 12 cells with 1,000 cohorts each.

**New simulation.** One fresh agent writes a generator from a brief, seeing none of the repository, the old generator, any results, or the candidate, and is told neither the candidate nor this rule. It runs only its own checks that the targets are right. Its code enters the repository as written, before the study. It must provide at least 24 cells with known targets, at least one third with the audited value as a predictor, at 25% and 50% whole-cluster audits, with 500 cohorts each.

Each cell runs all the methods above at m = 20 with a 60-second limit per request, and one over it counts as no interval.

## Data that guards

These have been seen, so they can only block adoption. Candidate results are run on the same cohorts and compared with the saved results for the option at the same settings.

- The 11 real Rwanda items at 25% and 50%.
- The existing 36 outcome cells and 6 predictor cells.
- All four Zambia audits on nine items.

## The rule

Within each deciding source, cells split into two classes: audited value as outcome, and as predictor. For a class, against the option, using the earlier scorer's measures:

- `ratio`: geometric mean across cells of the paired RMSE ratio.
- `distance`: mean across cells of the absolute distance of coverage from 95%.
- `width`: geometric mean across cells of the ratio of mean interval widths.
- `bias`: mean across cells of the absolute relative bias.
- `available`: share of requests giving an interval.

A class is **no worse** if `ratio` is at most 1.02, `distance` is at most the option's plus 0.25, and `available` is at least 99%. It is **better** if it is no worse, `width` is at most 1.25, and either `distance` is at most the option's minus 0.25 or `ratio` is at most 0.98. A predictor class must also have `bias` no larger than the option's to be better.

A source **passes** if both classes are no worse and at least one is better.

The guard **holds** if all of these are met, against the option:

- G1: real items at 25%. `ratio` at most 1.02, at least as many items at 92% or more, at least 10 of 11 at 98% or less, and `distance` at most the option's plus 0.25.
- G2: real items at 50%. The same four checks.
- G3: existing outcome cells. `ratio` at most 1.02, and at least as many cells within limits as the option, with limits as in the earlier scorer.
- G4: existing predictor cells. `ratio` at most 1.02, and `bias` no larger than the option's.
- G5: Zambia, each audit separately. `ratio` at most 1.02, `distance` at most the option's plus 0.5, and `available` at least the option's minus 1 point.

**ADOPT** if both deciding sources pass and the guard holds. The model replaces the option's in the next version, through a new commit whose results must replay this study's candidate requests. The option stays experimental until someone independent reviews it.

**NOT ADOPTED** otherwise. The option stays as it is, and the documents report what happened.

Thresholds do not move after results are read. The scorer checks itself against saved values. As a control on the cohorts, the runner reruns `ind_mice` on 110 real-item requests and they must match the saved results exactly.

**Description only:** each cell for every method, requests that fell back to the option, seconds per request, the error's intraclass correlation in SNAP, and cohort-bootstrap intervals for the differences in `distance` and `ratio`. Comparisons with `ind_mice` have no rule, since it is not the method being tested.

## Predictions fixed in advance

1. The outcome is NOT ADOPTED.
2. On SNAP, the outcome class has a `ratio` between 0.98 and 1.02, since most reviews find no error.
3. On the new simulation, the candidate's predictor class has `bias` below the option's.
4. On the new simulation, `ind_mice` has a smaller predictor-class `bias` than the candidate.

## Limits known in advance

SNAP is one paired item with three analyses, in dollars and not in health, and it leaves out cases found ineligible and overissuances that reach the full benefit, which removes some of the largest errors. Its covariates were recorded during the same review. A slope there may not transfer to facility counts. FNS re-reviews a subsample, but this file holds the state reviewer's finding and not a gold standard, and the error amount may reflect the reviewer's own mistakes. The new generator's agent shares a model with the author, and the author writes its brief. No human has reviewed this. The same assistant writes the plan, runner, scorer and write-up.

## Amendments

Made before the study ran and before any result was read.

**The generator's independence is weaker than planned.** The agent reports that it did not open the repository and fit only ordinary least squares on full true values. But its session loaded an index of the author's notes without being asked, and the one-line summaries name catmech results and methods, including `ind_mice`, linear MI and a PARTIAL verdict on the cluster option. They do not state the candidate or this rule. It reports that it did not design toward them and cannot show that they had no influence. So this source is not as clean as this plan says. The generator stands, because a replacement would load the same index, and the condition is disclosed in place of a rerun.

**Size.** The generator has 62 cells, 33 outcome and 29 predictor, against the required minimum of 24, so the new simulation makes 372,000 requests with both audit sizes. Six cells at 25% audit too few clusters for `ind_mice`, which then refuses. That has no effect on the rule, which judges the candidate and the option.

