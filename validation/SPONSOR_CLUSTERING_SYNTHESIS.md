# Sponsor clustering screen: results and interpretation

**Verdict: WEAK.** ClinicalTrials.gov discrepancies cluster within sponsors well beyond chance, and a sponsor's revealed trials help predict its other trials. Adaptive auditing did not find more than plain risk ranking, so the third condition failed. The script computed this verdict on D in the single evaluation run, under the [plan](SPONSOR_CLUSTERING_PLAN.md) frozen beforehand.

| Condition for CLUSTERING-SUPPORTED | Required | Result | Met |
|---|---|---|---|
| Test 1 p-value | below 0.01 | 0.0005 | Yes |
| AUC difference, history minus baseline | at least 0.02 | +0.088 | Yes |
| AUC difference, 95% lower limit | above 0 | 0.049 | Yes |
| Adaptive-to-risk yield ratio | at least 1.2 | 0.967 | **No** |
| Yield ratio, 95% lower limit | above 1 | 0.854 | **No** |

Source for all tables: [`sponsor-clustering/evaluation/summary.csv`](sponsor-clustering/evaluation/summary.csv), by its `section`, `outcome`, `subgroup` and `metric` columns, unless another file is named.

This tests clustering only. It cannot show errors spreading, because shared definitions, reporting habits, CROs, SOPs or training could all produce the same pattern with no spread. D is a difference between two sponsor-entered registry fields, not an error rate, and this does not cover clinical source-data errors.

## Data

- All 34,084 pool rows matched a sponsor name, a 100% match rate (27,208 of 27,208 evaluation rows). The names came from the same API filter on 2026-11-01 at 20:02 UTC, about four hours after the pool was fetched, in 36 pages and 5.4 MB. The query, page hashes and match rates are in [`fetch-manifest.json`](sponsor-clustering/fetch-manifest.json).
- Trimming and case-folding left 5,488 distinct names in the pool, the same number as before, so they merged nothing. Evaluation rows hold 4,817 sponsors.
- Variants of one organization remain split. In [`sponsor-map.csv`](sponsor-clustering/sponsor-map.csv), 5 distinct names contain "novartis", 3 "pfizer" and 12 "merck", though some of these may be separate organizations. One name can also cover many independent groups.
- The evaluation pool has 3,859 discrepancies on D (14.2%, matching the earlier audit) and 1,446 on D2 (5.3%).

| Set | Sponsors | Trials | D | D2 |
|---|---:|---:|---:|---:|
| Evaluation pool | 4,817 | 27,208 | 3,859 | 1,446 |
| Test 1, at least 3 trials | 1,254 | 22,912 | 3,315 | 1,219 |
| Test 2 fit: unrevealed trials of training sponsors | 466 | 5,327 | 776 | 283 |
| Test 2 scored: unrevealed trials of test sponsors | 465 | 5,864 | 870 | 335 |
| Test 3 population: all trials of test sponsors | 465 | 11,512 | 1,725 | 628 |

Test 1 has 127 nonempty strata. Across Test 2, 10,752 trials were revealed. Budget in Test 3 is 1,151 trials.

## Test 1: excess clustering

| Outcome | Statistic | Null mean (middle 95% of null) | Ratio | p |
|---|---:|---:|---:|---:|
| D | 2,110.4 | 1,048.4 (964.0 to 1,133.6) | 2.013 | 0.0005 |
| D2, descriptive | 1,561.7 | 1,171.5 (1,051.4 to 1,304.7) | 1.333 | 0.0005 |

Both p-values are at the floor of 1/2,001, so no permutation of 2,000 reached either statistic. Per-permutation values are in [`test1-null.csv`](sponsor-clustering/evaluation/test1-null.csv).

Post hoc description from the `term_D` column of [`sponsors.csv`](sponsor-clustering/evaluation/sponsors.csv): the 10 largest contributors make up 26.5% of the D statistic, and the 47 sponsors with at least 100 evaluation trials make up 25.8%. So a few large names do not account for most of it. This was not planned and is not a test.

## Test 2: prediction within held-out sponsors

| Outcome | Model | AUC (95% interval) | Brier (95% interval) |
|---|---|---:|---:|
| D | Baseline | 0.585 (0.548 to 0.617) | 0.1253 (0.1095 to 0.1422) |
| D | History | 0.674 (0.626 to 0.710) | 0.1194 (0.1062 to 0.1334) |
| D | **Difference** | **+0.088 (0.049 to 0.126)** | **-0.0059 (-0.0116 to -0.0013)** |
| D2, descriptive | Difference | +0.012 (-0.026 to 0.046) | -0.0001 (-0.0005 to 0.0003) |

History uses the logit of the shrunken revealed rate; its coefficient on D is 0.608 ([`coefficients.csv`](sponsor-clustering/evaluation/coefficients.csv)). All fits, D and D2, converged in at most 8 Newton steps. No bootstrap replicate was undefined. Intervals resample sponsors but hold one fixed reveal, so they leave out reveal-to-reveal variation.

## Test 3: audit yield

Discrepancies found per 100 audited, D:

| Strategy | Yield (95% interval) |
|---|---:|
| Random | 14.9 (12.6 to 17.7) |
| Risk | 22.1 (16.4 to 29.3) |
| Adaptive | 21.3 (15.5 to 27.7) |

| Ratio | Value (95% interval) |
|---|---:|
| **Adaptive to risk** | **0.967 (0.854 to 1.046)** |
| Adaptive to random | 1.431 (1.175 to 1.638) |
| Risk to random | 1.480 (1.216 to 1.743) |

Both strategies beat random, but adaptive did not beat risk. The interval's upper limit of 1.046 is well below the required 1.2. Per-seed counts are in [`test3-seeds.csv`](sponsor-clustering/evaluation/test3-seeds.csv) and replicate values in [`bootstrap.csv`](sponsor-clustering/evaluation/bootstrap.csv). Across seeds, the adaptive yield has a standard error of 0.08.

Adaptive's random stage flagged a mean of 54.6 of 465 sponsors, and its second stage never ran short of trials in flagged sponsors. That second stage alone found 27.7 per 100, 1.26 times the risk yield, which the random half then pulled down. This is post hoc arithmetic. The second stage depends on the random stage for its flags, no other split was run, and Test 3 does not use the history model. Any change to these would be a new plan and a new run, not a reading of this one.

## Descriptive results outside the verdict

**D2.** Excess clustering is present on Test 1 but smaller, with an excess over the null of 390 against 1,062 for D. History added no detectable prediction. Yields were 5.4 (4.6 to 6.3) for random, 8.5 (5.4 to 11.6) for risk and 8.9 (5.3 to 12.6) for adaptive, with an adaptive-to-risk ratio of 1.047 (0.803 to 1.281).

**By sponsor class.** Test 3 by class has point estimates only, as planned.

| Outcome | Measure | Industry | Other |
|---|---|---:|---:|
| D | Test 1 ratio, p | 2.374, 0.0005 | 1.587, 0.0005 |
| D | Sponsors, trials in Test 1 | 689, 9,838 | 565, 13,074 |
| D | AUC difference (95% interval) | +0.155 (0.056 to 0.245) | +0.037 (0.011 to 0.061) |
| D | Brier difference (95% interval) | -0.0141 (-0.0280 to -0.0023) | -0.0012 (-0.0033 to 0.0009) |
| D | Yields: random, risk, adaptive | 16.2, 27.8, 26.4 | 14.3, 19.8, 18.3 |
| D | Adaptive to risk | 0.951 | 0.928 |
| D2 | Test 1 ratio, p | 1.569, 0.0005 | 1.047, 0.274 |
| D2 | AUC difference (95% interval) | +0.029 (-0.062 to 0.109) | -0.001 (-0.034 to 0.029) |
| D2 | Adaptive to risk | 1.216 | 0.913 |

Clustering and the history gain are both stronger among industry sponsors. On D2, the other class shows no excess clustering at all. The D2 industry ratio of 1.216 has no interval, rests on 418 audits and is not evidence for adaptive auditing.

## What failed

- **The audit-yield condition failed on both parts.** Adaptive found slightly fewer discrepancies than risk ranking, and the interval excludes the required gain.
- **The AUC gain did not turn into an audit gain.** A real sponsor signal exists on D, but this adaptive strategy did not outperform ranking on trial features alone.
- **D2 did not support prediction.** Larger differences do not show the prediction gain that any difference does, which fits a pattern driven in part by small reporting differences. This data cannot establish that cause.
- **Part of the prediction missed**, as scored below.

## Prediction, scored

| Written before running | Result | Held |
|---|---|---|
| Verdict WEAK, about 75% | WEAK | Yes |
| Test 1 p at or near the floor, ratio about 2 to 4 | 0.0005, 2.013 | Yes, at the low end |
| Baseline AUC about 0.60 to 0.68 | 0.585 | **No**, lower |
| History adds about 0.03 to 0.08 | +0.088 | **No**, slightly higher |
| AUC lower limit above 0, small negative Brier difference | 0.049, -0.0059 | Yes |
| Yields: random about 14, risk 20 to 28, adaptive 18 to 25 | 14.9, 22.1, 21.3 | Yes |
| Adaptive to risk about 0.8 to 1.1, failing 1.2 | 0.967 | Yes |
| D2 clustering of the same kind, weaker in absolute count | Excess 390 against 1,062 | Yes, but the lost prediction gain was not predicted |
| Industry ratio above other class, low confidence | 2.374 against 1.587 | Yes |

## Interpretation limits

- WEAK means clustering beyond chance without enough audit gain to justify adaptive targeting in this registry setting. Here the prediction condition passed and the audit-yield condition failed.
- CLUSTERING-SUPPORTED would only have justified looking for record-level audit data with shared-resource metadata, and would have said nothing about direction or spreading. This result does not provide that justification. NO-CLUSTERING did not occur either, so the registry route has not been closed by a null.
- The result does not cover clinical source-data errors. Registry entry is not source-data verification, and enrolled and started can legitimately differ, so D is not an error rate.
- Clustering here can come from shared definitions, reporting habits, CROs, SOPs or training with no spread at all.
- Lead sponsor is the only shared factor tested. It is not a monitor, clerk, lab or system, and unmerged names can misstate the real unit.
- The pool holds only completed interventional trials with posted results.
- Adaptive here has one fixed form: half the budget at random, flags set once and a second stage ranked by the baseline model. Other forms are untested.

## Run and checks

- **Evaluation run:** once, 36.9 seconds ([`run-info.json`](sponsor-clustering/evaluation/run-info.json)). Plan, pool, sponsor map and script hashes were unchanged before and after ([`input-hashes.csv`](sponsor-clustering/evaluation/input-hashes.csv)).
- **Development run:** engineering only, 22 seconds, in [`sponsor-clustering/development/`](sponsor-clustering/development/). It also returned WEAK, with an adaptive-to-risk ratio of 1.156 (0.906 to 1.708). It was seen before the evaluation run. No code changed after the first development run, and both saved runs carry the same script hash.
- **Independent recomputation:** `tools/verify-artifacts.py` rebuilds in plain Python, from the pool and sponsor map, all 200 summary values of each run. That covers Test 1's strata and statistic, the train-test split and reveal, the logistic fits, AUC and Brier, every per-seed audit count, the interval limits from saved replicates, the first 25 bootstrap replicates' AUC, Brier and risk yield, and the verdict with a table of fixed cases.
- **Negative controls in the verifier:** a known-wrong Test 1 statistic is detected. In addition, 200 within-strata shuffles of evaluation D give a mean ratio of 0.999 to the saved null, with 3 of 200 below p = 0.01. The output is in [`verifier.log`](sponsor-clustering/verifier.log), and the whole verifier exits 0.
- **Not recomputed independently:** two things use NumPy streams that plain Python cannot repeat. For the saved Test 1 null draws, the verifier recomputes the p-value and ratio from them and compares their mean with its own 200 shuffles, to within 4 standard errors. For the adaptive and random yields inside bootstrap replicates, it checks only that their means sit within one standard deviation of the point estimates. Those replicates' audit seeds use the zero-based replicate number, while `bootstrap.csv` numbers replicates from 1. Passing these checks shows correspondence and repeatability, not inferential validity.
