# Test on real audit data: results and interpretation

On two real audits, no method won everywhere, and the contamination model did not earn a clear advantage over linear MI. Under the [rule fixed in advance](REAL_AUDIT_PLAN.md) the result is **inconclusive**: at a 25% audit of a Rwanda facility recount, the model had lower error than linear MI on two of three slopes, but it failed the coverage condition on all three and the return-rate condition on two. Three independent AI reviewers then tried to break this write-up. They found no error in the rule or the tables, and they did change how the results should be read. [Their corrections are listed below](#what-the-review-changed).

What holds after that review:

- **Point error favors staying near the reported values.** Across 20 cells the model had the lowest recovery error in 11, trusting reported values in 5 and linear MI in 4, all on ANC 4th visits.
- **Intervals favor linear MI.** Its clustered coverage was closest to 95% in 10 of 20 cells, with a worst case of 79.5% at facility-months at 25%, against 1 of 20 and 35.4% for the model. At that level neither reached the coverage of a recount.
- **Only one reported slope is clearly off.** The ANC 4th-visit gap excludes zero at the facility level under either definition of the target. The other gaps do not.
- **It is not a patient-level clinical source-data audit.** Rwanda holds facility counts, and the control comes from psychology.

There were **45,600 requests: 45,167 returned intervals, 426 convergence refusals, 3 refusals for unsupported designs and 4 other errors**. The run took 933 seconds at 14 workers, under the estimate of 21 minutes. [Scored output](real-audit-decisions.txt) comes from a separate script.

## What was checked, and how

A plain-Python check recomputes the pool targets, the published effect-size counts and the published Rwanda table: facility and exact-match counts for all eight elements, and within-range counts for six, with the other two off by one. It also recomputes the summaries and ratio limits from saved rows and re-derives Rule 2. A separate [replay](real-audit-replay/README.md) reproduces 600 requests from their seeds, across all six methods, with the runner's own functions. These checks show correspondence and repeatability. They do not show that a method is right. A reviewer's separate implementation, written from the plan, matched within Monte Carlo noise. Tables that were added after review are in [real-audit-review](real-audit-review/) and are post hoc.

## What the target is

**Every "full audit" here recounts one item and keeps the covariate as reported.** That is the estimand in the plan, because the package handles one audited target. The source recounted both sides, so the slope a real full recount would give is different, and which reported slopes look wrong depends on the choice.

| Slope | Reported | Recount of this item only (the target) | Recount of both sides | Reported minus target, facility-bootstrap 95% interval | Reported minus both-sides recount, interval | Rows disagreeing | Error SD / residual SD |
|---|---:|---:|---:|---:|---:|---:|---:|
| PNC per delivery | 1.035 | 0.973 | 1.016 | +0.062 [-0.009, 0.144] | +0.019 [-0.025, 0.069] | 101 of 200 | 0.81 |
| Syphilis tests per registration | 0.568 | 0.579 | 0.665 | -0.011 [-0.174, 0.112] | -0.097 [-0.253, 0.053] | 84 of 206 | 0.60 |
| Iron courses per registration | 1.020 | 0.894 | 0.954 | +0.126 [-0.010, 0.270] | +0.066 [-0.031, 0.192] | 56 of 186 | 0.78 |
| ANC 4th visits per registration | 0.435 | 0.256 | 0.284 | +0.179 [0.086, 0.260] | +0.151 [0.041, 0.243] | 179 of 207 | 0.78 |
| Effect sizes, SMD family, per log n | -0.094 | -0.087 | Not applicable | -0.006 | Not applicable | 41 of 147 | 0.30 |
| Effect sizes, r family, per log n | -0.055 | -0.057 | Not applicable | +0.003 | Not applicable | 34 of 202 | 0.18 |

Rwanda rows are facility-months at 69 health centers. Against the one-item target the reported slopes are 6% high on PNC, 2% low on syphilis, 14% high on iron and 70% high on ANC 4th visits. With both sides recounted they are 2% high, 15% low, 7% high and 53% high. So "biased" and "nearly right" below always mean relative to the one-item target.

The covariate was itself wrong on 15% to 23% of rows, often with the same overcount on both sides. At facility 64, reported registrations were 96, 112 and 78 against a recount of 73, 73 and 55. Reported iron copied the first set and recounted iron equaled the second. That single overcount is scored here as an iron error.

Three other features of the data matter below:

- **Copying.** The reported target equals the reported covariate on 150 of 186 iron rows and 144 of 206 syphilis rows, against 62 of 200 on PNC and 2 of 207 on ANC 4th visits.
- **Concentration.** Removing three facilities takes the gap from 0.126 to 0.001 on iron and from 0.062 to 0.004 on PNC, but only from 0.179 to 0.105 on ANC 4th visits. Removal must shrink a gap, so only the comparison across slopes means anything.
- **Recounts of zero.** Nine syphilis facility-months in five facilities have a recount of 0 with at least 10 reported, and iron has six in three facilities, such as 10, 0 and 0 recounted against 56, 57 and 54 reported at facility 40. These may be real, or may be registers left unfilled. Without them, the syphilis target is 0.684, the reported slope 0.617 and the gap -0.067.

The syphilis slope also depends heavily on a few facilities. Leaving one out moves it between 0.448 and 0.694, against 0.945 to 1.008 for PNC, which likely contributes to the +0.241 bias of audited rows alone.

## Method comparison

Each entry is **recovery RMSE, bias, clustered coverage**. Recovery RMSE is the distance from what a recount of that one item in the same cohort would give. Bias is against the pool slope. "Returned" counts the contamination model. Conditional on returned intervals.

**Rwanda facility-months** (500 cohorts per cell, coverage MCSE about 1 to 2 points)

| Analysis, audit | Returned | Recorded | Audited only | Linear MI | Raking | Contamination | Oracle bias, coverage |
|---|---:|---:|---:|---:|---:|---:|---:|
| PNC, 25% | 468/500 | 0.078, +0.059, 50.8% | 0.090, -0.015, 83.2% | 0.072, -0.007, 85.4% | 0.087, -0.018, 80.4% | 0.053, +0.022, 71.6% | -0.008, 95.2% |
| PNC, 50% | 500/500 | 0.078, +0.059, 50.8% | 0.055, -0.013, 89.8% | 0.044, -0.007, 86.0% | 0.050, -0.015, 88.6% | 0.042, +0.013, 86.0% | -0.008, 95.2% |
| Syphilis, 25% | 404/500 | 0.076, +0.003, 92.6% | 0.378, +0.241, 61.0% | 0.182, +0.059, 84.6% | 0.308, +0.197, 69.6% | 0.195, +0.085, 66.6% | +0.012, 93.2% |
| Syphilis, 50% | 494/500 | 0.076, +0.003, 92.6% | 0.231, +0.132, 75.2% | 0.106, +0.048, 85.2% | 0.151, +0.075, 82.0% | 0.089, +0.036, 83.2% | +0.012, 93.2% |
| Iron, 25% | 424/500 | 0.139, +0.128, 54.6% | 0.185, +0.028, 77.2% | 0.156, +0.014, 79.5% | 0.183, +0.027, 70.5% | 0.112, +0.093, 35.4% | +0.008, 90.2% |
| Iron, 50% | 500/500 | 0.139, +0.128, 54.6% | 0.105, +0.019, 85.6% | 0.085, +0.017, 84.6% | 0.093, +0.018, 84.4% | 0.080, +0.063, 70.4% | +0.008, 90.2% |
| ANC 4th visit, 25% | 469/500 | 0.186, +0.180, 27.6% | 0.133, +0.011, 78.8% | 0.085, +0.003, 87.2% | 0.115, +0.003, 80.4% | 0.124, +0.113, 49.7% | -0.000, 96.6% |
| ANC 4th visit, 50% | 499/500 | 0.186, +0.180, 27.6% | 0.067, +0.009, 90.0% | 0.041, -0.001, 95.2% | 0.048, +0.002, 93.4% | 0.080, +0.068, 78.0% | -0.000, 96.6% |

**Rwanda quarter totals** (independent rows, 200 cohorts per cell, coverage MCSE about 2 to 4 points)

| Analysis, audit | Returned | Recorded | Audited only | Linear MI | Raking | Contamination | Oracle bias, coverage |
|---|---:|---:|---:|---:|---:|---:|---:|
| PNC, 25% | 176/200 | 0.073, +0.057, 72.5% | 0.157, -0.002, 78.0% | 0.136, +0.022, 91.5% | 0.154, +0.007, 71.0% | 0.067, +0.036, 79.0% | -0.001, 96.0% |
| PNC, 50% | 200/200 | 0.073, +0.057, 72.5% | 0.082, -0.013, 95.5% | 0.072, -0.001, 90.5% | 0.073, -0.011, 92.5% | 0.050, +0.022, 88.5% | -0.001, 96.0% |
| Syphilis, 25% | 147/200 | 0.054, +0.037, 94.5% | 0.544, +0.387, 61.5% | 0.285, +0.156, 90.5% | 0.483, +0.357, 53.5% | 0.133, +0.078, 81.6% | +0.044, 94.0% |
| Syphilis, 50% | 198/200 | 0.054, +0.037, 94.5% | 0.280, +0.202, 80.5% | 0.110, +0.088, 90.5% | 0.189, +0.143, 84.5% | 0.080, +0.060, 90.9% | +0.044, 94.0% |
| Iron, 25% | 155/200 | 0.168, +0.158, 64.5% | 0.316, +0.097, 69.5% | 0.258, +0.033, 84.9% | 0.280, +0.064, 55.1% | 0.130, +0.111, 53.5% | +0.012, 95.0% |
| Iron, 50% | 199/200 | 0.168, +0.158, 64.5% | 0.141, +0.019, 86.0% | 0.097, +0.014, 92.0% | 0.110, +0.014, 80.5% | 0.089, +0.071, 86.4% | +0.012, 95.0% |
| ANC 4th visit, 25% | 190/200 | 0.213, +0.213, 37.0% | 0.220, +0.085, 91.5% | 0.109, -0.001, 97.0% | 0.193, +0.067, 79.0% | 0.144, +0.137, 65.3% | +0.006, 95.0% |
| ANC 4th visit, 50% | 199/200 | 0.213, +0.213, 37.0% | 0.097, +0.040, 90.0% | 0.053, -0.002, 94.0% | 0.067, +0.021, 90.0% | 0.101, +0.095, 83.9% | +0.006, 95.0% |

**Effect sizes, the control** (500 cohorts per cell)

| Analysis, audit | Returned | Recorded | Audited only | Linear MI | Raking | Contamination | Oracle bias, coverage |
|---|---:|---:|---:|---:|---:|---:|---:|
| SMD family, 25% | 463/500 | 0.026, -0.006, 95.8% | 0.249, +0.010, 96.4% | 0.069, +0.000, 96.4% | 0.172, +0.027, 90.8% | 0.034, +0.005, 94.6% | +0.001, 96.8% |
| SMD family, 50% | 500/500 | 0.026, -0.006, 95.8% | 0.108, +0.007, 97.0% | 0.031, +0.003, 96.0% | 0.053, +0.009, 94.8% | 0.024, +0.002, 96.0% | +0.001, 96.8% |
| r family, 25% | 489/500 | 0.004, +0.001, 97.2% | 0.037, -0.002, 96.6% | 0.006, -0.001, 97.0% | 0.022, -0.001, 94.0% | 0.004, -0.000, 96.3% | -0.002, 96.0% |
| r family, 50% | 500/500 | 0.004, +0.001, 97.2% | 0.020, -0.005, 93.8% | 0.003, -0.002, 95.8% | 0.007, -0.003, 93.8% | 0.002, -0.001, 96.0% | -0.002, 96.0% |

**RMSE against the pool slope**, which the plan listed and the first version left out. It does not charge a stable estimator for failing to follow the oracle's own sampling noise, which recovery RMSE does.

| Analysis, audit | Oracle | Recorded | Audited only | Linear MI | Raking | Contamination |
|---|---:|---:|---:|---:|---:|---:|
| PNC, 25% | 0.046 | 0.067 | 0.095 | 0.078 | 0.093 | 0.046 |
| Syphilis, 25% | 0.164 | 0.169 | 0.391 | 0.234 | 0.318 | 0.254 |
| Iron, 25% | 0.075 | 0.146 | 0.199 | 0.169 | 0.197 | 0.110 |
| ANC 4th visit, 25% | 0.050 | 0.192 | 0.142 | 0.096 | 0.124 | 0.133 |
| SMD family, 25% | 0.086 | 0.084 | 0.272 | 0.108 | 0.197 | 0.082 |
| r family, 25% | 0.018 | 0.018 | 0.040 | 0.018 | 0.028 | 0.018 |

Refusals are not random, and they flatter the model slightly. At 25%, other methods tend to do worse on the cohorts it refused: refused over returned recovery RMSE for audited rows alone was 0.861, 1.136, 1.258 and 1.305 across PNC, syphilis, iron and ANC 4th visits. Falling back to linear MI when the model refuses gives recovery RMSE 0.055, 0.198, 0.120 and 0.123, with clustered coverage 72.2%, 70.0%, 41.4% and 51.2%, so no conclusion changes.

## What limits the intervals

Default intervals covered only 65.4% to 87.2% even for the oracle at facility-months. Unequal variance and clustering both contribute. In the order PNC, syphilis, iron, ANC 4th visits, oracle coverage goes from 81.0%, 65.4%, 84.0%, 87.2% under the default, to 86.8%, 83.4%, 88.6%, 91.6% allowing unequal variance, to 95.2%, 93.2%, 90.2%, 96.6% clustered.

At 25%, linear MI's clustered coverage ran 8.6 to 10.7 points under the oracle. A 25% audit reaches about 17 facilities, but size alone is not the cause. Post hoc, a scratch version that resampled audited facilities in each imputation, which is not a package feature, covered 91.5%, 91.5%, 90.5% and 90.5% on the first 200 cohorts, against 84.0%, 80.5%, 81.0% and 86.5% for the package and 96.0%, 94.5%, 89.0% and 96.0% for the oracle. Its recovery RMSE was also lower (0.064, 0.162, 0.122, 0.083 against 0.069, 0.190, 0.155, 0.088). Most of linear MI's shortfall therefore seems to come from how the imputations handle uncertainty when rows are not independent; a [later planned test](CLUSTERED_SYNTHESIS.md) could not separate clustering from unequal variance. At quarter totals, with independent rows and the same audit size, linear MI covered 84.9% to 97.0%.

Audited rows alone and raking also use clustered variances on only about 17 clusters. Raking's clustered variance omits the small-sample factor used for every other method, which understates it by up to 17 points of coverage. With the factor, raking covers 85.0%, 73.8%, 74.7% and 83.4% at facility-months at 25%, and 96.6% and 98.4% on the controls. See [the table](real-audit-review/raking-like-with-like.csv).

The contamination model differs. On iron and ANC 4th visits its shortfall is mainly bias, which a variance fix cannot repair.

## Decision rules

**Rule 1: does a partial audit pay?** `linear_mi` against `recorded`, facility-months. The study's PNC verdict at 25% is what the rule returned, but it did not survive fresh seeds, because the study's 500 cohorts were a low draw. The [extension](real-audit-review/rule1-extension.csv) reran the comparison on 2,000 new cohorts per cell. It needs no MCMC and is post hoc.

| Analysis | Study, 25% | Extension, 25% | Against the pool slope, 25% | Reading at 25% | Extension, 50% |
|---|---:|---:|---:|---|---:|
| PNC | 0.923 [0.861, 0.991] | 1.009 [0.969, 1.048] | 1.136 [1.090, 1.179] | Rule said pays; a tie on fresh seeds, a cost against the pool | 0.622 [0.598, 0.646] |
| Syphilis | 2.393 [2.172, 2.626] | 2.577 [2.446, 2.715] | 1.469 [1.422, 1.520] | Costs | 1.385 [1.310, 1.462] |
| Iron | 1.125 [1.040, 1.216] | 1.153 [1.089, 1.234] | 1.205 [1.144, 1.278] | Costs | 0.557 [0.538, 0.579] |
| ANC 4th visit | 0.458 [0.430, 0.487] | 0.443 [0.428, 0.458] | 0.503 [0.486, 0.521] | Pays | 0.234 [0.225, 0.242] |

The percentile intervals reflect cohort resampling only. How much worse linear MI is on syphilis depends on what is held fixed: 2.4 in the study, 2.6 on fresh seeds, about 1.4 to 1.5 against the pool slope, 1.8 without the zero recounts, and 17.7 if the pool is fixed and only the audit is redrawn, because reported error is then just the gap of 0.011.

**Rule 2: does the contamination model add anything?** Outcome: **inconclusive**. Each condition must hold for an analysis to count. At least two had to count for "supported", and at least two had to show a lower limit above 1 for "contradicted".

| Analysis, 25% | Ratio against linear MI | Returned, need 90% | Clustered coverage, oracle | Quarter ratio | Counts? |
|---|---:|---:|---:|---:|---|
| PNC | 0.745 [0.690, 0.802] | 93.6% | 71.6%, 95.2% | 0.485 | No, coverage |
| Syphilis | 1.111 [0.991, 1.242] | 80.8% | 66.6%, 93.2% | 0.529 | No, error, returns and coverage |
| Iron | 0.717 [0.656, 0.792] | 84.8% | 35.4%, 90.2% | 0.499 | No, returns and coverage |

The stress case, ANC 4th visit, went against the model: 1.488 [1.400, 1.579] at 25% and 1.948 [1.835, 2.077] at 50%. The direction is stable. Against the pool slope the ratios are 0.584, 1.088, 0.646 and 1.405, and with the pool fixed they are 0.588, 0.973, 0.638 and 1.464, in the order PNC, syphilis, iron, ANC 4th visits. No audit-using method as implemented met the coverage condition on any of the three slopes, so the verdict says less than the coverage results above.

**Rule 3: what does auditing cost when there is little to correct?** Ratio against `recorded` at 25%, recovery RMSE:

| Analysis | Linear MI | Raking | Audited only | Contamination |
|---|---:|---:|---:|---:|
| SMD family | 2.705 [2.393, 3.040] | 6.706 [6.014, 7.441] | 9.742 [8.902, 10.741] | 1.306 [1.177, 1.444] |
| r family | 1.578 [1.435, 1.728] | 5.853 [5.233, 6.499] | 9.950 [9.036, 10.908] | 0.950 [0.888, 1.022] |

Against the pool slope, linear MI's cost is much smaller: 1.291 [1.211, 1.376] and 1.025 [0.992, 1.061]. The model stays near the reported values here, with coverage close to the oracle's. That means it did least harm where there was little to correct, not that it helped.

## Why the model behaved this way: post hoc

None of this was planned, so these are explanations tested after viewing, not findings of the study. They do not support one single cause.

**ANC 4th visits: errors grow with facility size.** The correlation of error with the covariate is +0.477, and the gap stays large without any three facilities. An additive error with one mean cannot remove that. Shuffling errors within catchment before each cohort, which removes the link, restores the model's clustered coverage to 93.9% (oracle 97.5%) with bias -0.003. Allowing separate error means by facility-size group with `error_strata` cut recovery RMSE from 0.124 to 0.071, bias from +0.111 to +0.031, and raised coverage from 51.9% to 86.0% on 129 paired cohorts, close to linear MI (0.086, +0.010, 88.4%). So this result is partly about configuration. The documentation did not suggest that setting before this study.

**PNC and iron: that explanation does not fit.** Nearly all of each gap sits in three facilities, and most wrong iron records copy the covariate, so they sit near the regression line and look like correct ones. With errors shuffled, the model's coverage stays at 74.9% and 63.1%, against 92.0% and 92.5% for the oracle, with standard errors about 0.6 of the empirical SD. Size groups changed almost nothing (recovery RMSE 0.051 against 0.052, and 0.110 against 0.110).

**Syphilis: no tracking at all.** The correlation within catchment is -0.017. The model still added bias (+0.085 against +0.003 for recorded), and with errors shuffled it covered 76.3% against 97.5%. Size groups made it worse (0.251 against 0.179).

So the model's intervals are too narrow even when errors do not track a covariate, and the cause on these three slopes is not identified. It would be wrong to read these results as saying errors unrelated to the other variables are safe for the model.

**Linear MI's low average bias hides single-audit misses.** Its mean error on iron is +0.111 when the audit misses all three key facilities and -0.091 when it reaches at least one, and on PNC +0.036 and -0.045. The average is near zero, but a single audit under-corrects or overshoots. On ANC 4th visits the split is narrower (+0.041 and -0.025). For the model the same figures are +0.119 and +0.067 on iron, +0.041 and +0.006 on PNC, and +0.131 and +0.100 on ANC 4th visits. See the [table](real-audit-review/conditional-on-key-facilities.csv).

**The reported bias the model keeps** is 38% on PNC, 73% on iron and 63% on ANC 4th visits at 25%, measured against the pool slope, or 45%, 72% and 62% measured against the same cohort's one-item recount. At 50% the figures are 23%, 50% and 38%, or 32%, 46% and 38%.

**Per-record output carried little information.** At 25%, AUC was 0.501, 0.525, 0.528 and 0.582. Expected error counts ran low on PNC (43.7 against 76.8), iron (26.6 against 42.0) and syphilis (54.1 against 64.0), and were close on ANC 4th visits and the controls. Posterior means were worse than the recorded values on three of four analyses (5.376 against 4.713 on PNC) and better on ANC 4th visits (7.553 against 9.787).

**The package's warnings do not cover this.** The warning that nearly all calibration rows disagree fired on 189 of 500 ANC 4th-visit cohorts at 25%, and never on iron at facility-months, where coverage was worst.

**Checking errors against the covariate on audited rows has low power.** The check is a natural idea, and it can help, but at a 25% audit a clustered test reached significance on only 17.8% of cohorts on PNC, 18.8% on iron and 56.6% on ANC 4th visits, with the wrong sign on 19.0%, 28.0% and 5.0%. A clear positive is informative. A null is not reassurance.

## Post hoc: two simple comparators

**Swap-in** corrects the audited records and leaves the rest. It is what an audit team would naturally do, and the plan omitted it. Rebuilding the cohorts from their seeds matched 15,200 saved oracle and recorded estimates to within 1e-7. See the [summary](real-audit-posthoc/swap-in-summary.csv).

| Facility-months | Recorded | Swap-in | Linear MI | Contamination, on its returned cohorts | Swap-in, same cohorts | Reported bias kept by contamination, swap-in |
|---|---:|---:|---:|---:|---:|---:|
| PNC, 25% | 0.0775 | 0.0634 | 0.0716 | 0.0533 | 0.0626 | 38%, 74% |
| Syphilis, 25% | 0.0762 | 0.0676 | 0.1824 | 0.1949 | 0.0694 | Reported bias near zero |
| Iron, 25% | 0.1387 | 0.1141 | 0.1562 | 0.1117 | 0.1133 | 73%, 79% |
| ANC 4th visit, 25% | 0.1856 | 0.1447 | 0.0849 | 0.1241 | 0.1446 | 63%, 75% |
| PNC, 50% | 0.0775 | 0.0484 | 0.0436 | 0.0419 | 0.0484 | 23%, 48% |
| Syphilis, 50% | 0.0762 | 0.0510 | 0.1056 | 0.0887 | 0.0510 | Reported bias near zero |
| Iron, 50% | 0.1387 | 0.0824 | 0.0853 | 0.0798 | 0.0824 | 50%, 55% |
| ANC 4th visit, 50% | 0.1856 | 0.1039 | 0.0413 | 0.0805 | 0.1037 | 38%, 49% |

Swap-in had lower recovery RMSE than recorded in all 20 cells, but that is largely arithmetic: its error is the cohort's reported gap minus the audited rows' share of it. With the pool fixed and only the audit redrawn, it beat recorded in 14 of 20 cells and lost on all four syphilis cells and both SMD cells. It is a baseline any method should have to beat where a gap is material, not a correction: it keeps most of the bias at 25%, and it gives no interval.

**A design-weighted difference estimator** uses the same proxy as raking and the known probabilities. Its recovery RMSE was lower than linear MI's in 15 of 20 cells and never more than 5% higher, while raking's was higher in all 20. So the raking column reflects this implementation, whose auxiliaries were narrowed after the pilot. It is not a verdict on design-based correction. *Correction, added later:* this estimator's correction term is zero by the normal equations, so it is least squares on fitted values and not a difference estimator. See the [later comparison](OUTSIDE_METHODS_SYNTHESIS.md#a-correction-to-an-earlier-write-up).

## Predictions, scored

| Prediction | Result |
|---|---|
| 1. Linear MI beats recorded at 25% on ANC 4th visits, not on PNC | **Half held as scored.** ANC 4th visit held at 0.458. PNC failed in the study's cohorts at 0.923 [0.861, 0.991], then held on 2,000 fresh cohorts at 1.009 [0.969, 1.048]. |
| 2. Contamination beats linear MI at 25% on PNC, syphilis and iron, not on ANC 4th visits | **Three of four held.** It won on PNC and iron and lost on ANC 4th visits. Syphilis failed, 1.111 [0.991, 1.242]. |
| 3. Default intervals undercover at facility-months; clustered are closer; contamination within 3 points of the oracle | **First two held, third failed badly.** Contamination's gap was 23.6 to 54.8 points. |
| 4. Recorded has the lowest recovery RMSE at 25% on both controls | **Held for SMD.** For r, contamination was marginally lower (0.0035 against 0.0037), a statistical tie at 0.950 [0.888, 1.022]. |
| 5. Raking has wider intervals than linear MI at 25% | **Held at facility-months** on clustered intervals (4 of 4; default 3 of 4) and on both controls. **Failed at quarter totals** (0 of 4), though raking's missing small-sample factor makes that comparison uneven: with it, 2 of 4. |
| 6. Per-record AUC stays below 0.75 at facility-months | **Held**, 0.501 to 0.619. |
| 7. Audit-using methods improve from 25% to 50% | **Held** in every cell on recovery RMSE. |

## What the review changed

Three AI reviewers worked independently: one rebuilt the results from raw data before reading any conclusions, one attacked the statistics, and one audited claims against the sources. No blocker was found, and every table cell reproduced. The first version of this document had these problems, each checked again before being changed here:

1. It called the target "a full audit" and sorted slopes into biased and nearly right without saying only one item is recounted.
2. It stated as fact that errors tracking facility size explained the model's failure. That fits ANC 4th visits and does not fit the other slopes.
3. It said linear MI had the lowest error among methods that removed the bias, which is false on PNC, and it left out linear MI's cost on iron and the cancellation behind its low average bias.
4. It reported Rule 1 on PNC as "pays" without testing whether that survived fresh seeds, which it does not, and it omitted RMSE against the pool, which the plan listed.
5. It said about 17 audited facilities could not support nominal intervals, when imputation design explains much of the shortfall.
6. It presented swap-in's 20 of 20 as a finding, when it is largely arithmetic.
7. It printed one slope as 0.569 instead of 0.568, said the calibrated error count failed everywhere, and described its checks as more independent than they were. A reviewer showed they would pass relabeled results, which is why the replay now exists.

## What this changes

- **The claim that the model is strongest when errors are large stays synthetic-only.** Rwanda's errors are moderate, with few extreme rows, and not separable, so gross errors remain untested on real data. On these moderate real errors, the model often had lower point error, but its intervals failed.
- **Linear MI stays the recommended start in a narrower form.** When a slope gap is plausible and an interval is needed, its average bias was near zero on the slopes with a gap and it had the best worst-case coverage. These results do not support starting there for accuracy: at a 25% audit its point error beat reported values clearly only on ANC 4th visits. In a single audit on PNC or iron it under-corrected or overshot.
- **Clustered imputation was the most useful lead, and a later test qualified it.** Resampling audited facilities in each imputation came much closer to the oracle's coverage on these slopes, in a post-hoc scratch result. Its [planned test](CLUSTERED_SYNTHESIS.md) on held-out items came out PARTIAL: it helped where linear MI undercovered and cost accuracy elsewhere.
- **A larger audit matters at least as much as the choice of method.** At 50%, linear MI beats reported values on three of four slopes, the methods move closer together and error falls in every cell.

## Limits

One health system, with 69 health centers that were not randomly selected. The recount was done in consultation with facility reporting staff, and a register can itself be wrong, which matters most for the rows with recounts of zero. Part of the ANC 4th-visit gap may reflect how gestational age was handled rather than recording slips. Resampling 69 facilities says little about transport to other settings. The analyses were chosen after seeing descriptive slopes, as the plan discloses.

The control has its own limits. The source oversampled outliers on purpose (42% to 43% of these rows), which it corrected for and this test does not, and that inflates the residual spread that makes errors look small: the error ratio without outliers is 0.36 and 0.42, against 0.30 and 0.18. The source's kappa of 0.71 describes assigning a row to a category, not the reference values. One g row just over the plan's tolerance of 0.049 counts as wrong although the source classed it as reproducible, which changes nothing material.

On process: the plan's first line says it was written before any model ran, and its later section reports pilots that ran those models on the same pools. Only statuses and timing are said to have been read from them. The plan was not under version control before this test, so "fixed in advance" rests on file hashes and timestamps taken at the run. `catmech_holdout_check()` was not part of the study. AI review is not independent human peer review. No package code changed.
