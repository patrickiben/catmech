# Real-data audit: results and interpretation

On real recording differences, the contamination model did **not** beat simple alternatives. It kept coverage close to the full-reference
oracle and cut interval width by about 48% relative to the v0.3 default, but on ClinicalTrials.gov both linear imputation using the recorded value and generalized raking
had lower bias and lower error. On NHANES it roughly tied them. Its clear advantage in the synthetic study depended on errors that are large and separable, which these two datasets mostly
do not have. The [plan was fixed first](REAL_DATA_PLAN.md), and the only earlier look was the disclosed 300-record probe.

This is not a clinical audit. No public record-level clinical source-data audit was found. Both datasets are public, have both a recorded and a reference value on every row, and contain real differences.
Nothing was injected. Cohorts of 600 were drawn with replacement from the evaluation split, so the pool's own coefficient, using references on every row, is the known target. There were **33,800 requests:
33,613 returned intervals, 76 refusals for unsupported designs, 83 convergence refusals and 28 other errors**, the latter mainly singular raking fits at low audit fractions. The run took 2,543 seconds at 14 workers,
well over my revised estimate of 25 to 30 minutes and my original estimate of 15 to 20. An independent plain-Python check recomputed both target coefficients and the pool audit counts.

See the [evidence figure](real-data-figures/real-data-evidence.pdf).

## What the whole-pool audit shows

| | ClinicalTrials.gov | NHANES 2017-18 |
|---|---|---|
| Recorded against reference | Enrollment against participant-flow Started, log scale | Self-reported against measured weight, kg |
| Evaluation rows | 27,208 | 3,908, of which 64 lack a self-report |
| Disagreement | 14.2% (Jeffreys interval 13.8% to 14.6%) | 67.3% beyond 1 kg (65.8% to 68.8%) |
| Mean and SD of error when disagreeing | +0.244 and 0.665 | -1.01 kg and 5.67 kg |
| Rates differ across | Phase (p = 0.0001), sponsor (p < 0.0001) | Sex (p < 0.0001), race (p = 0.0003) |
| Target coefficient (pool SE) | 0.0357 (0.0036) log months per log participant | 0.127 (0.014) mmHg per kg |

ClinicalTrials.gov enrollment is usually identical to the Started total, and when it differs it is usually larger, which fits enrolled participants who never started. Early-phase trials have the largest mean difference (+0.572).
This is a finding about registry fields, not evidence that either is wrong, because both are entered by sponsors and measure slightly different things.

## Method comparison

Each entry is **bias, coverage, mean interval width, RMSE**, conditional on returned intervals. Audit 25% has 1,000 cohorts (coverage MCSE about 0.9 points). Other cells have 200 (about 1.8 points). "Returned" counts the contamination
model. `contamination_strata` was nearly identical to `contamination` in every cell and is in the [full table](real-data-study/summary.csv).

**ClinicalTrials.gov**

| Cell | Returned | Contamination | Reference normal | Linear MI | Raking | Recorded as truth | Oracle |
|---|---:|---:|---:|---:|---:|---:|---:|
| Audit 5% | 164/200 | +0.0159, 90.9%, 0.137, 0.0423 | -0.0200, 93.5%, 0.529, 0.1259 | +0.0062, 94.5%, 0.133, 0.0356 | +0.0039, 90.9%, 0.133, 0.0337 | +0.0054, 93.0%, 0.097, 0.0272 | +0.0007, 93.0%, 0.096, 0.0268 |
| Audit 10% | 194/200 | +0.0145, 90.2%, 0.125, 0.0376 | -0.0022, 94.5%, 0.347, 0.0829 | +0.0031, 94.0%, 0.113, 0.0294 | +0.0022, 88.5%, 0.106, 0.0304 | +0.0054, 93.0%, 0.097, 0.0272 | +0.0007, 93.0%, 0.096, 0.0268 |
| Audit 25% | 996/1000 | +0.0068, 91.2%, 0.106, 0.0308 | -0.0003, 92.3%, 0.205, 0.0539 | +0.0001, 91.8%, 0.102, 0.0286 | +0.0001, 91.1%, 0.103, 0.0287 | +0.0043, 91.4%, 0.097, 0.0282 | -0.0003, 91.2%, 0.096, 0.0275 |
| Audit 50% | 200/200 | +0.0048, 92.5%, 0.100, 0.0281 | -0.0043, 94.0%, 0.141, 0.0368 | +0.0011, 94.5%, 0.098, 0.0269 | +0.0010, 93.0%, 0.103, 0.0269 | +0.0054, 93.0%, 0.097, 0.0272 | +0.0007, 93.0%, 0.096, 0.0268 |
| Targeted audit | 200/200 | +0.0071, 91.5%, 0.100, 0.0288 | +0.0019, 92.0%, 0.162, 0.0414 | +0.0034, 92.0%, 0.102, 0.0285 | +0.0023, 95.5%, 0.108, 0.0284 | +0.0054, 93.0%, 0.097, 0.0272 | +0.0007, 93.0%, 0.096, 0.0268 |

**NHANES**

| Cell | Returned | Contamination | Reference normal | Linear MI | Raking | Recorded as truth | Oracle |
|---|---:|---:|---:|---:|---:|---:|---:|
| Audit 10% | 186/200 | -0.0019, 94.1%, 0.165, 0.0435 | -0.0114, 91.9%, 0.499, 0.1334 | +0.0001, 94.5%, 0.180, 0.0490 | +0.0021, 89.4%, 0.168, 0.0510 | -0.0151, 93.0%, 0.149, 0.0422 | +0.0039, 94.5%, 0.140, 0.0365 |
| Audit 25% | 1000/1000 | -0.0066, 94.5%, 0.155, 0.0408 | -0.0063, 94.2%, 0.297, 0.0785 | -0.0025, 93.9%, 0.154, 0.0409 | -0.0022, 92.6%, 0.149, 0.0415 | -0.0180, 92.7%, 0.149, 0.0422 | +0.0001, 94.9%, 0.140, 0.0357 |
| Audit 50% | 200/200 | -0.0007, 94.5%, 0.149, 0.0393 | +0.0044, 91.5%, 0.199, 0.0538 | +0.0025, 94.0%, 0.145, 0.0382 | +0.0033, 93.0%, 0.145, 0.0380 | -0.0151, 93.0%, 0.149, 0.0422 | +0.0039, 94.5%, 0.140, 0.0365 |
| Targeted audit | 200/200 | -0.0027, 94.5%, 0.154, 0.0404 | +0.0027, 95.5%, 0.230, 0.0596 | +0.0057, 95.5%, 0.156, 0.0405 | +0.0041, 93.0%, 0.151, 0.0430 | -0.0151, 93.0%, 0.149, 0.0422 | +0.0039, 94.5%, 0.140, 0.0365 |

The oracle covers only 91.2% on ClinicalTrials.gov. Its reported SE (0.0244) is below its empirical SD (0.0275), consistent with non-constant residual variance, so every method here inherits that gap from the homoskedastic interval. As
planned, coverage is judged against the oracle. On ClinicalTrials.gov, 25%, paired against the contamination model:

- Linear MI has lower MSE by 1.31e-4 (MCSE 2.2e-5) and narrower intervals by 0.0036 (MCSE 0.0003).
- Raking has lower MSE by 1.26e-4 (MCSE 2.2e-5) and narrower intervals by 0.0025 (MCSE 0.0006).
- The contamination model has bias +0.0068 (MCSE 0.0010), which is more than treating the recorded value as truth (+0.0043, MCSE 0.0009). Linear MI and raking have essentially none.

Against the v0.3 default, intervals were 48% narrower on ClinicalTrials.gov and 48% narrower on NHANES, with MSE lower by 1.96e-3 and 4.50e-3. Coverage differences (-1.1 and +0.3 points) are within Monte Carlo noise. So the main practical
message survives: discarding the recorded value is expensive. Linear MI gets the same gain with fewer assumptions.

At low audit fractions on ClinicalTrials.gov, the contamination model is weakest. At 5%, it returned 164 of 200 intervals and had bias +0.0159, against 200 of 200 and +0.0062 for linear MI. The ordering partly reverses on NHANES at 10%, where its
RMSE (0.0435) beat linear MI (0.0490) and raking (0.0510), though 14 of 200 requests returned no interval.

## Why the model did not separate errors here

The model finds errors by asking whether a recorded value looks unusual given the other variables and given the calibrated error size. Real differences here are small compared with the spread of the target given the predictors. Among unreferenced records at 25%:

| | ClinicalTrials.gov | NHANES |
|---|---:|---:|
| Actual and expected errors per cohort | 63.7 and 65.7 | 297.6 and 296.4 |
| AUC for identifying errors | 0.543 | 0.500 |
| Brier score, model against prevalence alone | 0.127 against 0.121 | 0.222 against 0.220 |
| RMSE of posterior mean against recorded value | 0.277 against 0.252 | 4.81 against 4.73 |

The error *count* is calibrated, but the model cannot tell *which* record is wrong. Its posterior means are worse than the recorded values. In the synthetic missed-flags cell, by contrast, AUC was 0.962. **Per-record error probabilities should not be
used on data like these.** The interval is less affected because most of the benefit comes from correct records either way.

The hold-out check, using no simulated truth, shows this on the first 40 cohorts at 25%. Its alarm rule stayed quiet: alarms in 2 of 38 completed checks on ClinicalTrials.gov (2 more unresolved) and 1 of 40 on NHANES. Its scores tell the real story:
Brier 0.128 against 0.118 for the stratum rate and RMSE 0.395 against 0.248 for the recorded value on ClinicalTrials.gov, and 0.220 against 0.220 and 5.45 against 4.79 on NHANES. A model with no skill on held-out references is a signal to prefer linear MI or raking.
That reading comes from this study and has not been tested prospectively.

## Bootstrap-donor matching on real data

| Dataset, audit | Standard PMM | Bootstrap-donor PMM | Reference normal |
|---|---:|---:|---:|
| ClinicalTrials.gov, 10% | +0.0004, 92.5%, 0.334, 0.0828 | +0.0021, 97.0%, 0.504, 0.0889 | -0.0022, 94.5%, 0.347, 0.0829 |
| ClinicalTrials.gov, 25% | -0.0023, 91.9%, 0.202, 0.0545 | -0.0011, 97.9%, 0.271, 0.0557 | -0.0003, 92.3%, 0.205, 0.0539 |
| ClinicalTrials.gov, 50% | -0.0037, 93.0%, 0.143, 0.0375 | -0.0056, 96.5%, 0.165, 0.0371 | -0.0043, 94.0%, 0.141, 0.0368 |
| NHANES, 10% | -0.0115, 90.9%, 0.463, 0.1277 | -0.0046, 98.0%, 0.756, 0.1387 | -0.0114, 91.9%, 0.499, 0.1334 |
| NHANES, 25% | -0.0055, 92.6%, 0.293, 0.0781 | -0.0034, 96.9%, 0.380, 0.0795 | -0.0063, 94.2%, 0.297, 0.0785 |
| NHANES, 50% | +0.0041, 94.0%, 0.202, 0.0527 | +0.0031, 95.0%, 0.228, 0.0550 | +0.0044, 91.5%, 0.199, 0.0538 |

Standard PMM did not repeat its synthetic failure here: its coverage is within about 2 points of the normal method, with similar width. The bootstrap version overcovered (96.5% to 98.0% in five of six cells) with intervals up to 52% wider. The synthetic repair is therefore
conservative on these data, and the reason to prefer it over the normal method is weak.

## Predictions, scored

| Prediction | Result |
|---|---|
| 1. ClinicalTrials.gov 25%: coverage within 2 points of oracle, narrower than reference normal, linear MI and raking | **Failed in part.** Coverage matched (91.2% for both), and it was narrower than reference normal, but wider than linear MI and raking, with higher error. |
| 2. NHANES 25%: no clear advantage over linear MI or raking | **Held.** MSE differences are within noise, and raking covered 1.9 points lower with intervals 0.006 narrower. |
| 3. Recording errors are small, so recorded-as-truth may be close to oracle | **Held for ClinicalTrials.gov** (RMSE 0.0282 against 0.0275, coverage equal within noise). **Failed for NHANES**: bias -0.0180, about 14% of the target, RMSE 18% higher, coverage 2.2 points lower. |
| 4. Targeted audit: reference normal loses coverage, contamination and raking hold | **Failed for reference normal**, which held (92.0% and 95.5%). The synthetic loss did not appear when selection was symmetric in the recorded value. The other two held. |
| 5. Bootstrap donors cover better than standard PMM | **Held**, but as overcoverage rather than a repair. |

## What this changes

- Linear MI with the recorded value should be tried before the contamination model, now available as `catmech_reference(use_recorded = TRUE)`, and raking is well grounded when selection probabilities are known. The synthetic gross-error cell shows when this ordering can reverse.
- The contamination model's niche is narrower than the synthetic study implied: errors that are large relative to the target's conditional spread, gross errors that would corrupt a linear proxy, or fallible references with a credible declared SD.
- On ClinicalTrials.gov the audit finds little to correct for this slope. That is a legitimate audit outcome.

Limits: two public datasets, neither a trial's source-data verification. The ClinicalTrials.gov reference is not independent truth. NHANES survey weights are ignored because the pool is treated as a finite population. Resampling a pool says nothing about transport to other
data. Results were not retuned after viewing. [Pools and manifest](real-data/fetch-manifest.json) are frozen with their query, retrieval time and hashes.
