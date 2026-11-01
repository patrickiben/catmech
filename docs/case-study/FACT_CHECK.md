# Fact check

Written by [`fact_check.py`](fact_check.py). Last run: 2026-11-01. Each row recomputes the number from its source, and most sources are frozen results tables that the script only reads. A document fails if it contains a digit or number word that no row covers, or if a row's text no longer appears in it. Rerun with `python3 docs/case-study/fact_check.py`.

**Result: all 83 claims match their sources.** PDF: 1 page, generated from the current `CASE_STUDY.md`. Relative links and README tool paths checked: 72, 0 broken. External DOI links were not checked, because this check stays offline. Negative controls: a changed number, an added number, an added number word, a wrong claim value, an unnegated publication claim, a third-party reference and broken links were all caught.

## README.md front page

| As shown | What it says | Source file | Result |
|---|---|---|---|
| `0.4.1` | Package version | `DESCRIPTION` | Match |
| `25%` | Audit fraction: mean referenced rows over cohort size, ClinicalTrials.gov cell | `validation/real-data-study/summary.csv` | Match |
| `25%` | Audit fraction: mean referenced rows over cohort size, NHANES cell | `validation/real-data-study/summary.csv` | Match |
| `0.0286` | RMSE, linear MI, ClinicalTrials.gov, 25% | `validation/real-data-study/summary.csv` | Match |
| `0.0308` | RMSE, contamination model, ClinicalTrials.gov, 25% | `validation/real-data-study/summary.csv` | Match |
| `+0.0001` | Bias, linear MI, ClinicalTrials.gov, 25% | `validation/real-data-study/summary.csv` | Match |
| `+0.0068` | Bias, contamination model, ClinicalTrials.gov, 25% | `validation/real-data-study/summary.csv` | Match |
| `0.0409` | RMSE, linear MI, NHANES, 25% | `validation/real-data-study/summary.csv` | Match |
| `0.0408` | RMSE, contamination model, NHANES, 25% | `validation/real-data-study/summary.csv` | Match |
| `10%` | Audit fraction: mean referenced rows over cohort size, NHANES 10% cell | `validation/real-data-study/summary.csv` | Match |
| `0.0435` | RMSE, contamination model, NHANES, 10% | `validation/real-data-study/summary.csv` | Match |
| `0.0490` | RMSE, linear MI, NHANES, 10% | `validation/real-data-study/summary.csv` | Match |
| `14 of 200` | Model requests with no interval, NHANES, 10%: requested minus returned | `validation/real-data-study/summary.csv` | Match |
| `0.54` | Mean AUC for finding wrong records, ClinicalTrials.gov, 25% | `validation/real-data-study/corruption-scores.csv` | Match |
| `0.50` | Mean AUC for finding wrong records, NHANES, 25% | `validation/real-data-study/corruption-scores.csv` | Match |
| `about half` | Linear MI width over discard width, 25%: must be 0.45 to 0.55 in both datasets | `validation/real-data-study/summary.csv` | Match |
| `about twice` | Discard width over linear MI width, 25%: must be 1.8 to 2.2 in both datasets | `validation/real-data-study/summary.csv` | Match |
| `0.102` | Mean width, linear MI, ClinicalTrials.gov, 25% | `validation/real-data-study/summary.csv` | Match |
| `0.205` | Mean width, discarding recorded values, ClinicalTrials.gov, 25% | `validation/real-data-study/summary.csv` | Match |
| `0.154` | Mean width, linear MI, NHANES, 25% | `validation/real-data-study/summary.csv` | Match |
| `0.297` | Mean width, discarding recorded values, NHANES, 25% | `validation/real-data-study/summary.csv` | Match |
| `91.8%` | Coverage, linear MI, ClinicalTrials.gov, 25% | `validation/real-data-study/summary.csv` | Match |
| `91.2%` | Coverage, all references, ClinicalTrials.gov, 25% | `validation/real-data-study/summary.csv` | Match |
| `93.9%` | Coverage, linear MI, NHANES, 25% | `validation/real-data-study/summary.csv` | Match |
| `94.9%` | Coverage, all references, NHANES, 25% | `validation/real-data-study/summary.csv` | Match |
| `64.7%` | Coverage, linear MI, synthetic gross-and-small errors | `validation/contamination-study/summary.csv` | Match |
| `94.0%` | Coverage, contamination model, synthetic gross-and-small errors | `validation/contamination-study/summary.csv` | Match |
| `91.8%` | Coverage, contamination model, synthetic noisy reference with declared SD | `validation/contamination-study/summary.csv` | Match |
| `87.3%` | Coverage, contamination model, synthetic heavy-tailed target | `validation/contamination-study/summary.csv` | Match |
| `79.3%` | Coverage, contamination model, synthetic selection on true value with errors | `validation/contamination-study/summary.csv` | Match |
| `7.3%` | Coverage, reference normal, synthetic selective references | `validation/contamination-study/summary.csv` | Match |
| `52.7%` | Coverage, reference normal, synthetic noisy references | `validation/contamination-study/summary.csv` | Match |
| `10%` | Planner's limit on median width error | `validation/AUDIT_PLAN_SYNTHESIS.md` | Match |
| `5%` | Smallest target audit in the planner test | `validation/audit-plan-study/summary.csv` | Match |
| `11.3%` | Lowest planner median error at the 5% target, of four real-data cells | `validation/audit-plan-study/summary.csv` | Match |
| `16.8%` | Highest planner median error at the 5% target, of four real-data cells | `validation/audit-plan-study/summary.csv` | Match |
| `7.8%` | Largest planner median error at any larger target, real-data cells | `validation/audit-plan-study/summary.csv` | Match |
| `NO-GO` | Planner verdict | `validation/audit-plan-study/verdict.csv` | Match |
| `NO-GO` | Merged-idea verdict | `validation/merged-audit-pilot/study/verdict.csv` | Match |
| `set.seed(1)` | Example parameter, not a finding. Run by this check. | `README.md` | Match |
| `nrow(d) / 4` | Example parameter: a 25% audit. Run by this check. | `README.md` | Match |
| `agreement_tol = 1` | Example parameter: the 1 kg tolerance the frozen NHANES audit used | `validation/REAL_DATA_SYNTHESIS.md` | Match |
| `m = 20` | Imputations per fit in the real-data study | `validation/real-data-study/settings.json` | Match |
| `seed = 123` | Example parameter, not a finding. Run by this check. | `README.md` | Match |
| `0.126` | Example output: recorded_kg estimate, one random audit | `README.md, by running it` | Match |
| `95% interval` | Default level of catmech_pool_lm | `R/analysis.R` | Match |
| `0.097` | Example output: lower limit | `README.md, by running it` | Match |
| `0.155` | Example output: upper limit | `README.md, by running it` | Match |
| `0.127` | NHANES evaluation pool's coefficient using every measured weight | `validation/real-data-study/pool-targets.csv` | Match |
| `1,000` | Cohorts at a 25% audit | `validation/real-data-study/summary.csv` | Match |
| `600` | Records per cohort | `validation/real-data-study/settings.json` | Match |
| `v0.3` | Figure label for the old default, which discards recorded values | `docs/ASSESSMENT.md` | Match |
| `two public datasets` | Real-data pools in the study | `validation/audit-plan-study/pool-k.csv` | Match |
| `within Monte Carlo noise` | Paired MSE z for the model, NHANES 10%: must be between -2 and 0, and the most favorable cell | `validation/real-data-study/paired-comparisons.csv` | Match |

## docs/case-study/CASE_STUDY.md and PDF

| As shown | What it says | Source file | Result |
|---|---|---|---|
| `0.4.1` | Package version | `DESCRIPTION` | Match |
| `November 1, 2026` | Status date | `docs/ASSESSMENT.md` | Match |
| `Version 0.4` | Release that added the recording-error model | `NEWS.md` | Match |
| `25%` | Audit fraction: mean referenced rows over cohort size, ClinicalTrials.gov cell | `validation/real-data-study/summary.csv` | Match |
| `25%` | Audit fraction: mean referenced rows over cohort size, NHANES cell | `validation/real-data-study/summary.csv` | Match |
| `0.0286` | RMSE, linear MI, ClinicalTrials.gov, 25% | `validation/real-data-study/summary.csv` | Match |
| `0.0308` | RMSE, contamination model, ClinicalTrials.gov, 25% | `validation/real-data-study/summary.csv` | Match |
| `0.0539` | RMSE, discarding recorded values, ClinicalTrials.gov, 25% | `validation/real-data-study/summary.csv` | Match |
| `0.0275` | RMSE, all references, ClinicalTrials.gov, 25% | `validation/real-data-study/summary.csv` | Match |
| `0.0409` | RMSE, linear MI, NHANES, 25% | `validation/real-data-study/summary.csv` | Match |
| `0.0408` | RMSE, contamination model, NHANES, 25% | `validation/real-data-study/summary.csv` | Match |
| `0.0785` | RMSE, discarding recorded values, NHANES, 25% | `validation/real-data-study/summary.csv` | Match |
| `0.0357` | RMSE, all references, NHANES, 25% | `validation/real-data-study/summary.csv` | Match |
| `0.54` | Mean AUC for finding wrong records, ClinicalTrials.gov, 25% | `validation/real-data-study/corruption-scores.csv` | Match |
| `0.50` | Mean AUC for finding wrong records, NHANES, 25% | `validation/real-data-study/corruption-scores.csv` | Match |
| `91.8%` | Coverage, linear MI, ClinicalTrials.gov, 25% | `validation/real-data-study/summary.csv` | Match |
| `91.2%` | Coverage, all references, ClinicalTrials.gov, 25% | `validation/real-data-study/summary.csv` | Match |
| `91.2%` | Coverage, contamination model, ClinicalTrials.gov, 25% | `validation/real-data-study/summary.csv` | Match |
| `93.9%` | Coverage, linear MI, NHANES, 25% | `validation/real-data-study/summary.csv` | Match |
| `94.5%` | Coverage, contamination model, NHANES, 25% | `validation/real-data-study/summary.csv` | Match |
| `94.9%` | Coverage, all references, NHANES, 25% | `validation/real-data-study/summary.csv` | Match |
| `64.7%` | Coverage, linear MI, synthetic gross-and-small errors | `validation/contamination-study/summary.csv` | Match |
| `94.0%` | Coverage, contamination model, synthetic gross-and-small errors | `validation/contamination-study/summary.csv` | Match |
| `10%` | Planner's limit on median width error | `validation/AUDIT_PLAN_SYNTHESIS.md` | Match |
| `11.3%` | Lowest planner median error at the 5% target, of four real-data cells | `validation/audit-plan-study/summary.csv` | Match |
| `16.8%` | Highest planner median error at the 5% target, of four real-data cells | `validation/audit-plan-study/summary.csv` | Match |
| `7.8%` | Largest planner median error at any larger target, real-data cells | `validation/audit-plan-study/summary.csv` | Match |
| `30 of 600 records` | Records at the 5% target, cohort size | `validation/audit-plan-study/settings.json` | Match |
| `NO-GO` | Planner verdict | `validation/audit-plan-study/verdict.csv` | Match |
| `NO-GO` | Merged-idea verdict | `validation/merged-audit-pilot/study/verdict.csv` | Match |
| `m = 20` | Imputations per fit in the real-data study | `validation/real-data-study/settings.json` | Match |
| `1,000` | Cohorts at a 25% audit | `validation/real-data-study/summary.csv` | Match |
| `600` | Records per cohort | `validation/real-data-study/settings.json` | Match |
| `two public datasets` | Real-data pools in the study | `validation/audit-plan-study/pool-k.csv` | Match |
| `30%` | Smallest interval reduction, model against old default, four primary synthetic cells | `validation/contamination-study/summary.csv` | Match |
| `38%` | Largest interval reduction, model against old default, four primary synthetic cells | `validation/contamination-study/summary.csv` | Match |
| `three of four` | Primary synthetic cells whose Wilson interval is inside the screen 0.925 to 0.975 | `validation/contamination-study/summary.csv` | Match |
| `in the fourth` | Site-shift cell is the one that failed the screen | `validation/contamination-study/summary.csv` | Match |
| `nine real-data scenarios` | Paired contrasts of the model against linear MI | `validation/real-data-study/paired-comparisons.csv` | Match |
| `never clearly better` | Most favorable paired MSE z for the model across the nine cells must be above -2 | `validation/real-data-study/paired-comparisons.csv` | Match |
| `nearly halved` | RMSE reduction, linear MI against discarding, 25%: must be 40% to 50% in both datasets | `validation/real-data-study/summary.csv` | Match |
| `95% ranges` | Nominal level of the intervals scored for coverage | `validation/CONTAMINATION_PLAN.md` | Match |
| `two of five` | Advance real-data predictions that held | `validation/VALIDATION.md` | Match |
| `Two more ideas` | Merged idea and planner, both NO-GO | `validation/merged-audit-pilot/study/verdict.csv` | Match |
| `27,208` | ClinicalTrials.gov evaluation pool rows | `validation/audit-plan-study/pool-k.csv` | Match |
| `3,908` | NHANES evaluation pool rows | `validation/audit-plan-study/pool-k.csv` | Match |
| `1.31e-4` | Paired MSE difference, model minus linear MI, ClinicalTrials.gov, 25% | `validation/real-data-study/paired-comparisons.csv` | Match |
| `2.2e-5` | MCSE of that difference | `validation/real-data-study/paired-comparisons.csv` | Match |

## docs/case-study/BLURBS.md

| As shown | What it says | Source file | Result |
|---|---|---|---|
| `November 1, 2026` | Status date | `docs/ASSESSMENT.md` | Match |
| `25%` | Audit fraction: mean referenced rows over cohort size, ClinicalTrials.gov cell | `validation/real-data-study/summary.csv` | Match |
| `NO-GO` | Merged-idea verdict | `validation/merged-audit-pilot/study/verdict.csv` | Match |
| `two public datasets` | Real-data pools in the study | `validation/audit-plan-study/pool-k.csv` | Match |
| `Two other ideas` | Merged idea and planner, both NO-GO | `validation/merged-audit-pilot/study/verdict.csv` | Match |

## README.md lower sections, numbers added or moved in this rewrite

| As shown | What it says | Source file | Result |
|---|---|---|---|
| `86.7%` | Moved text: synthetic standard reference PMM coverage, 25% references | `validation/contamination-study/summary.csv` | Match |
| `95.9%` | Moved text: synthetic bootstrap-donor PMM coverage, 25% references | `validation/contamination-study/summary.csv` | Match |
| `51%` | Moved text: bootstrap-donor widths against reference normal on real data, largest increase. REAL_DATA_SYNTHESIS.md says 52%, from rounded table values | `validation/real-data-study/summary.csv` | Match |
| `about 380 of 600` | New text: mean records by formula when discarding, both pools (383 and 378) | `validation/audit-plan-study/pool-k.csv` | Match |
| `about 40 to 70` | New text: records by formula when keeping recorded values (41 and 66) | `validation/audit-plan-study/pool-k.csv` | Match |

## Statements checked by grep

| Statement | Source file | Result |
|---|---|---|
| Neither real dataset is a clinical source-data audit | `validation/VALIDATION.md` | Supported |
| On ClinicalTrials.gov both values are sponsor-entered, so the reference is not independent truth | `validation/REAL_DATA_SYNTHESIS.md` | Supported |
| Simulation data are synthetic, and the primary cells favor the model | `validation/CONTAMINATION_SYNTHESIS.md` | Supported |
| No new statistical method is claimed | `docs/ASSESSMENT.md` | Supported |
| Not on CRAN, with no outside peer review or clinical validation | `docs/GITHUB.md` | Supported |
| Linear MI is the recommended start and a supported option | `docs/ASSESSMENT.md` | Supported |
| The model could not tell which records were wrong | `validation/REAL_DATA_SYNTHESIS.md` | Supported |
| Every study had a written plan first, with earlier looks disclosed | `validation/VALIDATION.md` | Supported |
| Results were not retuned after viewing | `validation/AUDIT_PLAN_SYNTHESIS.md` | Supported |
| Coverage was judged against the oracle | `validation/REAL_DATA_SYNTHESIS.md` | Supported |
| Do not write the merged paper on that thesis | `docs/ASSESSMENT.md` | Supported |
| Merged idea: audits helped, targeting did not beat the most suspicious records, direction added nothing | `docs/ASSESSMENT.md` | Supported |
| The planner is not validated | `validation/AUDIT_PLAN_SYNTHESIS.md` | Supported |
| The planner test was not fully blind | `validation/AUDIT_PLAN_SYNTHESIS.md` | Supported |
| A bookkeeping bug and missed run-time estimates are on the record | `validation/VALIDATION.md` | Supported |
| Plain Python recomputes the reported summaries | `validation/VALIDATION.md` | Supported |
| AI use and its limits | `CONTRIBUTING.md` | Supported |
| Not published: this repository has no git remote | `git remote` | Supported |

## Excluded as not findings

| Document | Text | Why |
|---|---|---|
| README.md front page | This one audit | Describes the single illustrative run in the example |
| README.md front page | # 1. | Step number in a code comment |
| README.md front page | # 2. | Step number in a code comment |
| docs/case-study/CASE_STUDY.md and PDF | two poor options | Framing of the problem, not a measurement |
| docs/case-study/CASE_STUDY.md and PDF | until one passed | Ordinary phrase, not a count |
| docs/case-study/BLURBS.md | a checked one | Ordinary phrase, not a count |
| docs/case-study/BLURBS.md | One-line | Heading from the brief |
| docs/case-study/BLURBS.md | Three-line | Heading from the brief |
| docs/case-study/BLURBS.md | Two-sentence | Heading from the brief |

## The README example

The second R block on the front page ran as written, loading the package from this tree, with no warnings. `recorded_kg` estimate 0.1261, interval 0.0973 to 0.1549. The seeds were chosen before the first run and were not changed. The numbers are one illustrative audit, so another seed gives other values, and exact cross-platform equality is not promised.

## Paired check behind "never clearly better"

Paired MSE z, contamination model minus linear MI, from `validation/real-data-study/paired-comparisons.csv`. Negative favors the model. Conditional on returned intervals. This reading of the frozen table is descriptive and was not a planned test.

| Dataset | Audit | z |
|---|---|---:|
| ctgov | audit_05pct | +3.99 |
| ctgov | audit_10pct | +4.83 |
| ctgov | audit_25pct | +5.96 |
| ctgov | audit_50pct | +2.41 |
| ctgov | targeted_audit | +0.52 |
| nhanes | audit_10pct | -1.13 |
| nhanes | audit_25pct | -0.03 |
| nhanes | audit_50pct | +1.65 |
| nhanes | targeted_audit | -0.07 |

## Corrections from the fact check

1. Round 1. README and case study said the model's "one clear win" was the 10% NHANES audit. The paired table shows that edge is within Monte Carlo noise (z = -1.13), so both now say it was not clear. The marginal RMSE values were right, and the reading of them was too strong.
2. Round 1. The case study said the planner "predicted precision within 7.8%". That number is a median error, and the 90th percentile reached 20.6%. It now says median error.
3. Round 1. The first draft of the README example used `simulate_catmech_lab()`, whose default errors are gross, which is the documented failure case for linear MI. It now uses the public NHANES pool, with the study's 1 kg tolerance in the audit summary. Its comment on `fit$audit` was also wrong and was fixed: that table has one row per record and imputation.
4. Round 1. README text moved from the old version said bootstrap-donor intervals were up to 52% wider. Raw widths give 51.3% (0.7557 over 0.4994, NHANES 10%). The README now says 51% and names the comparator. `validation/REAL_DATA_SYNTHESIS.md` still says 52%, from rounded table values, because this task does not own it.
5. Round 2. Direct grep of every distinct numeric claim in the written records (the syntheses, validation record, assessment and planning guide), a path independent of the table recomputation. Of 65, 55 are printed there. The other 10 are the date format, example parameters and output, two phrases containing 95%, the 30-of-600 phrase and the corrected 51%. Nothing further needed fixing.
