# Validation record: catmech 0.4.1

Date: 2026-11-01. This is an experimental software release with independent implementation checks, prespecified synthetic studies and a public real-data
audit. It is not clinical validation, a general coverage guarantee, or a claim to new imputation methodology. The [v0.3.1 record](v0.3.1-VALIDATION.md) and earlier ones are preserved.

## Main findings

- **[Synthetic study](CONTAMINATION_SYNTHESIS.md), 52,950 requests.** The contamination model passed the coverage screen in three of four primary cells and failed narrowly in the fourth (93.7%, Wilson lower limit 0.920). Intervals were 30% to 38% narrower than the v0.3 default.
  It failed under heavy-tailed targets (87.3%), selection on the true value with errors present (79.3%) and undeclared reference noise (88.7%). These primary cells use the error form the model assumes. See the [figure](contamination-figures/contamination-evidence.pdf).
- **[Real-data audit](REAL_DATA_SYNTHESIS.md), 33,800 requests.** On ClinicalTrials.gov, linear MI and raking beat the model on bias and error. On NHANES they roughly tied. Against the v0.3 default, intervals were about 48% narrower on both. The model could not identify which records were wrong. Two of five advance predictions held, two failed in part or for one dataset, and one
  failed for its main claim. See the [figure](real-data-figures/real-data-evidence.pdf). Neither dataset is a clinical source-data audit.
- **[Go/no-go simulation](merged-audit-pilot/README.md), 28,800 records.** NO-GO for transmitter-targeted auditing under spillover.
- **[Audit planner test](AUDIT_PLAN_SYNTHESIS.md), 14,400 target fits and 7,200 plans.** NO-GO under a rule fixed in advance. Median error in predicting width was 11.3% to 16.8% at a 5% audit, against a 10% limit, and at most 7.8% elsewhere. The bootstrap interval for `k` undercovered. This is new in v0.4.1 and did not change any earlier result.
- **[Test on real audit data](REAL_AUDIT_SYNTHESIS.md), 45,600 requests.** A Rwanda health-facility recount, with a psychology effect-size audit as a control. The model's advantage over linear MI was inconclusive under a rule fixed in advance: lower error on two of three slopes at a 25% audit, but it failed the coverage condition on all three. Where a reported slope differed from a recount of that one item it kept 38% to 73% of the gap, with clustered coverage of 35% to 72% against 90% to 97% for the oracle. Linear MI's coverage held up better, though at a 25% audit its point error clearly beat trusting reported values on only one of four slopes. Neither audit holds patient values. No package code changed.

## Software and independent checks

| Check | Result | Evidence |
|---|---|---|
| Unit suite | 598 assertions passed, no warnings or skips | `v0.4.1-checks/tests.log` |
| Source build and package check | 0 errors, 0 warnings, 0 notes | `v0.4.1-checks/build.log`, `v0.4.1-checks/R-CMD-check.log` |
| Fresh archive installation | Namespace-only APIs and all five bundled examples ran from an isolated library | `v0.4.1-checks/installed-smoke.log` |
| Contamination posterior oracle | Independent Metropolis sampler on the marginal posterior, sharing no package code: 38/38 moments within five MCSEs (largest 1.96), 36/36 row quantities within 0.015, three wrong likelihoods detected, 15/15 contract controls | `contamination-oracle/` |
| Older workflows unchanged | 18 of 18 runs identical to frozen v0.3.1 in completions, audits, hashes and pooled results, rerun with v0.4.1 | `v0.4.1-checks/v0.3-unchanged-check.txt` |
| Post-study additions leave frozen runs alone | 12 of 12 runs identical to frozen v0.4.0 sources, rerun with v0.4.1 | `v0.4.1-checks/v0.4.0-postfreeze-check.txt` |
| New `use_recorded` option | Matches the study's linear MI on the same real-data cohorts in RMSE, coverage and width | `recorded-option-check/comparison.csv` |
| Serial and parallel pilots | Output files byte-identical between 1 and 14 workers | `contamination-pilot-check.txt` |
| Supplement for a runner bug | 1,731 score records and 120 hold-out records recovered; every estimate matched the main study to within 1e-9 | `contamination-supplement/` |
| Independent evidence recomputation | Plain Python recomputes all summaries, Wilson intervals and paired MCSEs for every study, both real-data target coefficients and pool audit counts, and the go/no-go verdicts, including the audit planner's summaries, pool `k` values and browser test vectors; deleting a planned cell is detected | `v0.4.1-checks/artifact-verification.log`, `tools/verify-artifacts.py` |
| Real-audit replay | 600 requests across all six methods, replayed from seeds with the runner's own functions, match the saved study; relabeled or rescaled results are caught | `real-audit-replay/`, `tools/verify-artifacts.py` |
| Review of the real-audit test | Three independent AI reviewers found no error in the rule or tables and corrected several readings, listed in the write-up. This is not human peer review. | `REAL_AUDIT_SYNTHESIS.md`, `real-audit-review/` |
| CRT port | Reproduces deposited naive results closely; adjusted RMSE 0.158 against 0.149 and coverage 0.946 against 0.963 | `merged-audit-pilot/reproduce/port-check.csv` |
| GitHub Actions | Three-platform workflow supplied; not run remotely | `.github/workflows/R-CMD-check.yaml` |

The local environment is R 4.6.0, mice 3.19.0 and posterior 1.7.0 on Apple silicon macOS. No new package dependency was added. Exact session records accompany each study. Minimum-version R, other platforms and a fresh-container restore have not been tested.

## Frozen studies and post-study changes

Plans preceded their runs: [synthetic](CONTAMINATION_PLAN.md), [real-data](REAL_DATA_PLAN.md) and [go/no-go](merged-audit-pilot/PLAN.md). Development runs and the changes they caused are disclosed in those plans. All three studies ran on the sources preserved in `v0.4.0-source/`, with hashes checked
before and after each run.

Two scoped changes came after those runs, as they did for v0.3.1: `catmech_reference(use_recorded = TRUE)` and the `advice` field of `catmech_holdout_check()`. The studies did not run either one directly. The first was checked against the study's linear MI on the same cohorts. The second is a reading of the real-data study and has not been tested prospectively.
The study runner also had a bookkeeping bug that dropped auxiliary records in six cells, recovered by the supplement above. Main estimates were unaffected.

Run times exceeded my estimates: the synthetic study took 40 minutes (estimated about 30), the real-data study 42 (estimated 15 to 20, later 25 to 30), the go/no-go 37 (estimated about 15) and the audit planner test 15 (estimated about 12).

Version 0.4.1 adds the audit planner and changes no existing function. Its [plan](AUDIT_PLAN_PLAN.md) preceded its evaluation run, which ran on the sources preserved in `v0.4.1-source/`. Development choices used the development split and are listed in the plan, along with the frozen tables already seen by hand.

The [real-audit test](REAL_AUDIT_PLAN.md) came after that release and ran on package sources identical to `v0.4.1-source/`. Its plan discloses the descriptive looks that preceded it and two changes made after its pilots, before any estimate was seen. It took 16 minutes against an estimate of 21. Its swap-in comparison and the tables in `real-audit-review/` are post hoc, are labeled as such and are not part of any rule. The plan was not under version control before that run, so its timing rests on file hashes and timestamps.

## Evidence limits

Real-data pools are public, and the ClinicalTrials.gov reference is sponsor-entered. The two public audits are facility counts and psychology effect sizes. No patient-level clinical source-data audit was available. Simulated error rates have not been checked against published audits, and the main fixture's 8% may be high. Hold-out checks had low power. Clustered records, non-constant variance, survival and interaction models remain limits.
Adversarial or AI-assisted review is not independent human refereeing. AI use is disclosed in `CONTRIBUTING.md`.

## Reproduce

```sh
Rscript -e 'testthat::test_local()'
R CMD build .
R CMD check --no-manual catmech_0.4.1.tar.gz
Rscript tools/check-contamination-oracle.R validation/my-contamination-oracle
Rscript tools/contamination-benchmark.R study validation/my-contamination-study
python3 -I tools/real-data-fetch.py my-raw-directory validation/my-real-data
Rscript tools/real-data-prepare.R my-raw-directory validation/my-real-data
Rscript tools/real-data-benchmark.R study validation/my-real-data-study
Rscript tools/merged-audit-pilot.R study validation/merged-audit-pilot/my-study 400
Rscript tools/audit-plan-benchmark.R study validation/my-audit-plan-study
Rscript tools/real-audit-prepare.R validation/real-audit/raw validation/my-real-audit
Rscript tools/real-audit-benchmark.R study validation/my-real-audit-study
python3 -I tools/real-audit-score.py my-decisions.txt validation/my-real-audit-study
Rscript tools/real-audit-posthoc.R validation/my-real-audit-posthoc
Rscript tools/real-audit-review-checks.R validation/my-real-audit-review
Rscript tools/real-audit-replay.R validation/my-real-audit-replay
python3 tools/verify-artifacts.py
```

Use new output directories. The real-audit benchmark always reads the pools in `validation/real-audit/`, and the post-hoc, review and replay scripts read the saved study. To reproduce the exact fingerprints, run from the frozen sources in `v0.4.0-source/`, or `v0.4.1-source/` for the audit planner test. A fresh ClinicalTrials.gov download will differ from the dated pool. Set `CATMECH_BENCHMARK_WORKERS` to change parallelism. Exact cross-platform numerical equality is not promised.
