# Validation record: catmech 0.4.0

Date: 2026-11-01. This is an experimental software release with independent implementation checks, prespecified synthetic studies and a public real-data
audit. It is not clinical validation, a general coverage guarantee, or a claim to new imputation methodology. The [v0.3.1 record](v0.3.1-VALIDATION.md) and earlier ones are preserved.

## Main findings

- **[Synthetic study](CONTAMINATION_SYNTHESIS.md), 52,950 requests.** The contamination model passed the coverage screen in three of four primary cells and failed narrowly in the fourth (93.7%, Wilson lower limit 0.920). Intervals were 30% to 38% narrower than the v0.3 default.
  It failed under heavy-tailed targets (87.3%), selection on the true value with errors present (79.3%) and undeclared reference noise (88.7%). These primary cells use the error form the model assumes. See the [figure](contamination-figures/contamination-evidence.pdf).
- **[Real-data audit](REAL_DATA_SYNTHESIS.md), 33,800 requests.** On ClinicalTrials.gov, linear MI and raking beat the model on bias and error. On NHANES they roughly tied. Against the v0.3 default, intervals were about 48% narrower on both. The model could not identify which records were wrong. Two of five advance predictions held, two failed in part or for one dataset, and one
  failed for its main claim. See the [figure](real-data-figures/real-data-evidence.pdf). Neither dataset is a clinical source-data audit.
- **[Go/no-go simulation](merged-audit-pilot/README.md), 28,800 records.** NO-GO for transmitter-targeted auditing under spillover.

## Software and independent checks

| Check | Result | Evidence |
|---|---|---|
| Unit suite | 515 assertions passed, no warnings or skips | `v0.4.0-checks/tests.log` |
| Source build and package check | 0 errors, 0 warnings, 0 notes | `v0.4.0-checks/build.log`, `R-CMD-check.log` |
| Fresh archive installation | Namespace-only APIs and all four bundled examples ran from an isolated library | `v0.4.0-checks/installed-smoke.log` |
| Contamination posterior oracle | Independent Metropolis sampler on the marginal posterior, sharing no package code: 38/38 moments within five MCSEs (largest 1.96), 36/36 row quantities within 0.015, three wrong likelihoods detected, 15/15 contract controls | `contamination-oracle/` |
| Older workflows unchanged | 18 of 18 runs identical to frozen v0.3.1 in completions, audits, hashes and pooled results | `v0.3-unchanged-check.txt` |
| Post-study additions leave frozen runs alone | 12 of 12 runs identical to frozen v0.4.0 sources | `v0.4.0-postfreeze-check.txt` |
| New `use_recorded` option | Matches the study's linear MI on the same real-data cohorts in RMSE, coverage and width | `recorded-option-check/comparison.csv` |
| Serial and parallel pilots | Output files byte-identical between 1 and 14 workers | `contamination-pilot-check.txt` |
| Supplement for a runner bug | 1,731 score records and 120 hold-out records recovered; every estimate matched the main study to within 1e-9 | `contamination-supplement/` |
| Independent evidence recomputation | Plain Python recomputes all summaries, Wilson intervals and paired MCSEs for every study, both real-data target coefficients and pool audit counts, and the go/no-go verdict; deleting a planned cell is detected | `v0.4.0-checks/artifact-verification.log`, `tools/verify-artifacts.py` |
| CRT port | Reproduces deposited naive results closely; adjusted RMSE 0.158 against 0.149 and coverage 0.946 against 0.963 | `merged-audit-pilot/reproduce/port-check.csv` |
| GitHub Actions | Three-platform workflow supplied; not run remotely | `.github/workflows/R-CMD-check.yaml` |

The local environment is R 4.6.0, mice 3.19.0 and posterior 1.7.0 on Apple silicon macOS. No new package dependency was added. Exact session records accompany each study. Minimum-version R, other platforms and a fresh-container restore have not been tested.

## Frozen studies and post-study changes

Plans preceded their runs: [synthetic](CONTAMINATION_PLAN.md), [real-data](REAL_DATA_PLAN.md) and [go/no-go](merged-audit-pilot/PLAN.md). Development runs and the changes they caused are disclosed in those plans. All three studies ran on the sources preserved in `v0.4.0-source/`, with hashes checked
before and after each run.

Two scoped changes came after those runs, as they did for v0.3.1: `catmech_reference(use_recorded = TRUE)` and the `advice` field of `catmech_holdout_check()`. The studies did not run either one directly. The first was checked against the study's linear MI on the same cohorts. The second is a reading of the real-data study and has not been tested prospectively.
The study runner also had a bookkeeping bug that dropped auxiliary records in six cells, recovered by the supplement above. Main estimates were unaffected.

Run times exceeded my estimates: the synthetic study took 40 minutes (estimated about 30), the real-data study 42 (estimated 15 to 20, later 25 to 30) and the go/no-go 37 (estimated about 15).

## Evidence limits

Real-data pools are two public datasets, and the ClinicalTrials.gov reference is sponsor-entered. No clinical source-data audit was available. Simulated error rates have not been checked against published audits, and the main fixture's 8% may be high. Hold-out checks had low power. Clustered records, non-constant variance, survival and interaction models remain limits.
Adversarial or AI-assisted review is not independent human refereeing. AI use is disclosed in `CONTRIBUTING.md`.

## Reproduce

```sh
Rscript -e 'testthat::test_local()'
R CMD build .
R CMD check --no-manual catmech_0.4.0.tar.gz
Rscript tools/check-contamination-oracle.R validation/my-contamination-oracle
Rscript tools/contamination-benchmark.R study validation/my-contamination-study
python3 -I tools/real-data-fetch.py my-raw-directory validation/my-real-data
Rscript tools/real-data-prepare.R my-raw-directory validation/my-real-data
Rscript tools/real-data-benchmark.R study validation/my-real-data-study
Rscript tools/merged-audit-pilot.R study validation/merged-audit-pilot/my-study 400
python3 tools/verify-artifacts.py
```

Use new output directories. To reproduce the exact fingerprints, run from the frozen sources in `v0.4.0-source/`. A fresh ClinicalTrials.gov download will differ from the dated pool. Set `CATMECH_BENCHMARK_WORKERS` to change parallelism. Exact cross-platform numerical equality is not promised.
