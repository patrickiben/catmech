# Validation record: catmech 0.7.1

Date: 2026-11-02. This is an experimental software release with independent implementation checks, prespecified synthetic studies and a public real-data
audit. It is not clinical validation, a general coverage guarantee, or a claim to new imputation methodology. The [v0.7.0 record](v0.7.0-VALIDATION.md), [v0.6.0 record](v0.6.0-VALIDATION.md), [v0.5.0 record](v0.5.0-VALIDATION.md), [v0.4.1 record](v0.4.1-VALIDATION.md), [v0.3.1 record](v0.3.1-VALIDATION.md) and earlier ones are preserved.

## Main findings

- **[Synthetic study](CONTAMINATION_SYNTHESIS.md), 52,950 requests.** The contamination model passed the coverage screen in three of four primary cells and failed narrowly in the fourth (93.7%, Wilson lower limit 0.920). Intervals were 30% to 38% narrower than the v0.3 default.
  It failed under heavy-tailed targets (87.3%), selection on the true value with errors present (79.3%) and undeclared reference noise (88.7%). These primary cells use the error form the model assumes. See the [figure](contamination-figures/contamination-evidence.pdf).
- **[Real-data audit](REAL_DATA_SYNTHESIS.md), 33,800 requests.** On ClinicalTrials.gov, linear MI and raking beat the model on bias and error. On NHANES they roughly tied. Against the v0.3 default, intervals were about 48% narrower on both. The model could not identify which records were wrong. Two of five advance predictions held, two failed in part or for one dataset, and one
  failed for its main claim. See the [figure](real-data-figures/real-data-evidence.pdf). Neither dataset is a clinical source-data audit.
- **[Go/no-go simulation](merged-audit-pilot/README.md), 28,800 records.** NO-GO for transmitter-targeted auditing under spillover.
- **[Audit planner test](AUDIT_PLAN_SYNTHESIS.md), 14,400 target fits and 7,200 plans.** NO-GO under a rule fixed in advance. Median error in predicting width was 11.3% to 16.8% at a 5% audit, against a 10% limit, and at most 7.8% elsewhere. The bootstrap interval for `k` undercovered. This is new in v0.4.1 and did not change any earlier result.
- **[Test on real audit data](REAL_AUDIT_SYNTHESIS.md), 45,600 requests.** A Rwanda health-facility recount, with a psychology effect-size audit as a control. The model's advantage over linear MI was inconclusive under a rule fixed in advance: lower error on two of three slopes at a 25% audit, but it failed the coverage condition on all three. Where a reported slope differed from a recount of that one item it kept 38% to 73% of the gap, with clustered coverage of 35% to 72% against 90% to 97% for the oracle. Linear MI's coverage held up better, though at a 25% audit its point error clearly beat trusting reported values on only one of four slopes. Neither audit holds patient values. No package code changed.
- **[Test of imputation for clustered audits](CLUSTERED_SYNTHESIS.md), 148,000 requests.** **PARTIAL** under a rule fixed in advance. On five held-out items of the same Rwanda recount at a 25% audit, `catmech_reference(cluster = )` moved coverage closer to the oracle's on all five and within 5 points on four, but only two of five met the accuracy condition. It helped where the existing option fell short and raised error elsewhere. It passed the simulation condition in 33 of 36 cells, but those simulations barely contain the clustering the option is for. One independent AI reviewer confirmed the verdict, found no bug and corrected several readings; its leads on rare events, within-facility correlation and stratified resampling are post hoc. This version of the option was later replaced.
- **[Test of replacements for the clustered resampling](CLUSTERED_FIX_SYNTHESIS.md), 210,000 requests.** **PARTIAL** under a new rule fixed in advance. Simulations with clustered residuals selected random-intercept imputation over a stratified jackknife draw, by a tie-break on error. On nine items of the same recount at a 25% audit, coverage was within its limits on eight and 89.0% on the ninth. Error was lower than the option without `cluster` on six common items and 4% to 5% higher on three rare ones, and both accuracy conditions failed. It replaces the resampling, and the option stays experimental. Against the resampling, error and interval stability improved and average coverage barely moved. It would have failed the first test's coverage condition, which the plan replaced before any real result with these candidates. The verdict sits within simulation noise of NO-GO, because `deliveries` cleared its coverage limit by less than one Monte Carlo standard error. The real items had all been seen under other methods, and this round has had no outside review.
- **[Corrected variances through imputation](VARIANCE_MI_SYNTHESIS.md), 148,500 requests.** Under a rule fixed in advance, CV3 replaces the default variance in `catmech_pool_lm(cluster = )`, narrowly, and CR2 fails on width and over-coverage. New simulated cells show `catmech_reference(cluster = )` biasing slopes toward zero when the audited value is a predictor, so that use now draws a warning. No outside review.
- **[Overall check before release](OVERALL_QC.md).** Five independent AI reviewers, no human. The formulas and recomputed numbers held. Three cases gave badly wrong intervals with no signal, and they now warn or refuse in version 0.7.1. Several claims were corrected, and some items remain open.
- **[Comparison with outside methods](OUTSIDE_METHODS_SYNTHESIS.md), 297,000 requests, with a post hoc rerun of 58,356.** **KEEP** under a rule fixed in advance that could have retired `catmech_reference(cluster = )`. The margin is narrow: with a `mice` default turned off, `mice` with `2l.lmer` matched the option's accuracy, and the option's lead was slightly wider intervals with higher coverage on real items. Its recipes for raking and the difference estimator gave intervals that were too narrow, and those methods also failed on accuracy. No package code changed, and there has been no outside review.

## Software and independent checks

| Check | Result | Evidence |
|---|---|---|
| Unit suite | 701 assertions passed, no warnings or skips | `v0.7.1-checks/tests.log` |
| Source build and package check | 0 errors, 0 warnings, 2 notes: the new-submission notice, and one example just over the time limit on a machine busy with other checks | `v0.7.1-checks/build.log`, `v0.7.1-checks/R-CMD-check.log` |
| Fresh archive installation | Namespace-only APIs, all five bundled examples and the clustered option ran from an isolated library | `v0.7.1-checks/installed-smoke.log` |
| Contamination posterior oracle | Independent Metropolis sampler on the marginal posterior, sharing no package code: 38/38 moments within five MCSEs (largest 1.96), 36/36 row quantities within 0.015, three wrong likelihoods detected, 15/15 contract controls | `contamination-oracle/` |
| Older workflows unchanged | 18 of 18 runs identical to frozen v0.3.1 in completions, audits, hashes and pooled results, rerun with v0.7.1 | `v0.7.1-checks/v0.3-unchanged-check.txt` |
| Post-study additions leave frozen runs alone | 12 of 12 runs identical to frozen v0.4.0 sources, rerun with v0.7.1 | `v0.7.1-checks/v0.4.0-postfreeze-check.txt` |
| Clustered-audit study | Plain Python recomputes its 296 summary groups, held-out pool targets and verdict from raw records, and catches a corrupted summary; 11,998 development requests match the earlier study's estimates and clustered standard errors | `v0.7.1-checks/artifact-verification.log`, `clustered-pilot-check.txt` |
| Replacement study | Plain Python recomputes its 420 summary groups, 11 pool targets, seeds, selection and verdict from raw records, and catches a corrupted summary. As controls, 27,000 requests reproduce the first clustered study and 18,000 reproduce the review's reruns | `v0.7.1-checks/artifact-verification.log`, `clustered-fix-pilot-check.txt` |
| Package replay | All 30,000 study requests for the selected candidate replay through the package's own functions to within 1e-10 | `clustered-fix-replay/` |
| Independent numerical checks | 50 of 50 scored checks pass against `lme4`, `mitools` and hand-written matrix code on fixed datasets. A coding check, not a validation of the method | `cluster-algebra-check.txt`, `VARIANCE_CHECK_RESULTS.md` |
| Recount under corrected variances | 24,500 requests reproduce the study's oracle estimates. The current variance is not a candidate under the reading fixed in advance, and `cr2`, `cv3` and `wcr` are. Rare counts are not repaired | `VARIANCE_CHECK_PLAN.md`, `VARIANCE_CHECK_RESULTS.md`, `variance-check-decisions.txt` |
| Corrected-variance study | Plain Python re-derives the rule's outcome, regenerates the scored output and post-hoc table, and traces the write-up's numbers. All 148,500 requests replay through the package under `cv3`. 90,000 reproduce the replacement study under the earlier variance, and 24,500 the earlier recount scoring | `VARIANCE_MI_PLAN.md`, `variance-mi-decisions.txt`, `variance-mi-replay/`, `v0.7.0-checks/artifact-verification.log` |
| Comparison with outside methods | Plain Python checks hashes and request counts, regenerates the scored output, the rerun's score and the post-hoc table, re-derives the accuracy failures behind the outcome, and traces the write-up's numbers. As controls, 660 requests reproduce the option's saved results, and 390 reproduce `mice` under the rerun's setting | `OUTSIDE_METHODS_PLAN.md`, `OUTSIDE_METHODS_SENSITIVITY_PLAN.md`, `outside-methods-decisions.txt`, `outside-methods-sensitivity-decisions.txt`, `outside-methods-check.txt` |
| Random-intercept draw | Variance components match a one-way analysis of variance by hand, the moment trace matches a matrix calculation, generalized least squares matches explicit matrix algebra, and imputed rows carry the cluster correlation | `../tests/testthat/test-clustered.R` |
| Cluster-robust variance | Matches a hand-written CR1 calculation, the unequal-variance form with one row per cluster, and Rubin's rules with clusters minus one | `../tests/testthat/test-clustered.R` |
| Review of the clustered test | One independent AI reviewer re-derived the results and found no error in the verdict, code or comparisons, and corrected several readings, listed in the write-up. A rerun of the study's own option reproduces all 8,998 of its returned requests exactly. This is not human peer review. | `CLUSTERED_SYNTHESIS.md`, `clustered-review/`, `clustered-posthoc.txt` |
| Package change after the clustered study | Documentation comments in `R/reference.R`, stating the verdict and its limits; its code matches the frozen copy | `v0.5.0-checks/post-study-diff.txt` |
| Package change after the replacement study | The selected candidate moved into `R/reference.R`, `cluster_strata` was removed, and the hold-out refit stopped passing it. Every other source file is byte-identical to the frozen copy | `v0.7.0-checks/post-study-diff.txt` |
| New `use_recorded` option | Matches the study's linear MI on the same real-data cohorts in RMSE, coverage and width | `recorded-option-check/comparison.csv` |
| Serial and parallel pilots | Output files byte-identical between 1 and 14 workers | `contamination-pilot-check.txt` |
| Supplement for a runner bug | 1,731 score records and 120 hold-out records recovered; every estimate matched the main study to within 1e-9 | `contamination-supplement/` |
| Independent evidence recomputation | Plain Python recomputes all summaries, Wilson intervals and paired MCSEs for every study, both real-data target coefficients and pool audit counts, and the go/no-go verdicts, including the audit planner's summaries, pool `k` values and browser test vectors; deleting a planned cell is detected | `v0.7.0-checks/artifact-verification.log`, `tools/verify-artifacts.py` |
| Real-audit replay | 600 requests across all six methods, replayed from seeds with the runner's own functions, match the saved study; relabeled or rescaled results are caught | `real-audit-replay/`, `tools/verify-artifacts.py` |
| Review of the real-audit test | Three independent AI reviewers found no error in the rule or tables and corrected several readings, listed in the write-up. This is not human peer review. | `REAL_AUDIT_SYNTHESIS.md`, `real-audit-review/` |
| CRT port | Reproduces deposited naive results closely; adjusted RMSE 0.158 against 0.149 and coverage 0.946 against 0.963 | `merged-audit-pilot/reproduce/port-check.csv` |
| GitHub Actions | Passed on Ubuntu, macOS and Windows with R 4.6.2, `--as-cran` and warnings treated as errors: 679 assertions with no failures, warnings or skips, and status OK with no notes, at commit `6ff5da6`. The one commit after it records this result, repairs a write-up sentence that the local verifier had flagged, and refreshes that verifier's log. It changes no code, and it skipped the workflow | `v0.7.0-checks/github-actions.txt`, `../.github/workflows/R-CMD-check.yaml` |

The local environment is R 4.6.0, mice 3.19.0 and posterior 1.7.0 on Apple silicon macOS. No new package dependency was added. Exact session records accompany each study. Minimum-version R and a fresh-container restore have not been tested. The GitHub workflow checks the package on the current R release on three platforms, but every study ran only on this machine.

## Frozen studies and post-study changes

Plans preceded their runs: [synthetic](CONTAMINATION_PLAN.md), [real-data](REAL_DATA_PLAN.md) and [go/no-go](merged-audit-pilot/PLAN.md). Development runs and the changes they caused are disclosed in those plans. All three studies ran on the sources preserved in `v0.4.0-source/`, with hashes checked
before and after each run.

Two scoped changes came after those runs, as they did for v0.3.1: `catmech_reference(use_recorded = TRUE)` and the `advice` field of `catmech_holdout_check()`. The studies did not run either one directly. The first was checked against the study's linear MI on the same cohorts. The second is a reading of the real-data study and has not been tested prospectively.
The study runner also had a bookkeeping bug that dropped auxiliary records in six cells, recovered by the supplement above. Main estimates were unaffected.

Run times exceeded the estimates: the synthetic study took 40 minutes (estimated about 30), the real-data study 42 (estimated 15 to 20, later 25 to 30), the go/no-go 37 (estimated about 15) and the audit planner test 15 (estimated about 12).

Version 0.4.1 adds the audit planner and changes no existing function. Its [plan](AUDIT_PLAN_PLAN.md) preceded its evaluation run, which ran on the sources preserved in `v0.4.1-source/`. Development choices used the development split and are listed in the plan, along with the frozen tables already seen by hand.

The [real-audit test](REAL_AUDIT_PLAN.md) came after that release and ran on package sources identical to `v0.4.1-source/`. Its plan discloses the descriptive looks that preceded it and two changes made after its pilots, before any estimate was seen. It took 16 minutes against an estimate of 21. Its swap-in comparison and the tables in `real-audit-review/` are post hoc, are labeled as such and are not part of any rule. The plan was not under version control before that run, so its timing rests on file hashes and timestamps.

The [clustered-audit test](CLUSTERED_PLAN.md) added package code after v0.4.1. Its plan was committed before the option existed, and no amendment was needed. It ran on the sources preserved in `v0.5.0-source/`, with hashes checked before and after, and took about 13 minutes. The unplanned runs in `clustered-review/` and the table in `clustered-posthoc.txt` are post hoc and part of no rule. One independent AI reviewer read this round; no human has.

The [replacement test](CLUSTERED_FIX_PLAN.md) changes `catmech_reference(cluster = )` for version 0.6.0. Its plan was committed before either candidate was in the package or had run on a real item, and no amendment was needed. It discloses a scratch run on simulated audits in `clustered-fix-scratch/`, and the earlier results on the same real items. The study ran at commit `8c81044`, on the sources preserved in `clustered-fix-source/`, with hashes checked before and after, and took about 11 minutes. The table in `clustered-fix-posthoc.txt` is post hoc and part of no rule. The first clustered study's runner and the review's rerun need the option as it was at that commit. No reviewer, AI or human, has read this round.

The [corrected-variance test](VARIANCE_MI_PLAN.md) changes the default in `catmech_pool_lm(cluster = )` for version 0.7.0. Its plan was committed before the runner was frozen and before the study ran, though the runner had already been drafted. It amended a small runner fault and a status label after the run and before reading results, with the originals preserved in `variance-mi-source/`. The study ran at commit `33a6bee` and took about 11 minutes. The table in `variance-mi-posthoc.txt` is post hoc and part of no rule. No reviewer, AI or human, has read this round either.

The [comparison with outside methods](OUTSIDE_METHODS_PLAN.md) changes no package code. Its plan, runner and scorer were committed before it ran, and it took about 49 minutes. A [rerun of two comparator settings](OUTSIDE_METHODS_SENSITIVITY_PLAN.md) was planned after results were read, and it cannot change the outcome. It took three attempts, and its plan records why. 

Version 0.7.1 adds warnings and refusals from the [overall check](OVERALL_QC.md). The package sources that the comparison ran on are kept in `outside-methods-source/`. Replaying all 148,500 requests of the corrected-variance study through it reproduces the saved results (`v0.7.1-checks/variance-mi-replay.txt`). GitHub Actions also passed on Ubuntu, macOS and Windows at version 0.7.1 (`v0.7.1-checks/github-actions.txt`), at its commit before the history rewrite.

## Evidence limits

Real-data pools are public, and the ClinicalTrials.gov reference is sponsor-entered. The two public audits are facility counts and psychology effect sizes. No patient-level clinical source-data audit was available. Simulated error rates have not been checked against published audits, and the main fixture's 8% may be high. Hold-out checks had low power. Non-constant variance, survival and interaction models remain limits, and clustered records are handled only experimentally.
Adversarial or AI-assisted review is not independent human refereeing. AI use is disclosed in `CONTRIBUTING.md`.

## Reproduce

```sh
Rscript -e 'testthat::test_local()'
R CMD build .
R CMD check --no-manual catmech_0.7.0.tar.gz
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
Rscript tools/clustered-benchmark.R study validation/my-clustered-study
python3 -I tools/clustered-score.py my-clustered-decisions.txt validation/my-clustered-study
Rscript tools/clustered-review-checks.R validation/my-clustered-review
python3 -I tools/clustered-posthoc.py
Rscript tools/clustered-fix-benchmark.R study validation/my-clustered-fix-study
python3 -I tools/clustered-fix-score.py my-clustered-fix-decisions.txt validation/my-clustered-fix-study
python3 -I tools/clustered-fix-posthoc.py
Rscript tools/clustered-fix-replay.R validation/my-clustered-fix-replay
Rscript tools/variance-mi-benchmark.R study validation/my-variance-mi-study
python3 -I tools/variance-mi-score.py my-variance-mi-decisions.txt validation/my-variance-mi-study
Rscript tools/variance-mi-replay.R validation/my-variance-mi-replay
python3 tools/verify-artifacts.py
```

Use new output directories. The real-audit and clustered benchmarks always read the pools in `validation/real-audit/`, and the post-hoc, review and replay scripts read the saved study. To reproduce the exact fingerprints, run from the frozen sources in `v0.4.0-source/`, `v0.4.1-source/` for the audit planner test, `v0.5.0-source/` for the clustered-audit test, or commit `8c81044` for the replacement test and the review's reruns, since version 0.6.0 replaced the resampling they call. Its replay needs commit `8982952`. A fresh ClinicalTrials.gov download will differ from the dated pool. Set `CATMECH_BENCHMARK_WORKERS` to change parallelism. Exact cross-platform numerical equality is not promised.
