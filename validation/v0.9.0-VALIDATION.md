# Validation record: catmech 0.9.0

Date: 2026-11-04. This is an experimental research-software release. It adds explicit posterior coefficient inference and a larger, faster sampling budget for the joint clustered model. The [v0.8.0 record](v0.8.0-VALIDATION.md) preserves the earlier release and its studies. The model and sampler are established methods; no novel-method or clinical-validation claim is made.

## What changed

`catmech_joint_summary(fit)` summarizes the population target coefficients in the exact imputation design, using every retained posterior draw. It reports original-unit means, posterior standard deviations, equal-tailed credible intervals, R-hat, bulk/tail ESS and Monte Carlo errors for means and endpoints. The intercept transformation includes predictor centering. This is an explicit model-based analysis under independence of site effects and selected covariates. It accepts no downstream formula and refuses an all-reference result with no fitted posterior. `catmech_pool_lm()` retains its CV3/Rubin calculations and separate interpretation.

The joint sampler now starts with 4,000 retained iterations per chain and extends to 32,000. Four chains, 1,000 burn-in iterations and all prior distributions are unchanged. R-hat must remain at most 1.05; bulk and tail ESS must each be at least 100 times the chain count. The original-unit intercept is also monitored. Trace construction uses recipient-pattern design averages and site weights instead of repeated row-level calculations. The predictive law, fresh-process isolation and exact-restart streams remain intact. See the [method and limitations](../docs/JOINT_CLUSTER_METHOD.md).

## Prospective comparison

The staged [protocol](JOINT_CALIBRATION_PLAN.md) was committed at `1355c74` before any candidate accuracy run. The implementation was committed at `719bd79` and the runner at `b0d4211`. A masked engineering pilot returned all 40 analyses. Development used 210 cohorts and 1,200 analysis requests; 1,198 returned, with both refusals in the old sampler profile. Its unfavorable null-slope coverage and contextual-covariate bias remain documented. The [development decision](JOINT_CALIBRATION_DECISION.md), committed at `a9034d8`, froze the unchanged candidate before confirmation. No interval multiplier, degrees-of-freedom adjustment or alternative prior was selected to improve observed coverage.

Confirmation uses 400 fresh cohorts in each of six original scenarios and 200 in the contextual-covariate stress: 2,600 cohorts, 6,400 MCMC fits and 12,800 analysis requests, including 2,600 complete-truth oracle analyses. Prior alternatives use the same first 200 cohorts in each of three fixed cells. Candidate CV3 and direct posterior summaries share each posterior fit, isolating the analysis choice. Calibration, relative interval-width sensitivity and fit availability have separate prespecified screens. Direct posterior standard deviations are not labeled robust standard errors. Results and Monte Carlo uncertainty belong in the [scientific synthesis](JOINT_CALIBRATION_SYNTHESIS.md).

## Confirmatory findings

All three prespecified screens passed. Of 12,800 analysis requests, 12,755 returned intervals; all 45 refusals arose in the old sampler profile. Every one of the 3,800 candidate fits and its two analyses returned. In the six original scenarios, fit refusals fell from 33/2,400 to 0/2,400; the paired refusal difference was -1.375 percentage points (95% Monte Carlo interval -1.840 to -0.910). The contextual stress contributed the remaining 12 old-profile refusals.

Across the four primary Gaussian scenarios, direct coverage was 96.75%, 96.25%, 96.50% and 93.00% (balanced, unequal-size, observed-MAR and null). Mean absolute distance from 95% was 1.625 percentage points, compared with 3.168 for old CV3 and 3.125 for candidate CV3. Both paired bootstrap improvement intervals excluded zero. These finite-grid screens do not establish exact nominal coverage; the null result remains a caution. Near-zero-ICC and t5 coverage was 94% and 95%, respectively.

On 600 common prior triples, the mean absolute log-width sensitivity statistic was 0.09757 for direct inference versus 0.12026 for candidate CV3, a descriptive 18.9% relative reduction. The paired difference was -0.02269 (95% Monte Carlo interval -0.02553 to -0.01985). The likelihood and covariance priors are unchanged. Point-estimate sensitivity was not uniformly reduced, so the claim concerns reported interval widths, not weaker prior influence on the fitted posterior.

The contextual-covariate stress exposed a substantial limit: direct coverage of the marginal target was 63.5%, with mean bias -0.1279; candidate CV3 covered 96.5% there. The fitted independent-site-effect assumption fails in that scenario. The new summary is an explicit model-based coefficient analysis, not a robust marginal-effect replacement. Complete tables, paired comparisons, interval scores and Monte Carlo uncertainty are in the [synthesis](JOINT_CALIBRATION_SYNTHESIS.md).

The independent [confirmation verifier](joint-calibration-confirmation/verification.json) checked all requests, 80,600 diagnostic rows, 27 unchanged source files, 1,520 summary values, 198 paired contrasts, every 2,000-replicate paired-bootstrap result and all adoption screens. Nine deliberate corruptions were rejected. Passing this audit establishes reporting correspondence; mathematical validity is addressed separately by the model, assumptions, oracle tests and simulation evidence.

## Mathematical and numerical checks

The new coefficient tests construct an independent Gaussian observed-data covariance matrix with known site and row covariances. They derive the normal coefficient posterior by generalized least squares and check means, variances and quantiles, including factor contrasts and nonzero predictor centers. A near-cancellation case detects failure of the original-scale intercept despite passing standardized-coefficient diagnostics. This checks extraction and transformations, not the unknown-covariance sampler's entire mathematics.

An end-to-end test also confirms identical full posterior traces and coefficient summaries for `m = 20` and `m = 40` at fixed data and seed. Direct coefficient inference uses all retained parameter draws; the CV3 validation remains specific to 40 completions.

The trace tests compare vectorized pattern means with individual-row conditioning across batch boundaries, repeated/reversed draw indices and one-pattern/intercept-only cases. Under explicit old controls, the engineering fixture produced bit-identical completions and raw parameter traces; previous conditional-mean trace columns differed by at most 3.8e-16. Trace processing was about 18 times faster on a single fixed-array benchmark, an illustrative timing rather than a controlled speed study.

All 22 selected historical convergence refusals passed at 16,000 retained draws under the new schedule. This challenge retained original inputs, seeds and prior scales, and computed no coefficient accuracy. It is a selected numerical check, not an estimate of a population failure rate. The six original data generators also reproduced 18 frozen cases byte-for-byte. [Engineering evidence](v0.9.0-development/engineering-check.txt), [challenge](v0.9.0-development/challenge/VERIFICATION.txt), [generator check](v0.9.0-development/dgm-compatibility.log).

## Release verification

All 978 package assertions passed without warnings or skips. `R CMD check --as-cran --no-manual` on the built `catmech_0.9.0.tar.gz` finished with zero errors, zero warnings and one routine incoming-check note, “New submission.” The archive contains exactly the tested namespace, R code, examples, new tests, README, NEWS and help. R's generated DESCRIPTION metadata was checked separately from its original fields.

| Check | Evidence |
|---|---|
| CRAN-style source check and all package assertions | `v0.9.0-checks/R-CMD-check.log`, `v0.9.0-checks/00check.log`, `v0.9.0-checks/tests.log` |
| Fresh source installation, all six examples, namespace/adapter/export checks | `v0.9.0-checks/install.log`, `v0.9.0-checks/installed-smoke.log` |
| Identical posterior traces and summaries at 20 versus 40 imputations | `v0.9.0-checks/joint-count-invariance.log` |
| Exact legacy compatibility: 18 v0.7.2, 18 v0.3.1 and 12 v0.4.0 workflows | `v0.9.0-checks/v0.7.2-unchanged-check.txt`, `v0.9.0-checks/v0.3-unchanged-check.txt`, `v0.9.0-checks/v0.4.0-postfreeze-check.txt` |
| Historical evidence and current learning-page checks | `v0.9.0-checks/artifact-verification.log`, `v0.9.0-checks/learning-page-check.log` |
| README facts and dated historical case study | `v0.9.0-checks/fact-check-final.log`, `../docs/case-study/FACT_CHECK.md` |
| Source archive contents, tested-code hashes and checksum | `v0.9.0-checks/build-integrity.txt`, `v0.9.0-checks/tested-source-hashes.json`, `v0.9.0-checks/source-package-sha256.txt` |
| Self-checked scientific figure and visual QA | `joint-calibration-figures/figure-provenance.json`, `joint-calibration-figures/visual-QA.txt` |

The dated v0.7.1 case-study HTML/PDF remains historical. Current GitHub Actions runs the package checks on Ubuntu, macOS and Windows; see the [workflow results](https://github.com/patrickiben/catmech/actions).

## Reproduction

Use a new output directory. From this release, run:

```sh
Rscript validation/run-joint-calibration.R confirmation validation/my-joint-calibration
python3 -I tools/verify-joint-calibration.py validation/my-joint-calibration --tamper-tests
Rscript tools/plot-joint-calibration.R validation/my-joint-calibration validation/my-joint-figures
```

The study runner refuses existing evidence directories, records every failure and warning, and checks source hashes before and after execution. Its default is eight outer workers; the saved confirmation uses sixteen, an execution choice that does not alter seed mapping. Each sampler uses at most two simultaneous fresh R workers. The independent verifier reconstructs the request grid, raw interval arithmetic, all scalar summaries, paired differences, prior triples and every paired bootstrap replicate. Only bootstrap index generation delegates to base R for its declared RNG stream; all reporting arithmetic is recomputed in Python. Deliberately damaged temporary copies must be rejected.

Historical v0.8.0 evidence now verifies against byte-identical sources in `v0.8.0-source/`, copied from public commit `b236293`. Its original runner is also frozen there. Earlier Rwanda, rare-count, predictor-role and other studies retain their original sources and findings; they do not validate the new posterior summary. AI agents performed implementation and review. There has been no independent human statistical review.
