# Review of 1.0.0 and fixes, 2026-11-07

This review covers the 28 commits that took catmech from 0.10.1 to 1.0.0 (662252a to
600462f). An AI assistant carried out the review and wrote the fixes at the maintainer's
request, and three further AI reviewers then checked the fixes (see the last section).
No independent statistician has reviewed them. The numbers in the tables, and the
equivalence with `mice::pool()`, come from
[`tools/review-fixes-2026-11-07.R`](../tools/review-fixes-2026-11-07.R), whose output is
[`review-fixes-2026-11-07/results.json`](review-fixes-2026-11-07/results.json). The sampler
equivalence comes from
[`tools/review-sampler-equivalence-2026-11-07.R`](../tools/review-sampler-equivalence-2026-11-07.R).
The check results, test counts and the one timing come from the runs described where
they appear.

## What held up

- `R CMD check --as-cran` on 600462f: no errors or warnings, and two notes (new submission;
  a local HTML Tidy version). All 1,754 test expectations passed, and GitHub CI passed on
  Ubuntu, macOS and Windows.
- The 1.0.0 rewrite of the marginal Gibbs sampler (`src/marginal.cpp`, and continuing chains
  instead of rerunning them in `R/marginal-sampler.R`) reproduces 0.10.1 exactly. With the
  same seed, the two builds gave identical completions, posterior summaries and per-attempt
  diagnostics: from 1,000 iterations (converged after doubling to 8,000) and from 100
  (refused at the 800-iteration cap). In one such fit on the reviewer's Mac the rewrite ran
  in 3.5 s against 6.0 s.
- Unclustered `catmech_pool_glm()` and `catmech_pool_coxph()` match `mice::pool()`, which
  also uses events minus coefficients as the complete-data degrees of freedom of a Cox fit.

## Problems found and fixed

### 1. `catmech_compare_strategies()` left the outcome out of the imputation model

When it chose the imputation predictors itself, it removed the analysis outcome. With the
audited value as a predictor, that biases the value's coefficient toward zero. In 100
simulated data sets (n = 600, 40% audit, recorded = true + N(0, 0.7²), true slope 0.80):

| Strategy | Mean slope (MC SE) |
| --- | --- |
| Validation subsample | 0.798 (0.008) |
| Naive, trusting recorded values | 0.520 (0.004) |
| Linear imputation, 1.0.0 default (outcome left out) | 0.633 (0.005) |
| Linear imputation, new default (outcome kept) | 0.800 (0.006) |
| Contamination model, 1.0.0 default | 0.617 (0.004) |
| Contamination model, new default | 0.694 (0.005) |

Every record carries error here, so the contamination model, which assumes some records
are exactly right, is misspecified; its remaining bias is expected.

While fixing this, a second problem turned up. The function passed `cluster` into
`catmech_reference()`, which switched on the experimental clustered imputation that biased
slopes when the audited value was a predictor (`VARIANCE_MI_SYNTHESIS.md`), and it gave the
validation and naive strategies model-based variances while the imputed strategies got
cluster-robust ones. `cluster` now only sets the variance, the same way in all four
strategies, over the coefficients each strategy can estimate. A contamination fit that fails
now warns with its reason instead of silently returning `NA` rows; the old tests had been
passing on such silent failures.

### 2. The GLM and Cox examples left the outcome out, and pooling did not warn

The examples for `catmech_pool_glm()`, `catmech_pool_coxph()`, `catmech_pool_lm()`, the plot,
`tidy()`, `glance()`, sensitivity and audit-trail help pages, and the clinical audit
vignette, imputed the lab value without the analysis outcome. In 200 simulated data sets
(n = 600, 40% audit, true log odds ratio 0.80, m = 10):

| Imputation model | Mean estimate (MC SE) | Coverage of nominal 95% (MC SE) |
| --- | --- | --- |
| Outcome left out | 0.611 (0.007) | 65.0% (3.4%) |
| Outcome included | 0.819 (0.011) | 95.0% (1.5%) |
| Reference: true values | 0.808 (0.008) | |
| Reference: recorded values | 0.497 (0.006) | |

The three pooling functions now warn when the audited value is an analysis predictor and
the outcome was not an imputation predictor (for a Cox model, the event indicator, read from
the `Surv()` call by argument). The check covers every engine, including `catmech_model()`,
and is skipped when nothing was imputed. The examples, the vignette and the tests now
include the outcome, and for Cox models the event indicator and the Nelson-Aalen cumulative
hazard (White and Royston, 2009, <doi:10.1002/sim.3618>). `catmech_pool_glm()` and
`catmech_pool_coxph()` also give the clustered-imputation warnings that `catmech_pool_lm()`
already gave.

### 3. The cluster-robust GLM variance was scaled by the dispersion

`catmech_glm_cluster_vcov()` used `vcov(fit)` as the bread of the sandwich, which carries the
estimated dispersion into the variance twice. Standard errors were therefore multiplied by
the dispersion for gaussian, Gamma and quasi families: too wide when it was above 1, and too
narrow when it was below 1, as in most Gamma fits. Cluster-robust standard error of the slope,
600 observations in 30 clusters:

| Outcome and family | Dispersion | 1.0.0 | Fixed |
| --- | --- | --- | --- |
| Overdispersed counts, Poisson | 1 | 0.0325 | 0.0325 |
| Overdispersed counts, quasi-Poisson | 3.06 | 0.0993 | 0.0325 |
| Gamma, log link | 0.41 | 0.0083 | 0.0205 |

The bread is now the unscaled inverse information, so the sandwich no longer depends on the
dispersion. Binomial and Poisson results do not change.

### 4. `catmech_audit_power()` had no validation and no label saying so

Without `plan` its factor is k = q / (1 - q), q = error rate (1 - correlation²)(1 - R²), with
`correlation` the recorded-true correlation over all records. Multiplying by the error rate
treats the unexplained variation as if it sat in erroneous records the analysis can identify.
At its defaults (N = 1,000, error rate 0.15, correlation 0.85, R² 0.30) it recommends an 8%
audit, with an expected width ratio of 1.18 to a full audit. For linear imputation given the
recorded value, the package's recommended start, the double-sampling regression estimator
at the same correlation has k = (1 - 0.85²)(1 - 0.30) = 0.194. At 8% that means a ratio of
about 1.80, and reaching 1.25 needs about 26%. The function is now labeled a rough,
unvalidated planning aid, its help page states the assumption, and without `plan` it also
reports this linear-imputation figure. That figure ignores the uncertainty in k, so it is not
an upper bound: at very high error rates it falls below the main recommendation. The function
still has not been validated.

### 5. Wording that claimed more than the evidence

- `catmech_audit_trail()` called its record tamper-evident and said the pooled intervals were
  tied to the hashes. The data and completion hashes are of serialized R data frames, the
  specification hash is of canonical JSON text, the pooled results are not hashed, and anyone
  who edits the data can regenerate the record. The wording now says so, an unknown package
  version is recorded as "unknown" rather than "1.0.0", and the HTML title is escaped.
- The vignette called the package "mathematically sound", said the audit-size curve
  "guarantees" a width, called the provenance record "Tamper-Evident" and "cryptographic"
  regulatory documentation, recommended `cluster_model = "joint"` for clustered sites against
  the README's own evidence, and gave an unsourced 10% to 15% audit rule. The README's method
  guide also recommended the joint cluster model. Both now follow the evidence, and the help
  pages, vignette and pkgdown index no longer call the audit-size output statistical power.
- `validation/VALIDATION.md`, `validation/v1.0.0-VALIDATION.md` and `docs/ASSESSMENT.md` called
  1.0.0 production-ready and a clinical toolkit; each now opens with a dated correction, with
  the original text kept. The 1.0.0 NEWS entry said "First CRAN release"; it now says the
  release was not submitted to CRAN and points to the corrections. `HANDOFF.md` said the
  package was "100% audited, verified, and ready for official CRAN submission" and that the
  documentation site was deployed, which it is not (GitHub Pages is not enabled); both are
  corrected, and its local paths are removed.
- `education/verify_educational_companion.py` named a private folder as a fallback path and
  silently skipped four gates when that folder was missing. It now reads the harness location
  only from `RESEARCH_VALIDATION_HARNESS_DIR`, reports skipped gates, and fails on them with
  `--require-harness`.

## Independent review of the fixes

Three AI reviewers, each with one lens (statistics and code, documentation and claims, tests
and reproducibility), read the changes and ran their own checks against builds of 0.10.1,
1.0.0 and the fixes. They raised 26 points. All were checked, and all except the version
number were fixed before commit; the version was then bumped to 1.0.1 (below). The ones that changed
results or claims:

- With `cluster =`, a validation-subsample fit with a coefficient the audited rows cannot
  estimate (a covariate constant among them) gave meaningless standard errors for linear
  models and an error for logistic ones. This came from the new cluster path in these fixes,
  not from 1.0.0. The cluster variances now drop aliased coefficients, and tests compare
  them with direct calculations on the estimable columns.
- The first draft of this record and NEWS said the GLM variance bug "inflated" standard
  errors. It scaled them by the dispersion, so they were too narrow below 1; section 3 now
  says so and shows a Gamma example.
- An unescaped percent sign cut a sentence short in the `catmech_pool_glm()` help page; the
  Cox help page wrongly said `mice::pool()` uses observations minus coefficients; the
  audit-trail help said every hash was of an R object; the audit-size help read `correlation`
  two ways; the README's new instruction sat inside a code block; and two README-linked files
  still said production-ready. All are corrected.
- Tests were added for each fix, including a contamination strategy that succeeds (the
  existing test's fits now fail to converge with the outcome included, which the new warning
  reports), the new GLM and Cox clustered-imputation warnings, the binomial cluster path, and
  the outcome check's handling of unimputed fits, `Surv()` arguments and `catmech_model()`.

## Not changed here

- `paper/` and the educational companion's text are untouched; the companion's files are
  locked by hashes and a human sign-off. The manuscript's workflow example pools a binary and
  a survival endpoint from an imputation that included only the continuous outcome, so its
  GLM and Cox fits now draw the new warning, as do the companion's executable labs C17 and
  C18 (they still pass). The manuscript also describes `catmech_audit_power()` as computing
  statistical power and the hashes as cryptographic checksums, and
  `paper/VALIDATION_LEDGER_2026-11-06.md` still contains local paths.
- The fixes are released as version 1.0.1, so provenance records can tell a fixed build from
  1.0.0; the README front page, the case study and the citation files were updated to match.
  `results.json` records 1.0.0 because it was generated before the version change, from the
  commit it names; the code it ran is the same.
