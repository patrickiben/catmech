# catmech

**An experimental R package for auditable imputation of continuous measurements
with quality flags, reference subsets and optional interval categories. Version 0.4.0.**

`catmech_model()` masks every flagged measurement and fits a conditional Gaussian
model using trusted observations. It can ignore categories, treat them as exact,
or use an externally specified category-error matrix. It preserves original data,
records assumptions and replacements, checks sampler diagnostics, and supports
linear-model pooling and comparison of declared sensitivity scenarios.
`catmech_reference()` uses genuinely measured reference values on an internal
subset, and `catmech_pmm()` provides an audited standard predictive-matching option.
`catmech_contamination()` models recording errors using those references, and
`catmech_recording_audit()` and `catmech_holdout_check()` report what the references
show on your own data.

The Gaussian model and data augmentation are established methods. The contribution
is a small, inspectable workflow with explicit failure handling. The earlier
retain/replace PMM engine remains available for **sensitivity analysis**; its
arbitrary retention of contaminated values caused severe bias and undercoverage.
See the [assessment](docs/ASSESSMENT.md), [model assumptions](docs/METHOD.md),
and [validation record](validation/VALIDATION.md).

## Install

```r
install.packages(c("digest", "jsonlite", "mice", "posterior", "truncnorm", "yaml"))
install.packages(".", repos = NULL, type = "source")
```

Run the second command from the repository directory, or install the source archive
with `R CMD INSTALL catmech_0.4.0.tar.gz`. R 4.1 or later is required. GitHub Actions
is configured for Linux, macOS and Windows; local checks do not establish that
the remote workflow has passed.

## Start with a declared trust policy

```r
library(catmech)

dat <- simulate_catmech_lab(n = 300, seed = 1)
# Hidden truth is for evaluation only, never an imputation predictor.
dat <- dat[setdiff(names(dat), c("lab_true", "corrupted"))]

fit <- catmech_model(
  dat, value = "lab_value",
  predictors = c("age", "sex", "inflammation", "outcome"),
  flag = "outlier_flag", category_mode = "ignore",
  m = 20, seed = 123
)

catmech_pool_lm(fit, outcome ~ lab_value + age + sex + inflammation)
head(fit$audit)
fit$model_diagnostics$summary
completed_one <- fit$completed[[1]]

# Use a new or empty directory; existing files are not overwritten.
write_catmech(fit, "output/example")
```

Including the downstream outcome in the imputation model matters. All selected
predictors must be complete. Every flagged target is replaced, including false
positive flags. Unflagged observed values are trusted and preserved exactly;
missing flags count as unflagged. This policy requires an ignorable trusted-data
selection mechanism given the predictors. It does not detect missed errors or
estimate whether any particular measurement was wrong.

## When flags miss errors: use actual reference information

Masking flagged values cannot identify errors among measurements still treated
as trusted. If an independently measured reference subset is available, use it
as the observed target and impute every unvalidated target:

```r
# reference_value contains actual reference measurements on selected rows,
# with NA elsewhere. Do not manufacture reference values from flags or W.
fit_reference <- catmech_reference(
  dat_with_reference, value = "lab_value", reference = "reference_value",
  predictors = c("age", "sex", "inflammation", "outcome"),
  method = "normal", m = 20, seed = 123
)
catmech_pool_lm(fit_reference,
  outcome ~ lab_value + age + sex + inflammation)
```

A random internal reference sample is a transparent design. Targeted verification
requires a suitable selection model and overlap; noisy references cannot be treated
as exact truth. The default omits recorded W as a predictor, sacrificing information
rather than silently retaining unvalidated values. This requires additional reference measurements. Read the
[reference method](docs/REFERENCE_METHOD.md) and
[synthetic example](inst/examples/reference.R).

For a donor-based alternative, `catmech_pmm(dat, "lab_value", predictors,
flag = "outlier_flag")` masks all flags and delegates standard type-1 PMM to `mice`.
`catmech_reference(..., method = "pmm")` restricts donors to reference rows, but
showed substantial undercoverage in the synthetic held-out study. Its pooled analysis
requires `allow_sensitivity = TRUE` and remains labeled `sensitivity_only`. `method = "pmm_abb"`
resamples donors for each imputation. That raised synthetic coverage from 86.7% to 95.9%, but it
overcovered on real data with intervals up to 52% wider, and standard PMM did not repeat its failure
there. The reference default remains normal. Predictors must be complete. An observed category can be selected explicitly as
an ordinary factor predictor; no interval constraint is then imposed. PMM cannot
extrapolate beyond donor values and does not guarantee robustness or coverage.
The new evaluation compares methods given the **same predictor information**.

## With references and recorded values: which workflow

Read what the references show before choosing a model. `catmech_recording_audit()` fits no
model. It reports disagreement rates and error sizes, tests whether rates differ across groups
such as site or flag, and compares referenced with unreferenced rows.

```r
catmech_recording_audit(dat_with_reference, "lab_value", "reference_value",
                        by = c("outlier_flag", "site"))
```

| Situation | Start with | Evidence and limit |
| --- | --- | --- |
| Most cases, especially errors small relative to the spread of the target | `catmech_reference(..., use_recorded = TRUE)`: linear imputation given the recorded value. Consider generalized raking when selection probabilities are known; the package does not implement raking. | Best or tied among reference-using methods on both real datasets at 25% audits. It failed under synthetic gross errors. |
| Large or gross errors, or a fallible reference with a credible declared SD | `catmech_contamination()` | Strongest in synthetic cells with these errors. Worse than linear MI and raking on ClinicalTrials.gov. It failed with heavy-tailed targets and with selection on the true value when errors were present. |
| Target residuals doubtful and precision not critical | `catmech_reference()` (default) | Assumes least about recorded values, with intervals about twice as wide on the real data. It failed with selective or noisy references. |

```r
fit <- catmech_contamination(
  dat_with_reference, "lab_value", "reference_value",
  predictors = c("age", "sex", "inflammation", "outcome"),
  corruption_strata = "outlier_flag", m = 20, seed = 125
)
catmech_pool_lm(fit, outcome ~ lab_value + age + sex + inflammation)
fit$contamination$corruption    # per-record posterior error probabilities
catmech_holdout_check(fit)      # scores unseen references, using no simulation truth
```

Here "contamination" means recording error, the usage in the selective-editing literature. It is
not treatment spillover between trial arms.

Selection into the audit may depend on recorded values, so a targeted audit is supported. Selection
on the true value beyond what the recorded value shows is not. `selection_delta` and `reference_sd`
are declared sensitivity assumptions to compare with `catmech_sensitivity()`. If the hold-out check shows
no skill against the recorded value or the stratum error rate, prefer linear MI or raking. On both real
datasets the model could not tell which records were wrong (AUC 0.54 and 0.50), so do not read its per-record
probabilities as labels on data like those. See the [method](docs/CONTAMINATION_METHOD.md),
[synthetic evidence](validation/CONTAMINATION_SYNTHESIS.md), [real-data evidence](validation/REAL_DATA_SYNTHESIS.md)
and [example](inst/examples/contamination.R).

## Use categories only under an explicit observation model

```r
labels <- c("Low", "Medium", "High")

# Example sensitivity assumption, not an estimated or clinically validated matrix.
H <- matrix(.05, 3, 3, dimnames = list(labels, labels))
diag(H) <- .90                  # Rows: true category. Columns: observed category.

noisy <- catmech_model(
  dat, "lab_value", c("age", "sex", "inflammation", "outcome"),
  flag = "outlier_flag", category = "lab_category",
  thresholds = c(-Inf, 5.5, 7, Inf), levels = labels,
  category_mode = "noisy", confusion = H, m = 20, seed = 124
)

catmech_sensitivity(
  list(ignore = fit, noisy_assumption = noisy),
  outcome ~ lab_value + age + sex + inflammation,
  term = "lab_value", expected_scenarios = c("ignore", "noisy_assumption")
)
```

`category_mode = "exact"` uses exact bins and rejects contradictions with trusted
numeric values. Noisy categories use the supplied likelihood instead of forcing
values into potentially incorrect bins. Category labels stay unchanged. The matrix
is fixed: uncertainty in its entries is not propagated. A category derived from
the same corrupted measurement is not independent corroboration.

Category models use four Gibbs chains by default. Failed R-hat or effective-sample
size checks raise a `catmech_convergence_error` retaining diagnostic traces. The
[complete example](inst/examples/model.R) shows how to preserve failed scenarios
in the sensitivity report. Its envelope is descriptive, not a confidence interval
with a new coverage guarantee.

## Scope and evidence

- One continuous target, complete selected predictors, independent records and
  a linear, homoskedastic Gaussian conditional model.
- Informative categories cover the real line, with right-closed intervals.
  Clinical boundary and rounding conventions must be reconciled first.
- Bayesian parameter uncertainty, posterior predictive completion, per-cell audit
  records, assumption/data hashes, and an externally generated `mids` adapter.
  Do not continue the adapter with `mice.mids()`.
- A prespecified held-out study compares the new workflow with ordinary trusted
  normal MI, trusted PMM, naive and oracle analyses, and the legacy method. It
  reports bias, coverage, MCSE, failures, and deliberate wrong-model stresses.
- Version 0.4 adds a synthetic study with linear MI and generalized raking as strong
  comparators, and a [real-data audit](validation/REAL_DATA_SYNTHESIS.md) on ClinicalTrials.gov and
  NHANES with real recording differences. Neither dataset is a clinical source-data audit.

The [v0.2 coverage protocol](validation/COVERAGE_PLAN.md) transfers useful controls from
the author's carryover, False Floor, calibration and assumption-audit research.
Independent numerical oracles check the category sampler and the reference
model's joint predictive distribution, including shared parameter uncertainty. Clinical-format
synthetic fixtures include missed flags, site shifts and rounded records.
[Transfer details](docs/RESEARCH_TRANSFER_PLAN.md) state what does and does not
carry over. The [new extension study](validation/EXTENSION_SYNTHESIS.md)
reports matched-information PMM comparisons, the value of actual reference
measurements, and failures under selective or inaccurate references. These are simulation studies, not clinical validation or a general
coverage theorem. Nonlinearity, missed contamination, informative selection,
incompatible analysis models and dependent records remain limitations.

## Legacy API

`catmech_spec()` and `catmech_impute()` retain the v0.1 PMM sensitivity semantics.
Their `corruption_probability` is an imposed replacement probability, not a
posterior corruption estimate. `catmech_pool_lm()` refuses this engine unless
`allow_sensitivity = TRUE`; that override labels the result as sensitivity-only.
See the [legacy method](docs/LEGACY_METHOD.md), [example](inst/examples/lab.R),
and [historical benchmark](validation/benchmark/RESULTS.md). Existing specifications
are not silently converted to the new model.

## Reproduce

```sh
Rscript -e 'install.packages(c("testthat", "pkgload", "roxygen2"))'
Rscript -e 'testthat::test_local()'
R CMD build .
R CMD check --no-manual catmech_0.4.0.tar.gz
Rscript tools/check-model-oracle.R validation/my-model-oracle
Rscript tools/coverage-benchmark.R confirmatory validation/my-heldout
Rscript tools/coverage-benchmark.R clinical validation/my-clinical
Rscript tools/check-reference-oracle.R validation/my-reference-oracle
Rscript tools/check-pooling-guards.R validation/my-pooling-guards
Rscript tools/extension-benchmark.R distribution validation/my-distribution
Rscript tools/extension-benchmark.R reference validation/my-reference-study
Rscript tools/extension-benchmark.R stress validation/my-reference-stress
Rscript tools/extension-benchmark.R controls validation/my-reference-controls
Rscript tools/check-contamination-oracle.R validation/my-contamination-oracle
Rscript tools/contamination-benchmark.R study validation/my-contamination-study
python3 -I tools/real-data-fetch.py my-raw-directory validation/my-real-data
Rscript tools/real-data-prepare.R my-raw-directory validation/my-real-data
Rscript tools/real-data-benchmark.R study validation/my-real-data-study
python3 tools/verify-artifacts.py
```

Use new output directories. Archived v0.1 to v0.4 sources retain the historical
benchmark inputs; their hashes are checked against those frozen versions. A fresh
ClinicalTrials.gov download will differ from the dated pool because the registry changes.
The v0.3.1 guards reject unsupported analyses and gate experimental reference-PMM
pooling; successful scalar estimates and intervals are unchanged. See
[validation](validation/VALIDATION.md) for evidence, environment and remaining gaps.

MIT license. Built on [mice](https://doi.org/10.18637/jss.v045.i03),
[data augmentation](https://doi.org/10.1080/01621459.1987.10478458), interval imputation
such as [hmi](https://doi.org/10.18637/jss.v095.i09), and
[multiple overimputation](https://doi.org/10.1177/0049124115585360).
AI-assisted work is disclosed in [CONTRIBUTING.md](CONTRIBUTING.md).
[GitHub preparation](docs/GITHUB.md) describes upload steps; no publication is
performed by preparing this repository.
