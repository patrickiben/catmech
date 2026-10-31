# catmech

**Experimental R package for auditable sensitivity completion of continuous
measurements with trusted interval categories and quality flags.**

catmech makes assumptions explicit, checks conflicts, restricts donors to compatible
intervals, and records how each requested value was handled. It combines established
imputation components with a small, inspectable workflow and `mice` analysis support.
It does not identify latent truth or estimate corruption probabilities.

**Status: research prototype, version 0.1.0.** The exploratory benchmark found good
performance in some exact-category scenarios and serious bias and undercoverage
when uncertain-flag retention or incorrect categories were used. It found no clear
general advantage over simply masking flagged measurements and using `mice` PMM.
Read the [assessment](docs/ASSESSMENT.md),
[method and limitations](docs/METHOD.md), and
[benchmark results](validation/benchmark/RESULTS.md) before interpreting pooled results.

## Install from this repository

```r
install.packages(c("digest", "jsonlite", "mice", "truncnorm", "yaml"))
install.packages(".", repos = NULL, type = "source")
```

Run the second command with the repository as the working directory. Alternatively,
install the release source archive with `R CMD INSTALL catmech_0.1.0.tar.gz`.
R 4.1 or later is required. The checked local environment is recorded in
[validation](validation/VALIDATION.md); the CI workflow checks three operating systems
when the repository is hosted on GitHub.

## A reproducible synthetic example

```r
library(catmech)

dat <- simulate_catmech_lab(n = 300, seed = 1)
# Truth is for validation only. Never use it as an imputation predictor.
dat <- dat[setdiff(names(dat), c("lab_true", "corrupted"))]

spec <- catmech_spec(
  value = "lab_value", category = "lab_category",
  levels = c("Low", "Medium", "High"),
  thresholds = c(-Inf, 5.5, 7, Inf),
  flag = "outlier_flag",
  corruption_probability = 1,
  conflict_policy = "reconstruct"
)

fit <- catmech_impute(
  dat, spec, predictors = c("age", "sex", "inflammation", "outcome"),
  m = 20, seed = 123
)
fit$diagnostics
head(fit$audit)
completed_one <- fit$completed[[1]]

# An analysis adapter, not a guarantee that the intervals are calibrated:
models <- with(fit$mids, lm(outcome ~ lab_value + age + sex + inflammation))
summary(mice::pool(models))

# Use a new or empty directory; existing files will not be overwritten.
write_catmech(fit, "output/example")
```

The simulator's default flags all indicate actual corruption, so this example
explicitly reconstructs every flag (`q = 1`). A flag on real data does not establish
that assumption. When a flag column is configured, `corruption_probability` must
be supplied; there is no default inferred probability. Explore several plausible
values and inspect whether conclusions change. A value equal to the marginal
fraction of corrupt flags need not be appropriate after conditioning on the
observed value and category.

The same specification can be read from YAML:

```r
spec <- read_catmech_spec(system.file("specs", "lab.yaml", package = "catmech"))
```

See [the example script](inst/examples/lab.R) for a complete run and a q sensitivity
grid. Its output uses synthetic data only.

## Supported scope

- One continuous target and explicitly selected, complete predictors.
- Ordered exact categories: `[a,b]` for the first interval and `(a,b]` thereafter.
- Missing target values, flagged observations, and explicit category-authority
  decisions for conflicting observations.
- Original trusted observations as donors, with type-1 predictive mean matching.
- Independent draws, per-cell audit records, assumption/data hashes and a `mids`
  adapter for analysis. The adapter is not a fitted chained-equations process and
  must not be resumed with `mice.mids()`.

Category errors, unflagged contamination, informative selection of trusted donors,
and poor donor support can invalidate the analysis. Missing flags are treated as
non-suspect, which is an assumption. A missing category is derived only from an
unflagged observed value. A missing or suspect value with no category is rejected.
Other incomplete predictors are rejected. No soft-category model, joint missingness
model, MNAR adjustment, topology or causal intervention is implemented.

## Policy and audit

`conflict_policy = "error"` is the default. Choosing `"reconstruct"` declares the
category authoritative, including conflicts on unflagged rows. It does not prove
which measurement was wrong.

`fallback = "error"` refuses insufficient within-category donor support. An opt-in
`"normal"` fallback uses a truncated regression draw and is separately audited.
Its stronger extrapolation assumption has not been benchmarked for coverage.
`min_donors = 1` permits sparse support; a single donor cannot express within-pool
uncertainty. Inspect `eligible_donors` and raise the minimum for your application.

`audit$reconstructed` records the decision to draw a replacement;
`audit$changed` records whether an observed numeric value actually changed.
The latter is missing for an originally missing measurement. These quantities
are different, and neither realized fraction is guaranteed to equal or exceed q.
Only requested cells appear in the audit; unchanged nonrequested rows stay in the
original and completed data. Derived-category row indices are in the manifest.

## Check and reproduce

```sh
Rscript -e 'install.packages(c("testthat", "pkgload", "roxygen2"))'
Rscript -e 'testthat::test_local()'
R CMD build .
R CMD check --no-manual catmech_0.1.0.tar.gz
Rscript tools/benchmark.R 100 10 validation/my-reproduction
```

The benchmark contains five scenarios, six methods, repeated datasets, failures,
warnings, Monte Carlo uncertainty and paired comparisons. Its plan precedes the
run and its source hashes identify the implementation tested. It is exploratory,
not a certification of proper multiple imputation.

## Related work and credit

The package builds on [mice](https://doi.org/10.18637/jss.v045.i03), predictive mean
matching, interval-data imputation such as
[hmi](https://doi.org/10.18637/jss.v095.i09), and the literature on
[multiple overimputation](https://doi.org/10.1177/0049124115585360).
The [prior-art assessment](docs/ASSESSMENT.md) explains the overlap and the narrow
software contribution. No methodological priority or superiority is claimed.
The [research connections](docs/RESEARCH_CONNECTIONS.md) identify useful evaluation
and application ideas from the author's other projects, with their transfer limits.

MIT license. AI-assisted review and implementation are disclosed in
[CONTRIBUTING.md](CONTRIBUTING.md). See [GitHub preparation](docs/GITHUB.md) for upload
instructions and [the integration map](docs/INTEGRATION_MAP.md) for changes from the
supplied local prototype.
