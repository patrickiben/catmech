# Pool a prespecified linear analysis of completed datasets

Uses Rubin's scalar rules with the complete-data residual degrees of
freedom. A result is conditional on the imputation/analysis models and
their assumptions. Failed convergence is rejected. The legacy
retain/replace sensitivity engine requires explicit opt-in and is
labeled as sensitivity output.

## Usage

``` r
catmech_pool_lm(
  result,
  formula,
  conf_level = 0.95,
  allow_sensitivity = FALSE,
  cluster = NULL,
  cluster_variance = c("cv3", "cr1s")
)
```

## Arguments

- result:

  A catmech_result from \[catmech_model()\], \[catmech_pmm()\],
  \[catmech_reference()\], \[catmech_contamination()\] or
  \[catmech_impute()\].

- formula:

  A prespecified formula for \[stats::lm()\]. Missing analysis values
  are rejected rather than silently dropping records.

- conf_level:

  Confidence level between zero and one.

- allow_sensitivity:

  Explicitly allow pooled calculations for the legacy sensitivity engine
  or experimental partial-reference PMM. Results retain a
  sensitivity-only label; this does not establish their inferential
  validity.

- cluster:

  Optional name of a complete column in the completed data identifying
  clusters. Each completed analysis then uses a cluster-robust variance,
  pooled by Rubin's rules with complete-data degrees of freedom of
  clusters minus one. A fit from \`catmech_reference(cluster = )\` uses
  its own column by default. With few clusters, or a few influential
  ones, these intervals can be too narrow. A coefficient determined by a
  single cluster, such as a site type found at one site, draws a
  warning, because its interval can be far too narrow. This variance
  repairs the analysis step only: linear imputation without \`cluster\`
  still undercovered on real items after it.

- cluster_variance:

  With \`cluster\`, \`"cv3"\` (the default) for the delete-one-cluster
  jackknife, or \`"cr1s"\` for the form with factor (G / (G - 1))
  ((n - 1) / (n - p)), the default through version 0.6.0. A planned test
  chose \`"cv3"\`, narrowly, as documented in
  \`validation/VARIANCE_MI_SYNTHESIS.md\`. It gives wider intervals than
  \`"cr1s"\`.

## Value

A data.frame with term, estimate, std.error, df, lower, upper,
within_variance, between_variance, m, n and interpretation, plus
clusters when a cluster-robust variance was used.

## Details

When every audited row agreed and recorded values were kept, the
interpretation ends in \`\_recorded_kept\`, and a warning repeats that
those records are treated as correct. When the workflow declares that no
target cells were imputed, all completed analysis frames must agree and
ordinary complete-data t intervals are returned. A zero estimated
between-imputation variance alone does not activate this rule.
Predictive matching and reference-subset results retain distinct
interpretation labels; none guarantees coverage under misspecification.
Partial-reference PMM requires explicit sensitivity opt-in because
held-out studies found material undercoverage even with random accurate
references. This package policy does not assert that every PMM
application is invalid. Only univariate linear models with at least one
named coefficient are supported.

When the audited value is an analysis predictor, the analysis outcome
must be among the imputation \`predictors\`; otherwise its coefficient
is biased toward zero and a warning says so. When the audited value is
the outcome, the analysis predictors belong in the imputation model
instead.

## Examples

``` r
set.seed(123)
n <- 40
df <- data.frame(
  lab = stats::rnorm(n, 10, 2),
  age = stats::rnorm(n, 50, 10),
  outcome = stats::rnorm(n, 5, 1)
)
df$lab_ref <- df$lab
df$lab_ref[11:n] <- NA
# lab is a predictor below, so the outcome joins the imputation model.
res <- catmech_reference(df, value = "lab", reference = "lab_ref",
                         predictors = c("age", "outcome"), m = 2, seed = 1)
pooled <- catmech_pool_lm(res, outcome ~ lab + age)
print(pooled)
#> -- Catmech Pooled Linear Model --
#> Completions (m): 2 | Observations (n): 40
#> 
#>         term estimate std.error      df   lower  upper
#>  (Intercept)   4.9555    1.3799 10.3797  1.8960 8.0150
#>          lab   0.0499    0.0801 15.4733 -0.1204 0.2203
#>          age  -0.0082    0.0173 15.7352 -0.0449 0.0286
#>               interpretation
#>  conditional_reference_model
#>  conditional_reference_model
#>  conditional_reference_model
```
