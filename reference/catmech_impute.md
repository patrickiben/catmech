# Generate auditable sensitivity completions

Supports one continuous target with complete predictors. Each completion
uses an independent parameter and donor draw; no iterative chain is
needed for this restricted setting. Categories must be observed, except
that a missing category can be derived from an unflagged observed value.
Predictors must be selected explicitly to avoid accidental use of
identifiers or simulation truth. Include downstream outcomes when
appropriate for the intended analysis. Completed data and Rubin pooling
are experimental and do not establish identification, proper MI, or
interval coverage.

## Usage

``` r
catmech_impute(data, spec, predictors, m = 20L, seed = 1L)
```

## Arguments

- data:

  A data.frame. The original is stored without modification.

- spec:

  A \[catmech_spec()\] or equivalent validated list.

- predictors:

  Character vector of complete predictor column names. Use character(0)
  for an intercept-only model. The category, target and flag cannot be
  predictors. Category information restricts the donor support.

- m:

  Number of completed datasets, at least 2.

- seed:

  Nonnegative integer random seed. The caller's RNG is restored.

## Value

A catmech_result containing original, completed, mids, audit,
diagnostics, spec and manifest. mids is an analysis adapter created from
the supplied completions, not evidence of mice chain convergence.

## Examples

``` r
d <- simulate_catmech_lab(150, seed = 12)
d <- d[setdiff(names(d), c("lab_true", "corrupted"))]
s <- catmech_spec("lab_value", "lab_category", c("Low", "Medium", "High"),
  c(-Inf, 5.5, 7, Inf), flag = "outlier_flag", corruption_probability = 1,
  conflict_policy = "reconstruct")
fit <- catmech_impute(d, s, c("age", "sex", "inflammation", "outcome"), m = 3)
fit$diagnostics
#>   imputation requested reconstructed changed threshold_conflicts fallback_draws
#> 1          1        42            42      10                   8              0
#> 2          2        42            42      10                   8              0
#> 3          3        42            42      10                   8              0
#>   threshold_violations
#> 1                    0
#> 2                    0
#> 3                    0
```
