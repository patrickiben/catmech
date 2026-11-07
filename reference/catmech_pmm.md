# Audited predictive mean matching after masking quality flags

Masks every flagged target and uses standard type-1 predictive mean
matching from \`mice\`. Categories may be selected as ordinary complete
predictors; they never force interval membership or cause trusted
observations to change.

## Usage

``` r
catmech_pmm(
  data,
  value,
  predictors,
  flag = NULL,
  flag_yes = "Yes",
  flag_no = "No",
  donors = 5L,
  m = 20L,
  seed = 1L
)
```

## Arguments

- data:

  A data frame with unique column names. Original data are retained.

- value:

  Name of the numeric measurement column.

- predictors:

  Explicit complete predictor names, excluding the target and flag.
  Include downstream outcomes when appropriate. A category can be an
  ordinary factor predictor; no category likelihood is imposed.

- flag:

  Optional quality-flag column name.

- flag_yes, flag_no:

  Distinct suspect and non-suspect labels. Missing flags are treated as
  non-suspect, a substantive trust assumption.

- donors:

  Positive number of nearest donors, default five. When any imputation
  is needed, fewer trusted observations than \`donors\` is an error.

- m:

  Number of completions, at least two.

- seed:

  Nonnegative integer seed. The caller's RNG state is restored.

## Value

A \`catmech_result\`. The \`mids\` adapter stores a masked target as its
observed data; the unchanged raw measurements remain in \`original\`.

## Details

Calls \`mice::mice.impute.pmm\` with \`matchtype = 1\` and \`ridge =
1e-5\`, using only originally trusted observations as donors. The
trusted design must have full rank and positive residual degrees of
freedom and variance whenever imputation is required. Completed values
lie in the observed donor support; missing extremes, nonignorable donor
selection and misspecification can bias inference. PMM is not a general
robustness or coverage guarantee.

The delegated method returns donor values, not donor row identifiers.
The audit records \`donor_row\` only when a value matches exactly one
trusted row. For ties, \`donor_row\` is \`NA\`, and \`donor_value_rows\`
lists all value-compatible trusted rows. This list does not identify the
actual nearest-neighbor pool. \`eligible_donors\` counts all trusted
rows; \`candidate_pool_size\` records the requested nearest-donor count
used by the delegated matching algorithm.

## Examples

``` r
set.seed(123)
df <- data.frame(
  x = stats::rnorm(30, 50, 10),
  y = stats::rnorm(30, 10, 2),
  flag = c(rep("No", 25), rep("Yes", 5))
)
res <- catmech_pmm(df, value = "y", predictors = "x", flag = "flag",
                   flag_yes = "Yes", flag_no = "No", donors = 3, m = 2, seed = 1)
summary(res)
#> <summary.catmech_result: trusted_pmm>
#>   Target value:    y
#>   Predictors:      x
#>   Observations:    30
#>   Completions (m): 2
#>   Imputed rows:    5
#>   Convergence:     passed
#> 
#> Parameter diagnostics (first 8 monitors):
#> [1] parameter rhat      ess_bulk  ess_tail 
#> <0 rows> (or 0-length row.names)
```
