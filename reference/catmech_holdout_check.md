# Score an imputation model on reference rows it did not see

Hides folds of reference values, refits the same workflow and compares
its predictions with the hidden measurements. This is the one check of
imputation calibration available on real data. It uses no simulation
truth.

## Usage

``` r
catmech_holdout_check(
  fit,
  folds = 5L,
  m = 40L,
  level = 0.9,
  seed = 1L,
  fold_by = c("row", "cluster")
)
```

## Arguments

- fit:

  A \`catmech_result\` from \[catmech_contamination()\] or
  \[catmech_reference()\].

- folds:

  Number of folds, at least two.

- m:

  Number of imputations per fold refit, at least twenty.

- level:

  Nominal central coverage level, between zero and one.

- seed:

  Nonnegative integer seed. Caller RNG is preserved.

- fold_by:

  \`"row"\` (default) for stratified row-level cross-validation, or
  \`"cluster"\` for cluster-level holdout cross-validation where all
  reference rows in held-out clusters are hidden together. Requires a
  clustered fit with \`spec\$cluster\`.

## Value

A list with \`status\`, \`summary\` (one row of metrics), \`rows\` (one
per held-out reference), \`flags\`, \`failures\`, \`advice\` and
\`interpretation\`. \`advice\` is nonempty when predictions are no
better than the recorded value or the stratum error rate. On two real
datasets that happened while \`status\` stayed clean, and simpler
methods did as well or better. It is a reading from that study, not a
tested rule.

A list with \`status\`, \`summary\`, \`rows\`, \`flags\`, \`failures\`,
\`advice\` and \`interpretation\`.

## Details

Folds are stratified by whether the reference agreed with the recorded
value. Reported measures include predictive interval coverage,
randomized probability integral transforms, root mean squared error of
the model and of the recorded value, and for the contamination model the
Brier score against the stratum error rate as well as expected against
observed errors. \`status\` is \`miscalibration_detected\` when coverage
is significantly below nominal, the expected and observed error counts
differ at the 1 level. A clean result is not proof of validity: held-out
rows share the reference selection process, the check has limited power
with few errors, and it does not address nonignorable selection or a
reference that is itself wrong.

## Examples

``` r
set.seed(123)
df <- data.frame(val = stats::rnorm(30, 10, 2), age = stats::rnorm(30, 50, 10))
df$ref <- df$val + stats::rnorm(30, 0, 0.2)
res <- catmech_reference(df, value = "val", reference = "ref", predictors = "age", m = 2, seed = 1)
chk <- catmech_holdout_check(res, folds = 2L, m = 20L, seed = 1L)
print(chk)
#> -- Catmech Holdout Cross-Validation Check --
#> Status:              no_miscalibration_detected
#> Folds completed:     2 (holdout rows: 30, partitioned by: row)
#> Predictive coverage: 0.8333 (nominal: 0.9) (uniform PIT KS p-value: 0.3983)
#> Predictive RMSE:     model = 2.2008 vs recorded = 0.1711
#> Mean interval width: 6.8334
#> Flags:               none
#> 
#> Advice:
#>   * Predictions are no closer to held-out references than the recorded values are.
#>   * Consider catmech_reference(use_recorded = TRUE) or generalized raking.
```
