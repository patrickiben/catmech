# Compare audit analysis strategies

Fits a user-specified regression model across multiple audit analysis
strategies on a dataset with recorded values and an audit reference
subsample:

1.  **Validation subsample**: Audited records only (discards un-audited
    recorded values).

2.  **Naive analysis**: Trusts un-audited recorded values as correct
    without audit adjustment.

3.  **Linear Multiple Imputation**: Calibrates un-audited values given
    recorded values (\[catmech_reference(use_recorded = TRUE)\]).

4.  **Contamination Model**: Imputes un-audited values under a
    reference-calibrated error model (\[catmech_contamination()\]).

## Usage

``` r
catmech_compare_strategies(
  data,
  value,
  reference,
  formula,
  predictors = NULL,
  family = stats::gaussian(),
  m = 10L,
  seed = 123L,
  conf_level = 0.95,
  agreement_tol = 1,
  cluster = NULL,
  ...
)
```

## Arguments

- data:

  A data frame containing analysis variables, recorded values, and audit
  references.

- value:

  Character string naming the recorded value column.

- reference:

  Character string naming the audit reference column (with \`NA\` for
  un-audited rows).

- formula:

  Model formula for the downstream analysis (e.g. \`outcome ~
  lab_value + age\`).

- predictors:

  Character vector of predictor column names to include in the
  imputation models. If \`NULL\`, non-target predictors from \`formula\`
  present in \`data\` are automatically detected.

- family:

  Error distribution and link function for GLM (default
  \[stats::gaussian()\] for linear models; use \[stats::binomial()\] for
  logistic regression).

- m:

  Number of completed datasets to impute (default 10L).

- seed:

  Random seed for reproducible imputation (default 123L).

- conf_level:

  Confidence level (default 0.95).

- agreement_tol:

  Tolerance for agreement audit in contamination model (default 1.0).

- cluster:

  Optional cluster column name for cluster-robust pooling.

- ...:

  Additional arguments passed to model estimation.

## Value

An object of class \`c("catmech_strategy_comparison", "data.frame")\`
containing:

- strategy:

  Strategy label (\`validation_subsample\`, \`naive_trust_recorded\`,
  \`linear_mi_recorded\`, \`contamination_model\`).

- term:

  Model coefficient term.

- estimate:

  Point estimate.

- std.error:

  Standard error.

- df:

  Degrees of freedom.

- lower, upper:

  Confidence interval bounds.

- ci_width:

  Confidence interval width (upper - lower).

- rel_efficiency:

  Relative precision (inverse variance ratio) relative to the validation
  subsample.

- n_total, n_audited:

  Total observations and audited records count.

## Examples

``` r
set.seed(123)
n <- 50
df <- data.frame(
  outcome = stats::rnorm(n, 10, 2),
  lab = stats::rnorm(n, 5, 1),
  age = stats::rnorm(n, 50, 10)
)
df$ref <- df$lab + stats::rnorm(n, 0, 0.2)
df$ref[16:n] <- NA
comp <- catmech_compare_strategies(df, value = "lab", reference = "ref",
                                   formula = outcome ~ lab + age, m = 2, seed = 1)
print(comp)
#> -- Catmech Audit Strategy Comparison --
#> Cohort Size (N): 50 | Audited Records (n): 15 (30.0% audit fraction)
#> 
#>              strategy        term estimate std.error ci_width rel_efficiency
#>  validation_subsample (Intercept)   7.5921    4.4095  19.2150         1.0000
#>  validation_subsample         lab   0.3596    0.6044   2.6339         1.0000
#>  validation_subsample         age   0.0195    0.0639   0.2785         1.0000
#>  naive_trust_recorded (Intercept)  10.1927    2.1990   8.8477         4.0209
#>  naive_trust_recorded         lab  -0.0657    0.3017   1.2140         4.0128
#>  naive_trust_recorded         age   0.0045    0.0276   0.1111         5.3562
#>    linear_mi_recorded (Intercept)  10.4097    2.2997   9.2673         3.6766
#>    linear_mi_recorded         lab  -0.1013    0.3163   1.2758         3.6523
#>    linear_mi_recorded         age   0.0037    0.0279   0.1125         5.2390
#>   contamination_model (Intercept)  10.4012    2.2959   9.2647         3.6887
#>   contamination_model         lab  -0.1002    0.3192   1.2898         3.5863
#>   contamination_model         age   0.0039    0.0277   0.1115         5.3264
```
