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
  imputation models. If \`NULL\`, every variable in \`formula\` other
  than \`value\` and \`reference\` is used, including the outcome. When
  the audited value is an analysis predictor, leaving the outcome out of
  the imputation model biases its coefficient toward zero; a supplied
  vector without the outcome draws a warning from the pooling step.

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

  Optional complete cluster column. Every strategy then uses a
  cluster-robust variance with complete-data degrees of freedom of
  clusters minus one: the delete-one-cluster jackknife for linear
  models, as in \[catmech_pool_lm()\], and the CR1S sandwich for other
  families, as in \[catmech_pool_glm()\]. It is not passed to
  \`catmech_reference(cluster = )\`, whose clustered imputation biased
  slopes toward zero in simulations when the audited value was a
  predictor (\`validation/VARIANCE_MI_SYNTHESIS.md\`); a formula term
  naming the same column is still used as an ordinary predictor.

- ...:

  Currently unused.

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

The imputation predictors used are stored in the
\`imputation_predictors\` attribute. If the contamination model fails,
its rows are \`NA\` and a warning gives the reason.

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
#> Imputation predictors: outcome, age
#> 
#>              strategy        term estimate std.error ci_width rel_efficiency
#>  validation_subsample (Intercept)   7.5921    4.4095  19.2150         1.0000
#>  validation_subsample         lab   0.3596    0.6044   2.6339         1.0000
#>  validation_subsample         age   0.0195    0.0639   0.2785         1.0000
#>  naive_trust_recorded (Intercept)  10.1927    2.1990   8.8477         4.0209
#>  naive_trust_recorded         lab  -0.0657    0.3017   1.2140         4.0128
#>  naive_trust_recorded         age   0.0045    0.0276   0.1111         5.3562
#>    linear_mi_recorded (Intercept)   9.8453    2.8575  15.3433         2.3813
#>    linear_mi_recorded         lab  -0.0121    0.4368   2.8823         1.9145
#>    linear_mi_recorded         age   0.0062    0.0278   0.1123         5.2680
#>   contamination_model (Intercept)  10.1765    2.2185   8.9360         3.9505
#>   contamination_model         lab  -0.0629    0.3059   1.2323         3.9030
#>   contamination_model         age   0.0046    0.0276   0.1112         5.3567
```
