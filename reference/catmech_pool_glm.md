# Pool a prespecified generalized linear analysis of completed datasets

Fits a generalized linear model (\[stats::glm()\]) across each completed
dataset and pools parameter estimates and standard errors using Rubin's
rules with Barnard & Rubin (1999) small-sample degrees of freedom.

## Usage

``` r
catmech_pool_glm(
  result,
  formula,
  family = stats::binomial(),
  conf_level = 0.95,
  exponentiate = FALSE,
  allow_sensitivity = FALSE,
  cluster = NULL,
  ...
)
```

## Arguments

- result:

  A \`catmech_result\` from \[catmech_model()\], \[catmech_pmm()\],
  \[catmech_reference()\], \[catmech_contamination()\] or
  \[catmech_impute()\].

- formula:

  A prespecified formula for \[stats::glm()\].

- family:

  A description of the error distribution and link function to be used
  in the model. Defaults to \`stats::binomial()\` (for logistic
  regression). Can be a character string naming a family function or the
  result of a call to a family function.

- conf_level:

  Confidence level between zero and one (default 0.95).

- exponentiate:

  Logical: if \`TRUE\`, adds exponentiated estimates and confidence
  intervals (\`estimate_exp\`, \`lower_exp\`, \`upper_exp\`), commonly
  interpreted as odds ratios (for logistic regression) or incidence rate
  ratios (for Poisson regression).

- allow_sensitivity:

  Explicitly allow pooled calculations for sensitivity output.

- cluster:

  Optional name of a column identifying clusters for cluster-robust
  inference. Each completed analysis then uses the CR1S sandwich, factor
  (G / (G - 1)) ((n - 1) / (n - p)), built from the unscaled inverse
  information so that it does not depend on an estimated dispersion, and
  pooling uses complete-data degrees of freedom of clusters minus one.

- ...:

  Additional arguments passed to \[stats::glm()\] (e.g. \`control =
  list(maxit = 50)\`).

## Value

A data frame of class \`c("catmech_pooled_glm", "data.frame")\`
containing:

- term:

  Model coefficient names.

- estimate:

  Pooled log-scale parameter estimate.

- std.error:

  Total pooled standard error.

- df:

  Barnard & Rubin adjusted degrees of freedom.

- lower, upper:

  Confidence interval bounds on the link scale.

- estimate_exp, lower_exp, upper_exp:

  Exponentiated parameters and bounds (if \`exponentiate = TRUE\`).

- within_variance:

  Average within-imputation variance.

- between_variance:

  Between-imputation variance.

- m, n:

  Number of completions and total observations.

- family, link:

  The GLM family and link function used.

- interpretation:

  Inference label.

## Details

When the audited value is an analysis predictor, include the analysis
outcome among the imputation \`predictors\`. Imputing a predictor
without the outcome biases its coefficient toward zero; a warning says
so. In a check recorded in \`validation/REVIEW_FIXES_2026-11-07.md\`, a
logistic log odds ratio of 0.80 came out at 0.61, with 65% coverage of
nominal 95% intervals, when the outcome was left out.

## Examples

``` r
set.seed(123)
n <- 50
df <- data.frame(
  lab = stats::rnorm(n, 10, 2),
  age = stats::rnorm(n, 50, 10),
  y = stats::rbinom(n, 1, 0.4)
)
df$lab_ref <- df$lab
df$lab_ref[11:n] <- NA
# lab is a predictor below, so the outcome y joins the imputation model.
res <- catmech_reference(df, value = "lab", reference = "lab_ref",
                         predictors = c("age", "y"), m = 2, seed = 1)
pooled_glm <- catmech_pool_glm(res, y ~ lab + age, family = stats::binomial(),
                               exponentiate = TRUE)
print(pooled_glm)
#> -- Catmech Pooled GLM: binomial (logit link) --
#> Completions (m): 2 | Observations (n): 50
#> 
#>         term estimate std.error     df    lower   upper estimate_exp lower_exp
#>  (Intercept)  -8.5357    6.8331 1.6939 -43.6640 26.5926       0.0002    0.0000
#>          lab   1.2050    0.8170 1.3122  -4.8278  7.2379       3.3369    0.0080
#>          age  -0.0776    0.0624 4.0811  -0.2494  0.0942       0.9253    0.7792
#>     upper_exp
#>  3.540139e+11
#>  1.391153e+03
#>  1.098800e+00
```
