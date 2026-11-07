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
  inference.

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
res <- catmech_reference(df, value = "lab", reference = "lab_ref",
                         predictors = "age", m = 2, seed = 1)
pooled_glm <- catmech_pool_glm(res, y ~ lab + age, family = stats::binomial(),
                               exponentiate = TRUE)
print(pooled_glm)
#> -- Catmech Pooled GLM: binomial (logit link) --
#> Completions (m): 2 | Observations (n): 50
#> 
#>         term estimate std.error      df   lower  upper estimate_exp lower_exp
#>  (Intercept)  -1.2630    2.6105  3.8111 -8.6550 6.1289       0.2828    0.0002
#>          lab   0.1023    0.1897  2.4596 -0.5839 0.7886       1.1078    0.5577
#>          age  -0.0035    0.0367 44.8465 -0.0775 0.0705       0.9965    0.9254
#>  upper_exp
#>   458.9272
#>     2.2004
#>     1.0731
```
