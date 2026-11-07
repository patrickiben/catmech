# Pool a prespecified Cox proportional hazards analysis of completed datasets

Fits a Cox proportional hazards model (\[survival::coxph()\]) across
each completed dataset and pools parameter estimates and standard errors
using Rubin's rules with Barnard & Rubin (1999) small-sample degrees of
freedom based on event counts.

## Usage

``` r
catmech_pool_coxph(
  result,
  formula,
  conf_level = 0.95,
  exponentiate = TRUE,
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

  A prespecified formula for \[survival::coxph()\] (e.g. \`Surv(time,
  status) ~ lab_value + age\`).

- conf_level:

  Confidence level between zero and one (default 0.95).

- exponentiate:

  Logical: if \`TRUE\` (default), adds exponentiated estimates and
  confidence intervals (\`estimate_exp\`, \`lower_exp\`, \`upper_exp\`),
  representing Hazard Ratios (HR).

- allow_sensitivity:

  Explicitly allow pooled calculations for sensitivity output.

- cluster:

  Optional name of a column identifying clusters for cluster-robust
  inference.

- ...:

  Additional arguments passed to \[survival::coxph()\] (e.g. \`ties =
  "efron"\`).

## Value

A data frame of class \`c("catmech_pooled_coxph", "data.frame")\`
containing:

- term:

  Model coefficient names.

- estimate:

  Pooled log hazard ratio parameter estimate.

- std.error:

  Total pooled standard error.

- df:

  Barnard & Rubin adjusted degrees of freedom.

- lower, upper:

  Confidence interval bounds on the log hazard scale.

- estimate_exp, lower_exp, upper_exp:

  Exponentiated parameters and bounds (Hazard Ratios) if \`exponentiate
  = TRUE\`.

- within_variance:

  Average within-imputation variance.

- between_variance:

  Between-imputation variance.

- m, n, nevent:

  Number of completions, total observations, and events.

- interpretation:

  Inference label.

## Examples

``` r
set.seed(123)
n <- 50
df <- data.frame(
  time = stats::rexp(n, 0.1),
  status = stats::rbinom(n, 1, 0.5),
  lab = stats::rnorm(n, 10, 2),
  age = stats::rnorm(n, 50, 10)
)
df$lab_ref <- df$lab
df$lab_ref[11:n] <- NA
res <- catmech_reference(df, value = "lab", reference = "lab_ref",
                         predictors = "age", m = 2, seed = 1)
pooled_cox <- catmech_pool_coxph(res, survival::Surv(time, status) ~ lab + age,
                                 exponentiate = TRUE)
print(pooled_cox)
#> -- Catmech Pooled Cox PH Model --
#> Completions (m): 2 | Observations (n): 50 | Events: 27
#> 
#>  term estimate std.error      df   lower  upper estimate_exp lower_exp
#>   lab   0.0701    0.0675 21.3470 -0.0701 0.2103       1.0726    0.9323
#>   age   0.0203    0.0237 20.8258 -0.0290 0.0696       1.0205    0.9714
#>  upper_exp
#>     1.2341
#>     1.0721
```
