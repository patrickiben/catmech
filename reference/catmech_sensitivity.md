# Inspect a finite set of prespecified imputation assumptions

Inspired by assumption-path audits: keep failed and unresolved scenarios
visible, and do not select the most favorable result. The interval
envelope is descriptive; it is not a confidence interval over arbitrary
misspecification. It covers only the supplied finite set and inherits
each member's limitations.

## Usage

``` r
catmech_sensitivity(
  results,
  formula = NULL,
  term,
  reference = 0,
  conf_level = 0.95,
  expected_scenarios = NULL,
  cluster = NULL,
  cluster_variance = c("cv3", "cr1s"),
  method = c("pool_lm", "pool_glm", "pool_coxph", "marginal_summary", "joint_summary"),
  family = stats::binomial()
)
```

## Arguments

- results:

  A uniquely named list of catmech_result objects or error conditions
  captured with tryCatch(..., error = identity). Supply every
  prespecified scenario, including failures.

- formula:

  A prespecified linear analysis formula, identical across scenarios.
  When \`method = "marginal_summary"\`, \`formula\` is optional and
  defaults to the population regression of the true measurement on the
  continuous predictors; it can also specify a custom response and
  covariate subset (e.g. \`outcome ~ recorded + x\`). When \`method =
  "joint_summary"\`, \`formula\` is not used.

- term:

  Name of the coefficient to compare.

- reference:

  Reference value, usually zero, for interval-direction reporting.

- conf_level:

  Conditional confidence level passed to \[catmech_pool_lm()\] or
  \[catmech_marginal_summary()\] / \[catmech_joint_summary()\].

- expected_scenarios:

  Optional character vector naming the complete prespecified set.
  Omitted members are recorded as failures; unplanned extras are
  rejected. Without this argument only supplied scenarios can be
  checked.

- cluster, cluster_variance:

  Passed to \[catmech_pool_lm()\], so that fits made with and without
  \`catmech_reference(cluster = )\` can be compared under the same
  variance. Without \`cluster\`, a fit with its own cluster column uses
  a cluster-robust variance and a fit without one a model-based
  variance, and their widths are not comparable.

- method:

  Estimation method: \`"pool_lm"\` (default) for multiple-imputation
  pooling of linear fits, \`"pool_glm"\` for generalized linear model
  pooling, \`"pool_coxph"\` for Cox proportional hazards model pooling
  via \[catmech_pool_coxph()\], \`"marginal_summary"\` for population
  marginal inference from \[catmech_marginal_summary()\], or
  \`"joint_summary"\` for fixed-effect summaries from
  \[catmech_joint_summary()\].

- family:

  Error distribution and link function for \`method = "pool_glm"\`.
  Defaults to \`stats::binomial()\`.

## Value

A list with table, status, complete, successful_envelope and
interpretation. An incomplete scenario set always has an unresolved
status.

## Examples

``` r
set.seed(123)
df <- data.frame(
  lab = stats::rnorm(30, 10, 2),
  age = stats::rnorm(30, 50, 10),
  y = stats::rnorm(30, 5, 1)
)
df$ref <- df$lab
df$ref[11:30] <- NA
res1 <- catmech_reference(df, value = "lab", reference = "ref",
                          predictors = "age", m = 2, seed = 1)
res2 <- catmech_reference(df, value = "lab", reference = "ref",
                          predictors = "age", use_recorded = TRUE, m = 2, seed = 2)
#> Warning: All 10 audited rows agree with their recorded values, so 20 unaudited recorded values were kept with no imputation uncertainty. An audit of this size cannot rule out an error rate of about 30.0%, and intervals here treat those records as correct.
sens <- catmech_sensitivity(
  list(standard = res1, recorded = res2),
  formula = y ~ lab + age, term = "lab"
)
#> Warning: Every audited row agreed, so unaudited recorded values were kept, and these intervals treat them as correct.
print(sens)
#> -- Catmech Sensitivity Analysis --
#> Status:              unresolved_in_all_scenarios (complete: TRUE)
#> Successful envelope: [-0.114, 0.5184]
#> Scenarios:           2 evaluated
#> 
#>  scenario status estimate std.error   lower  upper          direction
#>  standard     ok   0.2246    0.1256 -0.0692 0.5184 includes_reference
#>  recorded     ok   0.0564    0.0828 -0.1140 0.2268 includes_reference
```
