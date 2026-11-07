# Canonical prior sensitivity grid for joint marginal cluster models

Automatically computes the prespecified 2-axis prior sensitivity
protocol from the Joint Marginal validation specification:

1.  Covariance scale axis: \\\text{prior\\scale} \in \\0.25, 1.0,
    4.0\\\\, holding regression scale fixed at 10.

2.  Regression scale axis: \\\text{regression\\scale} \in \\5, 10,
    20\\\\, holding covariance scale fixed at 1.0.

## Usage

``` r
catmech_marginal_sensitivity_grid(
  fit,
  term = NULL,
  formula = NULL,
  prior_scales = c(0.25, 1, 4),
  regression_scales = c(5, 10, 20),
  conf_level = 0.95,
  reference = 0
)
```

## Arguments

- fit:

  A fitted `catmech_result` from
  [`catmech_reference`](https://patrickiben.github.io/catmech/reference/catmech_reference.md)
  with `cluster_model = "joint_marginal"`, or a named list of
  precomputed scenarios.

- term:

  Name of the target coefficient to assess. If `NULL`, defaults to the
  first modeled predictor or formula covariate.

- formula:

  Optional formula passed to
  [`catmech_marginal_summary`](https://patrickiben.github.io/catmech/reference/catmech_marginal_summary.md)
  for custom outcome specifications.

- prior_scales:

  Numeric vector of covariance prior scales to evaluate (default
  `c(0.25, 1, 4)`).

- regression_scales:

  Numeric vector of regression prior scales to evaluate (default
  `c(5, 10, 20)`).

- conf_level:

  Credible level (default 0.95).

- reference:

  Reference value for direction checks (default 0).

## Value

An object of class `"catmech_marginal_sensitivity_grid"` containing:

- table:

  Data frame of sensitivity results per scenario.

- axis_summary:

  Summary comparing mean absolute shifts along each axis to MCSE bounds.

- statistic:

  Primary sensitivity statistic:
  \\\max(\text{mean\\abs\\change}\_{\text{cov}},
  \text{mean\\abs\\change}\_{\text{reg}})\\.

- mcse_mean:

  Monte Carlo standard error of the posterior mean under the baseline
  prior.

- term:

  The evaluated coefficient term.

- baseline:

  Baseline scenario result row.

- status:

  Sensitivity status from
  [`catmech_sensitivity`](https://patrickiben.github.io/catmech/reference/catmech_sensitivity.md).

- complete:

  Logical indicating whether all evaluated scenarios succeeded.

- envelope:

  Credible interval envelope across successful scenarios.

- raw_sensitivity:

  Underlying output from
  [`catmech_sensitivity`](https://patrickiben.github.io/catmech/reference/catmech_sensitivity.md).

## Details

Evaluates the average absolute change in the posterior mean under the
alternative hyperparameters along each axis, and returns the primary
sensitivity statistic (the maximum of the two axis-level mean changes)
along with Monte Carlo standard errors.
