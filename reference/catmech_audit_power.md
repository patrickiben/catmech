# Rough audit-size planning aid (not validated)

Evaluates, over a grid of audit fractions, how often the expected
interval width stays within a target ratio of a full audit's width
(\\\rho \le \rho\_{target}\\), using the two-phase relation \\\rho =
\sqrt{1 + k (1 - f) / f}\\. The factor \\k\\ comes from a fitted
\[catmech_audit_plan()\] or from assumed error parameters.

## Usage

``` r
catmech_audit_power(
  plan = NULL,
  N = NULL,
  target_ratio = 1.25,
  power_target = 0.8,
  error_rate = 0.15,
  correlation = 0.85,
  r_squared = 0.3,
  fractions = seq(0.02, 0.4, by = 0.02),
  n_sim = 2000L,
  seed = 42L
)
```

## Arguments

- plan:

  Optional fitted \`catmech_audit_plan\` object. If supplied, \\k\\ is
  drawn from a log-normal distribution fitted to the plan's estimate and
  interval for \\k\\.

- N:

  Total cohort size (number of patient records). Required if \`plan\` is
  \`NULL\`.

- target_ratio:

  Maximum acceptable confidence interval width relative to an ideal full
  audit (default 1.25, meaning within 25% of the full-audit width).

- power_target:

  Target share of plausible \\k\\ values whose width ratio meets
  \`target_ratio\` (default 0.80); see Details for why this is not the
  power of a test.

- error_rate:

  Anticipated error or corruption rate \\p \in (0, 1)\\ (default 0.15).
  Ignored if \`plan\` is supplied.

- correlation:

  Anticipated correlation between recorded and true values over all
  records (default 0.85). Ignored if \`plan\` is supplied.

- r_squared:

  Anticipated auxiliary predictor \\R^2\\ (default 0.30). Ignored if
  \`plan\` is supplied.

- fractions:

  Numeric vector of audit fractions to evaluate, strictly between 0
  and 1. Defaults to \`seq(0.02, 0.40, by = 0.02)\`.

- n_sim:

  Number of Monte Carlo draws of \\k\\ (default 2000).

- seed:

  Optional integer seed for reproducible draws.

## Value

An object of class \`catmech_audit_power\` containing:

- table:

  Data frame of evaluated audit fractions, audit sizes, expected width
  ratios with 95% ranges, and \`power\`, the share of \\k\\ draws
  meeting \`target_ratio\`.

- recommended:

  Smallest evaluated audit fraction and size whose \`power\` reaches
  \`power_target\`, or the largest fraction if none does.

- target_ratio:

  The target width ratio evaluated.

- power_target:

  The target share of \\k\\ draws.

- N:

  Total cohort size.

- k_summary:

  Point estimate and uncertainty bounds for the two-phase variance
  factor \\k\\.

- linear_mi:

  Without \`plan\`: the linear-imputation \\k\\ described above and the
  audit fraction and size it needs to meet \`target_ratio\`; \`NULL\`
  with \`plan\`.

- status:

  Always \`"not_validated"\`.

## Details

**Status: not validated.** No study has checked this function's
recommendations, and the planner it builds on, \[catmech_audit_plan()\],
returned NO-GO under a rule fixed in advance. Treat the output as a
rough guide, not a sample-size guarantee.

Without \`plan\`, \\k = q / (1 - q)\\ with \\q\\ = \`error_rate\` (1 -
\`correlation\`^2) (1 - \`r_squared\`), where \`correlation\` is the
correlation between recorded and true values over all records.
Multiplying by the error rate treats the unexplained variation as if it
sat in erroneous records the analysis can identify and correct. That is
optimistic: on ClinicalTrials.gov and NHANES the contamination model
could not tell which records were wrong (AUC 0.54 and 0.50). For linear
imputation given the recorded value, the package's recommended start,
the double-sampling regression estimator has \\k \approx (1 -
\rho^2)(1 - R^2)\\ at the same correlation \\\rho\\, about 1 /
\`error_rate\` times larger when errors are rare. The result reports the
audit fraction that this \\k\\ needs (\`linear_mi\`), and the printout
shows both. That figure ignores the uncertainty in \\k\\, so it can fall
below the main recommendation, for example at very high error rates.

"Power" here is the share of draws of \\k\\, from an assumed log-normal
spread around its estimate, whose width ratio meets the target. It is
not the power of a hypothesis test.

## Examples

``` r
pwr <- catmech_audit_power(N = 500, target_ratio = 1.3, n_sim = 200, seed = 1)
print(pwr)
#> -- Catmech Audit-Size Planning Aid (not validated) --
#> Cohort Size (N):    500
#> Target Ratio:       <= 1.30 (within 30% of full audit width)
#> Target share of k:  80%
#> Two-Phase Factor k: 0.030 [0.010, 0.117] (prior_specification)
#> 
#>  fraction n_audit expected_ratio ratio_lower ratio_upper power
#>      0.02      10          1.646       1.226       2.591 0.105
#>      0.04      20          1.360       1.116       1.949 0.480
#>      0.06      30          1.249       1.077       1.681 0.735
#>      0.08      40          1.189       1.057       1.530 0.860
#>      0.10      50          1.151       1.045       1.432 0.920
#>      0.12      61          1.125       1.037       1.362 0.955
#>      0.14      70          1.106       1.031       1.310 0.970
#>      0.16      80          1.092       1.027       1.270 0.990
#>      0.18      90          1.080       1.023       1.237 1.000
#>      0.20     100          1.071       1.020       1.211 1.000
#>      0.22     110          1.063       1.018       1.189 1.000
#>      0.24     120          1.056       1.016       1.170 1.000
#>      0.26     130          1.051       1.015       1.154 1.000
#>      0.28     140          1.046       1.013       1.140 1.000
#>      0.30     151          1.042       1.012       1.128 1.000
#>      0.32     160          1.038       1.011       1.117 1.000
#>      0.34     170          1.035       1.010       1.107 1.000
#>      0.36     181          1.032       1.009       1.099 1.000
#>      0.38     190          1.030       1.008       1.091 1.000
#>      0.40     200          1.027       1.008       1.084 1.000
#> 
#> Recommended Audit Size: n = 40 records (8.0% audit fraction)
#> Meets the target in 86.0% of k draws (expected ratio: 1.189)
#> Linear imputation (double-sampling regression, k = 0.194): about 22.0% (n = 110) for ratio <= 1.30, ignoring uncertainty in k.
#> Power here is the share of plausible k values meeting the target, not the power of a test.
#> Not validated; see ?catmech_audit_power for the assumptions.
```
