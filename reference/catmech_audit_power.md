# Prospective sample size and statistical power calculation for clinical data audits

Computes the statistical power and sample size required to achieve a
target confidence interval precision ratio (\\\rho \le \rho\_{target}\\)
in a two-phase clinical audit. Can evaluate power prospectively from
prior clinical error parameters or empirically from a fitted
\[catmech_audit_plan()\].

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

  Optional fitted \`catmech_audit_plan\` object. If supplied, the
  empirical distribution of \\k\\ is used to evaluate power.

- N:

  Total cohort size (number of patient records). Required if \`plan\` is
  \`NULL\`.

- target_ratio:

  Maximum acceptable confidence interval width relative to an ideal full
  audit (default 1.25, meaning within 25% of the full-audit width).

- power_target:

  Desired statistical power (default 0.80).

- error_rate:

  Anticipated error or corruption rate \\p \in (0, 1)\\ (default 0.15).
  Ignored if \`plan\` is supplied.

- correlation:

  Anticipated correlation between recorded and true values (default
  0.85). Ignored if \`plan\` is supplied.

- r_squared:

  Anticipated auxiliary predictor \\R^2\\ (default 0.30). Ignored if
  \`plan\` is supplied.

- fractions:

  Numeric vector of audit fractions to evaluate, strictly between 0
  and 1. Defaults to \`seq(0.02, 0.40, by = 0.02)\`.

- n_sim:

  Number of Monte Carlo draws for evaluating the sampling distribution
  of the precision ratio (default 2000).

- seed:

  Optional integer seed for reproducible power evaluation.

## Value

An object of class \`catmech_audit_power\` containing:

- table:

  Data frame of evaluated audit fractions, audit sample sizes, expected
  relative widths, and statistical power.

- recommended:

  Recommended minimum audit fraction and sample size reaching
  \`power_target\`.

- target_ratio:

  The target precision ratio evaluated.

- power_target:

  The target statistical power.

- N:

  Total cohort size.

- k_summary:

  Point estimate and uncertainty bounds for the two-phase variance
  factor \\k\\.

## Examples

``` r
pwr <- catmech_audit_power(N = 500, target_ratio = 1.3, n_sim = 200, seed = 1)
print(pwr)
#> -- Catmech Audit Power & Sample Size Calculator --
#> Cohort Size (N):    500
#> Target Ratio:       <= 1.30 (within 30% of full audit width)
#> Target Power:       80%
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
#> Achieves 86.0% power (expected ratio: 1.189)
```
