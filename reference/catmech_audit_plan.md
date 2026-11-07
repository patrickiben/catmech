# Plan how many records to audit from a random pilot

Answers "how many records should we audit?" for linear multiple
imputation with the recorded value, \[catmech_reference()\] with
\`use_recorded = TRUE\`. It predicts how much wider the interval for one
regression coefficient will be at each audit fraction than with a full
audit, using the classical two-phase form \`width(f) = width_full \*
sqrt(1 + k \* (1 - f) / f)\`. It predicts precision only. It says
nothing about bias or coverage, and it does not replace
\[catmech_recording_audit()\].

## Usage

``` r
catmech_audit_plan(
  data,
  value,
  reference,
  predictors,
  formula,
  term,
  fractions = c(0.05, 0.1, 0.25, 0.5, 0.75),
  target_ratio = 1.25,
  conf_level = 0.95,
  bootstrap = 500L,
  m = 20L,
  seed = 1L
)
```

## Arguments

- data:

  A data frame holding every record in the cohort, audited or not.

- value:

  Name of the numeric recorded value, with \`NA\` where nothing was
  recorded.

- reference:

  Name of the numeric reference value, observed on the pilot audit only.
  The pilot must be a random sample of the records.

- predictors:

  Complete auxiliary columns available on every record, as in
  \[catmech_reference()\]. Include the analysis outcome.

- formula:

  The planned linear analysis, as passed to \[catmech_pool_lm()\]. It
  must be linear in \`value\`, so transform the column first rather than
  writing \`log(value)\`.

- term:

  Name of the coefficient the audit is meant to protect.

- fractions:

  Audit fractions to report, each above zero and at most one.

- target_ratio:

  Width relative to a full audit that counts as good enough, above one.
  Default 1.25, meaning within 25% of the full-audit width.

- conf_level:

  Level of the analysis interval and of the plan's intervals.

- bootstrap:

  Bootstrap replicates of the pilot, at least 100.

- m:

  Imputations for the pilot's own fits, at least two.

- seed:

  Nonnegative integer seed. The caller's RNG state is restored.

## Value

A list with:

- \`status\`:

  \`"no_flags"\` or \`"review_flags"\`.

- \`flags\`:

  Short names of every warning raised.

- \`k\`:

  One row each for \`recorded\` (keep the recorded value, linear MI) and
  \`discard\` (the \[catmech_reference()\] default): \`estimate\`,
  \`lower\`, \`upper\`, \`classical\`, the empirical version described
  below, and \`rubin_check\`, the value implied by the pilot fit's own
  between-imputation variance.

- \`table\`:

  One row per fraction: \`fraction\`, \`records\`,
  \`additional_records\` beyond the pilot, \`small_audit\`, which marks
  a size the function would refuse as a pilot, \`relative_width\` with
  \`relative_lower\` and \`relative_upper\` from the bootstrap interval
  for \`k\`, \`predictive_lower\` and \`predictive_upper\` for the width
  a single audit of that size may realize, \`width\` with
  \`width_lower\` and \`width_upper\` anchored on the pilot fit, and the
  same for discarding recorded values: \`discard_relative_width\`,
  \`discard_relative_lower\`, \`discard_relative_upper\` and
  \`discard_width\`.

- \`answer\`:

  Records needed to reach \`target_ratio\`: \`fraction\` and \`records\`
  from the estimate, \`records_lower\` and \`records_upper\` from the
  interval for \`k\`, and \`smallest_grid_fraction\`, the smallest of
  \`fractions\` predicted to reach it, or \`NA\`.

- \`pilot\`:

  The pilot's size, disagreements, fits and diagnostics.

- \`assumptions\`:

  The assumptions the plan rests on, as text.

## Details

Classical theory for an efficient two-phase estimator gives \`Var(f) /
Var(1) = 1 + k \* (1 - f) / f\` with \`k = Var(psi - E\[psi \| aux\]) /
Var(psi)\`, where \`psi\` is the influence function of the coefficient
using true values and \`aux\` is what is known on every record. This is
not a new result.

The planner evaluates that ratio under the working models, because the
interval from Rubin's rules comes from them. The numerator is \`Var(psi
\| aux)\` under the normal linear imputation model fitted on the pilot's
reference rows, averaged over all records. The denominator is
\`Var(psi)\` under the homoskedastic analysis model. The result has a
closed form and needs no simulation. The \`classical\` column instead
uses empirical variances on the pilot rows. On unequal-variance data the
two differ: the first describes the reported interval and the second the
true sampling variance. Use \`estimate\` to plan interval width.

\`width\` is anchored on the pilot's own fit: its within-imputation
variance estimates the full-audit variance, and \`relative_width\`
scales the result. The pilot's total width is not the anchor, because
its between-imputation variance is too noisy with few imputations.

Assumptions: the pilot is a simple random sample of these records; the
analysis is linear MI with Rubin's rules through \[catmech_pool_lm()\];
the reference is exact and measures the same quantity; the later audit
is random and is imputed from the same model. The formula is a
large-sample one. It understates small audits, most of all when
discarding recorded values, because little data then drives a large
share of the variance.

The intervals answer different questions. \`relative_lower\` and
\`relative_upper\` come from a bias-corrected and accelerated bootstrap
for \`k\`, so they describe the expected width. A single small audit can
land far from that, because its own residual variance is noisy;
\`predictive_lower\` and \`predictive_upper\` try to include this. A
pilot can miss rare errors altogether, and then both ranges are too low
and too narrow, so neither is a guarantee.

The function refuses a pilot with fewer than 30 rows holding both values
or fewer than 20 residual degrees of freedom in the imputation model,
and a pilot design that is rank deficient. It warns below 100 rows or 10
disagreements, when no pilot row disagrees, when the pilot looks
selective on the recorded value (either z above 3.29, for its mean or
spread), and when the outcome is missing from \`predictors\`. It also
warns when gross errors could distort a linear fit: a pilot discrepancy
above 5 robust SDs of the reference, or a recorded value above 10 robust
SDs from the median. Linear MI failed in the synthetic gross-error cell
of the contamination study, so do not rely on a plan carrying that
warning.

## Examples

``` r
set.seed(123)
n <- 120
df <- data.frame(
  outcome = stats::rnorm(n, 10, 2),
  lab = stats::rnorm(n, 5, 1),
  age = stats::rnorm(n, 50, 10)
)
df$ref <- df$lab + stats::rnorm(n, 0, 0.2)
df$ref[41:n] <- NA
plan <- catmech_audit_plan(df, value = "lab", reference = "ref",
                           predictors = c("age", "outcome"), formula = outcome ~ lab + age,
                           term = "lab", bootstrap = 100L, m = 2L, seed = 1L)
#> Warning: Small pilot: k rests on 40 rows, so its interval is wide. Below 100 rows, plan on the upper limit.
summary(plan)
#> -- Catmech Audit Size Plan --
#> Target ratio:       <= 1.25x full-audit width
#> Recommended audit:  8 records (fraction: 0.067) [95% CI: 6, 13]
#> Pilot audit:        40 / 120 records (disagreements: 40)
#> Status:             review_flags
#> Flags:              small_pilot
#> 
#> Relative width across planned fractions:
#>  fraction records additional_records relative_width relative_lower
#>      0.05       6                  0          1.328          1.246
#>      0.10      12                  0          1.167          1.123
#>      0.25      30                  0          1.058          1.043
#>      0.50      60                 20          1.020          1.014
#>      0.75      90                 50          1.007          1.005
#>  relative_upper predictive_lower predictive_upper
#>           1.487            1.001            3.762
#>           1.254            1.000            2.331
#>           1.091            1.000            1.295
#>           1.031            1.000            1.126
#>           1.011            1.000            1.055
```
