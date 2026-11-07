# Complete a target using a reference-calibrated recording-error model

Uses the recorded measurement without trusting it. Each recorded value
is either an exact recording of the true target or a recording with
additive normal error. Rows with both a recorded value and a reference
calibrate that mixture, and every other recorded value receives a
posterior probability of being wrong. This keeps the information in
correct records, which \[catmech_reference()\] discards by default.

## Usage

``` r
catmech_contamination(
  data,
  value,
  reference,
  predictors,
  corruption_strata = NULL,
  error_strata = NULL,
  agreement_tol = 0,
  reference_sd = 0,
  selection_delta = 0,
  prior = list(),
  m = 20L,
  seed = 1L,
  burnin = 200L,
  iterations = 300L,
  max_iterations = 1200L,
  chains = 4L
)
```

## Arguments

- data:

  A data frame; raw target measurements remain in \`original\`.

- value:

  Name of the numeric recorded target to replace in completions.

- reference:

  Name of a numeric reference column, observed only where independently
  measured.

- predictors:

  Explicit complete predictor names, excluding \`value\` and
  \`reference\`. Include the scientific outcome and relevant design
  variables.

- corruption_strata:

  Optional categorical columns with separate corruption rates, usually a
  quality flag and any site or source effect. These must be complete on
  recorded rows and have calibration rows.

- error_strata:

  Optional categorical columns with a separate error mean and standard
  deviation, such as a site with a systematic recording shift.

- agreement_tol:

  Recorded values within this distance of the reference count as
  correct, for example half a rounding unit. The default requires
  numerical equality apart from a relative fuzz of \`1e-8\`.

- reference_sd:

  Declared standard deviation of independent reference measurement
  error. The default treats references as exact. A positive value is an
  assumption to vary in a sensitivity analysis; it is not estimated. It
  cannot be combined with \`agreement_tol\`.

- selection_delta:

  Sensitivity shift added to the target mean of unreferenced rows. Zero
  assumes ignorable reference selection. Compare several values with
  \[catmech_sensitivity()\]; the data cannot identify it.

- prior:

  Optional named list overriding \`pi_a\` and \`pi_b\` (Beta prior for
  each corruption rate, default 0.5 and 0.5), \`mu0\`, \`kappa0\`,
  \`nu0\` and \`tau0\` (normal-scaled-inverse-chi-squared prior for each
  error mean and variance, defaults 0, 0.1, 2 and the residual SD of the
  reference regression).

- m:

  Number of completions, at least two.

- seed:

  Nonnegative integer seed, with caller RNG restored.

- burnin:

  Discarded Gibbs iterations per chain.

- iterations:

  Retained iterations per chain, at least 100.

- max_iterations:

  If the convergence screen fails, the sampler restarts with twice as
  many iterations, up to this limit, and then refuses. The manifest
  records \`iterations_used\`. Set equal to \`iterations\` to disable
  restarts.

- chains:

  Number of Gibbs chains, at least four.

## Value

A \`catmech_result\`. Its \`contamination\` element holds \`corruption\`
(one row per record, with \`posterior_corruption\`), \`parameters\`
(posterior summaries of the recording and target models) and
\`recording_audit\` (from \[catmech_recording_audit()\]).

## Details

The model is \`T \| Z ~ Normal(Z beta, sigma^2)\`, \`A ~
Bernoulli(pi\[stratum\])\`, \`W = T\` when \`A = 0\`, and \`W \| T ~
Normal(T + mu\[stratum\], tau\[stratum\]^2)\` when \`A = 1\`. Each scan
updates the corruption rates and error parameters by slice sampling with
\`A\` and \`T\` integrated out, draws \`A\` and then \`T \| A\`, and
draws the regression parameters from their conjugate conditional. With
exact references, \`A\` is known on calibration rows. The reference
design must have full column rank, which keeps the regression posterior
proper under \`p(beta, sigma^2)\` proportional to \`1/sigma^2\`.

Reference selection and missing recorded values must be ignorable given
the predictors, recorded value and strata. Selection may therefore
depend on the recorded value, as in a targeted audit of unusual records;
\[catmech_reference()\] supports that only with \`use_recorded = TRUE\`.
Errors must be nondifferential given the strata: they may not depend on
the outcome or true value once the strata are fixed.

Recorded values that look unusual relative to \`Z beta\` are more likely
to be classified as errors, so a non-normal or nonlinear target model
can misclassify genuine extremes. Flags computed from the recorded value
itself double count that information. Run \[catmech_holdout_check()\]
and compare with \[catmech_reference()\] rather than assuming the model
holds. This is a contamination model used with validation data, not a
new identification result. Without calibration rows the call is refused.

## Examples

``` r
d <- simulate_catmech_lab(200, seed = 8, corruption_probability = .5)
d$gold <- NA_real_
d$gold[seq(1, 200, 4)] <- d$lab_true[seq(1, 200, 4)]
d[c("lab_true", "corrupted")] <- NULL
fit <- catmech_contamination(d, "lab_value", "gold",
  c("age", "inflammation", "outcome"), corruption_strata = "outlier_flag",
  m = 5, seed = 12, iterations = 200)
fit$contamination$parameters[1:3, ]
#>           parameter       mean          sd      lower      upper
#> 1  beta:(Intercept) 1.65684791 0.320879535 0.98405733 2.26290704
#> 2          beta:age 0.02657317 0.005778867 0.01568617 0.03789621
#> 3 beta:inflammation 0.39514122 0.057874641 0.28254894 0.51070676
```
