# Method, assumptions and interpretation

## What is being generated

For each participant i, let Y_i be an underlying continuous measurement, W_i its
possibly missing or corrupted recorded value, C_i an observed interval category,
and F_i a quality flag. X_i contains completely observed predictors selected by the
analyst. The package generates alternative completed measurements conditional on
declared assumptions. It does not estimate an identified individual truth.

The example downstream estimand is the coefficient of Y in a linear regression
for an outcome. Its units are outcome units per measurement unit. A completed
dataset is an input to that analysis, not a guarantee of its validity.

## Required assumptions

1. Each category is an exact coarsening of Y at the supplied thresholds. Rounded,
   stale, independently miscoded, or corrupted categories do not meet this premise.
2. Unflagged observed values that agree with their category are sufficiently
   trustworthy to serve as original donors. Missing flags are treated as unflagged.
3. The predictor model and trusted donor selection transport to the missing or
   reconstructed records. Informative selection, unmeasured effects and poor
   overlap may violate this premise. No MAR or MNAR mechanism is estimated.
4. An explicitly chosen probability q defines a sensitivity replacement policy for
   flagged, observed, category-compatible values. It is not learned from the data.
5. If conflicts are reconstructed, the analyst grants the category authority over
   the continuous value. Default behavior instead rejects a conflict.

The unrestricted first and last intervals do not bound extreme errors. A very
large corrupted value can still be compatible with the highest category. The
retention policy can preserve that error, as the uncertain-flag benchmark shows.

## Algorithm

Validate the complete specification and input labels. Resolve missing categories
only from unflagged, observed, in-domain values. Detect conflicts on every row,
including unflagged observations. Preserve original inputs separately.

For each independent completion:

- Missing values and authorized conflicts require reconstruction.
- A flagged, observed, in-bound value is reconstructed with Bernoulli probability
  q; otherwise it remains exactly unchanged.
- Fit the continuous target on X using only original trusted observations.
  Use QR reduction for globally redundant columns and reject nonestimable
  recipient predictions. With n trusted rows and design rank r, draw residual
  variance as SSE / chi-squared(n-r), then regression coefficients from the
  corresponding normal distribution. Positive residual variance and degrees of
  freedom are required.
- For each reconstruction, restrict original donors to the recipient's category.
  Compare fitted donor predictions to the drawn recipient prediction (type-1 PMM).
  Break distance ties randomly, choose the nearest available donors up to
  `n_donors`, and sample one. A singleton pool always returns that sole donor.
- If support is below `min_donors`, fail by default. The optional normal fallback
  draws from the regression predictive normal restricted to the interval. It is
  a stronger assumption, records its use, and is not a validated interval-likelihood
  or measurement-error model.
- Verify all completed values against the interval convention and record decisions.

The regression excludes the category and flag as covariates in the high-level
API. Categories restrict donor eligibility. For the direct low-level `mice`
method, the analyst controls the predictor matrix and must supply the immutable
original response in `original_y`. The high-level API is the supported audited
workflow; a low-level invocation does not automatically create its provenance.

Complete predictors and one incomplete target make repeated chain iterations
unnecessary in this implementation. An independent parameter draw is shared across
recipients in each completion. The `mids` adapter stores the completions for `with`,
`complete` and pooling. It does not provide MCMC convergence diagnostics.

## q is a policy, not a corruption model

For an observed flagged row that already agrees with its category, the completion
rule is a point mass at W with weight 1-q, and a donor-based reconstruction with
weight q. No corruption density p(W|Y,F), validation sample, or posterior update
is fitted. The value of q can be varied to explore assumptions, but it should not
be called the estimated probability that a record is wrong.

If fraction f of the observed flagged rows conflicts with its category and those
conflicts are authorized for reconstruction, the expected reconstruction fraction
is f + (1-f)q for a scalar q. The realized fraction fluctuates. The actual changed
fraction can be smaller because a replacement can equal the observed value.
Missing target values are not part of this observed-flagged denominator.

Even when half of the flags truly indicate corruption, setting q=0.5 for every
compatible flagged record does not recover the appropriate posterior. Observed
magnitude, category and predictors may carry additional information. The benchmark
deliberately exposes this failure rather than tuning q to remove it.

## Inference boundary

Gaussian coefficient draws do not, by themselves, make the entire PMM and
retain/replace scheme proper multiple imputation. The optional fallback has no
general compatibility proof with the downstream model. Standard Rubin pooling is
available through mice as an analysis calculation; its variance and coverage need
separate validation for the intended setting.

The exploratory benchmark includes clean missingness, actual corruption, false
positive flags, heavy tails and category miscoding. It establishes neither broad
robustness nor clinical suitability. In particular, visibly constraint-consistent
completions can yield badly biased inference when the category premise is false.
An audit trail makes those assumptions inspectable; it cannot make them true.

## Provenance and reproducibility

The assumption hash uses canonical typed JSON, preserving list order, signed
infinities, missing values and numeric precision while sorting named maps. Null
fields are preserved; integer and equivalent double values intentionally share
the same numeric encoding. Equivalent YAML and R specifications are normalized
before hashing. This hash covers declared assumptions, not every execution input.

The manifest separately records the original-data fingerprint, selected predictors,
seed, imputation count, RNG kinds and R/package versions. The data fingerprint is
based on R serialization version 2 and can depend on data classes and environment;
it is not a portable semantic-data identifier. A canonical assumption hash is not
a full execution hash. The benchmark adds source-file hashes and environment records.
Seeded runs are reproducible in the same software and numerical environment; exact
cross-platform equality is not promised.
