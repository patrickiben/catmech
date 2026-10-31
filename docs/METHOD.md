# Model, assumptions and interpretation

Version 0.2 introduced `catmech_model()`. The unchanged PMM retain/replace engine
is documented separately in [LEGACY_METHOD.md](LEGACY_METHOD.md). Version 0.3
adds ordinary flag-masking PMM and an explicit reference-subset workflow, described
in [REFERENCE_METHOD.md](REFERENCE_METHOD.md). They do not alter the Gaussian
category likelihood described here.

## Target and trust policy

Let Y be the latent continuous measurement, W its recorded value, F a quality
flag, C an optional category, and X complete selected predictors including an
intercept. The downstream estimand must be stated separately; the examples use
the coefficient of Y in a linear outcome regression.

All missing W and all flagged W are treated as unavailable. Observed unflagged W
are treated as exact Y and preserved. Missing flags are treated as unflagged.
False positive flags cost information but do not result in arbitrary retention
of suspect values. Missed flags leave residual contamination and can invalidate
inference. This is a masking policy, not a fitted probability of corruption.

The model assumes `Y | X, beta, sigma² ~ Normal(X beta, sigma²)` independently
across records, with reference prior `p(beta, sigma²) ∝ 1/sigma²`. Selection into
the trusted set must be ignorable conditional on X. The trusted design must have
full column rank, more observations than coefficients and positive residual
variance. Constants, aliases, incomplete predictors and unsupported levels are
rejected. The trusted posterior provides a proper anchor before categories are used.
For example, flags created by thresholding the measurement itself can remove
genuine extremes and violate ignorability. Masking is not automatically conservative;
the primary random-flag simulations do not validate that selection mechanism.

## Optional category likelihood

Let `b(Y)` denote the interval containing Y and let `H[k,c] = Pr(C=c | b(Y)=k)`.
Rows of H are true bins and columns observed bins. The supplied likelihood assumes
category errors depend on Y only through its bin, given the modeled information.
The endpoints must span `(-Inf, Inf)`; finite outer bounds would require additional
truncation normalization incompatible with this conjugate parameter update.
Intervals are right closed: `(a,b]`, with the first interval including its lower
endpoint. Finite clinical reference boundaries may follow a different convention.

- `ignore`: categories have no likelihood and are not used to repair numeric data.
- `exact`: H is the identity matrix. Trusted-value contradictions are errors.
- `noisy`: a named nonnegative row-stochastic H is supplied externally. A trusted
  disagreement is allowed only if its likelihood is positive.

Missing categories contribute no likelihood. Observed categories with zero
probability in every bin are rejected. Trusted numeric values and recorded category
labels are never changed to make them agree. The category likelihood for a trusted
Y is a constant with respect to beta and sigma, so it does not change their posterior.

H must come from external validation or explicitly declared sensitivity assumptions.
It is not fitted here; uncertainty in H is not integrated out. Categories measured
at another time, in different units, or derived from the same erroneous W need a
different observation model. A noisy H does not automatically solve those problems.

## Posterior computation

Ignoring categories gives ordinary Bayesian normal imputation. If n trusted rows
and p coefficients give least-squares estimate beta-hat and residual sum of squares
SSE, draw `sigma² = SSE / chi-squared(n-p)`, then
`beta ~ Normal(beta-hat, sigma² (X'X)^(-1))` using a QR factorization.
Draw each unavailable Y from its normal predictive conditional on this shared
parameter draw. Each completion has an independent parameter draw. These are
established calculations, also implemented in
[mice.impute.norm](https://amices.org/mice/reference/mice.impute.norm.html).

With informative categories on masked records, Gibbs data augmentation alternates:

1. Draw regression parameters using the full currently completed data and the same
   reference-prior conjugate update. The degrees of freedom are n_total - p, not the number of trusted records
   minus p.
2. For each masked row with observed category c, sample bin k with weight
   `H[k,c] * {Phi((upper[k]-X beta)/sigma) - Phi((lower[k]-X beta)/sigma)}`.
   Then sample Y from the corresponding truncated normal. Rows with missing
   categories have unrestricted normal draws.

The marginal category likelihood is thus included in the posterior. Merely fitting
trusted data and clipping draws to bins would generally omit category information
from the parameter posterior. Exactly uninformative likelihood columns reduce to
the independent engine. Log interval probabilities are evaluated in the nearer
normal tail to reduce cancellation. Unrepresentable draws cause explicit failure.
This is conventional [data augmentation](https://doi.org/10.1080/01621459.1987.10478458),
not a new sampling algorithm.

By default, four chains start with separate trusted-posterior draws. Defaults discard 200 steps
per chain and retain 400. Parameter traces and the mean latent value are screened
using rank-normalized split/folded R-hat and bulk/tail ESS from `posterior`. Every
screened variable must have R-hat <= 1.05 and both ESS values >= 100 times the number
of chains. Nonfinite diagnostics also fail. A `catmech_convergence_error` retains
traces and summaries; it does not return a successful analysis. These are finite
screens, not proof of convergence, exploration of every mode, or coverage.
[Diagnostic reference](https://mc-stan.org/rstan/reference/Rhat.html).

Requested completions are spaced through retained iterations and balanced across
chains. They are not claimed independent. More imputations can reduce imputation
Monte Carlo error, but cannot repair an incorrect model. Increase sampling effort
when needed based on diagnostics, and retain failed runs when evaluating a fixed
procedure rather than retrying them invisibly.

## Pooling and assumption comparisons

`catmech_pool_lm()` fits the requested model to each completion and uses
`mice::pool.scalar()` with complete-data sample size and coefficient count for
Rubin pooling and Barnard-Rubin degrees of freedom. Missing analysis data,
singular fits, nonfinite or degenerate uncertainty and failed sampler diagnostics
are rejected. Legacy sensitivity results require an explicit opt-in and retain
their sensitivity label. Standard PMM and reference workflows have separate
interpretation labels. Partial-reference PMM requires explicit sensitivity opt-in
after its adverse held-out coverage result; fully referenced PMM is an ordinary
complete-data endpoint. When explicit workflow metadata establishes that no target
was imputed, pooling returns the ordinary complete-data regression interval after
checking that analysis data agree across completions. A coincidentally zero
between-imputation variance alone does not trigger that rule. The `mids` object is an adapter for externally generated
completions, not a fitted chained-equations process; never resume it.

A coherent imputation likelihood does not establish compatibility with every
analysis. Linear transformations, interactions, nonlinearity, selection, residual
distribution and dependent observations need appropriate modeling and validation.
The tested linear-Gaussian settings do not cover all estimands or clinical designs.

`catmech_sensitivity()` takes a named list of results or caught error conditions
from the same original data. `expected_scenarios` catches omitted prespecified
cases and rejects extras. Errors and missing scenarios remain visible and make the
overall conclusion unresolved. Successful intervals are classified relative to a
declared reference; the smallest lower and largest upper bound form a descriptive
envelope. Agreement over a finite set is not robustness over all plausible models,
and the envelope has no claimed confidence level. A caller must supply its complete
prespecified scenario list for this omission check to be meaningful.

## What this change repairs, and what it cannot

The v0.1 uncertain-flag failure retained some severely corrupted values under an
arbitrary q. Its intervals were already wide, so the main remedy is removal of
that bias source, not interval inflation. Wrong categories caused a second failure
when treated as authoritative. Ignoring them avoids that premise; an appropriate
H uses them probabilistically. Incorrect H remains a deliberate stress condition.

Unflagged contamination, informative trusted-data selection, heteroskedastic or
nonlinear regressions, category-error dependence, uncertainty in H
and repeated/clustered outcomes remain outside the fitted model. Clinical BLQ/PK
censoring is outside the validated scope, even where a single censoring interval
is mathematically representable by the category likelihood. The laboratory site-shift fixture also includes random errors and reduced flag
detection, so its failure cannot be attributed to the shift alone. It diagnoses
residual error; it does not validate random-effects or cluster-robust inference. No individual latent truth is identified by the model.

## Provenance

Original data, target-cell changes, assumptions, seeds, RNG kinds, versions and
sampler traces are retained. Assumption hashes use canonical typed JSON, including
signed infinities and nulls. The model specification stores H as named rows/columns.
Data hashes use R serialization version 2 and are environment-dependent fingerprints,
not portable semantic identifiers. Benchmarks additionally record source hashes
before and after execution. Equal hashes establish correspondence, not correctness.
Exact cross-platform numerical equality is not promised.
