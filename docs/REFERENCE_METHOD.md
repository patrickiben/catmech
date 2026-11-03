# Completion from independently obtained reference measurements

`catmech_reference()` addresses a different data situation from flag-based
completion: a subset of records has a genuine, sufficiently accurate reference
measurement of the target. The remaining true values are unavailable even when
the original recorded values are present and unflagged. This can address errors
missed by a flagging rule, provided the reference sample and imputation model
support the intended analysis. A flag, a plausibility check, and agreement with
a category are not themselves reference measurements.

This is an implementation of established missing-data and measurement-error
methods. Cole, Chu and Greenland described multiple imputation for measurement
error using a validation substudy in 2006. The default below is the simpler
special case that does not use the error-prone measurement as an auxiliary
predictor. It does not estimate an error process or each record's probability
of corruption. [Cole et al., 2006](https://doi.org/10.1093/ije/dyl097)

## Data contract

```r
fit <- catmech_reference(
  data,
  value = "recorded_value",
  reference = "reference_value",
  predictors = c("outcome", "age", "site"),
  method = "normal",
  m = 20,
  seed = 123
)
```

The columns have distinct roles:

| Input | Meaning and treatment |
| --- | --- |
| `value` | Original recorded target. The returned completions replace this column with reference or imputed values. The original remains in `fit$original`. |
| `reference` | Finite numeric measurements only on genuinely referenced rows, with `NA` elsewhere. These values are copied exactly into the completed target. This column is preserved. |
| `predictors` | Explicit, completely observed predictors of the true target. Include the analysis outcome and relevant analysis/selection covariates when appropriate. The target and reference columns cannot be predictors. |
| Other columns | Preserved; neither flags nor category labels automatically affect imputation. |

Every row without a reference value is imputed. An unflagged recorded target is
not promoted to a reference. Conversely, an observed reference is used even
when it differs from the recorded target. `fit$mids$data` represents the
imputation problem with reference values observed and all other target values
missing. The original data, completed data, and imputation inputs are therefore
different, deliberately documented objects.

Reference values must come from an independent measurement, adjudication, or
verification process appropriate to the estimand. In a simulation, the analysis
input must expose frozen truth only on the prespecified reference rows; the
remaining truth belongs in the evaluator, not in the imputer's data frame.
The function cannot establish that a supplied reference is scientifically valid.

## Statistical model and assumptions

Write the true target as `T`, selected complete predictors as `Z`, and the
indicator of an available reference as `R`. The normal option fits
`T | Z ~ Normal(Z beta, sigma^2)` to the reference rows, using the reference prior
`p(beta, sigma^2) proportional to 1 / sigma^2`. Each completion draws regression
parameters and then draws the unavailable targets. It propagates fitted-model
and residual uncertainty; it does not substitute predicted means for every
unavailable target. This is the standard Bayesian normal-imputation construction
also documented by `mice`. [mice normal imputation](https://amices.org/mice/reference/mice.impute.norm.html)

For partially referenced data, the reference design must have full column rank,
positive residual degrees of freedom, and positive residual variance. These
checks establish a usable fitted model, not representative sampling or
adequate sample size. If all rows are referenced, completion is deterministic
and no imputation model is needed.

The key identifying assumption is `T independent of R given Z`: after
conditioning on the selected predictors, reference availability does not
depend on the unavailable true target. There must also be adequate reference
support across the predictor space in which values are imputed. A random
internal reference sample is a transparent design. A sample selected using
observed outcome, site, recorded measurement, or flags requires those selection
variables and suitable functional forms in the model, as well as overlap. The
general internal-validation MI argument and its model requirements are described
in the STRATOS guidance. [Shaw et al., 2020, section 2.4](https://pmc.ncbi.nlm.nih.gov/articles/PMC7272296/)

Examples that violate or strain these assumptions include verification limited
to one site, references collected only for extreme recorded values, and
verification triggered by an unrecorded judgment correlated with the true
target. Including a site variable does not create reference information for
unsupported sites. Selection depending on the unobserved truth after
conditioning on the predictors remains nonignorable.

A noisy reference is not a gold standard. Treating a second error-prone
measurement as exact changes the problem and can retain bias. Repeated-measure
models are an established alternative when a suitable error model is credible;
for example, `smcfcs` supports latent normal targets with classical-error
replicates. This function does not implement that error model.
[smcfcs measurement-error vignette](https://cran.r-project.org/web/packages/smcfcs/vignettes/coverror.html)

## Optional use of the recorded measurement

The default does not use the recorded target to predict its true value. This
choice makes the resulting completions invariant to changes in that recorded
column, conditional on fixed reference values, selected predictors, and seed.
It sacrifices potentially useful information and cannot be assumed to improve
precision over analysis restricted to the reference sample.

To use the recorded measurement as a linear predictor, set `use_recorded = TRUE`:

```r
fit_with_recorded <- catmech_reference(
  data, "recorded_value", "reference_value",
  predictors = c("outcome", "age", "site"),
  method = "normal", use_recorded = TRUE, m = 20, seed = 123
)
```

Records with a recorded value are imputed from a model fitted on references that also have one.
Records without one are imputed from a model on the predictors alone, fitted on all references. If every
calibration row agrees exactly, there is no residual variance to model and those records keep their recorded
value, with no added uncertainty, a warning and an interpretation ending in `_recorded_kept`. Intervals then treat those
records as correct, so a small audit that happened to miss large errors gives intervals that do not cover the truth. This was the strong comparator in the version 0.4 studies, and it was best or
tied among reference-using methods on ClinicalTrials.gov and NHANES at 25% audits. See [the real-data evidence](../validation/REAL_DATA_SYNTHESIS.md).

This fits the predictive relationship between reference truth, the auxiliary
measurement, and other predictors. That relationship may be nonlinear or
heteroscedastic when errors are a mixture of correct values and occasional
large defects. Merely adding the recorded measurement does not make a linear
normal model correct. A recorded value can also be informative about reference
selection, so omitting it is not automatically valid in a targeted audit.

If the recorded measurement is missing, it cannot be passed as an incomplete
predictor. The default can still run without it when the chosen predictors are
complete and support the required missingness assumption. It does not
automatically impute the auxiliary measurement, fill missing measurements with
zero, or borrow a model across unsupported missingness patterns. Those choices
need a prespecified joint or pattern-specific model and suitable reference
information. A validation sample containing no references in a relevant
missing-measurement subgroup cannot by itself verify extrapolation to that
subgroup.

### Clustered audits: experimental

With `use_recorded = TRUE`, `cluster = "site"` adds a random intercept for the cluster to the recorded-value
regression: `target = x'beta + c_i + e_ij`. If the audit is stratified, put the strata in `predictors`. In each imputation:

1. Variance components come from moments on the audited rows (Henderson, 1953, method 3). The fit with cluster
   indicators gives the within-cluster sum of squares, and the rest of the pooled fit's gives the between-cluster part.
   Each is divided by a chi-square draw, and the intercept variance is the moment solution, floored at zero.
2. Coefficients are drawn from the generalized least squares fit given those components.
3. An unaudited cluster gets a new intercept, an audited cluster with unaudited rows gets one conditional on its
   residuals, and every imputed row gets normal noise.

Records with no recorded value use the same steps without that predictor. If the audited rows cannot separate the two
variances, as with one row per cluster, the intercept variance is zero, with a warning. If the predictors absorb every
audited cluster, as with one audited cluster per stratum, the fit is refused when unaudited clusters need imputing. Fewer
than five audited clusters are refused.
`catmech_pool_lm()` applies a delete-one-cluster jackknife variance (`cv3`) to each completed dataset and pools by
Rubin's rules, with complete-data degrees of freedom of clusters minus one. A fit made without `cluster` can still use
that variance through `catmech_pool_lm(cluster = )`. A [planned test](../validation/VARIANCE_MI_SYNTHESIS.md) chose it
narrowly (MacKinnon, Nielsen and Webb, 2023, [doi:10.1016/j.jeconom.2022.04.001](https://doi.org/10.1016/j.jeconom.2022.04.001), describe
the variance; Barnard and Rubin, 1999, [doi:10.1093/biomet/86.4.948](https://doi.org/10.1093/biomet/86.4.948), the degrees of freedom).
`cluster_variance = "cr1s"` gives the earlier default. A coefficient that a single cluster determines draws a warning: a
pseudoinverse keeps the variance computable, but it leaves out that cluster's contribution, and in a simulated check such
intervals covered about a third of the time.

**Use `cluster` in the fit only when the audited value is the analysis outcome.** In simulations with the audited value
as a predictor, the option biased that slope toward zero, and most when errors shifted whole clusters. A new intercept
drawn independently of a cluster's outcomes appears to act like added noise in a predictor, though that explanation has
not been tested. `catmech_pool_lm()` warns about this use. `mice` with `2l.lmer` and `2l.pan` attenuated the slope in the
same cells, so the model accounts for a large part of it, though this implementation attenuated more than `2l.lmer` in
all of them, for reasons not established.

The same model is available in `mice`, and a [planned comparison](../validation/OUTSIDE_METHODS_SYNTHESIS.md) found the same
accuracy once the `mice` check that can drop the recorded value was off. This method draws its variance components in
closed form, without MCMC or `lme4`, and its intervals were slightly wider. An [independently written setup](../validation/INDEPENDENT_SETUPS_SYNTHESIS.md) with separate
within-cluster and between-cluster slopes looks better on most measures. [Adding those slopes alone](../validation/WITHIN_BETWEEN_SYNTHESIS.md) was tested and not adopted: it
raised error on small real audits when the audited value was the outcome, so what else separates that setup remains open.

This is the usual random-intercept imputation model (Schafer and Yucel, 2002,
[doi:10.1198/106186002760180608](https://doi.org/10.1198/106186002760180608); Grund, Lüdtke and Robitzsch, 2018,
[doi:10.1177/1094428117703686](https://doi.org/10.1177/1094428117703686)), and not a new method, with a moment-based draw of the
variance components (Henderson, 1953, [doi:10.2307/3001853](https://doi.org/10.2307/3001853)) in place of a full posterior. That draw is
exact in distribution only for balanced clusters, and it is weakest with one or two large sites. The model assumes one normal
intercept per cluster, so cluster-specific slopes and non-constant variance are not modeled.

Its [planned test](../validation/CLUSTERED_FIX_SYNTHESIS.md) came out PARTIAL: coverage held on eight of nine items of a
real recount at a 25% audit, with lower error on six common items, but error on three rare items ran 4% to 5% above the option
without `cluster`, and the accuracy conditions failed. It replaced an [earlier version](../validation/CLUSTERED_SYNTHESIS.md)
that resampled audited clusters in each imputation, within `cluster_strata` if given, which also came out PARTIAL. A bootstrap
mean within a stratum of n audited clusters has (n - 1)/n of the usual variance, and a stratum with one contributes none. A
rescaled bootstrap (Rao and Wu, 1988) overshot in the second test, and a stratified jackknife draw reached higher coverage with
much wider intervals.

## Predictive mean matching and analysis

`method = "pmm"` uses the public `mice` predictive mean matching implementation,
with donor values restricted to the reference rows and `donors = 5L` by
default. This preserves observed reference support and can be useful when a
normal residual model is doubtful. It does not guarantee calibrated inference
under heavy tails, nonlinear relationships, sparse reference samples, or poor
overlap. It cannot generate target values outside the reference donor support.
Changing the donor pool changes the procedure and should be prespecified.

`method = "pmm_abb"` resamples donors with replacement for each imputation before the same
matching. It addresses attenuated between-imputation variance and raised synthetic coverage from 86.7% to 95.9%
in version 0.4. On real data it overcovered with wider intervals, so it is not a default.

In the held-out random-25%-reference experiment of version 0.3, reference PMM covered only
85.1% of nominal 95% intervals and exactly matched the independent mice PMM
comparator. Version 0.3.1 therefore requires `allow_sensitivity = TRUE` when
pooling a partial-reference PMM result, and labels it `sensitivity_only`. This
is a conservative package policy based on adverse evidence in the tested regime,
not a theorem against every PMM application. The normal reference method remains
the default. Fully referenced data need no imputation and retain ordinary pooling.
See the [full extension evidence](../validation/EXTENSION_SYNTHESIS.md).
[mice PMM documentation](https://amices.org/mice/reference/mice.impute.pmm.html),
[Morris et al., 2014](https://doi.org/10.1186/1471-2288-14-75)

Use an analysis model compatible with the imputation model, and include the
outcome in covariate imputation. Linear normal imputation is not a general
solution for survival models, nonlinear effects, or interactions. Substantive
model compatible MI is an established alternative when the analysis requires
those structures. Its reference implementation includes an explicit internal
validation example in which the error-prone measurement is an auxiliary
predictor and the outcome enters through the substantive model.
[smcfcs vignette](https://cran.r-project.org/web/packages/smcfcs/vignettes/coverror.html)

Audit-based correction is not necessarily more efficient, or even necessary
for every estimand. For example, the effect of a randomized treatment can remain
unbiased in certain settings with covariate errors while the coefficient of
the mismeasured covariate is biased. The target parameter, error process,
reference sample size, and comparison analysis all matter.
[Shepherd, Shaw and Dodd, 2012](https://pmc.ncbi.nlm.nih.gov/articles/PMC3728661/)

## What the software and validation can establish

The useful contribution is an explicit reference-data contract, preservation
of original measurements, auditable completion decisions, and reproducible
comparison of assumptions. Normal MI, PMM, internal-validation error correction,
and Rubin pooling are established methods; this interface is not a new
identification result or a guarantee of improved coverage.

Implementation checks should verify reference preservation, independence from
unused recorded values, absence of hidden truth outside the reference sample,
and correct posterior predictive moments. For the normal model, with reference
design `A`, unavailable-row design `A0`, residual sum of squares `SSE`, and
`nu = n_reference - n_parameters`, the unavailable targets have a multivariate
Student distribution with location `A0 beta_hat` and scale matrix
`SSE / nu * (I + A0 %*% solve(t(A) %*% A) %*% t(A0))`. When `nu > 2`, its covariance is
that matrix multiplied by `nu / (nu - 2)`. Checking cross-record covariance
tests whether shared parameter uncertainty is propagated, beyond testing
individual marginal means.

The [independent reference oracle](../validation/reference-oracle/README.md)
implements this calculation without the package's fitting or sampling helpers.
Its 12,000-completion run passed all 14 mean/covariance checks and 12 data-contract
controls, with fixed five-MCSE moment tolerances and unchanged runtime source
hashes. Deliberately wrong zero cross-record covariance and plug-in variance
answers were detected. This is implementation evidence, not a coverage study.

Simulation comparisons should include reference-only analysis and standard
`mice` methods using the same reference rows and predictors. Reference fraction,
selection, reference noise, unsupported sites, and model misspecification
should be varied separately. A synthetic improvement over a method that was
given less information does not establish a superior imputation algorithm.
Finite simulation coverage describes the simulated setting; it is not evidence
from patients or a general clinical performance claim.

## Research boundary: uncertain category confusion probabilities

The existing noisy-category model takes a fixed matrix
`H[true_category, observed_category]`. External gold-standard classification
counts could support row-wise Dirichlet distributions with explicitly specified
positive prior concentrations. This would quantify uncertainty about the
category observation model, subject to transportability, rather than calibrate
the legacy retain/replace probability or identify unflagged numeric corruption.
Multinomial misclassification models with gold-standard information and
Dirichlet priors are established prior art.
[Wang, Wang and Swartz, 2020](https://doi.org/10.1002/cjs.11556)

Drawing matrices from their external-data posterior and running each as a
fixed scenario cuts feedback from the main study into the matrix. It therefore
defines a modular uncertainty analysis, not the full joint posterior. A full
joint model would update the matrix using the current trusted/latent true bins
and observed labels, and would require additional diagnostics and an independent
oracle. Cutting feedback is an explicit modeling choice, not a synonymous name
for ordinary joint Bayesian inference.
[Liu and Goudie, 2025](https://academic.oup.com/jrsssb/article/87/4/1171/8096422)

This extension is not implemented here. It should follow, rather than obscure,
validation of the reference workflow. If explored, retain distinct scenarios
and external-count provenance; do not pool repeated batches under different
matrices as if they were an ordinary single set of MI draws. Nested
multiple-model imputation has its own combining rules and assumptions.
[Siddique, Harel and Crespi, 2012](https://pmc.ncbi.nlm.nih.gov/articles/PMC3596844/)
