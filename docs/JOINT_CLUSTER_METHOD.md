# Joint clustered reference imputation

Version 0.8.0 adds the explicit `cluster_model = "joint"` option to
`catmech_reference(method = "normal", use_recorded = TRUE, cluster = "site")`.
It supports both recorded-present and recorded-missing unaudited rows within the
same site. The default `"moment"` option preserves the earlier behavior.

## Model and contribution

Let `T` be the true measurement, observed only at checked reference rows, and
`W` the recorded measurement, possibly missing. For complete selected covariates
`x`, the joint model is

```
(T_ij, W_ij)' = B' x_ij + u_j + e_ij
u_j ~ N_2(0, Psi),   e_ij ~ N_2(0, Sigma).
```

The two covariance matrices are unrestricted positive definite matrices.
Different sites are independent conditional on the population parameters;
rows are independent conditional on their site effect. Both measurements share
one site vector in each completion. In particular, the unobserved true values
in the two recording patterns can remain dependent after fitting the covariates.

This model class and its Gibbs sampler are established methods. `catmech`
delegates parameter sampling to [`pan`](https://search.r-project.org/CRAN/refmans/pan/html/pan.html),
which implements multivariate mixed-model imputation following Schafer and
Yucel (2002). The package contribution is the reference-data interface,
preservation of audited values, mixed-pattern conditional prediction,
diagnostics and audit trail. It is not a new estimator or a claim that an
existing undercoverage problem has generally been solved.

## What a completion conditions on

For each posterior parameter draw, the site effect is sampled conditional on
**all original observed components** in that site. If `O_i` indexes the observed
components at row `i`, its contribution to site precision is the embedding of
`Sigma[O_i,O_i]^-1` in the corresponding two coordinates. Summing these
contributions and adding `Psi^-1` gives the posterior precision of `u_j`.
The same precision-weighted residuals determine its posterior mean.

An unaudited row with observed `W` then receives the normal conditional draw
of `T` given `W`, `x`, and that shared `u_j`. If both components are absent,
the draw uses the marginal `T` distribution conditional on `x` and `u_j`.
Reference `T` values are copied exactly. An entirely unobserved site uses a
new population site draw. Rows with neither component observed do not add
information to the parameter likelihood.

The returned completed target contains true-value imputations. The original
recorded values, including missingness, remain unchanged in `original`;
the reference column is unchanged everywhere. Even complete agreement on a
small calibration subset does not turn this joint model into deterministic
acceptance of unaudited records.

## Priors, computation and support

Each response is centered and scaled using its observed mean and sample
standard deviation. The defaults are independent inverse-Wishart priors
`Sigma ~ IW(4, 0.5 I_2)` and `Psi ~ IW(4, 0.5 I_2)` on that standardized scale,
with the flat (improper) fixed-effect prior used by `pan`. Here `IW(df, S)` has mean
`S / (df - 3)` in two dimensions. Thus each default covariance prior has mean
`0.5 I_2`; these are substantive, empirically scaled priors, not
noninformative priors. `joint_control = list(prior_scale = k)` multiplies both
scale matrices by positive `k`. Inspect sensitivity to this choice, especially
with few informative sites or weak covariance identification.

Defaults use four chains, 1,000 burn-in iterations and 2,000 retained iterations
per chain, extending the retained run up to 8,000 when needed. The controls
`chains`, `burnin`, `iterations`, and `max_iterations` are available through
`joint_control`. Each chain runs in a fresh local R worker: `pan` 2.0 retains a
cached normal draw that its public seed argument does not reset, so in-process
calls can depend on earlier `pan` use. Isolation gives each declared chain a
fresh sampler state and is also applied on automatic extensions. This adds
process-startup cost. Rank-normalized R-hat must be at most 1.05, with bulk and tail
effective sample sizes each at least 100 times the chain count. A failed
diagnostic screen refuses completion and carries the diagnostics in a
`catmech_convergence_error`; earlier numerical failures can raise ordinary
errors before valid traces exist. Passing these checks is evidence about numerical
sampling, not proof of model fit or interval calibration.

The input must provide identifiable covariate designs, paired calibration
observations, multiple informative sites and replication within sites.
Proper covariance priors alone do not make an unsupported covariance split
identified by the data. Support checks are conservative engineering rules,
not thresholds that guarantee reliable estimation. The fitted specification
records the actual controls, transformations and priors.
Even with five gold sites, paired measurements can concentrate in very few
sites, leaving the cross-response site covariance strongly prior-dependent.
The initial simulations have broad paired support and do not validate that
concentrated-reference setting.

## Scientific limits and validation

References are assumed exact and comparable in quantity, units and time.
Missingness of either component must be ignorable given the observed data
under this joint model, with appropriate support. Selection depending on
unobserved true values after conditioning is not repaired. Random slopes,
nonlinear relationships, heteroskedasticity, non-Gaussian count distributions
and occasional gross recording errors are outside this model.
Including a site effect does not by itself repair audit selection driven by
unobserved site quality.

Use `catmech_pool_lm()` for a compatible prespecified linear analysis; its
default site variance remains the delete-one-cluster jackknife followed by
Rubin pooling. Few influential sites can still defeat that analysis.
The new tests concern the audited value as the **outcome**. Using it as a
predictor needs separate validation of the imputation/analysis relationship.
The earlier moment model's predictor-slope attenuation is neither attributed
to this implementation nor declared cured by it.

`catmech_holdout_check()` retains the joint model and its controls, but splits
reference rows rather than entire sites. Other references in the same site can
still inform a held-out row, so this does not test transport to wholly
unaudited sites. Its Wilson and uniform-PIT flags use row-level calculations
without a cluster correction. With dependent rows, treat those flags as
descriptive diagnostics, not calibrated tests of site-level predictive validity.

The [frozen validation plan](../validation/JOINT_CLUSTER_PLAN.md) and
[results](../validation/JOINT_CLUSTER_SYNTHESIS.md) distinguish mathematical
implementation checks, matched-model simulations, prior sensitivity and
misspecification stress. Historical Rwanda and rare-count studies concern
earlier models; their findings are preserved and do not validate this option.
There has been no independent human statistical review.

In the initial prespecified study, the default covariance prior with 40
completions returned 1,185 of 1,200 requested fits. Coverage among returned
intervals ranged from 96.5% to 99.5% across the six scenarios; the balanced
Gaussian case was conservative. Counting failed fits in the denominator gave
93.0% in the t5 stress case, where eight fits refused convergence. Covariance
prior changes materially altered widths and coverage. These results validate
neither exact 95% calibration nor a general solution to undercoverage. See the
synthesis for Monte Carlo intervals, paired comparisons and all failures.
