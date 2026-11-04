# Joint marginal inference with contextual site effects

The opt-in `cluster_model = "joint_marginal"` models the true measurement,
recorded value and selected continuous predictors together. It addresses the
specific failure of a conditional random-intercept regression when its site
effects are associated with the predictors. It does not make arbitrary
misspecification harmless.

```r
fit <- catmech_reference(data, value = "recorded", reference = "gold",
  predictors = "x", use_recorded = TRUE, cluster = "site",
  cluster_model = "joint_marginal", m = 40, seed = 412)
catmech_marginal_summary(fit)
```

## Model and target

For row `i` at site `g`, model the vector `Z = (T, W, X)` as

```
Z_gi = mu + u_g + e_gi
u_g ~ Normal(0, Psi)
e_gi ~ Normal(0, Sigma).
```

Sites are independent. Site and row vectors are independent of each other;
each covariance matrix is unrestricted and positive definite. Correlations
between the predictor's and target's site components are estimated. Thus the
model does not require `u_T` to be independent of observed `X`. In particular,
rows with neither T nor W still inform the model through their observed X.

Let `V = Psi + Sigma`, with the recorded value excluded from the set of analysis
predictors. The population marginal regression is

```
beta = solve(V_XX, V_XT)
alpha = mu_T - transpose(beta) %*% mu_X.
```

This is the regression for an individual under the model. Site sizes must be
noninformative for the stated population interpretation. Each posterior draw
of means and covariances gives one draw of the coefficients, including the
population predictor distribution's uncertainty. The intercept and slopes are
transformed jointly to original units. Summaries use all retained draws across
all chains, independently of the number of completed datasets.

The contextual validation example has `x = a_site + z`,
`T = 0.5 + 0.7*x + 0.4*a_site + v_T + e_T`, with independent unit-variance
`a_site` and `z`. Its within-site slope is 0.7, its between-site slope is 1.1,
and its marginal slope is 0.9. Adding site means to a regression and reporting
the remaining x coefficient would answer the within-site question. The new
summary explicitly answers the marginal question. See
[Mundlak (1978)](https://people.stern.nyu.edu/wgreene/Econometrics/Mundlak-1978.pdf)
for the established correlated-effects problem; this implementation is not
claimed as a new statistical method.

## Prior and coefficient sensitivity

Means have a flat prior. Residual and site covariance matrices have independent
[Huang--Wand hierarchical priors](https://www.matt-p-wand.net/publicns/Huang13.pdf).
For a covariance matrix C of dimension d, use

```
C | a ~ inverse-Wishart(d + 1, 4 * diag(1 / a))
a_k ~ inverse-gamma(1/2, 1 / A_k^2)
A_k = sqrt(joint_control$prior_scale).
```

The inverse-gamma density is proportional to
`a^(-shape-1) * exp(-rate/a)`. This gives marginal half-t standard deviations
with two degrees of freedom and scale A, and marginal uniform pairwise
correlations. These parameters are not jointly independent. The default A is
one in empirically standardized response units; the prior is not noninformative
in every respect. Observed means and standard deviations used for scaling are
retained in the result metadata.

The scale sensitivity comparison multiplies squared half-t scale by 0.25 or
four, thereby halving or doubling A. The released conditional joint model uses
a fixed inverse-Wishart prior with different shape and moments. Equal named
multipliers across these families do not make the priors equivalent. The new
model and prior change together, so a comparison with the old procedure does
not attribute improvements separately to either component.

The coefficient sensitivity statistic is the average absolute change in the
posterior mean under the two alternative scales, in original units. Interval
width sensitivity is reported separately. Finite-chain error can contribute to
small measured differences; posterior mean Monte Carlo errors are retained.

## Sampler and missing patterns

An original compiled Gibbs sampler performs these updates:

1. Draw missing vector components conditional on the observed components and
   the current shared site effect. Observed components stay fixed.
2. Draw the population mean after integrating site effects, then draw each site
   effect conditional on that mean. This is a joint mean/site-effect block.
3. Draw each covariance from its inverse-Wishart full conditional and update
   the residual-covariance inverse-gamma auxiliaries.
4. Interweave site-scale moves described below, marginalizing the site-covariance
   auxiliaries and then redrawing them from their full conditional.

Final completed T values are drawn afresh from the parameter draws conditional
on the original observations. Each site gets one shared multivariate site draw
across all its missing patterns. Existing augmentation values do not become
additional evidence. Exact reference values are copied back without numerical
rounding changes. Other columns remain unchanged.

Near a zero site variance, a centered Gibbs covariance update can mix slowly.
The numerical amendment adds two scale proposals per response per sweep. With
`eta ~ Normal(0,1)` and `c = exp(eta)`, scale one column of the site effects by c
and its covariance by `Psi' = D Psi D`, where D has c in that coordinate.
The site-effect quadratic form is unchanged. The full transformation Jacobian
is `c^(G+d+1)`; its `c^G` factor cancels the normal site-density determinant
change. The remaining log acceptance ratio is the completed-data log-likelihood
difference plus the marginal covariance log-prior difference plus
`(d+1)*eta`. The hierarchical auxiliaries are redrawn conditional on the resulting
covariance before any update uses them. This is a posterior-preserving
Metropolis move, not a modification of the prior or interval widths.

The first masked pilot and its source snapshot are preserved. Its four
convergence refusals prompted this numerical change before any candidate
coverage was examined. The same pilot cohorts are replayed for a paired
engineering check; the development and confirmation seeds remain fresh.

## Diagnostics, scope and evidence

The default four chains start with 1,000 burn-in and 4,000 retained iterations,
extending the same seeded chains to 8,000, 16,000 and at most 32,000 retained
iterations. Rank-normalized R-hat must be at most 1.05 and bulk and tail ESS at
least 100 times the number of chains. Monitors include means, covariance entries,
hierarchical auxiliaries, original-unit marginal coefficients and conditional
means for both missing-target recording patterns. Passing the diagnostic screen
does not prove convergence or coverage.

The public engine requires complete continuous numeric predictors, sufficient
paired measurement support, at least five observed sites for each measurement,
and within-site replication. Constant, globally aliased, site-constant or
within-site aliased predictors are refused. The all-reference shortcut still
copies exact observations without fitting a posterior, so the marginal summary
refuses that shortcut explicitly.

The independent implementation checks include dense observed-Gaussian
conditioning, a known-covariance posterior oracle, prior marginal distributions,
and posterior comparisons with `pan` under a matched likelihood and prior.
The prospective repeated-sampling study is specified in the
[validation plan](../validation/JOINT_MARGINAL_PLAN.md). Its results, failures and
Monte Carlo uncertainty determine the empirical claims; these mathematical
checks alone establish no coverage guarantee.

This model adds a Gaussian predictor-distribution assumption. Its validation
scope does not include categorical or missing predictors, informative site size,
nonlinear associations, imperfect reference measurements, MNAR, rare counts or
the repaired measurement used as an analysis predictor. The existing
`cluster_model = "joint"`, `catmech_joint_summary()` and CV3/Rubin analysis keep
their distinct interpretations.
