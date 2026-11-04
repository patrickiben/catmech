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

Sites are independent, and row-error vectors are mutually independent. Site and
row vectors are independent of each other;
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
With no predictors, the model is bivariate and the reported intercept is the
population target mean `mu_T`.

The contextual validation example has `x = a_site + z`,
`T = 0.5 + 0.7*x + 0.4*a_site + v_T + e_T`, with independent unit-variance
`a_site` and `z`. Its within-site slope is 0.7, its latent between-site component slope is 1.1,
and its marginal slope is 0.9. Adding site means to a regression and reporting
the remaining x coefficient would answer the within-site question. The new
summary explicitly answers the marginal question. See
[Mundlak (1978)](https://people.stern.nyu.edu/wgreene/Econometrics/Mundlak-1978.pdf)
for the established correlated-effects problem; this implementation is not
claimed as a new statistical method.

## Prior and coefficient sensitivity

The population mean has a flat prior. Each row/site covariance Q is expressed
through a predictor covariance S, measurement regression B, and conditional
measurement covariance Omega. For measurements R = (T,W), use

```
S = Q_XX
B = solve(S, Q_XR)
Omega = Q_RR - transpose(B) %*% S %*% B
Q_RR = Omega + transpose(B) %*% S %*% B
Q_XR = S %*% B.
```

S and Omega have independent
[Huang--Wand hierarchical priors](https://www.matt-p-wand.net/publicns/Huang13.pdf).
For either covariance block C of dimension q,

```
C | a ~ inverse-Wishart(q + 1, 4 * diag(1 / a))
a_k ~ inverse-gamma(1/2, 1 / A_k^2)
A_k = sqrt(joint_control$prior_scale).
B | Omega ~ matrix-normal(0, regression_scale^2 * I, Omega)
regression_scale = 10 by default.
```

The inverse-gamma density is proportional to
`a^(-shape-1) * exp(-rate/a)`. **Within S and Omega**, marginal standard deviations
have half-t distributions with two degrees of freedom and scale A, and pairwise
correlations are marginally uniform. Those properties do not apply to the
reconstructed full Q. Standard deviations and correlations are not jointly
independent. The matrix-normal regression prior is proper and broad in the
empirically standardized units; its coefficient covariance is Omega times the
squared regression scale. Observed means and standard deviations used for
scaling are retained in the result metadata. These scaling constants are held
fixed in posterior simulation: they define a data-dependent prior on the
original-scale parameters, rather than another uncertain level of the model.

This factorization makes the regression prior explicit. The first development
candidate instead put a single Huang--Wand prior on each full covariance.
It improved contextual coverage but retained material downward coefficient bias.
That result, its failed calibration screens and its exact sources remain in the
[first development review](../validation/JOINT_MARGINAL_DEVELOPMENT_REVIEW.md).
The replacement prior was specified before new development cohorts were drawn;
no interval multiplier or nominal-level adjustment was selected.

A flat prior for B would be unsafe at a zero predictor site variance. In the
single-predictor case, write its site standard deviation as t and C = t B.
As t approaches zero, the observed-data likelihood can stay positive for bounded
C while the flat-B integration volume grows like `t^(-2)`. A proper finite-scale
regression prior avoids that divergent boundary. The default scale ten remains
an assumption that must be checked, not a claim of prior-free inference.

Sensitivity is assessed along two separate axes: covariance `prior_scale` at
0.25, one and four; regression scale at five, ten and twenty with covariance
scale held at one. For each axis, measure the average absolute change in the
posterior coefficient mean under its two alternatives. The candidate's primary
sensitivity statistic is the **larger** of these two values. Compare it with
the old conditional model's covariance-scale statistic on the same cohorts,
requiring all fits to return. Width sensitivity is separate. This comparison
covers the declared perturbations; it does not establish insensitivity to every
prior family or combination of hyperparameters. Finite-chain error can
contribute to small measured changes, so posterior mean Monte Carlo errors
are also retained.

## Sampler and missing patterns

An original compiled Gibbs sampler performs these updates:

1. Draw missing vector components conditional on the observed components and
   the current shared site effect. Observed components stay fixed.
2. Draw the population mean after integrating site effects, then draw each site
   effect conditional on that mean. This is a joint mean/site-effect block.
3. Update S from predictor sums of squares. Update Omega after integrating B,
   then B from its matrix-normal conditional. This is a joint regression and
   conditional-covariance block. Reconstruct the full row and site covariances.
4. Interweave site-scale moves described below, integrating both sets of site
   covariance auxiliaries and then redrawing them from their full conditionals.

For step 3, let the completed, centered vectors be the row residuals when
updating Sigma, or the site effects when updating Psi. Partition their matrix
as `[Y, X]`, with n vectors, two measurement columns in Y and k predictor columns
in X. Given the current auxiliaries, set

```
K = transpose(X) %*% X + I / regression_scale^2
M = solve(K, transpose(X) %*% Y)
S | rest ~ inverse-Wishart(n + k + 1,
  transpose(X) %*% X + 4 * diag(1 / a_X))
Omega | rest, B integrated out ~ inverse-Wishart(n + 3,
  transpose(Y - X %*% M) %*% (Y - X %*% M) +
  transpose(M) %*% M / regression_scale^2 + 4 * diag(1 / a_R))
B | Omega, rest ~ matrix-normal(M, solve(K), Omega).
```

For a block C of dimension q, each auxiliary is then drawn independently as
`inverse-gamma((q + 2) / 2, 2 * solve(C)[j,j] + 1 / A_j^2)`.
The residual-covariance auxiliaries are updated in step 3; the site-covariance
auxiliaries are updated after step 4. With k = 0, omit S and B and use the
two-measurement covariance update alone.

Final completed T values are drawn afresh from the parameter draws conditional
on the original observations. Each site gets one shared multivariate site draw
across all its missing patterns. Existing augmentation values do not become
additional evidence. Exact reference values are copied back without numerical
rounding changes. Other columns remain unchanged.

Near a zero site variance, a centered Gibbs covariance update can mix slowly.
The numerical amendment adds two scale proposals per vector coordinate,
including each predictor, per sweep. With
`eta ~ Normal(0,1)` and `c = exp(eta)`, scale one column of the site effects by c
and its covariance by `Psi' = D Psi D`, where D has c in that coordinate.
The site-effect quadratic form is unchanged. The full transformation Jacobian
is `c^(G+d+1)`; its `c^G` factor cancels the normal site-density determinant
change. The remaining log acceptance ratio is the completed-data log-likelihood
difference plus the marginal covariance log-prior difference plus
`(d+1)*eta`. The full-covariance prior density includes the factorization
Jacobian: `p(Q) = p(S) p(Omega) p(B | Omega) / |S|^2`. The hierarchical
auxiliaries for both covariance blocks are redrawn conditional on the resulting
covariance before any update uses them. This is a posterior-preserving
Metropolis move, not a modification of the prior or interval widths.

The first masked pilot and its source snapshot are preserved. Its four
convergence refusals prompted this numerical change before any candidate
coverage was examined. The same pilot cohorts were replayed for a paired
engineering check. After the subsequent prior revision, a further masked pilot
checks the revised procedure; its new development seeds and the unused
confirmation seeds are declared in the protocol.

## Diagnostics, scope and evidence

The default four chains start with 1,000 burn-in and 4,000 retained iterations,
extending the same seeded chains to 8,000, 16,000 and at most 64,000 retained
iterations. Rank-normalized R-hat must be at most 1.05 and bulk and tail ESS at
least 100 times the number of chains. Monitors include means, covariance entries,
hierarchical auxiliaries, conditional covariance blocks, within/between
regression matrices, original-unit marginal coefficients and conditional
means for both missing-target recording patterns. Passing the diagnostic screen
does not prove convergence or coverage.

The public engine requires complete continuous numeric predictors, sufficient
paired measurement support, at least five observed sites for each measurement,
and within-site replication. After adjustment for predictors and paired-site
indicators, the paired measurements must retain at least two residual degrees
of freedom and two numerically independent residual directions. The smallest
singular value of the standardized residual matrix must exceed
`sqrt(number_of_pairs) * 1e-8`. This conservative support rule excludes exact
measurement agreement and singular residual relationships under the weak
covariance prior. It can also exclude otherwise identifiable sparse paired
designs; it is a supported-data restriction, not a necessary identifiability
criterion. Constant, globally aliased, site-constant or within-site aliased
predictors are refused. The all-reference shortcut still
copies exact observations without fitting a posterior, so the marginal summary
refuses that shortcut explicitly.

The independent implementation checks include dense observed-Gaussian
conditioning, a known-covariance posterior oracle, prior marginal distributions,
and posterior comparisons with `pan` under a matched likelihood and prior.
The prospective repeated-sampling study is specified in the
[validation plan](../validation/JOINT_MARGINAL_PLAN.md). Its results, failures and
Monte Carlo uncertainty determine the empirical claims; these mathematical
checks alone establish no coverage guarantee.
The [completed confirmation](../validation/JOINT_MARGINAL_SYNTHESIS.md) passed
the contextual, ordinary calibration and slope-sensitivity screens. Its strict
availability screen failed because two candidate fits refused to return; the
option remains experimental and the diagnostic thresholds remain unchanged.

The repeated-sampling study has one continuous predictor. Multiple-predictor
support has algebraic and numerical checks, rather than a separate coverage
study. This model adds a Gaussian predictor-distribution assumption. Its validation
scope does not include categorical or missing predictors, informative site size,
nonlinear associations, imperfect reference measurements, MNAR, rare counts or
the repaired measurement used as an analysis predictor. The existing
`cluster_model = "joint"`, `catmech_joint_summary()` and CV3/Rubin analysis keep
their distinct interpretations.

Version 0.10.1 changes only the marginal default maximum from 32k to 64k;
all convergence checks remain fixed. The original confirmation used 32k and
remains historical evidence. The [numerical replay](../validation/MARGINAL_EXTENSION_SYNTHESIS.md)
uses the same cohorts and must not be counted as independent confirmation.
