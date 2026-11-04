# Prospective study: contextual dependence and coefficient prior sensitivity

Written 2026-11-04 before candidate accuracy runs. The completed v0.9.0
study is exposed development evidence, including its contextual failure. This
protocol uses new cohorts for development and fresh cohorts for confirmation.
The candidate, controls and any amendments must be committed before scoring.
Amendment 2 specifies the current candidate; earlier sections preserve the
prospective design and history of the first scientific candidate.

## Scientific question and fixed target

Can an explicit joint model of the true measurement, recorded measurement and
continuous predictors recover the population marginal coefficient when site
effects and predictors are associated? Does a hierarchical covariance prior
reduce the remaining sensitivity of that coefficient to declared scale changes?
These are separate claims. A successful implementation or a low failure rate
does not establish either claim.

The candidate is requested explicitly with `cluster_model = "joint_marginal"`.
It models each row's vector `Y = (T, W, X)` as `mu + u_site + e_row`, with
independent multivariate Gaussian site and row effects having covariance `Psi`
and `Sigma`. All selected predictors must be complete, continuous numeric
variables. Predictors are modeled random variables, rather than fixed regressors.
With `V = Psi + Sigma`, the marginal regression target is

```
beta = solve(V_XX, V_XT)
alpha = mu_T - transpose(beta) %*% mu_X.
```

`catmech_marginal_summary()` summarizes this transformation of every retained
posterior draw, on original units. It uses the posterior mean, posterior SD and
type-7 equal-tailed 95% credible interval. Modeling the predictor distribution
propagates uncertainty in the population projection; plugging an observed
cluster-mean projection into an otherwise conditional interval would not.
This target is neither a within-site structural coefficient nor the least-squares
coefficient of the realized completed sample. The existing `joint` engine and
`catmech_joint_summary()` keep their meanings and defaults.

The candidate uses a flat prior for `mu` and separate Huang--Wand hierarchical
inverse-Wishart priors for `Sigma` and `Psi`. For dimension `d`, `nu = 2`,
`a_k ~ inverse-gamma(1/2, 1/A_k^2)` and
`C | a ~ inverse-Wishart(nu + d - 1, 2*nu*diag(1/a))`.
The inverse-gamma density convention is proportional to
`a^(-shape-1)*exp(-scale/a)`. Each marginal standard deviation has a half-t
distribution with two degrees of freedom and scale `A_k`; pairwise correlations
are marginally uniform but are not asserted independent of the scales.
After the recorded empirical standardization, `A_k = sqrt(prior_scale)` for
every component of both covariance matrices. The baseline scale of one is fixed
by a unit-scale prior specification, not selected for observed coverage.

Changing `prior_scale` from one to 0.25 or four therefore halves or doubles the
half-t scale. The old inverse-Wishart prior uses the same variance-scale
multipliers, but these are different prior families with different shapes and
moments; they are not equivalent priors. The study compares sensitivity over
declared multiplicative ranges, not a prior-free estimator. Model and prior
change together, so the primary comparison does not isolate their causal
contributions. No unsupported claim attributes all improvement to one component.

## Methods, information and numerical controls

For each public cohort fit the following workflows with identical observed data:

1. `baseline`: exact frozen v0.9.0 sources from commit
   `ae8f94086710ea1acc0df937f45ee847b6a0976b`, stored under
   `validation/v0.9.0-source`. Use `cluster_model = "joint"`. Summarize the same
   fit using the frozen direct summary (`baseline_direct`) and frozen
   CV3/Rubin helper (`baseline_cv3`).
2. `marginal`: the candidate `joint_marginal` engine and new marginal summary
   (`marginal_direct`); also analyze its completions with the unchanged
   CV3/Rubin helper (`marginal_cv3`) as a secondary check.
3. `oracle`: all true `T`, analyzed with the frozen complete-data CV3 helper
   (`oracle_cv3`). This uses additional information and is not a posterior oracle.

Use 40 completed datasets for each imputed fit. Direct summaries use every
retained parameter draw independently of this number. The oracle uses two
identical complete datasets only to access the no-imputation analysis helper.
No unselected true `T` enters an imputer's input or predictor list.

Initial controls for both numerical procedures are four chains, 1,000 burn-in
iterations, 4,000 retained iterations initially and a maximum of 32,000, using
the fixed extension sequence 4k, 8k, 16k, 32k. Rank-normalized R-hat must be at
most 1.05 and bulk and tail ESS at least `100 * chains` for all required finite
monitors. The marginal coefficient transformation is monitored directly.
Support refusals, sampler failures, nonfinite diagnostics and summary-only
failures are recorded, not replaced by successful new cohorts. No discretionary
retry or per-cohort extension is allowed. An engineering-only global control
amendment may follow the masked pilot and must precede development accuracy.

On the first prespecified prior-sensitivity cohorts, refit **both** imputation
models at scales 0.25 and four and retain both analyses. The prior cohorts are
balanced Gaussian, near-zero ICC, observed MAR and contextual X; the selection
does not depend on returned intervals or observed sensitivity.

## Nine data-generating cells

The first seven cells preserve the exact generating code from
[the v0.9.0 protocol](JOINT_CALIBRATION_PLAN.md). There are 24 sites and 240 rows.
Balanced sizes are ten; unbalanced sizes repeat 5, 8, 12, 15 six times. Site
sizes are predetermined and independent of latent effects. In noncontextual
cells generate independent standard-normal `x`, then bivariate independent
site effects `v` and row errors `e` with covariances

```
Omega = [[0.50, 0.45], [0.45, 0.60]]
Sigma = [[0.50, 0.55], [0.55, 0.80]]
T = 0.5 + structural*x + v_T + e_T
W = 0.8 + 1.05*x + v_W + e_W.
```

In contextual cells generate independent `a_site ~ N(0,1)` and
`z_row ~ N(0,1)`, set `x = a_site + z_row`, and add `0.4*a_site` to `T`.
The population marginal slope is `structural + 0.4/2`. This formula applies to
both predetermined site-size schedules. All methods are scored against the
same marginal target, including the misspecified old direct model.

| Index | Cell | Change | Marginal target |
| --- | --- | --- | ---: |
| 1 | `gaussian_balanced` | Ordinary Gaussian, structural 0.7 | 0.7 |
| 2 | `gaussian_unbalanced` | Unbalanced sizes, structural 0.7 | 0.7 |
| 3 | `near_zero_icc` | Multiply `Omega` by 0.02 | 0.7 |
| 4 | `observed_mar` | Observed-data masking below | 0.7 |
| 5 | `null_slope` | Structural coefficient zero | 0 |
| 6 | `student_t5` | Covariance-matched multivariate t5 site/row effects | 0.7 |
| 7 | `contextual_x` | Balanced contextual X, structural 0.7 | 0.9 |
| 8 | `contextual_unbalanced` | Unbalanced contextual X, structural 0.7 | 0.9 |
| 9 | `contextual_null` | Balanced contextual X, structural -0.2 | 0 |

For t5, divide a bivariate Gaussian draw by `sqrt(chisq(5)/3)` to preserve its
covariance. The t5 cell is an explicit distributional misspecification stress.
The ordinary independent-X cells put the predictor's site variance on the
zero boundary, a deliberate stress for the candidate's covariance model.

Randomly select six entirely unaudited sites in every cohort. Ordinarily,
select exactly 60 reference rows from the other 18 sites and mask
`ceiling(0.3 * site_size)` W values independently in each site. References are
exact. All unselected reference values are `NA`. The unbalanced rounding masks
35% of W versus 30% in balanced cells; that contrast does not isolate imbalance.

For observed MAR, W missingness is Bernoulli with log odds `a + 0.8*x`;
solve `a` so the cohort's mean probability is 0.30. Reference selection in the
18 eligible sites is Bernoulli with log odds `b + 0.6*x + 0.6*w_observed`, where
`w_observed` is standardized observed W and zero for missing W. Solve `b` so
expected total audits equal 60. Retain the Bernoulli outcomes without forcing
support or retrying. Neither unobserved W nor unobserved T affects selection.

## Stages and source freezes

The masked engineering pilot uses one cohort per cell, plus both alternative
priors in all four sensitivity cells. It exports timing, support and diagnostic
results but no coefficient estimates, coverage, bias or interval widths.
Development uses 30 cohorts in each of nine cells (270 total), with alternative
priors on all 30 in the four sensitivity cells. Its accuracy is explicitly
development evidence. A scientific change after examining it requires a dated
amendment and distinct development seeds; keep all failed attempts. It may not
select interval multipliers, nominal levels or degrees of freedom for coverage.

After the candidate and interpretation are committed, confirmation uses
**400 fresh cohorts per cell (3,600 total)**, with alternative priors on the
first 200 in each of the four fixed sensitivity cells. This is 10,400 MCMC fits
plus 3,600 oracle workflows and **24,400 analysis requests**. The development
grid is 1,020 MCMC fits plus 270 oracles and 2,310 analysis requests; the masked
pilot is 34 MCMC fits plus nine oracles and 77 analysis requests.

For cell index `j` and replicate `r`, use `base + j*100000 + r`:

| Stage | Data seed base | Imputation seed base | Cohort bootstrap seed |
| --- | ---: | ---: | ---: |
| Pilot | 1670000000 | 1680000000 | none |
| Development | 1690000000 | 1700000000 | 1730000000 |
| Confirmation | 1710000000 | 1720000000 | 1740000000 |

No stage overwrites another directory. At 95% coverage, 400 cohorts give
binomial MCSE about 1.09 percentage points. Prior comparisons with 200 cohorts
per cell use paired uncertainty, not an assumption of independent fitted models.
Worker count is an execution choice (at most 16); seed mapping is invariant.
Record session/package/compiler information and every diagnostic extension.

Hash every current R source, `src/*.c`, `src/*.cpp`, `src/*.h`, `src/Makevars*`, DESCRIPTION,
NAMESPACE, this plan and runner, and frozen baseline sources before and after
each stage. Workers verify those hashes before and after their assigned cohort.
No scored stage completes if source hashes differ or the full planned grid is
absent. Record interrupted executions as incomplete; never silently combine
source versions or remove difficult requests.

## Separate prospective adoption screens

These thresholds are research judgments for this synthetic scale, not universal
statistical laws. Every failed component remains visible.

1. **Implementation and target:** fixed-covariance multivariate-normal
   conditioning oracles must verify mixed-pattern means, variance and shared
   site covariance; independent algebra verifies population slope and intercept,
   scaling transformations and Huang--Wand full conditionals. Check prior
   marginals or equivalent distributional identities, row/site alignment,
   reference preservation, caller RNG, affine recoding, legacy results and
   diagnostic refusals. A material error blocks adoption regardless of coverage.
2. **Contextual correction:** for each of the three contextual cells, the new
   direct summary must have conditional coverage between 92.5% and 97.5% and
   absolute bias at most 0.02. The equal-cell paired coverage-per-request
   improvement over the old direct summary must have a lower two-sided 95%
   Monte Carlo bound above zero, with no negative cell point contrast. Its
   equal-cell mean absolute bias must also be lower. Report the paired
   bootstrap uncertainty for the absolute-bias difference. All comparisons use
   the unchanged marginal targets 0.9, 0.9 and zero.
3. **Ordinary calibration:** in balanced, unbalanced, MAR and null Gaussian
   cells, candidate direct conditional coverage must remain within 92.5--97.5%
   and absolute bias at most 0.02. The baseline already approached 95%; this
   screen does not require a statistically significant improvement over it.
   Report baseline/new coverage distances and paired uncertainty. Near-zero ICC
   and t5 results remain prominent stress results rather than being hidden by
   this four-cell screen.
4. **Coefficient prior sensitivity:** on common complete triples for both
   direct procedures, calculate
   `S_beta = (abs(beta_.25-beta_1) + abs(beta_4-beta_1))/2` on original units.
   The equal-cell mean of `S_beta_candidate - S_beta_baseline` must have its
   upper two-sided 95% Monte Carlo bound below zero, and no cell may have a
   positive point contrast. For `K=4` fixed cells the overall MCSE is
   `sqrt(sum(cell_MCSE^2))/K`. Report all triple denominators, component signed
   shifts, per-scale bias/coverage and failures. Relative width sensitivity
   `mean(abs(log(width_alt/width_1)))` is separate and cannot substitute for
   the coefficient gate. Changes in prior-dependent posterior SD are not used
   as denominators to manufacture reduced sensitivity. Finite-chain posterior
   means contribute Monte Carlo noise to absolute shifts; compare reported
   mean MCSEs with their magnitudes and do not interpret small shifts as exact
   prior independence.
5. **Availability:** each cell must return at least 99% of candidate fits and
   at least 99% of direct intervals, without weaker diagnostics. The candidate
   must have no more observed default-fit refusals than the frozen v0.9.0
   procedure in any cell. Report summary-only failures separately and paired
   failure-rate differences with MCSEs. This is a nonregression screen, not a
   demand for a statistically demonstrable reduction from a near-zero baseline.

For bootstrap uncertainty use 2,000 cohort resamples within each fixed cell,
preserving methods and prior fits together. Set
`RNGkind("Mersenne-Twister", "Inversion", "Rejection")`; loop bootstrap
replicates first and the nine cells in table order second, sampling 400 indices
with replacement per cell in confirmation (30 in development). This one fixed
stream supplies contextual absolute-bias and ordinary coverage-distance
contrasts. The independent verifier may call base R solely to reproduce these
indices and computes all reporting statistics independently.

## Required reporting and independent verification

Save raw public-cohort hashes, seeds, targets, statuses, error classes, warnings,
point estimates, posterior SD or pooled SE, limits, endpoint/mean simulation
errors, diagnostics and attempts. Fit failures propagate to both associated
analyses; a failed direct summary does not erase a valid CV3 result. Report
requested and returned counts, coverage per request (a refusal cannot cover),
coverage conditional on return, Wilson intervals and binomial MCSEs, bias and
MCSE, RMSE and delta-method MCSE, empirical SD, RMS reported spread, mean widths
and 95% interval scores with MCSEs. The interval score is
`width + 40*max(0, lower-truth) + 40*max(0, truth-upper)`.

Paired differences use exactly common successful cohorts for numerical
summaries and all common requested cohorts for availability/coverage-per-request.
Never compare a 200-cohort sensitivity subset to a 400-cohort default summary.
Unreturned intervals have no width or interval score and remain in failure
denominators. Monte Carlo precision does not measure model uncertainty.

An independent verifier must reconstruct the full planned grid, seeds,
denominators, summaries, paired metrics, prior triples, bootstrap contrasts and
separate screens from raw records, and check source hashes. Tamper checks must
reject deleted cells, changed estimates/coverage, altered summaries and changed
source provenance. Existing studies keep their own frozen source versions.

Passing these finite synthetic checks supports an explicit Gaussian marginal
analysis for the tested reference/missingness designs. It does not establish
validity for categorical/incomplete predictors, informative cluster size,
nonlinear predictor relations, rare counts, imperfect references, MNAR, or
using repaired T as an analysis predictor. No claim of a novel statistical
method or of independent human review follows from this implementation.

## Amendments

### Amendment 1: masked-pilot numerical interweaving, before development

The first completed masked pilot used source commit `8e65e4b`: nine cohorts,
34 MCMC fits and nine oracles, with 77 requested analyses. Four MCMC requests
failed the unchanged convergence screen at the fixed maximum; their eight
analysis refusals are retained. The pilot exported no estimates, interval
widths, bias or coverage. Its 35 source files are frozen under
`validation/joint-marginal-pilot-source`; the independent verifier checked all
646 diagnostic rows and the complete grid. This is numerical evidence only.

Before development, add an exact Metropolis--Hastings interweaving move to the
candidate Gibbs sampler. For each response component and sweep, make two
normal log-scale proposals with standard deviation one, jointly rescaling the
corresponding site-effect column and its site-covariance row and column. The
acceptance ratio uses the observation likelihood, the marginal Huang--Wand
site-covariance prior with its auxiliaries integrated out, the site-effect
density and the transformation Jacobian. Redraw those auxiliary variables from
their full conditionals after these moves. The target distribution, data,
covariance priors, transformations, numerical budgets and diagnostic thresholds
remain unchanged. Record the candidate sampler identifier in settings.

Retain the first masked pilot unchanged and run the same pilot cohorts and
seeds under the amended source in `joint-marginal-pilot-interwoven`. This paired
engineering rerun is permitted because no pilot accuracy was exported or used.
Review mixing, refusals and timing, and repeat independent sampler/conditional
oracles. Development and confirmation still use their original untouched seeds
and counts. No scored accuracy preceded this amendment.

### Amendment 2: conditional-regression covariance priors and fresh development

The completed first development stage is preserved under
`joint-marginal-development`, with its exact 35 source files under
`joint-marginal-stage1-source`. Its [development review](JOINT_MARGINAL_DEVELOPMENT_REVIEW.md)
records both failed calibration screens, improved contextual coverage and
remaining positive-context bias. The old direct model had contextual coverage
63.3%, 56.7% and 63.3%; the first candidate had 96.7%, 93.3% and 100%, but its
positive-context biases were -0.0491 (MCSE 0.0180) and -0.0593 (0.0203).
The finite-chain coefficient scale-sensitivity screen passed and all 1,020
MCMC requests returned. These results are development evidence, not confirmation.

The new scientific candidate changes the covariance prior factorization,
keeping the full joint Gaussian likelihood, marginal target and imputation law.
For each row/site covariance `Q`, partition predictors X from measurements
`Y = (T,W)` and use

```
S = Q_XX
B = solve(S, Q_XY)
Omega = Q_YY - transpose(B) %*% S %*% B
Q = [[Omega + B' S B, B' S], [S B, S]].
```

Independently for row and site components, put Huang--Wand nu=2 priors on `S`
and `Omega`, with half-t scales `sqrt(prior_scale)` on their own standard
deviations after the same empirical standardization. Given `Omega`, put a
proper matrix-normal prior on the predictor-by-two coefficient matrix:
`B | Omega ~ MN(0, row_covariance = regression_scale^2 * I,
column_covariance = Omega)`. The mean-vector prior remains
flat. The full covariance's marginal standard deviations are no longer claimed
to have the original full-matrix half-t marginals.

The default `regression_scale = 10` is a broad, proper conditional-regression
scale fixed before the next accuracy run, with separate sensitivity values five
and 20. This parameterization directly states beliefs about measurement-on-X
regressions, rather than using a full-covariance correlation prior as an implicit
regression prior. A flat prior for the site regression is deliberately avoided:
when predictor site variance approaches zero, the regression is weakly
identified and an improper prior need not produce a proper posterior. The new
prior remains an assumption; it is not selected by a grid search for 95%
coverage. Correlation shrinkage was a scientific diagnosis to investigate, not
an independently established causal explanation for all first-stage bias.

The exact Gibbs/full-conditional and Metropolis--Hastings updates must be checked
against independent conditional matrix-normal/inverse-Wishart calculations,
prior simulations and fixed-covariance conditioning oracles. The scale
interweaving moves retain two proposals per component per sweep with normal
log-scale standard deviation one, using the new prior density and its correct
Jacobian; redraw both conditional-block auxiliary sets after those moves.
Monitor the original full covariances and all previous parameters, plus both
conditional measurement covariances and both regression matrices. For one X
and both missing-W patterns this is 35 monitors. R-hat/ESS thresholds and the
4k/8k/16k/32k schedule are unchanged. Record the factorized sampler identifier
and public `regression_scale` in settings; the frozen v0.9.0 model does not
receive or acquire this new control.

Before the factorized-source freeze, add a conservative paired-residual support
screen. Among rows with both exact T and recorded W, residualize both
standardized measurements against X and paired-row site indicators. Require
at least two residual degrees of freedom and rank two for the resulting
two-column residual matrix. Its smaller singular value must exceed
`sqrt(number_of_paired_rows) * 1e-8` on the standardized scale. This restriction
protects the weak covariance-scale prior from an exactly singular measurement
residual likelihood. It is a conservative sufficient support restriction, not
a necessary identification condition or a general proof of posterior propriety;
it can refuse very sparse paired patterns that other information might identify.
Keep all such refusals in the planned fit and analysis denominators, without
resampling or replacing cohorts. The existing `joint` engine is unchanged.

Run a masked engineering pilot on the original nine pilot cohorts/seeds in a
new `joint-marginal-pilot-factorized` directory after freezing this source.
Then run **fresh development D2**, 30 cohorts in every cell, using data base
**1750000000**, imputation base **1760000000** and bootstrap seed
**1770000000**. Cell-index/replicate offsets remain unchanged. No D2 accuracy
precedes this amendment. The original unused confirmation seed bases
1710000000/1720000000 and bootstrap seed 1740000000 remain reserved; no
confirmation begins until the D2 decision and selected source are committed.

On every prespecified sensitivity cohort, run the same three old-model fits
at inverse-Wishart covariance variance-scale settings 0.25, one and four.
Run **five** candidate fits:
HW settings 0.25, one and four with regression scale ten, plus regression scales
five and 20 with HW setting one. Keep both direct and CV3 summaries for all fits.
The candidate's prior-scale setting affects its HW blocks; the regression-scale
setting affects its conditional matrix-normal prior. No averaging of extra
axes is permitted to dilute a larger sensitivity.

The strengthened primary coefficient metric is

```
S_HW = mean(abs(beta_HW=.25 - beta_default), abs(beta_HW=4 - beta_default))
S_B  = mean(abs(beta_B=5 - beta_default), abs(beta_B=20 - beta_default))
S_candidate = max(S_HW, S_B)
S_baseline = baseline S_IW (the same mean absolute shift for its inverse-Wishart scale).
```

Use only cohorts where the old three and candidate five direct analyses all
return. Compute the paired candidate-minus-baseline contrast on these common
sets. Keep the original equal-cell weighting, upper 95% Monte Carlo bound
below zero and no-positive-cell requirements. Report each axis separately,
all signed shifts and every failure. The secondary width statistic uses the
same maximum across the two candidate mean absolute log-width ratios; label it
explicitly as the maximum-axis statistic and also report both component axes.
The comparison spans declared ranges for different prior families, not
identical priors or unrestricted prior robustness. Keep finite-chain mean-MCSE
cautions and do not subtract simulation error from absolute shifts.

The extra two candidate fits and four analyses per sensitivity cohort produce:

- Factorized pilot: 42 MCMC fits plus nine oracle workflows, 93 analyses.
- D2: 1,260 MCMC fits plus 270 oracle workflows, **1,530 workflows and 2,790 analyses**.
- Confirmation: 12,000 MCMC fits plus 3,600 oracle workflows,
  **15,600 workflows and 27,600 analyses** from the same 3,600 fresh cohorts.

All other data-generating cells, targets, calibration/bias screens,
availability screen, diagnostic thresholds and reporting rules remain fixed.
Preserve the earlier failed scientific candidate and its reporting verifier as
`tools/verify-joint-marginal-stage1.py`; the new main verifier checks the amended
grid, both sensitivity axes, stronger gate and new diagnostic monitor set.
This amendment supersedes the initial candidate's covariance-prior description,
D seeds, sensitivity grid/metric and request counts for all future runs.

### Amendment 3: one-paired-site input correction before confirmation

After completed D2 and before confirmation, a read-only input audit identified
a support-check implementation error: constructing factor contrasts for paired
rows from only one site raises a generic R contrasts error. Such data can still
meet the separate observed-site, paired-residual-rank and degrees-of-freedom
requirements. For this case only, construct the paired-site design as a column
of ones. Leave the existing multi-site design path unchanged. This corrects
input handling without changing the likelihood, prior, posterior sampler,
controls, estimand, interval summary or scoring rules.

Preserve all 35 D2 source files under
`validation/joint-marginal-factorized-development-source`, archived at commit
`18a0542`. Its source fingerprint is
`3a6c65470a6198dfd9c85f99c71217446b6f17206113322d0a1c351746009c5c`.
The [equivalence audit](v0.10.0-development/single-paired-site-equivalence.json)
regenerated all 270 D2 public cohorts from their frozen seeds and matched every
recorded input hash. Each cohort used at least 12 paired sites, so the corrected
one-site branch was never entered. All prepared model objects and their
serialized bytes were identical before and after the correction. This audit
performed no fitting or accuracy evaluation; no D2 result was replaced.

The [D2 decision](JOINT_MARGINAL_D2_DECISION.md) freezes the unchanged scientific
candidate for the originally reserved confirmation seeds and amended 27,600-
analysis grid. Freeze the corrected source and this amendment before launching
confirmation. Independently verify D2 against its archived source root and
confirmation against its own before/after source manifest. Every declared
calibration, bias, strongest-axis sensitivity and availability screen remains
unchanged; any unfavorable confirmation result remains reportable.
