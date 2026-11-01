# Using recorded values without trusting them

"Contamination" here means recording error, the usage in the selective-editing literature. It does not mean treatment
spillover between trial arms, the usage in cluster randomized trials.

**Try linear imputation first when you need an interval.** On real recording differences, `catmech_reference(..., use_recorded = TRUE)` and generalized raking
matched or beat this model. In synthetic studies it was strongest when errors were large and separable, or when a reference had a credible
declared SD. A [test on real audit data](../validation/REAL_AUDIT_SYNTHESIS.md) with moderate errors did not reproduce the first of these, and gross, separable errors remain untested on real data.
At a 25% audit the model often had lower point error than linear imputation, but where a reported slope differed from a recount of that one item it kept 38% to 73% of the gap, and its intervals covered 35% to 72% against 90% to 97% for the recount.
Its intervals were too narrow on most slopes even with the link between errors and the covariate removed. Post hoc, separate error means by facility size (`error_strata`) largely repaired the one slope whose errors grew with size, and did not help the others.
See [synthetic evidence](../validation/CONTAMINATION_SYNTHESIS.md) and [real-data evidence](../validation/REAL_DATA_SYNTHESIS.md).

`catmech_contamination()` addresses a gap in `catmech_reference()`. The reference workflow
treats every unreferenced target as missing, so it throws away the recorded value W even
when most recorded values are exactly right. Trusting W, as flag masking does, fails when
flags miss errors. This model sits between the two: each recorded value is either an exact
recording or a recording with additive error, and the reference rows calibrate that mixture.

## Data contract

```r
audit <- catmech_recording_audit(data, "recorded_value", "reference_value", by = c("flag", "site"))

fit <- catmech_contamination(
  data,
  value = "recorded_value",
  reference = "reference_value",
  predictors = c("outcome", "age", "site"),
  corruption_strata = "flag",
  m = 20,
  seed = 123
)
catmech_holdout_check(fit)
```

| Input | Meaning and treatment |
| --- | --- |
| `value` | The recorded measurement W, which may be missing. It is used as evidence and never promoted to truth. Originals stay in `fit$original`. |
| `reference` | Independently obtained measurements on a subset, `NA` elsewhere. By default they are exact and copied into every completion. |
| `predictors` | Complete predictors of the true target, including the analysis outcome and design variables. |
| `corruption_strata` | Categorical columns with separate corruption rates, usually the flag and any site or source effect. Every stratum holding unreferenced records needs calibration rows. |
| `error_strata` | Categorical columns with a separate error mean and SD, for example a site with a systematic shift. |
| `agreement_tol` | Distance within which W and the reference count as agreeing, for example half a rounding unit. |
| `reference_sd` | Declared SD of independent reference error, an assumption to vary rather than an estimate. |
| `selection_delta` | Declared shift for unreferenced targets, a sensitivity analysis for nonignorable selection. |

At least one row must have both W and a reference, and the reference design must have full
column rank. Otherwise the call is refused.

## Model

Write T for the true target, Z for the predictors and A for an unobserved indicator that W
is wrong:

- `T | Z ~ Normal(Z beta, sigma^2)`, with `selection_delta` added for unreferenced rows.
- `A ~ Bernoulli(pi[stratum])`.
- `W = T` when `A = 0`.
- `W | T ~ Normal(T + mu[stratum], tau[stratum]^2)` when `A = 1`.

Exact references reveal A and the error size on calibration rows. For an unreferenced row with
W, integrating out T gives the mixture

```
(1 - pi) * Normal(W; Z beta, sigma^2) + pi * Normal(W; Z beta + mu, sigma^2 + tau^2)
```

so each record receives a posterior probability of being wrong. If `A = 0` its target is W. If
`A = 1` its target is a precision-weighted blend of the regression prediction and `W - mu`.
With a declared `reference_sd`, A is also latent on referenced rows and the reference enters
through a bivariate normal term.

Priors are `p(beta, sigma^2)` proportional to `1 / sigma^2`, `pi ~ Beta(0.5, 0.5)`,
`mu | tau^2 ~ Normal(0, tau^2 / 0.1)`, and `tau^2` scaled inverse chi-squared with 2 degrees of
freedom and scale equal to the residual SD of the reference regression. The error prior matters
when few errors are observed, and it can be overridden through `prior`. With a declared noisy
reference, `sigma^2` receives a proper prior because no row then pins the target exactly.

## Assumptions, and which ones are checkable

| Assumption | If it fails | Evidence available on real data |
| --- | --- | --- |
| Reference selection is ignorable given Z, W and the strata | Bias | Nothing fully tests this. `selection_delta` shows how far conclusions move. |
| References are exact, or their SD is declared correctly | Bias | Near-universal disagreement in the audit is a warning sign. |
| The target is linear and normal given Z | Genuine extremes can be classified as errors | Hold-out check, comparison with `catmech_reference()` |
| Errors are additive normal within strata | Misclassification, though gross errors are usually still caught | Audit of observed error sizes, hold-out check |
| Errors do not depend on outcome, value or other predictors given the strata | Bias. This failed on one slope of a [real facility audit](../validation/REAL_AUDIT_SYNTHESIS.md), where errors grew with size and `error_strata` by size largely repaired it, post hoc. | Audit by candidate columns, the slope of the error on analysis covariates in audited rows (low power, so a null is no reassurance), hold-out check, limited power |
| Missing W is ignorable given Z and the strata | Bias, because W rows inform the target model | Compare W-missing and W-observed reference rows |
| Flags are not computed from W | The same evidence counts twice | Knowledge of how flags were created |

Reference selection may depend on W. A targeted audit of flagged or unusual records is therefore
supported, as long as selection does not depend on the true value beyond what W and Z say. That is
a weaker requirement than `catmech_reference()` has, because that workflow omits W from its model.

## Computation

Each scan updates the corruption rates and error parameters by slice sampling (Neal, 2003) with A
and T integrated out, then draws A and T jointly, then draws `beta` and `sigma^2` from their
conjugate conditional. This two-block structure avoids the slow mixing of the plain data-augmentation
sampler, which conditions the error parameters on imputed errors. When a stratum has fewer than five
observed errors, a scale move (Liu and Sabatti, 2000) also crosses the funnel in the error prior.

Four chains are screened with the same rule as `catmech_model()`: R-hat at most 1.05 and bulk and tail
ESS at least 100 per chain. A failed screen restarts with twice as many iterations, up to
`max_iterations`, and then raises a `catmech_convergence_error`. `manifest$iterations_used` records the
length used. The [independent oracle](../validation/contamination-oracle/README.md) compares this sampler
with a separate Metropolis sampler on the marginal posterior.

## Checks on real data

`catmech_recording_audit()` fits no model. It reports disagreement rates and error sizes overall and by
stratum, tests whether rates differ across candidate columns, and compares recorded values on referenced
and unreferenced rows. Read it before choosing strata.

`catmech_holdout_check()` hides folds of references, refits and scores predictions against the hidden
measurements. It reports predictive coverage, randomized PITs, RMSE and, for this model, the Brier score
and expected against observed error counts. It uses no simulation truth and is the one check of imputation
calibration available on real data. It shares the reference selection process, has limited power with few
errors and cannot detect a reference that is itself wrong. A clean result is not proof of validity. In the synthetic study it
caught a strong undeclared site shift 18 times in 40, but only once in 40 for heavy-tailed targets, where coverage failed. Its `advice`
field reports when the model shows no skill against the recorded value or the stratum error rate, as happened on ClinicalTrials.gov and NHANES. That
is a reason to prefer linear MI or raking. This reading comes from the real-data study and has not been tested prospectively.

## Bootstrap donors for reference matching

`catmech_reference(method = "pmm_abb")` addresses a separate failure in the v0.3 evidence. Standard PMM
covered 85.1% with random 25% references because it conditions on one fixed donor set, which attenuates
between-imputation variance (Gaffert, Meinfelder and Bosch). The new option resamples donors with replacement
for each imputation, an approximate Bayesian bootstrap (Rubin and Schenker, 1986; Siddique and Belin, 2008),
before the same `mice` type-1 matching. `method = "pmm"` is unchanged and stays gated. Coverage rose from 86.7% to 95.9% in the synthetic
25% cell. On real data standard PMM did not repeat its failure and the bootstrap version overcovered, so this is a limited repair.

## What this does not establish

This combines established pieces: contamination models from selective editing (Di Zio and Guarnera, 2013; the
SeleMix R package fits them without validation data), multiple imputation from a validation subsample (Cole, Chu
and Greenland, 2006; Giganti et al., 2020), and data augmentation. No claim of a new identification result or
novel estimator is made. Design-based alternatives such as generalized raking (Lumley, Shaw and Dai, 2011; Oh et
al., 2021) need no recording model and remain better grounded when the recording model is doubtful and selection
probabilities are known. The posterior corruption probability is a model output and does not identify which
record is actually wrong. The workflow has been evaluated on real recording differences and on two real audits, one of health-facility counts and one of
psychology effect sizes. No patient-level clinical source-data audit has been used.

## Sources

- Di Zio and Guarnera (2013). A contamination model for selective editing. *J Off Stat* 29(4), 539-555. [doi:10.2478/jos-2013-0039](https://doi.org/10.2478/jos-2013-0039)
- Guarnera and Buglielli. [SeleMix: Selective Editing via Mixture Models](https://cran.r-project.org/package=SeleMix). R package.
- Cole, Chu and Greenland (2006). Multiple-imputation for measurement-error correction. *Int J Epidemiol* 35(4). [doi:10.1093/ije/dyl097](https://doi.org/10.1093/ije/dyl097)
- Giganti et al. (2020). Accounting for dependent errors in predictors and time-to-event outcomes using electronic health records, validation samples, and multiple imputation. *Ann Appl Stat*. [Project Euclid](https://projecteuclid.org/euclid.aoas/1593449337)
- Lumley, Shaw and Dai (2011). Connections between survey calibration estimators and semiparametric models for incomplete data. *Int Stat Rev* 79(2), 200-220. [doi:10.1111/j.1751-5823.2011.00138.x](https://doi.org/10.1111/j.1751-5823.2011.00138.x)
- Oh, Shepherd, Lumley and Shaw. Raking and regression calibration: methods to address bias from correlated covariate and time-to-event error. [arXiv:1905.08330](https://arxiv.org/abs/1905.08330)
- Gaffert, Meinfelder and Bosch. [Towards an MI-proper predictive mean matching](https://core.ac.uk/download/357228032.pdf). Discussion paper.
- Rubin and Schenker (1986). Multiple imputation for interval estimation from simple random samples with ignorable nonresponse. *JASA* 81(394), 366-374. [doi:10.1080/01621459.1986.10478280](https://doi.org/10.1080/01621459.1986.10478280)
- Siddique and Belin (2008). Multiple imputation using an iterative hot-deck with distance-based donor selection. *Stat Med*. [doi:10.1002/sim.3001](https://doi.org/10.1002/sim.3001)
- Neal (2003). Slice sampling. *Ann Stat* 31(3), 705-767. [doi:10.1214/aos/1056562461](https://doi.org/10.1214/aos/1056562461)
- Liu and Sabatti (2000). Generalised Gibbs sampler and multigrid Monte Carlo for Bayesian computation. *Biometrika* 87(2), 353-369. [doi:10.1093/biomet/87.2.353](https://doi.org/10.1093/biomet/87.2.353)
