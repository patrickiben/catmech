# Validation record: catmech 0.10.0

Date: 2026-11-04. This experimental research release adds an explicit joint
Gaussian model for population marginal coefficients when measurements and
continuous predictors share site effects. Its likelihood and prior are declared
model choices built from established methods. No novel-method or clinical-
validation claim is made. The [v0.9.0 record](v0.9.0-VALIDATION.md) preserves the
previous release, including its contextual stress failure.

## Model and target

`catmech_reference(cluster_model = "joint_marginal", use_recorded = TRUE,
cluster = "site")` models the true measurement, recorded value and complete
continuous predictors jointly. The population marginal regression uses their
total site-plus-row covariance. `catmech_marginal_summary()` transforms every
retained posterior draw to original-unit coefficients and reports posterior
means, standard deviations, equal-tailed credible intervals, convergence
measures and Monte Carlo errors. The summary is independent of the number of
completed datasets. It is an association under the declared population model,
not a causal effect, a within-site slope or an arbitrary downstream regression.

The new covariance prior factors each full covariance into predictor covariance,
conditional measurement covariance and a regression matrix. The first two
blocks have Huang–Wand priors; the regression matrix has a proper matrix-normal
prior. `prior_scale` and `regression_scale` control separate parts. Half-t
standard-deviation and uniform pairwise-correlation properties apply to the
covariance blocks, not to the reconstructed full covariance. Empirical scaling
constants are fixed when defining these priors.

The compiled Gibbs sampler includes posterior-preserving site-scale moves to
improve mixing near zero variance. Four chains, 1,000 burn-in iterations and
4,000–32,000 retained iterations per chain are used. R-hat must be at most 1.05;
bulk and tail ESS must each be at least 100 times the chain count. Failures
remain explicit. The conditional `joint` engine, its priors and its summary
remain available, and CV3/Rubin pooling retains its separate interpretation.
See the [model, equations and restrictions](../docs/JOINT_MARGINAL_METHOD.md).

## Prospective development and confirmation

The [initial protocol and amendments](JOINT_MARGINAL_PLAN.md) preserve the
sequence of decisions. The first masked numerical pilot had four MCMC refusals.
An exact scale-interweaving update then returned every fit in the repeated
masked pilot without changing that candidate's posterior target.

The first accuracy development used 270 cohorts. It substantially improved
contextual coverage but retained contextual slope attenuation; the
[failed candidate and decision](JOINT_MARGINAL_DEVELOPMENT_REVIEW.md) remain
available with exact archived sources. The factorized proper regression prior
was specified before a new masked pilot and a fresh 270-cohort development.
That second development returned all 1,260 MCMC fits and reduced the stronger
coefficient-sensitivity statistic by 68.3%, but its contextual and ordinary
calibration screens still failed. The [D2 decision](JOINT_MARGINAL_D2_DECISION.md)
records those failures and freezes the scientific candidate for confirmation.
D1 and D2 use different cohorts; their difference does not identify a causal
prior-factorization effect.

A one-paired-site input construction error was corrected before confirmation.
All 270 D2 inputs and prepared objects were verified byte-identical because
none entered that branch. The corrected scientific sources and plan were
frozen at `295fea8`, before any confirmation outcomes were evaluated. The exact
35-file tree is preserved under `joint-marginal-confirmation-source/`.

Confirmation uses 400 fresh cohorts in each of nine scenarios, including three
contextual scenarios. Prior comparisons use the first 200 cohorts in each of
four fixed scenarios: three old-model fits and five new-model fits per cohort.
This is 3,600 cohorts, 12,000 MCMC fits plus 3,600 complete-truth oracle workflows,
and 27,600 analysis requests. The new sensitivity score takes the larger
coefficient shift across its covariance-scale and regression-scale axes within
each cohort. Comparisons use common complete old-three/new-five sets, with
Monte Carlo uncertainty. The different prior families and finite-chain error
remain part of the interpretation.

Contextual calibration, ordinary Gaussian calibration, coefficient sensitivity
and availability have separate prespecified screens. Near-zero site variance
and Student-t5 errors are reported as stress scenarios. No result can be
replaced, and no interval multiplier, prior scale or diagnostic threshold is
chosen from confirmation to make a screen pass.

## Confirmation findings

The contextual, ordinary Gaussian and slope prior-sensitivity screens passed.
The availability nonregression screen **failed**. All 3,600 cohorts and 12,000
MCMC requests are retained: 27,596 of 27,600 analysis requests returned intervals.
Two candidate default fits missed the unchanged bulk-ESS threshold at 32,000
retained iterations, causing all four analysis refusals. The baseline had no
refusals. Full findings and Monte Carlo uncertainty are in the
[scientific synthesis](JOINT_MARGINAL_SYNTHESIS.md).

| Contextual scenario | Old direct coverage | New marginal coverage | New slope bias (MCSE) |
| --- | ---: | ---: | ---: |
| Balanced | 58.50% | 95.75% (383/400) | -0.00558 (0.00514) |
| Unequal site sizes | 59.25% | 94.00% (376/400) | 0.00029 (0.00537) |
| Marginal null | 67.00% | 94.74% (378/399) | 0.00788 (0.00510) |

Coverage in the table is conditional on a returned interval; the contextual
null covers 94.50% of all 400 requests. The equal-cell coverage-per-request
gain is 33.17 percentage points (95% Monte Carlo interval 30.18–36.16 points).
Mean absolute contextual bias decreases from 0.13059 to 0.00458; the paired
bootstrap difference is -0.12600 (95% interval -0.13012 to -0.11758).

Ordinary Gaussian coverage is 95.75%, 96.74%, 96.00% and 95.75% for balanced,
unequal-size, observed-MAR and null scenarios. The mean absolute coverage
distance from 95% is slightly greater than the baseline's, and its paired
bootstrap difference includes zero. Intervals in those four ordinary calibration cells are 12.2%–16.3% wider.
Near-zero-variance and Student-t5 stress coverage is 95.25% and 96.50%.
Small negative biases remain, including observed MAR (-0.01224, MCSE 0.00454)
and t5 (-0.01106, MCSE 0.00420). Passing the finite-grid screens does not
establish zero bias, exact nominal coverage or general superiority.

All 800 common old-three/new-five prior sets are available. Mean measured slope
sensitivity falls from 0.0101795 to 0.0030397, a 70.1% reduction. The paired
difference is -0.0071398 (MCSE 0.0002845; 95% Monte Carlo interval -0.0076974 to
-0.0065821); each scenario's contrast is negative. The score uses the larger
new sensitivity axis within each cohort. These remaining shifts are near
finite-chain posterior-mean errors, so no prior-independence claim follows.

Default-fit availability is 3,598/3,600 for the candidate versus 3,600/3,600 for
the baseline; including prior alternatives it is 6,798/6,800 versus 5,200/5,200.
Each candidate cell returns at least 99% of fits and intervals. The stricter
no-observed-increase condition nevertheless fails in unequal-size Gaussian and
contextual-null cells. No failed fit was replaced, no threshold was weakened,
and no confirmation outcome was used for statistical retuning.

The [independent audit](joint-marginal-confirmation/verification.json) verified
all 27,600 requests, 305,600 diagnostic rows, 35 frozen sources, 3,534 summary
values, 504 paired contrasts, all 2,000 bootstrap replicates and all four screen
calculations. Nine deliberate corruptions were rejected. This establishes
reporting correspondence; the mathematical tests and explicit assumptions
address different questions.

## Mathematical and numerical evidence

Independent tests check dense observed-Gaussian conditional means and covariance,
known-covariance parameter posteriors, matched-posterior comparison with `pan`,
Huang–Wand prior marginal distributions, proper matrix-normal/IW updates,
transformed prior densities and scale-move Jacobians. Interweaving on/off
comparisons check posterior means and quantiles, including covariance
auxiliaries. Multiple-predictor affine transformations, original-unit intercepts,
original-observation-only prediction and shared site draws are tested separately.

Public tests cover reference anchors, `mice` adapters, exports, caller RNG
preservation, support refusals, holdout forwarding, the one-paired-site design,
and identical posterior traces/summaries at different imputation counts. The
final source package passed all 1,523 assertions with no warnings or skips.
The release checks and namespace-correction evidence are recorded below.

Exact legacy comparisons retained 48 historical workflow outputs. The historical
artifact verifier independently regenerated saved study summaries and checked
the learning-page numbers and code examples. Separate AI reviews checked the
mathematics/API and native implementation; [review scope](v0.10.0-checks/AI_REVIEW.md)
is documented. These checks establish implementation correspondence, not
universal inferential validity. There has been no independent human statistical
review.

## Release verification

`R CMD check --as-cran --no-manual` on the built `catmech_0.10.0.tar.gz`
finished with zero errors, zero warnings and one routine incoming-check note,
“New submission.” All 1,523 assertions passed without skips. A fresh source
installation passed all seven bundled examples, namespace-only calls, native
sampling, posterior summaries, reference anchors, adapters and JSON/RDS exports.

The pre-release check identified one unqualified `ave()` call. After confirmation
completed, it was changed to `stats::ave`. The
[engineering correction](v0.10.0-checks/NAMESPACE_CORRECTION.md) preserves the
original study and exact source archive. All 3,600 public inputs and full model
preparations were regenerated and verified byte-identical, and a stats-detached
namespace check confirmed the correction. No fit or accuracy calculation was
rerun. Release source differs from confirmation by this one explicit
qualification; it is not mislabeled as an identical source tree.

The archive audit checks all 69 members against the checkout, including the
20 scientific source files and explicit qualification receipt. R's generated
DESCRIPTION fields are checked semantically, and its empty partial-Rd cache
is independently serialized and compared byte-for-byte. Separate controls
accept the intended cases and reject altered sources, receipts and unexpected
archive contents. The source package SHA-256 is recorded with the build evidence.

| Check | Evidence |
| --- | --- |
| Final source check and all 1,523 assertions | `v0.10.0-checks/R-CMD-check.log`, `00check.log`, `tests.log` |
| Fresh installation and all examples | `v0.10.0-checks/install.log`, `installed-smoke.log` |
| Source archive integrity and SHA-256 | `v0.10.0-checks/archive-verification.json`, `source-package-sha256.txt` |
| Exact post-confirmation namespace equivalence | `v0.10.0-checks/release-equivalence.json` |
| Archive acceptance and rejection controls | `v0.10.0-checks/synthetic-release-archive-final-controls.json` |
| Historical compatibility and evidence | `v0.10.0-checks/v0.7.2-unchanged-check.txt`, `v0.3-unchanged-check.txt`, `v0.4.0-postfreeze-check.txt`, `artifact-verification.log` |
| Independent confirmation and corrupted-record checks | `joint-marginal-confirmation/verification.json` |
| Figure arithmetic, provenance and visual inspection | `joint-marginal-figures/verification.json`, `figure-provenance.json`, `visual-QA.txt` |
| New README numbers and historical case study | `../docs/case-study/FACT_CHECK.md` |

The dated case-study HTML/PDF remains historical. GitHub Actions runs package
checks on Ubuntu, macOS and Windows; [workflow results](https://github.com/patrickiben/catmech/actions)
are separate from the local build evidence. The package now needs a C++17 build
toolchain; [building instructions](../docs/BUILDING.md) describe its dependencies.

## Reproduction and scope

Use new output directories:

```sh
Rscript validation/run-joint-marginal.R confirmation validation/my-joint-marginal
python3 -I tools/verify-joint-marginal.py validation/my-joint-marginal --tamper-tests
Rscript tools/plot-joint-marginal.R validation/my-joint-marginal validation/my-marginal-figures
python3 -I tools/verify-joint-marginal-figure.py validation/my-joint-marginal validation/my-marginal-figures
```

The runner records every planned request and source hashes before and after the
stage. The saved study uses sixteen outer workers; execution order does not
change the declared seed mapping. Independent Python verification reconstructs
the request grid, diagnostic screens, interval arithmetic, paired contrasts,
prior comparisons, summaries and every bootstrap result. Only generation of
the declared bootstrap indices delegates to base R.

The repeated-sampling evidence concerns one continuous predictor, independent
sites, noninformative site sizes, exact reference measurements and ignorable
observation patterns. Multiple-predictor support has mathematical/numerical
checks, not a separate coverage study. Categorical, missing, site-constant and
aliased predictors are unsupported. The conservative paired residual-rank
restriction can reject otherwise identifiable sparse designs. This release
does not establish robustness to arbitrary contextual effects, informative site
sizes, nonlinear relationships, imperfect references, MNAR, rare counts or the
repaired measurement used as an analysis predictor. Historical real-audit
findings remain valuable evidence for their own workflows; they do not validate
this new joint predictor model.
