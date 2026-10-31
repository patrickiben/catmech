# Coverage evidence and remaining limits

Evaluated 2026-10-31 against the frozen v0.2 source. This is synthetic validation,
not clinical validation or proof for arbitrary analysis models.

## The two original failures improved through model changes

On the **same new held-out data sets**, masking every flagged measurement and
using ordinary normal posterior-predictive imputation removed the large biases
caused by the legacy retention rule and forced incorrect categories. The
uncertain-flag intervals became substantially narrower while coverage improved;
this was not a repair by widening intervals.

Each row below uses 1,000 independent data sets and 20 imputations. Parentheses
contain Monte Carlo standard errors (MCSEs). Coverage is for nominal 95% intervals.

| Scenario | Workflow | Bias (MCSE) | Coverage (MCSE) | Mean interval width (MCSE) |
|---|---|---:|---:|---:|
| false_positive_flags | Legacy retain/replace + forced categories | -0.296 (0.003) | 0.397 (0.015) | 0.513 (0.004) |
| false_positive_flags | New flag-masking normal model, categories ignored | -0.003 (0.002) | 0.952 (0.007) | 0.255 (0.001) |
| false_positive_flags | Independent standard normal MI | -0.003 (0.002) | 0.953 (0.007) | 0.255 (0.001) |
| category_error | Legacy retain/replace + forced categories | -0.154 (0.002) | 0.346 (0.015) | 0.251 (0.001) |
| category_error | New flag-masking normal model, categories ignored | -0.002 (0.002) | 0.952 (0.007) | 0.255 (0.001) |
| category_error | Independent standard normal MI | -0.002 (0.002) | 0.952 (0.007) | 0.255 (0.001) |

The independent `mice.impute.norm` comparator performs essentially the same as
the new default. This supports the implementation and the changed treatment of
suspect values; it does not establish a new statistical method or superiority
over standard normal MI. The strongest PMM comparator also performs well in
these settings. See the [full summary](coverage-confirmatory/summary.csv) and
[paired comparisons](coverage-confirmatory/paired-comparisons.csv).

## Prespecified primary criterion

The four primary Gaussian settings met the development criterion that the
coverage Wilson interval lie wholly within [.925, .975]:

| New default, categories ignored | Coverage (MCSE) | 95% Wilson interval |
|---|---:|---:|
| missing_only | 0.950 (0.007) | [0.935, 0.962] |
| perfect_flags | 0.963 (0.006) | [0.949, 0.973] |
| false_positive_flags | 0.952 (0.007) | [0.937, 0.964] |
| category_error | 0.952 (0.007) | [0.937, 0.964] |

This is evidence for these four settings under ignorable flag selection and
trustworthy unflagged measurements. It is not a tolerance for clinical error.
The criterion applies to the primary 1,000-replicate settings, not to the smaller
Gibbs subsets or deliberate misspecification stresses.

## Category likelihoods, refusals and uncertainty

With 10% category miscoding, the noisy-category Gibbs model given the correct
fixed confusion matrix returned 250/250 intervals: bias -0.00084 (MCSE .00413),
coverage .932 (MCSE .01592; Wilson interval [.894, .957]), and mean width .24723.
This smaller study is consistent with nominal coverage but does not establish
it with the precision of the primary study. The confusion matrix was supplied
from the synthetic design; it was not estimated or validated on clinical data.

Declaring only 2% category error on those same data produced bias -0.02037
(MCSE .00418) and coverage .920 (MCSE .01716; Wilson interval [.880, .948]).
A coherent computation cannot remove bias from a wrong category likelihood.
The exact-category model correctly refused all 1,000 requests with contradictory
trusted numeric/category observations. Its conditional coverage is **undefined**,
not zero. It did not silently rewrite trusted numeric observations.

The exact-category Gibbs subsets returned 250 intervals each for missing-only,
perfect-flag and uncertain-flag settings, with coverage .976, .984 and .932,
respectively. These differences and their wider Monte Carlo intervals should be
read as smaller-sample evidence. Paired interval widths were smaller than the
standard normal comparator, but this does not establish broad superiority.
See [request outcomes](coverage-confirmatory/request-status.csv) and
[all replicate records](coverage-confirmatory/replicates.csv).

## Negative controls and model misspecification

For the null-effect control, new-default coverage was .944 (MCSE .00727),
corresponding to a false-positive rate .056. Legacy coverage was .983
(MCSE .00409): a low false-positive rate under the null did not reveal its severe
attenuation of the nonzero effect in the uncertain-flag scenario. This implements
the False Floor lesson that one-sided verification can miss another failure.

Under heavy-tailed latent residuals, new-default coverage fell to .925
(MCSE .00833; Wilson interval [.907, .940]); standard normal MI gave .921,
while trusted PMM gave .947. The Gaussian workflow therefore does **not** pass
a broad robustness claim. The exact-category heavy-tail subset returned 249/250
intervals with conditional coverage .9317. One request failed the bulk-ESS screen
(225.3 below 400; R-hat 1.003, tail ESS 434.7), so no interval was returned and no
retry was substituted. This also shows why R-hat alone was not the convergence gate.

## Clinical-format transfer checks

These separate fixtures used 250 data sets per scenario, 600 records per data
set and the same declared scientific model. Category-constrained engines were
excluded because the clinical NORMAL interval includes both reference endpoints;
the current engine uses a different explicit interval convention.

| New default, categories ignored | Bias (MCSE) | Coverage (MCSE) |
|---|---:|---:|
| Clean | 0.006 (0.003) | 0.948 (0.014) |
| Errors detected by flags | 0.005 (0.003) | 0.952 (0.014) |
| Some errors missed by flags | -0.073 (0.003) | 0.644 (0.030) |
| Site shift + poorly detected random errors | -0.134 (0.003) | 0.148 (0.022) |
| Records rounded to two decimals | 0.005 (0.003) | 0.948 (0.014) |

Independent standard normal MI behaves nearly identically. The remaining failures
come primarily from residual recording error in the allegedly trusted data. The
paired oracle decomposition gave selection differences of .00184 (MCSE .00153)
and .00334 (MCSE .00153), versus residual-error differences of -.07848
(MCSE .00166) and -.14126 (MCSE .00205), for missed flags and the compound site
scenario. See the [clinical summary](coverage-clinical/summary.csv) and
[bias decomposition](coverage-clinical/bias-decomposition-summary.csv).

The site scenario combines a +.20 fixed shift at one site, 8% random recording
errors and flag sensitivity of 20%, compared with 60% in the missed-flag scenario.
The analysis includes SITEID, which can absorb a uniform site offset. These
results do not isolate the incremental effect of the site shift or validate
inference under random cluster-correlated outcomes.

Retrospective bias-centering raised the new-default coverage to .916 in both
failed clinical stresses. That diagnostic uses hidden simulation truth and is
not an available correction for real data. Ignoring missed corruption, selecting
trusted records nonignorably, or supplying inaccurate category probabilities
remains outside the validated claim. Removing every flagged observation is not
universally conservative, particularly if flags select genuine extreme values.

## Execution and evidence boundary

The synthetic study recorded 44,750 requests: 43,749 returned intervals, 1,000
expected exact-category refusals and one convergence failure. The clinical study
returned all 7,500 intervals, including the scientifically poor results above.
All failures remain in the denominators. Source hashes were checked before and
after both runs. One-worker and four-worker engineering pilots gave identical
scientific outputs before held-out execution; optimized analysis calculations
agreed with formula-based `lm` to 2.11e-15 in an independent construction check.

The work did not tune methods, intervals or seeds using these held-out outcomes.
The [protocol](COVERAGE_PLAN.md), [pilot comparison](coverage-pilot-check.txt),
[analysis check](coverage-analysis-check.txt), and source/environment records in
each result directory define the scope. The bias-centered diagnostic and its
bootstrap uncertainty are descriptive only. Further validation is needed for
nonlinear relations, informative selection, estimated confusion matrices and
other downstream models.
