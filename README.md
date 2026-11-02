# catmech

**An experimental R package for regression when every record has a recorded value and only an audited subset has a checked reference value. Version 0.7.1.**

**The question:** can you keep error-prone recorded values for precision without trusting them, and which method should you try first?

**The finding:** on public data with real recording differences, simple linear multiple imputation (MI) given the recorded value matched or beat `catmech_contamination()`, the recording-error model this package added. On a real audit the model often had lower point error, but its intervals failed. Start with linear MI when you need an interval.

- **Real data, 25% audits.** Linear MI had lower RMSE than the model on ClinicalTrials.gov (0.0286 against 0.0308) and less bias (+0.0001 against +0.0068). They tied on NHANES (0.0409 against 0.0408). The model's largest edge, at the 10% NHANES audit (0.0435 against 0.0490), was within Monte Carlo noise, and 14 of 200 requests returned no interval there. It could not tell which records were wrong (AUC 0.54 and 0.50).
- **What held up.** Keeping the recorded value pays. Linear MI's intervals were about half as wide as those from discarding it (0.102 against 0.205, and 0.154 against 0.297), with coverage near the all-references oracle (91.8% against 91.2%, and 93.9% against 94.9%). Against simply trusting the recorded value, a small audit often did not reduce error, and nothing tells you in advance when it will.
- **Where linear MI fails.** Under synthetic gross errors it covered 64.7%, against 94.0% for the model. It also fell to 38.0% with noisy references and 74.0% with selection on the true value. If every audited row happens to agree, it keeps the recorded values with no added uncertainty and warns you, and in that case its intervals can be badly wrong.
- **A real audit was inconclusive.** On a [health-facility recount in Rwanda](validation/REAL_AUDIT_SYNTHESIS.md), the rule fixed in advance came out inconclusive at a 25% audit. The model had lower error than linear MI on two of three slopes, but where a reported slope differed from a recount of that item its intervals covered 35% to 72%, against 90% to 97% for the recount. Linear MI's coverage held up better (79.5% at worst), and it cost accuracy where reported values were close: about 1.4 to 2.6 times the error of trusting them on the syphilis slope, on the measures used here, and more if the pool is held fixed. These are facility counts, not patient values, and gross errors remain untested on real data.
- **Clustered audits: experimental.** `catmech_reference(cluster = )` adds a random intercept for each audited site. It replaced a [first version](validation/CLUSTERED_SYNTHESIS.md) that resampled sites, in which resampling within strata understated variance. Under a [new rule fixed in advance](validation/CLUSTERED_FIX_SYNTHESIS.md) the replacement, like the first version, came out **PARTIAL**, on a Rwanda recount at a 25% audit. On six common items where linear MI fell short of a recount, it raised coverage under the variance then in use (from 79.6% to 91.0% on early antenatal visits) with lower error, though `syphilis` reached only 89.0%. On three rare items its error ran 4% to 5% higher, which failed the accuracy condition. That verdict sits within simulation noise of NO-GO. In [simulations where the audited value was a predictor](validation/VARIANCE_MI_SYNTHESIS.md), the option biased the slope, so use it only when the audited value is the outcome. In a [planned comparison with outside methods](validation/OUTSIDE_METHODS_SYNTHESIS.md), it was kept, narrowly: `mice` with the same model matched its accuracy once a default that can drop the recorded value was off, and its lead was somewhat wider intervals. A [test on another health system](validation/SECOND_SYSTEM_SYNTHESIS.md) came out mixed: errors there barely clustered, and the option added little. [Setups written by fresh agents without seeing results](validation/INDEPENDENT_SETUPS_SYNTHESIS.md) did not meet the rule either, but the best of them looks better than the option on most measures, so do not treat the option as the best available. An AI assistant wrote the plans, rules, scoring and write-ups of these rounds, and no independent statistician has reviewed them.
- **Audit planning: not validated.** `catmech_audit_plan()` returned **NO-GO** under a rule fixed in advance. Its median error in predicting interval width missed the 10% limit at a 5% audit (11.3% to 16.8%). Larger targets stayed within 7.8%, but that pattern was not the test. It is a documented calculation, not a validated tool, so it is left out of the example below.

## Which method to start with

| Situation | Start with | Evidence and limit |
| --- | --- | --- |
| Most cases | `catmech_reference(..., use_recorded = TRUE)`: linear imputation given the recorded value. The package does not implement raking. In the recipes tried here, which are not its best possible use, raking tied on recording-difference datasets and had higher error on the real facility audit and in the later comparison. | Best or tied among reference-using methods on ClinicalTrials.gov and NHANES at 25% audits. On a real facility audit its average bias was near zero and its worst-case coverage the best, but it rarely had the lowest point error, and on the slopes tested first its clustered coverage, under the earlier variance, ran about 9 to 11 points below a recount's at 25%. Failed under synthetic gross errors (64.7% coverage), noisy references (38.0%) and selection on the true value (74.0%). |
| Gross, separable errors, or a fallible reference with a credible declared SD | `catmech_contamination()` | Strongest in synthetic cells with these errors (coverage 94.0% under gross errors, 91.8% with a declared reference SD), which remain untested on real data. Worse than linear MI and raking on ClinicalTrials.gov. On a real facility audit with moderate errors it often had lower point error, but its intervals were too narrow. Failed with heavy-tailed targets (87.3%) and with selection on the true value when errors were present (79.3%). |
| Target residuals doubtful and precision not critical | `catmech_reference()` (default) | Assumes least about recorded values, with intervals about twice as wide on the real data. Failed with selective references (7.3%) and noisy ones (52.7%). |
| Audit clustered by site or facility | Linear imputation with `catmech_pool_lm(cluster = )`. The other methods have not been tested with this variance. Add `catmech_reference(cluster = )` only as an experimental check. | The pooled variance is a cluster jackknife, chosen narrowly in a planned test. It repairs the analysis step only: with the audited value as the outcome, linear imputation with it still covered 81% to 86% on the worst items of a real recount at a 25% audit, against 94% to 97% for a recount. That shortfall did not appear in the same way on [another health system](validation/SECOND_SYSTEM_SYNTHESIS.md), where errors barely clustered. The imputation option came out PARTIAL in both of its planned tests, and should not be used when the audited value is a predictor. `mice` offers the same model, with a default to watch. |

"Contamination" here means recording error, not treatment spillover between trial arms.

## Try it

```r
install.packages(c("digest", "jsonlite", "mice", "posterior", "truncnorm", "yaml"))
install.packages(".", repos = NULL, type = "source")   # from the repository directory
```

```r
library(catmech)

# Public NHANES data: self-reported and measured weight. Keep measured weight
# only on a random 25% of rows, as an audit would.
d <- read.csv("validation/real-data/nhanes-pool.csv")
d <- d[d$split == "evaluation", ]
set.seed(1)
audited <- sample(nrow(d), nrow(d) / 4)
d$audit_kg <- NA_real_
d$audit_kg[audited] <- d$reference_kg[audited]

# 1. Read what the audit shows. This fits no model.
catmech_recording_audit(d, "recorded_kg", "audit_kg", by = "sex", agreement_tol = 1)

# 2. Start with linear MI given the recorded value.
fit <- catmech_reference(
  d, "recorded_kg", "audit_kg",
  predictors = c("age", "sex", "race", "height", "sbp"),
  use_recorded = TRUE, m = 20, seed = 123
)
catmech_pool_lm(fit, sbp ~ recorded_kg + age + sex + race + height)
# recorded_kg: 0.126, 95% interval 0.097 to 0.155. With every measured weight: 0.127.
head(fit$audit)   # for each record and imputation: what was replaced, and why
```

This one audit of the whole evaluation pool is an illustration. It is not part of the frozen study.

## The evidence

[![Dot plots for two public datasets at a 25% audit. Linear MI and the contamination model have similar coverage and interval width, close to the all-references oracle. The v0.3 default, which discards recorded values, has intervals about twice as wide. On ClinicalTrials.gov the model has higher RMSE than linear MI.](validation/real-data-figures/real-data-evidence.png)](validation/real-data-figures/real-data-evidence.pdf)

Random 25% audit, 1,000 cohorts of 600 from each evaluation pool. "Reference normal (v0.3)" discards recorded values. Select the figure for its PDF.

## Limits

- **Not clinical.** No real dataset here is a patient-level clinical source-data audit. The Rwanda audit holds facility counts, the effect-size audit comes from psychology, and on ClinicalTrials.gov both values are sponsor-entered, so the reference is not independent truth. A few public datasets say nothing about transport.
- **Synthetic simulations.** All simulation data are synthetic, and the primary cells use the error form the model assumes, which favors it.
- **No novelty or guarantee.** The methods are established, and no new statistical method is claimed. Even the oracle covered only 91.2% on ClinicalTrials.gov.
- **No outside review.** Not on CRAN, with no outside peer review or clinical validation. AI-assisted work is [disclosed](CONTRIBUTING.md).

## Read next

[Assessment](docs/ASSESSMENT.md) | [Real-audit test](validation/REAL_AUDIT_SYNTHESIS.md) | [Clustered option against outside methods](validation/OUTSIDE_METHODS_SYNTHESIS.md) | [Setups by fresh agents](validation/INDEPENDENT_SETUPS_SYNTHESIS.md) | [Another health system](validation/SECOND_SYSTEM_SYNTHESIS.md) | [Rare counts](validation/RARE_COUNTS_SYNTHESIS.md) | [Its variance](validation/VARIANCE_MI_SYNTHESIS.md) | [Its replacement test](validation/CLUSTERED_FIX_SYNTHESIS.md) | [Real-data evidence](validation/REAL_DATA_SYNTHESIS.md) | [Synthetic evidence](validation/CONTAMINATION_SYNTHESIS.md) | [Audit planner test](validation/AUDIT_PLAN_SYNTHESIS.md) | [Merged-paper go/no-go](validation/merged-audit-pilot/README.md) | [Validation record](validation/VALIDATION.md) | [Case study](docs/case-study/CASE_STUDY.md)

<!-- front-page-end: docs/case-study/fact_check.py checks every number above this line. -->

---

**Details below:** [package overview](#package-overview) | [install](#install) | [flags with no references](#flags-with-no-references-declare-a-trust-policy) | [reference subsets](#when-flags-miss-errors-use-actual-reference-information) | [recording audit and contamination model](#with-references-and-recorded-values-workflow-details) | [clustered audits](#clustered-audits-experimental) | [audit planning](#planning-an-audit-not-validated) | [categories](#use-categories-only-under-an-explicit-observation-model) | [scope and evidence](#scope-and-evidence) | [legacy API](#legacy-api) | [reproduce](#reproduce) | [license, credits and AI use](#license-credits-and-ai-use)

## Package overview

catmech does auditable imputation of continuous measurements with quality flags, reference subsets and optional interval categories.

`catmech_model()` masks every flagged measurement and fits a conditional Gaussian
model using trusted observations. It can ignore categories, treat them as exact,
or use an externally specified category-error matrix. It preserves original data,
records assumptions and replacements, checks sampler diagnostics, and supports
linear-model pooling and comparison of declared sensitivity scenarios.
`catmech_reference()` uses genuinely measured reference values on an internal
subset, and `catmech_pmm()` provides an audited standard predictive-matching option.
`catmech_contamination()` models recording errors using those references, and
`catmech_recording_audit()` and `catmech_holdout_check()` report what the references
show on your own data.

The Gaussian model and data augmentation are established methods. The contribution
is a small, inspectable workflow with explicit failure handling. The earlier
retain/replace PMM engine remains available for **sensitivity analysis**; its
arbitrary retention of contaminated values caused severe bias and undercoverage.
See the [assessment](docs/ASSESSMENT.md), [model assumptions](docs/METHOD.md),
and [validation record](validation/VALIDATION.md).

## Install

The [commands above](#try-it) install the dependencies and then the package. Run the second
command from the repository directory, or build the source archive with `R CMD build .` and
install it with `R CMD INSTALL catmech_0.7.0.tar.gz`. R 4.1 or later is required, although only R 4.6
has been tested. GitHub Actions checks Linux, macOS and Windows, and all three passed. See the
[validation record](validation/VALIDATION.md).

## Flags with no references: declare a trust policy

```r
library(catmech)

dat <- simulate_catmech_lab(n = 300, seed = 1)
# Hidden truth is for evaluation only, never an imputation predictor.
dat <- dat[setdiff(names(dat), c("lab_true", "corrupted"))]

fit <- catmech_model(
  dat, value = "lab_value",
  predictors = c("age", "sex", "inflammation", "outcome"),
  flag = "outlier_flag", category_mode = "ignore",
  m = 20, seed = 123
)

catmech_pool_lm(fit, outcome ~ lab_value + age + sex + inflammation)
head(fit$audit)
fit$model_diagnostics$summary
completed_one <- fit$completed[[1]]

# Use a new or empty directory; existing files are not overwritten.
write_catmech(fit, "output/example")
```

Including the downstream outcome in the imputation model matters. All selected
predictors must be complete. Every flagged target is replaced, including false
positive flags. Unflagged observed values are trusted and preserved exactly;
missing flags count as unflagged. This policy requires an ignorable trusted-data
selection mechanism given the predictors. It does not detect missed errors or
estimate whether any particular measurement was wrong.

## When flags miss errors: use actual reference information

Masking flagged values cannot identify errors among measurements still treated
as trusted. If an independently measured reference subset is available, use it
as the observed target and impute every unvalidated target:

```r
# A synthetic example with an audit of whole sites. In real use, reference_value contains
# actual reference measurements on selected rows, with NA elsewhere. Do not manufacture
# reference values from flags or W.
dat_with_reference <- simulate_catmech_lab(n = 600, seed = 1)
dat_with_reference$site <- rep(seq_len(30), each = 20)
audited_sites <- dat_with_reference$site %in% seq(1, 30, by = 3)
dat_with_reference$reference_value <- ifelse(audited_sites, dat_with_reference$lab_true, NA)
dat_with_reference[c("lab_true", "corrupted")] <- NULL

fit_reference <- catmech_reference(
  dat_with_reference, value = "lab_value", reference = "reference_value",
  predictors = c("age", "sex", "inflammation", "outcome"),
  method = "normal", m = 20, seed = 123
)
catmech_pool_lm(fit_reference,
  outcome ~ lab_value + age + sex + inflammation)
```

A random internal reference sample is a transparent design. Targeted verification
requires a suitable selection model and overlap; noisy references cannot be treated
as exact truth. The default omits recorded W as a predictor, sacrificing information
rather than silently retaining unvalidated values. This requires additional reference measurements. Read the
[reference method](docs/REFERENCE_METHOD.md) and
[synthetic example](inst/examples/reference.R).

For a donor-based alternative, `catmech_pmm(dat, "lab_value", predictors,
flag = "outlier_flag")` masks all flags and delegates standard type-1 PMM to `mice`.
`catmech_reference(..., method = "pmm")` restricts donors to reference rows, but
showed substantial undercoverage in the synthetic held-out study. Its pooled analysis
requires `allow_sensitivity = TRUE` and remains labeled `sensitivity_only`. `method = "pmm_abb"`
resamples donors for each imputation. That raised synthetic coverage from 86.7% to 95.9%, but it
overcovered on real data with intervals up to 51% wider than the normal method, and standard PMM did not repeat its failure
there. The reference default remains normal. Predictors must be complete. An observed category can be selected explicitly as
an ordinary factor predictor; no interval constraint is then imposed. PMM cannot
extrapolate beyond donor values and does not guarantee robustness or coverage.
The extension study compares methods given the **same predictor information**.

## With references and recorded values: workflow details

Here "contamination" means recording error, the usage in the selective-editing literature. It is
not treatment spillover between trial arms.

Read what the references show before choosing a model, then pick a method from
[the table above](#which-method-to-start-with). `catmech_recording_audit()` fits no
model. It reports disagreement rates and error sizes, tests whether rates differ across groups
such as site or flag, and compares referenced with unreferenced rows.

```r
catmech_recording_audit(dat_with_reference, "lab_value", "reference_value",
                        by = c("outlier_flag", "site"))
```

When the table points to the contamination model:

```r
fit_contamination <- catmech_contamination(
  dat_with_reference, "lab_value", "reference_value",
  predictors = c("age", "sex", "inflammation", "outcome"),
  corruption_strata = "outlier_flag", m = 20, seed = 125
)
catmech_pool_lm(fit_contamination, outcome ~ lab_value + age + sex + inflammation)
fit_contamination$contamination$corruption    # per-record posterior error probabilities
catmech_holdout_check(fit_contamination)      # scores unseen references, using no simulation truth
```

Selection into the audit may depend on recorded values, so a targeted audit is supported. Selection
on the true value beyond what the recorded value shows is not. `selection_delta` and `reference_sd`
are declared sensitivity assumptions to compare with `catmech_sensitivity()`. If the hold-out check shows
no skill against the recorded value or the stratum error rate, prefer linear MI. On ClinicalTrials.gov
and NHANES the model could not tell which records were wrong (AUC 0.54 and 0.50), so do not read its per-record
probabilities as labels on data like those. The model also assumes errors do not track the other variables. You can check
that on the audited rows, but a null is no reassurance: at a 25% audit of a [real facility recount](validation/REAL_AUDIT_SYNTHESIS.md), a clustered test
reached significance on 18% to 57% of audits where a gap existed. There the model's intervals were too narrow with or
without such tracking. See the [method](docs/CONTAMINATION_METHOD.md),
[synthetic evidence](validation/CONTAMINATION_SYNTHESIS.md), [real-data evidence](validation/REAL_DATA_SYNTHESIS.md)
and [example](inst/examples/contamination.R).

## Clustered audits: experimental

Records from one site are not independent, so pool with a cluster-robust variance:

```r
fit_recorded <- catmech_reference(
  dat_with_reference, "lab_value", "reference_value",
  predictors = c("age", "sex", "inflammation", "outcome"),
  use_recorded = TRUE, m = 20, seed = 126
)
catmech_pool_lm(fit_recorded, outcome ~ lab_value + age + sex + inflammation, cluster = "site")
```

That applies a delete-one-cluster jackknife variance (`cv3`) to each completed dataset, then Rubin's rules with
complete-data degrees of freedom of clusters minus one. A [planned test](validation/VARIANCE_MI_SYNTHESIS.md) chose it,
narrowly, over the earlier default, which undercovered [even with every row checked](validation/VARIANCE_CHECK_RESULTS.md).
Intervals can still be too narrow with few clusters or a few influential ones, and this variance does not repair a shortfall in the imputation itself: linear imputation still undercovered on real items with it. A coefficient that a single cluster determines, such as a site type found at one site, draws a warning, because its interval can be far too narrow. Use `cluster_variance = "cr1s"` to reproduce
earlier results.

If the audited value is the analysis outcome, you can also add a random intercept for the cluster to the imputation
model by passing `cluster` to the fit. `catmech_pool_lm()` then uses that column by default:

```r
fit_clustered <- catmech_reference(
  dat_with_reference, "lab_value", "reference_value",
  predictors = c("age", "sex", "inflammation"),
  use_recorded = TRUE, cluster = "site", m = 20, seed = 126
)
catmech_pool_lm(fit_clustered, lab_value ~ age + sex + inflammation)
```

**Do not use `cluster` in the fit when the audited value is a predictor.** In simulations, it biased that slope toward
zero, and `catmech_pool_lm()` warns about this use. The option's [planned test](validation/CLUSTERED_FIX_SYNTHESIS.md),
with the audited value as the outcome, came out PARTIAL: it helped where the first example's imputation fell short of a
recount and cost some accuracy on rare items. It replaced an [earlier version](validation/CLUSTERED_SYNTHESIS.md) that
resampled audited clusters. If the audit is stratified, put the strata in `predictors`; `cluster_strata` no longer exists. With one audited cluster per stratum, nothing can estimate between-cluster variance, so the fit is refused.
It assumes one normal intercept per cluster. Treat the option as a check, not a default. Fewer than five audited clusters
are refused. Very wide intervals remain possible with few audited clusters, so read the `df` column. See the
[method](docs/REFERENCE_METHOD.md#clustered-audits-experimental).

The same model is available elsewhere. In a [planned comparison](validation/OUTSIDE_METHODS_SYNTHESIS.md), `mice` with `2l.lmer` under
its defaults did much worse, because by default it drops a predictor correlated at 0.99 or more with the imputed variable,
which can remove the recorded value when records are good. Check `loggedEvents` if you use `mice` for audited values. With
that check off, in a rerun made after seeing the results, it matched the option's accuracy on the same cohorts, with slightly
narrower intervals and slightly lower coverage on the real items. Turning it off takes a small wrapper, because `eps = 0` is
passed on to `lmer`, which rejects it; see [the rerun's runner](tools/outside-methods-sensitivity.R). Raking and a difference
estimator had higher error on these small audits. At a 25% audit, the recorded value alone had lower error than every
audit-based method on 6 of 11 real items (five clearly, one tied), and far higher error on others, so a small audit does not
always improve a good record, and nothing here tells you in advance which case you have.

[Comparator setups written by fresh AI agents](validation/INDEPENDENT_SETUPS_SYNTHESIS.md), who saw no results, are in `tools/independent-setups/`.
None met the rule, so the option stays, but `ind_mice`, with separate within-cluster and between-cluster slopes, reached 92% on
all 11 real items at a 25% audit against the option's nine, at about the same accuracy, and largely removed the bias when the
audited value was a predictor. It missed the rule by covering above 98% on two items. These are unreviewed research code and not package features.

Short intervals on count data come from [heavy-tailed cluster effects and not from rarity](validation/RARE_COUNTS_SYNTHESIS.md). Rare counts
with no cluster effect cover well, and common counts with heavy-tailed effects fall just as short. No fix was found in the
package, and no warning worked, because the failures come when an influential cluster is missing from the sample, which the
sample cannot show. Results can also turn on one or two clusters: on [another health system](validation/SECOND_SYSTEM_SYNTHESIS.md), one district
of four clinics reversed the sign of a slope.

## Planning an audit: not validated

`catmech_audit_plan()` estimates from a random pilot how much narrower one coefficient's interval would be at each
audit fraction, keeping or discarding the recorded value, and how many records reach a chosen width. It applies classical
two-phase theory, so this is not a new result. Its prospective test on two public datasets returned **NO-GO** under a rule
fixed in advance, and it should not be described as validated. Small pilots estimated its key quantity `k` poorly, and the
interval for `k` undercovered. By formula, using each whole pool's `k`, keeping the recorded value cuts the audit needed to come within 25%
of the full-audit width from about 380 of 600 records to about 40 to 70 on these data. Read the
[planning guide](docs/AUDIT_PLANNING.md), the [test results](validation/AUDIT_PLAN_SYNTHESIS.md) and the
[example](inst/examples/audit-plan.R) before using it.

## Use categories only under an explicit observation model

```r
labels <- c("Low", "Medium", "High")

# Example sensitivity assumption, not an estimated or clinically validated matrix.
H <- matrix(.05, 3, 3, dimnames = list(labels, labels))
diag(H) <- .90                  # Rows: true category. Columns: observed category.

noisy <- catmech_model(
  dat, "lab_value", c("age", "sex", "inflammation", "outcome"),
  flag = "outlier_flag", category = "lab_category",
  thresholds = c(-Inf, 5.5, 7, Inf), levels = labels,
  category_mode = "noisy", confusion = H, m = 20, seed = 124
)

catmech_sensitivity(
  list(ignore = fit, noisy_assumption = noisy),
  outcome ~ lab_value + age + sex + inflammation,
  term = "lab_value", expected_scenarios = c("ignore", "noisy_assumption")
)
```

`category_mode = "exact"` uses exact bins and rejects contradictions with trusted
numeric values. Noisy categories use the supplied likelihood instead of forcing
values into potentially incorrect bins. Category labels stay unchanged. The matrix
is fixed: uncertainty in its entries is not propagated. A category derived from
the same corrupted measurement is not independent corroboration.

Category models use four Gibbs chains by default. Failed R-hat or effective-sample
size checks raise a `catmech_convergence_error` retaining diagnostic traces. The
[complete example](inst/examples/model.R) shows how to preserve failed scenarios
in the sensitivity report. Its envelope is descriptive, not a confidence interval
with a new coverage guarantee.

## Scope and evidence

- One continuous target, complete selected predictors, independent records (apart
  from the experimental `cluster` option) and a linear, homoskedastic Gaussian
  conditional model.
- Informative categories cover the real line, with right-closed intervals.
  Clinical boundary and rounding conventions must be reconciled first.
- Bayesian parameter uncertainty, posterior predictive completion, per-cell audit
  records, assumption/data hashes, and an externally generated `mids` adapter.
  Do not continue the adapter with `mice.mids()`.
- A prespecified held-out study compares the new workflow with ordinary trusted
  normal MI, trusted PMM, naive and oracle analyses, and the legacy method. It
  reports bias, coverage, MCSE, failures, and deliberate wrong-model stresses.
- Version 0.4 adds a synthetic study with linear MI and generalized raking as strong
  comparators, and a [real-data audit](validation/REAL_DATA_SYNTHESIS.md) on ClinicalTrials.gov and
  NHANES with real recording differences. Neither dataset is a clinical source-data audit.
- Version 0.4.1 adds the [audit planner test](validation/AUDIT_PLAN_SYNTHESIS.md), which returned NO-GO, and
  changes no existing function.
- A later [test on real audit data](validation/REAL_AUDIT_SYNTHESIS.md) uses a Rwanda health-facility recount and a
  psychology effect-size audit. It changes no package code.
- Version 0.5.0 adds the `cluster` options and their [planned test](validation/CLUSTERED_SYNTHESIS.md), which came
  out PARTIAL. Calls without `cluster` are numerically unchanged.
- Version 0.6.0 replaces the resampling in `cluster` with a random intercept after a
  [second planned test](validation/CLUSTERED_FIX_SYNTHESIS.md), which also came out PARTIAL, and removes `cluster_strata`.
  Calls without `cluster` are numerically unchanged.
- Version 0.7.0 changes the default variance in `catmech_pool_lm(cluster = )` to a delete-one-cluster jackknife after a
  [third planned test](validation/VARIANCE_MI_SYNTHESIS.md), and adds a warning. Calls without `cluster` are numerically unchanged.
- A later [comparison with outside methods](validation/OUTSIDE_METHODS_SYNTHESIS.md) kept the `cluster` option under a rule
  that could have retired it. It changes no package code.
- Version 0.7.1 adds warnings and refusals found in an [overall check](validation/OVERALL_QC.md). Numbers from earlier
  studies are unchanged.
- Later planned tests on [another health system](validation/SECOND_SYSTEM_SYNTHESIS.md), of [rare counts](validation/RARE_COUNTS_SYNTHESIS.md) and of
  [independently written setups](validation/INDEPENDENT_SETUPS_SYNTHESIS.md) change no package code.

The [v0.2 coverage protocol](validation/COVERAGE_PLAN.md) transfers useful controls from
the author's carryover, False Floor, calibration and assumption-audit research.
Independent numerical oracles check the category sampler and the reference
model's joint predictive distribution, including shared parameter uncertainty. Clinical-format
synthetic fixtures include missed flags, site shifts and rounded records.
[Transfer details](docs/RESEARCH_TRANSFER_PLAN.md) state what does and does not
carry over. The [extension study](validation/EXTENSION_SYNTHESIS.md)
reports matched-information PMM comparisons, the value of actual reference
measurements, and failures under selective or inaccurate references. These are simulation studies, not clinical validation or a general
coverage theorem. Nonlinearity, missed contamination, informative selection,
incompatible analysis models and dependent records remain limitations.

## Legacy API

`catmech_spec()` and `catmech_impute()` retain the v0.1 PMM sensitivity semantics.
Their `corruption_probability` is an imposed replacement probability, not a
posterior corruption estimate. `catmech_pool_lm()` refuses this engine unless
`allow_sensitivity = TRUE`; that override labels the result as sensitivity-only.
See the [legacy method](docs/LEGACY_METHOD.md), [example](inst/examples/lab.R),
and [historical benchmark](validation/benchmark/RESULTS.md). Existing specifications
are not silently converted to the new model.

## Reproduce

```sh
Rscript -e 'install.packages(c("testthat", "pkgload", "roxygen2"))'
Rscript -e 'testthat::test_local()'
R CMD build .
R CMD check --no-manual catmech_0.7.0.tar.gz
Rscript tools/check-model-oracle.R validation/my-model-oracle
Rscript tools/coverage-benchmark.R confirmatory validation/my-heldout
Rscript tools/coverage-benchmark.R clinical validation/my-clinical
Rscript tools/check-reference-oracle.R validation/my-reference-oracle
Rscript tools/check-pooling-guards.R validation/my-pooling-guards
Rscript tools/extension-benchmark.R distribution validation/my-distribution
Rscript tools/extension-benchmark.R reference validation/my-reference-study
Rscript tools/extension-benchmark.R stress validation/my-reference-stress
Rscript tools/extension-benchmark.R controls validation/my-reference-controls
Rscript tools/check-contamination-oracle.R validation/my-contamination-oracle
Rscript tools/contamination-benchmark.R study validation/my-contamination-study
python3 -I tools/real-data-fetch.py my-raw-directory validation/my-real-data
Rscript tools/real-data-prepare.R my-raw-directory validation/my-real-data
Rscript tools/real-data-benchmark.R study validation/my-real-data-study
Rscript tools/merged-audit-pilot.R study validation/merged-audit-pilot/my-study 400
Rscript tools/audit-plan-benchmark.R study validation/my-audit-plan-study
Rscript tools/clustered-benchmark.R study validation/my-clustered-study
python3 -I tools/clustered-score.py my-clustered-decisions.txt validation/my-clustered-study
Rscript tools/clustered-review-checks.R validation/my-clustered-review
python3 -I tools/clustered-posthoc.py
Rscript tools/clustered-fix-benchmark.R study validation/my-clustered-fix-study
python3 -I tools/clustered-fix-score.py my-clustered-fix-decisions.txt validation/my-clustered-fix-study
python3 -I tools/clustered-fix-posthoc.py
Rscript tools/clustered-fix-replay.R validation/my-clustered-fix-replay
Rscript tools/oracle-variance-check.R validation/my-variance-check
python3 -I tools/oracle-variance-score.py my-variance-decisions.txt validation/my-variance-check
Rscript tools/check-cluster-algebra.R my-cluster-algebra-check.txt
Rscript tools/variance-mi-benchmark.R study validation/my-variance-mi-study
python3 -I tools/variance-mi-score.py my-variance-mi-decisions.txt validation/my-variance-mi-study
python3 -I tools/variance-mi-posthoc.py
Rscript tools/variance-mi-replay.R validation/my-variance-mi-replay
Rscript tools/outside-methods-benchmark.R study validation/my-outside-methods-study
python3 -I tools/outside-methods-score.py validation/my-outside-methods-study my-outside-methods-decisions.txt
Rscript tools/outside-methods-sensitivity.R validation/my-outside-methods-sensitivity
python3 -I tools/outside-methods-sensitivity-score.py
python3 -I tools/outside-methods-posthoc.py
Rscript tools/outside-methods-twophase-check.R my-twophase-check.txt
Rscript tools/second-system-benchmark.R study validation/my-second-system-study
python3 -I tools/second-system-score.py validation/my-second-system-study my-second-system-decisions.txt
Rscript tools/second-system-posthoc.R
Rscript tools/rare-counts-benchmark.R study validation/my-rare-counts-study
python3 -I tools/rare-counts-score.py validation/my-rare-counts-study my-rare-counts-decisions.txt
Rscript tools/independent-setups-benchmark.R study validation/my-independent-setups-study
python3 -I tools/independent-setups-score.py validation/my-independent-setups-study my-independent-setups-decisions.txt
python3 -I tools/independent-setups-posthoc.py
python3 tools/verify-artifacts.py
```

Use new output directories. The `outside-methods` and `independent-setups` runners also need `survey`, `lme4` and `pan`, `tools/second-system-prepare.R` needs `readxl`, and `tools/check-cluster-algebra.R` needs `lme4` and `mitools`. [Validation](validation/VALIDATION.md) lists further commands, including the real-audit study. `tools/clustered-benchmark.R`, `tools/clustered-review-checks.R` and
`tools/clustered-fix-benchmark.R` ran on the resampling that version 0.6.0 replaced, so run those three from commit
`8c81044`. `tools/clustered-fix-replay.R` needs the variance that was the default through 0.6.0, so run it from `8982952`. Archived v0.1 to v0.5.0 sources retain the historical
benchmark inputs; their hashes are checked against those frozen versions. A fresh
ClinicalTrials.gov download will differ from the dated pool because the registry changes.
The v0.3.1 guards reject unsupported analyses and gate experimental reference-PMM
pooling; successful scalar estimates and intervals are unchanged. See
[validation](validation/VALIDATION.md) for evidence, environment and remaining gaps.

## License, credits and AI use

MIT license for code and documents. Data files keep their own terms, described below. Built on [mice](https://doi.org/10.18637/jss.v045.i03),
[data augmentation](https://doi.org/10.1080/01621459.1987.10478458), interval imputation
such as [hmi](https://doi.org/10.18637/jss.v095.i09), and
[multiple overimputation](https://doi.org/10.1177/0049124115585360).
AI-assisted work is disclosed in [CONTRIBUTING.md](CONTRIBUTING.md).
[GitHub notes](docs/GITHUB.md) record the repository's history and metadata.

Real-audit data come from [Nshimyiryo et al. (2020)](https://doi.org/10.1371/journal.pone.0235823),
*PLOS ONE*, under CC BY, and [Maassen et al. (2020)](https://doi.org/10.1371/journal.pone.0233107),
*PLOS ONE*, with data on [OSF](https://osf.io/7nsmd/) under CC0. See
[the data note](validation/real-audit/README.md) for what was copied and changed. Real-data pools come from
NHANES (National Center for Health Statistics) and ClinicalTrials.gov (U.S. National Library of Medicine); see
[their note](validation/real-data/README.md) for sources, changes and terms of use.
