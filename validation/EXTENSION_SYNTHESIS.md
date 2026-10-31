# Frozen v0.3 extension evidence and release interpretation

The substantive improvement is a default normal-MI workflow anchored by actual
reference measurements. In this synthetic experiment, a random, exact 25% reference
sample removed the large bias caused by trusting missed recording errors. It does
not identify truth without those measurements. PMM improved one heavy-tail comparison
but failed other settings, especially the partial-reference experiment; it is not a
universal robustness fix.

This report uses the locally prespecified [plan](EXTENSION_PLAN.md) and frozen v0.3.0
numerical implementation. There were **132,580 requests: 132,540 returned intervals,
40 expected no-reference refusals, and zero other errors**. Every run verified the
same 12 source/plan hashes before and after execution. See the
[pilot check](extension-pilot-check.txt), [completed integrity check](extension-evidence-check.txt),
and the [versioned v0.3.0 source snapshot](v0.3.0-source/SOURCE.txt). The later v0.3.1 release guards reject
unsupported analysis inputs and require an explicit sensitivity override for pooled
partial-reference PMM. They do not alter the imputation algorithms or retune these
results. Complete-reference PMM still reduces to ordinary complete-data analysis.

The planned estimand is the adjusted outcome coefficient of the error-free continuous
measurement: 0.8, except the prespecified null experiment (0). Tables show estimates
of repeated-sampling performance with Monte Carlo standard errors (MCSE), not standard
errors for one scientific dataset. Conditional coverage and per-request covering yield
agree in all performance settings because all requested intervals were returned.
No-reference refusals remain visible in the control denominator.

## Equal-information normal/PMM comparisons

Each cell uses 1,000 independent synthetic datasets, n=300 and 20 imputations. The
2×2 comparison holds predictor information constant within normal/PMM contrasts.
Categories are ordinary factor covariates, not forced interval constraints. Their
inclusion is also a model change: even a Gaussian base DGM becomes truncated and
nonlinear after conditioning on an exact target-derived category.

| Setting | Method | Bias (MCSE) | Coverage (MCSE) | 95% Wilson interval | Width (MCSE) |
|---|---|---:|---:|---:|---:|
| Heavy tails | Normal; no category | 0.0058 (0.0029) | 0.9170 (0.0087) | [0.8983, 0.9325] | 0.2727 (0.0016) |
| Heavy tails | Normal + category | 0.0056 (0.0030) | 0.9110 (0.0090) | [0.8917, 0.9271] | 0.2646 (0.0014) |
| Heavy tails | PMM; no category | 0.0090 (0.0024) | 0.9400 (0.0075) | [0.9235, 0.9531] | 0.2784 (0.0019) |
| Heavy tails | PMM + category | 0.0138 (0.0024) | 0.9260 (0.0083) | [0.9081, 0.9406] | 0.2719 (0.0020) |
| Heavy tails + category errors | Normal; no category | 0.0001 (0.0026) | 0.9130 (0.0089) | [0.8939, 0.9289] | 0.2725 (0.0014) |
| Heavy tails + category errors | Normal + category | -0.0013 (0.0026) | 0.9120 (0.0090) | [0.8928, 0.9280] | 0.2698 (0.0014) |
| Heavy tails + category errors | PMM; no category | 0.0036 (0.0025) | 0.9280 (0.0082) | [0.9103, 0.9424] | 0.2782 (0.0015) |
| Heavy tails + category errors | PMM + category | 0.0060 (0.0025) | 0.9300 (0.0081) | [0.9125, 0.9442] | 0.2793 (0.0015) |
| Strong MAR | Normal; no category | -0.0049 (0.0024) | 0.9490 (0.0070) | [0.9336, 0.9610] | 0.2939 (0.0011) |
| Strong MAR | Normal + category | 0.0023 (0.0022) | 0.9430 (0.0073) | [0.9269, 0.9557] | 0.2723 (0.0007) |
| Strong MAR | PMM; no category | -0.0152 (0.0026) | 0.9110 (0.0090) | [0.8917, 0.9271] | 0.2803 (0.0011) |
| Strong MAR | PMM + category | -0.0122 (0.0024) | 0.9300 (0.0081) | [0.9125, 0.9442] | 0.2754 (0.0008) |

With heavy tails and no category, PMM increased coverage by **0.0230 (paired MCSE
0.00725)** relative to normal MI. The corresponding interval-width difference was
+0.00574 (MCSE 0.00103). Its coverage was nevertheless 0.940, with Wilson interval
[0.9235, 0.9531]; this is evidence of improvement in this experiment, not a general
coverage guarantee. Under strong MAR, PMM instead lost **0.0380 coverage (paired MCSE
0.00840)** and its bias was −0.01525 (MCSE 0.00265). Adding categories did not supply
a uniformly successful remedy. The heavy-tail/category-error and heavy-tail/no-error
scenarios use different held-out streams; differences between those rows alone do
not isolate a causal category-error effect.

The six Gaussian no-category normal coverage estimates were: missing-only 0.940,
perfect flags 0.958, false-positive flags 0.953, category errors 0.954, null effect
0.955, and strong MAR 0.949. **Five of six met the prespecified screen** requiring
the 95% Wilson interval to lie within [0.925, 0.975]. Missing-only narrowly failed:
MCSE 0.00751 and Wilson [0.9235, 0.9531]. The independent mice normal comparator also
returned 0.940 in that cell. A simulation screen with finite Monte Carlo uncertainty
is not a universal validation theorem, and the failed screen was not retuned away.

The new PMM outputs exactly match independent mice PMM estimates, SEs, degrees of
freedom and intervals in all 16,000 distribution pairs. The normal comparisons also
track the independent normal implementation. These are audited workflows around
established methods, not claims of a statistically superior new PMM algorithm.
[Raw summaries](extension-distribution/summary.csv) and
[paired contrasts](extension-distribution/paired-comparisons.csv) retain every cell.

## What actual reference measurements change

The clinical-format fixture uses n=600, a random Bernoulli(0.25) reference sample,
20 imputations, complete outcome/covariates/site, and no recorded-W or category
auxiliaries. The mean reference count was 150.141. Only selected reference cells
contain simulated truth in the supplied data; all other reference cells are NA.
Normal/reference methods use exactly the same gold cells as their independent
mice comparators and reference-only regression.

| Setting | Method | Bias (MCSE) | Coverage (MCSE) | 95% Wilson interval | Width (MCSE) |
|---|---|---:|---:|---:|---:|
| Same 1,000 cohorts | Full-truth oracle | 0.0016 (0.0013) | 0.9520 (0.0068) | [0.9369, 0.9636] | 0.1612 (0.0002) |
| Same 1,000 cohorts | Reference-only regression | 0.0048 (0.0027) | 0.9530 (0.0067) | [0.9381, 0.9645] | 0.3304 (0.0009) |
| Same 1,000 cohorts | Reference normal (default) | -0.0027 (0.0022) | 0.9520 (0.0068) | [0.9369, 0.9636] | 0.2782 (0.0013) |
| Same 1,000 cohorts | Independent mice normal | -0.0027 (0.0022) | 0.9470 (0.0071) | [0.9313, 0.9593] | 0.2790 (0.0013) |
| Same 1,000 cohorts | Reference PMM (experimental) | 0.0068 (0.0022) | 0.8510 (0.0113) | [0.8276, 0.8717] | 0.2122 (0.0010) |
| Same 1,000 cohorts | Independent mice PMM | 0.0068 (0.0022) | 0.8510 (0.0113) | [0.8276, 0.8717] | 0.2122 (0.0010) |
| Clean records | v0.2 mask-only normal | 0.0005 (0.0013) | 0.9640 (0.0059) | [0.9506, 0.9739] | 0.1723 (0.0003) |
| Missed errors | v0.2 mask-only normal | -0.0771 (0.0015) | 0.5810 (0.0156) | [0.5502, 0.6112] | 0.1747 (0.0003) |
| Site shift + missed errors | v0.2 mask-only normal | -0.1393 (0.0015) | 0.1270 (0.0105) | [0.1078, 0.1491] | 0.1727 (0.0004) |

The reference-normal coverage was **0.952 (MCSE 0.00676)**, Wilson [0.9369, 0.9636],
meeting the prespecified screen. Its bias was −0.00273 (MCSE 0.00217). Independent
normal MI gave 0.947 coverage; the paired difference was 0.0050 (MCSE 0.00436), with
no substantive advantage established for the package implementation.

Relative to reference-only regression, normal MI shortened mean intervals by
0.05221 (paired MCSE 0.00106), while the coverage difference was −0.0010 (MCSE
0.00755). It uses the observed outcome and covariate information from all records,
conditional on its model. It does not make imputed values equivalent to new assays.
The full-truth oracle remains substantially more precise.

Relative to the old mask-only workflow, adding references reduced missed-error bias
from −0.07707 to −0.00273 and increased coverage by 0.3710 (paired MCSE 0.01731).
In the compound site/recording-error stress, coverage increased by 0.8250 (MCSE
0.01321). These are **different-information comparisons**: the successful method
receives new measurements. Its intervals are wider than mask-only intervals, but
its bias also changes substantially; the improvement is not an arbitrary interval
inflation. In clean data, discarding reliable recorded W has a precision cost:
normal-reference width exceeds mask-only width by 0.10584 (MCSE 0.00130).

The site setting combines a +0.20 shift at one site, 8% random errors and only 20%
flag sensitivity; SITEID is included. It does not isolate the incremental effect
of a uniform site shift, which site adjustment can partly absorb.

All five recording scenarios deliberately reuse the same 1,000 cohorts, reference
subsets and imputation seeds. Reference-normal, reference-PMM, gold-only and oracle
numerical results are identical across them. This verifies that unused W does not
leak into reference fitting. It is **one paired reference experiment, not five
independent demonstrations of coverage**.
[Reference summaries](extension-reference/summary.csv) and
[paired contrasts](extension-reference/paired-comparisons.csv) provide the full evidence.

## Why partial-reference PMM is gated in the release

Under the ideal random, exact 25% reference design, reference PMM covered only
**0.851 (MCSE 0.01126)**, Wilson [0.8276, 0.8717]. This is a confirmed failure, not
an inconclusive small deviation from 95%. The independent mice PMM implementation
produced exactly the same fits. Its empirical sampling SD was 0.07099 (bootstrap
MCSE 0.00172), but its mean reported SE was 0.05346 (MCSE 0.00025). Normal-reference
MI had empirical SD 0.06864 and mean SE 0.06867. Retrospective bias centering only
raised PMM coverage to 0.863, so point-estimate bias alone does not explain the
undercoverage. That diagnostic uses hidden truth and is not a correction for users.

The v0.3.1 policy therefore keeps partial-reference PMM experimental and requires
`allow_sensitivity = TRUE` for `catmech_pool_lm()`, labeling its output
`sensitivity_only`. The default reference method remains normal. The PMM algorithm
and all unfavorable evidence are retained unchanged; the gate is an interpretation
and analysis safeguard, not an accuracy-tuned estimator.

## Reference quantity, selection and accuracy limits

The 10% and 50% sensitivity studies each use 250 cohorts; the 25% primary row uses
1,000. Reference uniforms are shared, so subsets are nested where applicable.
The missed-error and compound-site reference results are again identical, not
independent replications. At 10%, 25% and 50%, mean reference counts were 59.008,
150.141 and 298.472, respectively.

| Setting | Method | Bias (MCSE) | Coverage (MCSE) | 95% Wilson interval | Width (MCSE) |
|---|---|---:|---:|---:|---:|
| 10% references; R=250 | Reference normal | -0.0204 (0.0070) | 0.9600 (0.0124) | [0.9279, 0.9781] | 0.4714 (0.0080) |
| 10% references; R=250 | Reference PMM | 0.0033 (0.0082) | 0.7480 (0.0275) | [0.6907, 0.7978] | 0.3000 (0.0072) |
| 10% references; R=250 | Reference-only regression | 0.0057 (0.0085) | 0.9560 (0.0130) | [0.9229, 0.9753] | 0.5672 (0.0062) |
| 25% references; R=1000 | Reference normal | -0.0027 (0.0022) | 0.9520 (0.0068) | [0.9369, 0.9636] | 0.2782 (0.0013) |
| 25% references; R=1000 | Reference PMM | 0.0068 (0.0022) | 0.8510 (0.0113) | [0.8276, 0.8717] | 0.2122 (0.0010) |
| 25% references; R=1000 | Reference-only regression | 0.0048 (0.0027) | 0.9530 (0.0067) | [0.9381, 0.9645] | 0.3304 (0.0009) |
| 50% references; R=250 | Reference normal | -0.0021 (0.0032) | 0.9600 (0.0124) | [0.9279, 0.9781] | 0.2059 (0.0013) |
| 50% references; R=250 | Reference PMM | 0.0025 (0.0033) | 0.9080 (0.0183) | [0.8657, 0.9379] | 0.1876 (0.0010) |
| 50% references; R=250 | Reference-only regression | 0.0009 (0.0036) | 0.9600 (0.0124) | [0.9279, 0.9781] | 0.2313 (0.0009) |

Normal MI retained broad intervals at 10% references, with a small negative bias
of −0.02037 (MCSE 0.00695); successful coverage alone does not erase that bias. More
references improved precision. PMM undercoverage persisted across all three tested
fractions, including 0.748 at 10% and 0.908 at 50%. These smaller studies were
prespecified sensitivity analyses, not additional primary validation claims.

| Setting | Method | Bias (MCSE) | Coverage (MCSE) | 95% Wilson interval | Width (MCSE) |
|---|---|---:|---:|---:|---:|
| Nonignorable references; R=1,000 | Reference normal | 0.2639 (0.0024) | 0.0550 (0.0072) | [0.0425, 0.0709] | 0.2808 (0.0013) |
| Nonignorable references; R=1,000 | Independent mice normal | 0.2638 (0.0024) | 0.0560 (0.0073) | [0.0434, 0.0720] | 0.2812 (0.0014) |
| Nonignorable references; R=1,000 | Reference PMM | 0.2809 (0.0025) | 0.0140 (0.0037) | [0.0084, 0.0234] | 0.2106 (0.0008) |
| Nonignorable references; R=1,000 | Reference-only regression | 0.2682 (0.0029) | 0.1360 (0.0108) | [0.1161, 0.1586] | 0.3388 (0.0011) |
| Noisy references; R=500 | Reference normal | -0.1116 (0.0029) | 0.6780 (0.0209) | [0.6358, 0.7175] | 0.2749 (0.0019) |
| Noisy references; R=500 | Independent mice normal | -0.1116 (0.0029) | 0.6740 (0.0210) | [0.6317, 0.7136] | 0.2762 (0.0019) |
| Noisy references; R=500 | Reference PMM | -0.1063 (0.0030) | 0.5120 (0.0224) | [0.4683, 0.5556] | 0.2141 (0.0015) |
| Noisy references; R=500 | Reference-only regression | -0.1043 (0.0036) | 0.7720 (0.0188) | [0.7332, 0.8066] | 0.3213 (0.0012) |

In the nonignorable-selection stress, reference selection depends on a latent
residual beyond the supplied complete predictors and outcome. A **cohort-specific**
intercept is solved so the average selection probability is 0.25. Normal-reference
MI then has bias +0.26385 and coverage 0.055. Reference-only regression also fails;
its coverage is 0.136. The source of failure is not fixed by choosing a simpler
analysis of the same unrepresentative references.

The noisy-reference stress adds independent assay error with SD 0.08 to selected
references while retaining the error-free estimand. Normal-reference bias is
−0.11163 and coverage 0.678. Independent normal MI fails similarly. Calling an
inaccurate measurement “gold” does not make it exact.

A paired anchor diagnostic supports these mechanisms: for nonignorable selection,
the gold-only retained-truth slope differs from the full-truth oracle by +0.26658
(MCSE 0.00263), with zero reference measurement error. For noisy references, the
selection component is +0.00553 (MCSE 0.00322) and the supplied-reference-error
component is −0.11272 (MCSE 0.00174). These decompose the retained-record analysis,
not the MI estimator itself. They use simulation truth and cannot be computed as a
real-data fix. See [stress summaries](extension-stress/summary.csv),
[paired comparisons](extension-stress/paired-comparisons.csv), and
[anchor decomposition](extension-stress/bias-decomposition-summary.csv).

## Functional controls and claim boundary

All 40 complete-reference calls preserve every target exactly and agree with direct
linear regression in estimates, SEs, degrees of freedom and confidence limits;
maximum absolute numerical difference is 2.66×10^-15. All 40 no-reference calls
refuse despite recorded W being available. See [control records](extension-controls/replicates.csv).

This work supports an audited implementation of normal MI given a credible internal
reference subset and compatible models. It also supplies a standard PMM option for
explicit comparison, with a strong safeguard for the failed partial-reference use.
It does not establish a new identification result, guarantee inference under heavy
tails or arbitrary MAR, infer unmeasured truth without references, repair selective
validation automatically, or account for unknown reference-assay error. No private
or real clinical measurements were used in this evaluation.
