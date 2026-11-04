# Joint clustered reference model: initial validation

The explicit joint option supports mixed recorded-present and recorded-missing unaudited rows within one site using a shared conditional site effect. Mathematical and software checks support that implementation. The simulation results support a bounded working-model release, **not a claim of nominal 95% calibration or uniform improvement**: intervals were conservative in the balanced Gaussian cell, covariance-prior scale materially changed their width, and convergence refusals reduced availability in the t stress.

## Design and provenance

The [prospective plan](JOINT_CLUSTER_PLAN.md) fixed six cells with 200 independent cohorts each, 24 sites and 240 rows. There are **1,200 independent cohorts and 4,200 method requests**, with methods paired within cohorts; prior comparisons reuse the first 100 cohorts in three fixed cells. Pilot and aborted executions add no independent scientific evidence. All imputed analyses here use **40 completions** and the same cluster CV3 variance/Rubin pooling. Joint requests use the fixed chain controls in the plan, with the prespecified covariance-prior alternatives. The public API's smaller default number of completions was not assessed by this study.

The reference measurements were exact; six randomly chosen sites were wholly unaudited. Ordinary cells audited 25% of rows from the remaining sites. The MAR cell used observed `x`, available `W`, and recording-pattern information, with cohort-specific logistic intercepts setting expected audit and missingness rates; latent unobserved values did not enter selection. The imputer received gold only in selected cells. The oracle alone received all true values. Balanced cells hid 30% of recorded values; the prescribed rounding within unequal site sizes hid 35% in the unbalanced cell. Thus their contrast does not isolate imbalance alone. The t5 cell preserved the specified covariance matrices but changed conditional distributional shape as well as tails.

The completed study used frozen source commit `4467ca1` and `pan` 2.0 and took 1,904.9 seconds. [Raw requests](joint-cluster-study-final/records.csv), [summaries](joint-cluster-study-final/summary.csv), [paired comparisons](joint-cluster-study-final/paired.csv), [diagnostics](joint-cluster-study-final/diagnostics.csv), [errors](joint-cluster-study-final/errors.csv), [settings](joint-cluster-study-final/settings.json), and [source hashes](joint-cluster-study-final/source-hashes-before.csv) are retained. Before/after/current source hashes agree.

## Coverage and refusals

There were **4,178 returned intervals and 22 convergence refusals**: 15 among 1,200 default-prior joint requests, five among 300 prior-quarter requests, and two among 300 prior-four requests. There were no other request failures in the completed run. Baseline and oracle always returned intervals. Refusals were retained and not retried manually.

Coverage per request counts a refusal as no covering interval. Coverage among returns conditions on an interval being produced; it cannot describe availability. MCSEs below are in percentage points. The Wilson intervals quantify Monte Carlo uncertainty in the reported coverage proportion, not uncertainty for an individual fitted coefficient. Full Wilson intervals, including for failure rates and the comparators, are in the linked summary.

| Cell | Joint intervals / requests | Joint coverage among returns: % (MCSE; Wilson interval) | Joint coverage per request: % (MCSE) | Plain recorded: % (MCSE) | Oracle: % (MCSE) |
| --- | ---: | ---: | ---: | ---: | ---: |
| Gaussian balanced | 198 / 200 | 99.5 (0.5; 97.2–99.9) | 98.5 (0.9) | 96.0 (1.4) | 96.5 (1.3) |
| Gaussian unbalanced | 198 / 200 | 98.0 (1.0; 94.9–99.2) | 97.0 (1.2) | 95.0 (1.5) | 95.0 (1.5) |
| Near-zero ICC | 200 / 200 | 96.5 (1.3; 93.0–98.3) | 96.5 (1.3) | 95.5 (1.5) | 95.0 (1.5) |
| Observed-data MAR | 198 / 200 | 98.0 (1.0; 94.9–99.2) | 97.0 (1.2) | 97.0 (1.2) | 98.5 (0.9) |
| Null slope | 199 / 200 | 98.0 (1.0; 94.9–99.2) | 97.5 (1.1) | 95.5 (1.5) | 93.5 (1.7) |
| Student-t5 stress | 192 / 200 | 96.9 (1.3; 93.4–98.6) | 93.0 (1.8) | 96.0 (1.4) | 96.0 (1.4) |

The plain-recorded baseline itself had observed coverage of 95.0%–97.0% across these cells, so this study does not demonstrate repair of an undercoverage deficit in that comparator. The balanced joint result, 197 covers among 198 returned intervals, is conservative: its Monte Carlo Wilson interval is 97.2%–99.9%, above the nominal 95%. We did not recalibrate or shorten intervals after seeing this. The experiment does not isolate how the covariance prior, analysis variance and pooled degrees of freedom contribute to that conservatism. Other cells have only 200 cohorts, so a result near 95% does not establish one-percentage-point calibration.

The t5 stress returned 192 intervals from 200 requests. Its 96.9% coverage among returns becomes 93.0% per request. The 4.0% refusal rate has a Wilson interval of 2.0%–7.7%. These are different quantities; the latter result should not be hidden by reporting only successful fits. Even zero observed refusals, as in the near-zero-ICC cell, has a nonzero Wilson upper bound (1.9% for 200 requests).

## Error and precision

Bias, RMSE and width summaries are conditional on an interval being returned; denominators therefore differ slightly between methods. The [paired file](joint-cluster-study-final/paired.csv) reports comparisons on the same successful cohorts, and separately compares per-request coverage over every shared request.

| Cell | Joint bias (MCSE) | Joint RMSE | Plain RMSE | Joint mean width (MCSE) | Plain mean width |
| --- | ---: | ---: | ---: | ---: | ---: |
| Gaussian balanced | -0.0003 (0.0057) | 0.0800 | 0.0899 | 0.3950 (0.0039) | 0.4182 |
| Gaussian unbalanced | -0.0018 (0.0063) | 0.0889 | 0.1058 | 0.4082 (0.0046) | 0.4397 |
| Near-zero ICC | 0.0085 (0.0053) | 0.0754 | 0.0723 | 0.3373 (0.0034) | 0.3025 |
| Observed-data MAR | -0.0005 (0.0063) | 0.0883 | 0.1027 | 0.4407 (0.0052) | 0.4932 |
| Null slope | -0.0017 (0.0060) | 0.0850 | 0.1000 | 0.4078 (0.0043) | 0.4328 |
| Student-t5 stress | -0.0033 (0.0062) | 0.0854 | 0.1015 | 0.4079 (0.0055) | 0.4291 |

Joint point estimates had small observed bias. On complete pairs, the joint-minus-plain squared-error difference was negative in five cells and positive in near-zero ICC. This does not establish uniform accuracy superiority. For example, balanced-cell mean width was smaller by 0.0224 (paired MCSE 0.0040), while near-zero-ICC width was larger by 0.0348 (0.0027). In the MAR cell, width was smaller by 0.0518 (0.0059). The joint model's gains here are not obtained by simply widening the existing baseline intervals. The baseline is an existing supported recorded-value model with the same observed information, not a claim to have exhausted stronger competing joint models.

## Prespecified prior sensitivity

The following comparisons use the **same first 100 cohorts** for each alternative and the default scale of one. Coverage differences use all 100 requests, including refusals. Width differences use complete pairs only. These alternatives were sensitivity checks; no prior was selected from these results.

| Cell | Covariance prior scale | Refusals / 100 | Paired change in coverage per request, percentage points (MCSE) | Paired change in width (MCSE; complete pairs) |
| --- | ---: | ---: | ---: | ---: |
| Gaussian balanced | 0.25 | 2 | -2.0 (1.4) | -0.0213 (0.0015; 98) |
| Gaussian balanced | 4 | 1 | +0.0 (0.0) | +0.0711 (0.0030; 99) |
| Near-zero ICC | 0.25 | 0 | +0.0 (0.0) | -0.0235 (0.0022; 100) |
| Near-zero ICC | 4 | 0 | +2.0 (1.4) | +0.0689 (0.0033; 100) |
| Observed-data MAR | 0.25 | 3 | -3.0 (1.7) | -0.0208 (0.0020; 97) |
| Observed-data MAR | 4 | 1 | +1.0 (1.0) | +0.0801 (0.0037; 98) |

The scale-four prior increased width by about 0.069–0.080 in these paired comparisons. It produced 100% coverage among returned intervals in the balanced and near-zero-ICC subsets. That supports treating these empirically scaled covariance priors as substantive assumptions, not as neutral defaults. This limited scale check does not establish robustness to other covariance-prior families or degrees of freedom. Where both the observed paired coverage difference and its MCSE are zero, the indicators agreed in every sampled pair; this does not establish a zero population difference.

## Implementation checks and engineering restart

Independent conditional-Gaussian checks covered means and full covariance, including cross-pattern dependence within a site and entirely unobserved sites. Direct `pan` comparisons, exact gold preservation, original-data preservation, RNG reproducibility, controls and convergence refusal were also tested; see [the joint tests](../tests/testthat/test-joint-cluster.R) and [method specification](../docs/JOINT_CLUSTER_METHOD.md). Such checks establish implementation properties, not model fit in a new application.

The first held-out execution was interrupted because an external CRAN-style check disallowed four simultaneous child workers. Its [350 saved requests](joint-cluster-study/ABORTED.json) were retained, with no accuracy inspection before the execution-only fix. Chains were then batched into at most two fresh workers, keeping the same seeds and statistical settings. The final [72-request pilot](joint-cluster-pilot-final/PILOT_EQUIVALENCE.txt) reproduced the original pilot's status and diagnostics exactly; the original pilot remains tied to commit `554b0c9`.

The independent Python [verifier](../tools/verify-joint-cluster.py) checked the entire planned grid, every seed, source hashes, failure denominators and 966 scalar summaries. **All 350 saved interrupted-run requests exactly match their final counterparts**, including statistical output, status and diagnostics, with only elapsed time and source fingerprint excluded. Temporary-copy negative controls that removed an entire cell, changed a summary bias, or flipped raw coverage were rejected. [Verification log](joint-cluster-study-final/VERIFICATION.txt).

## Scope of the result

This is established multivariate Gaussian multilevel imputation integrated into an auditable reference workflow. It supplies a coherent option for the previously unsupported mixed pattern. It is not a new statistical method, a general repair for undercoverage, or validation of the earlier moment model. The reference truth was exact and selection ignorable; paired support was broad. The repaired value was the analysis **outcome**. Imperfect or nonignorable references, concentrated paired support, analysis-predictor use, contextual effects, random slopes, nonlinear means and clinical counts require separate work. Historical Rwanda and rare-count findings do not validate this new option. No independent human statistical review has been completed.
