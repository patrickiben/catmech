# Prospective staged study: joint-model inference, priors and computation

Written 2026-11-04 before candidate accuracy runs. Development and confirmation are separate. The completed v0.8.0 study is exposed development evidence; its numerical results will not be relabeled as fresh validation. Amendments and the selected candidate must be committed before confirmatory runs.

## Questions and fixed estimand

This study separates three questions: whether direct posterior inference improves the calibration of the current pooled intervals; whether its relative sensitivity to the covariance prior decreases; and whether a larger, faster numerical sampling budget reduces convergence refusals without weaker diagnostics. Success in one question is not success in the others.

The initial statistical candidate keeps the v0.8.0 joint Gaussian likelihood, flat fixed-effect prior and inverse-Wishart covariance priors unchanged. It summarizes the population mean coefficients `B_T` from all retained posterior draws, using the explicit `catmech_joint_summary(result, level = 0.95)` API. With independent zero-mean site effects, these coefficients also describe the marginal mean of `T` given the selected covariates. They are not the coefficients of a particular completed sample's OLS fit, and they are not cluster-robust inference under arbitrary misspecification.

The summary returns the coefficients for the exact selected model-matrix columns; it does not accept an arbitrary reduced analysis formula. On original units, slopes are `b_k = s_T * b*_k / s_Xk`, and the intercept is `c_T + s_T*b*_0 - sum(c_Xk*b_k)`. The reported estimate is the posterior mean, spread is the posterior SD, and the interval consists of type-7 equal-tail posterior quantiles. Report original-scale R-hat, bulk/tail ESS and Monte Carlo SEs for the mean and both endpoints. Failed convergence, a nonjoint workflow, and a result with no fitted coefficient traces are refused. In particular, the all-reference bypass is not silently turned into a posterior fit. `catmech_pool_lm()` retains its current interpretation and calculations.

## Diagnosis already exposed

The v0.8.0 balanced-cell joint result had mean reported variance / empirical estimate variance about 1.23 and a mean interval critical multiplier about 2.24, with median pooled degrees of freedom about 10, although there were 24 sites. The MAR and near-zero-ICC cells behaved differently. The complete-truth CV3 analysis itself showed variance inflation in some cells. These findings do not identify a single cause or justify multiplying widths or retuning degrees of freedom. They motivate separating the model posterior from the OLS/CV3/Rubin analysis.

The recorded 22 convergence refusals concerned ESS rather than final R-hat. They may be replayed as an explicitly selected numerical challenge set, examining diagnostics and runtime only. Such a selected set is not a coverage study and is not independent evidence about the frequency of future refusals.

## Procedures and contrasts

For every cohort, fit the following procedures on identical observed `W`, exact selected references and complete `x`:

1. `old_cv3`: frozen v0.8.0 joint imputation, prior scale one, four chains, burn-in 1,000, initial 2,000 retained iterations and maximum 8,000; 40 imputations analyzed by the frozen CV3/Rubin helper. The exact frozen sources are under `validation/v0.8.0-source`, from commit `b236293`. Load them in an isolated environment; do not substitute an unrecorded installed package.
2. `new_cv3`: the initial candidate joint imputer with the same likelihood and priors, four chains, burn-in 1,000, initial 4,000 retained iterations and maximum 32,000, followed by the current CV3/Rubin analysis of 40 imputations.
3. `direct`: posterior `B_T` summary from the **same fitted candidate object** as `new_cv3`, using every retained coefficient draw rather than just the 40 imputation draws. The primary term is `x`.
4. `oracle_cv3`: all true `T` values followed by the frozen complete-data CV3 analysis. This uses information unavailable to the other procedures and measures the analysis's own finite-sample behavior; it is not a posterior oracle.

The numerical candidate retains R-hat at most 1.05 and both ESS measures at least 100 times the chain count, with all entries finite. It additionally monitors the original-scale intercept. The fixed extension schedule is 4,000, 8,000, 16,000, 32,000. The initial candidate keeps the restart streams and priors unchanged; trace-vectorization is an execution improvement. No manual extension of an individual scored request or retry until success is allowed.

For each prespecified prior-sensitivity cohort, also fit covariance-prior scales 0.25 and 4 with the candidate numerical controls, returning both CV3 and direct summaries from each shared fit. The residual and site inverse-Wishart scales are changed together, exactly as in v0.8.0. The old-profile comparator is run only at prior scale one.

The contrasts are: old versus new CV3 for the numerical procedure; new CV3 versus direct posterior inference for the analysis choice; the three covariance-prior scales within each analysis; and complete-truth versus incomplete-data results to distinguish analysis limitations from missing/reference information. These do not isolate a single mechanism when a stress cell changes several distributional properties.

## Cohorts

Reuse the six fully specified generating mechanisms in [JOINT_CLUSTER_PLAN.md](JOINT_CLUSTER_PLAN.md): Gaussian balanced, Gaussian unbalanced, near-zero ICC, observed-data MAR, null slope, and variance-matched multivariate t5 stress. Preserve their means, covariance matrices, site counts/sizes, masking rules and exact reference measurement. In particular, unbalanced rounding hides 35% of `W`, compared with 30% in balanced cells; that contrast does not isolate cluster-size imbalance.

Add `contextual_x` as an explicit assumption-violation stress. There are 24 balanced sites of 10 rows. Generate independent `a_j ~ N(0,1)` and `z_ij ~ N(0,1)` and set `x_ij = a_j + z_ij`. Generate `v_j` with the ordinary bivariate site covariance and row errors with the ordinary covariance. Set

```
T_ij = 0.5 + 0.7*x_ij + 0.4*a_j + v_Tj + e_Tij
W_ij = 0.8 + 1.05*x_ij + v_Wj + e_Wij.
```

Use the ordinary random audit and recording masks, not the MAR mask. The marginal population slope in `T ~ x` is **0.9**, because `Cov(a,x)/Var(x) = 1/2`. The within-site structural coefficient is 0.7. Do not score the direct random-intercept coefficient against 0.7 while comparing the other analyses against 0.9. All reported marginal-mean coverage uses 0.9; explicitly label this as a misspecification/estimand stress, because the fitted model's independent-site-effect assumption fails. It is not part of the primary calibration adoption gate, and failure here must limit claims about marginal effects under contextual dependence.

Hidden true values are supplied only to the complete-truth oracle and scoring. Nonselected references remain `NA` in every imputer's input. The same public-cohort fingerprint, data seed and imputation seed identify paired requests; the oracle's additional truth is separately declared.

## Stage D: bounded development

Run 30 new cohorts per cell, including `contextual_x`: 210 cohorts. Each produces the four main analysis requests above. On all 30 development cohorts in balanced Gaussian, near-zero ICC and MAR, run the two alternative prior scales and both analyses. That is **1,200 analysis requests**, from 600 MCMC fits plus 210 complete-truth oracle analyses. CV3 and direct summaries reuse each candidate fit.

The development data seed is `1610000000 + cell_index*100000 + replicate`; imputation seed is `1620000000 + cell_index*100000 + replicate`. Cell indices follow the six-cell order in the earlier plan, then `contextual_x` seventh. Do not add successful-cohort replacements or remove difficult datasets.

Development may expose accuracy. It can reveal implementation defects, model/estimand incompatibility, and whether direct inference reduces relative prior sensitivity. It cannot be used to select a multiplicative interval correction, alter the nominal level, or choose priors by how closely a small sample happens to cover at 95%. All versions and failed attempts are retained with their source hashes. Any scientific candidate change requires an amendment before its accuracy run.

For runtime and indexing only, an engineering pilot may use one cohort per cell on distinct seed bases `1590000000` and `1600000000`. Its coefficient accuracy is not exported. The pilot includes the alternative-prior fits in the same three fixed cells. It therefore contains 40 analysis requests. No development or confirmation starts before the respective source freeze.

## Prior-sensitivity branch, decided before confirmation

For each cohort with all three prior fits and the relevant analysis returned, define

```
S = ( abs(log(width_0.25 / width_1))
    + abs(log(width_4 / width_1)) ) / 2.
```

Compare `S_direct - S_new_cv3` on the same complete triples, reporting denominators and a paired MCSE in each of the three fixed cells. The overall contrast is the equal-weight mean over those cells; its MCSE is `sqrt(sum(cell_MCSE^2)) / 3`, with the division outside the square root. This compares **relative** sensitivity and cannot be satisfied merely by making every interval narrower. Also report paired estimate shifts and failures by scale; width is not the whole prior-sensitivity story.

The direct-only candidate supports a claim of reduced relative width sensitivity only when the upper two-sided 95% Monte Carlo confidence bound for that overall contrast is below zero and no cell has a positive point contrast. In development, failure of this screen activates a **separate prospective covariance-prior development amendment before confirmation**. The next principled candidate to consider is a hierarchical inverse-Wishart/half-t construction such as Huang and Wand, with its actual joint prior, scale choice, sampler and validation specified in that amendment. It is not silently substituted here, and half-t marginal scales do not imply independent scales and correlations.

Direct inference may be retained as a calibration improvement even if this sensitivity screen fails, but it cannot be presented as a solution to prior sensitivity. Because the user requests improvement in all three areas, the prior-development branch must be investigated rather than declaring the entire task complete from calibration alone. No final confirmation data are used to choose that additional model.

## Stage C: fresh confirmation after candidate freeze

After development and any prior amendment, commit the selected candidate, controls and adoption interpretation. Then run **400 fresh cohorts per original cell** and **200 fresh contextual-X cohorts**: 2,600 cohorts. Use the first 200 cohorts in balanced Gaussian, near-zero ICC and MAR for the alternative prior scales. With the initial candidate set, this gives **12,800 analysis requests**, from 6,400 MCMC fits plus 2,600 complete-truth oracle analyses. If a later prior candidate changes that grid, amend the grid before any confirmation run.

The confirmation data seed is `1630000000 + cell_index*100000 + replicate`; imputation seed is `1640000000 + cell_index*100000 + replicate`. At 95% coverage, 400 cohorts give binomial MCSE about 1.09 percentage points; 200 give about 1.54. These are scenario-specific checks, not a guarantee of universal calibration.

## Separate adoption screens

The numerical thresholds below are prospective research judgments for this synthetic scale, not general statistical laws. Report their uncertainty and every failed screen. Do not change thresholds after seeing confirmation.

1. **Implementation and estimand:** independent original-unit coefficient/contrast tests, a known-covariance Gaussian coefficient oracle, transformation/equivariance checks, exact preservation contracts, source integrity and diagnostic refusal tests must pass. A material error blocks adoption regardless of favorable coverage.
2. **Calibration claim:** use the four primary Gaussian cells (balanced, unbalanced, MAR and null). Conditional coverage of the selected direct procedure must lie between 92.5% and 97.5% in each, and the equal-cell mean absolute distance from 95% must be smaller than both old CV3 and candidate CV3. Report paired Monte Carlo uncertainty for the distance improvements using a prespecified 2,000-replicate bootstrap of cohorts within each fixed cell, resampling each cohort's methods together. A claim of demonstrated improvement requires its 95% percentile interval to exclude zero in the improvement direction. Absolute bias must not exceed 0.02 in any primary cell; report bias MCSE rather than treating that screen as proof of zero bias. Near-zero ICC and t5 results remain prominently reported even though outside this primary gate.
3. **Prior sensitivity claim:** the complete-triple screen specified above must pass on fresh confirmation. Report per-scale returned/per-request coverage, relative and absolute width changes, estimate shifts and every refusal. If a later hierarchical prior is selected, its predeclared sensitivity perturbations and comparator belong in its amendment.
4. **Convergence/availability claim:** thresholds are unchanged. On each of the six original cells the candidate fit must return at least 99% of requests and have no greater observed refusal rate than the old profile. The equal-cell paired refusal-rate difference must favor the candidate, with its upper 95% normal Monte Carlo bound below zero, for a claim of demonstrated reduction. Report the contextual stress separately. A selected numerical challenge set cannot satisfy this gate by itself.

The mean 95% interval score, `width + 40*(lower-truth)` when truth is below the interval and `width + 40*(truth-upper)` when above, is a secondary precision/calibration summary on returned intervals and on complete pairs. It cannot override failed calibration, sensitivity or implementation screens. A refusal has no interval score and remains visible in the separate all-request coverage/availability summaries.

No all-purpose GO/NO-GO label is used to hide a failed component. A limited explicit model-based API can remain useful while a broader improvement claim fails. Any failure must narrow release wording; new algorithms or tuning after confirmation require a new plan and new seeds.

## Reporting, pairing and reproducibility

For each cell, procedure and prior scale, record requested and returned intervals, errors, all warnings, coverage per request with refusals counted as no covering interval, coverage conditional on return, Wilson intervals and binomial MCSEs, bias, RMSE, empirical estimate SD, mean reported spread, interval width, interval score and their MCSEs. A direct posterior SD is labeled as posterior spread, not as a robust frequentist SE. Keep df and within/between variance only where a CV3/Rubin procedure actually supplies them. Save direct endpoint MCSEs and all model diagnostics/extension attempts.

When a shared candidate fit fails, both analyses are recorded as refused for that reason. If posterior summarization alone fails its additional original-coefficient checks, retain a valid CV3 result and record the direct refusal separately. Do not silently force identical analysis availability. Pair analysis contrasts on common successful cohorts and also compare coverage/failure indicators on all shared requests. Do not compare a 200-cohort prior subset against the full 400-cohort default result.

The independent verifier checks the exact planned grid and seeds, pairing, raw coefficient/interval/coverage arithmetic, failure denominators, summary MCSEs, prior triples and source hashes. It must reject temporary-copy mutations that remove an entire cell, alter a summary, or flip coverage. Nonfinite diagnostics constitute a failed screen. Development and confirmation use different directories; no output directory is overwritten.

Hash the candidate R sources, `DESCRIPTION`, `NAMESPACE`, this plan and runner, and the exact frozen baseline sources before and after each stage. Record `pan`/R/package versions, controls and worker counts. Parallelism is an execution choice; it cannot alter the seed mapping. Save a completed-run marker only after the planned grid and hashes verify. Interrupted execution is retained and labeled; no accuracy-based restart is permitted.

## Methodological basis

- Xie and Meng (2017), [uncongenial imputation and analysis](https://www3.stat.sinica.edu.tw/sstest/j27n4/j27n41/j27n41-10.html), motivate separating model-posterior inference from an independently chosen completed-data analysis.
- Barnard, McCulloch and Meng (2000), [covariance separation](https://www.rob-mcculloch.org/some_papers_and_talks/papers/published/barnard-meng-mcculloch.pdf), and Gelman (2006), [variance priors](https://sites.stat.columbia.edu/gelman/bayescomputation/Gelman2006.pdf), motivate treating covariance assumptions explicitly rather than labeling a convenient prior noninformative.
- Huang and Wand (2013), [marginally weakly informative covariance priors](https://www.matt-p-wand.net/publicns/Huang13.pdf), is a possible separately specified development branch, not the initial candidate in this plan.
- Bartlett and Hughes (2020), [bootstrap inference with MI under uncongeniality](https://arxiv.org/abs/1911.09980), provides a possible future route for retaining a frequentist OLS target; cluster-level resampling and its substantial computation would need a separate protocol.

## Amendments

None yet. No candidate accuracy runs preceded this plan. Before accuracy runs, the cohort-bootstrap seeds were fixed at `1650000000` for development and `1660000000` for confirmation. The engineering pilot produces no accuracy bootstrap.

Before accuracy runs, bootstrap ordering was fixed explicitly: use R `RNGkind("Mersenne-Twister", "Inversion", "Rejection")`, the stage-specific seed above, then loop over bootstrap replicates first and `primary_cells` in the listed order second. Draw `sample.int(n, n, replace = TRUE)` for each cell and reuse those cohort indices across methods. An independent Python verifier may use base R solely to reproduce these indices, then recompute all statistics itself; it does not call the runner's scoring functions.
