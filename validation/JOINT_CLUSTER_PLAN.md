# Prospective validation: joint recorded/reference cluster model

Written 2026-11-03 before the engineering pilot or scored runs. This evaluates an explicit opt-in joint model; it does not replace the historical cluster studies or change their results. Commit this plan before running the pilot. Record any engineering amendment below before inspecting held-out accuracy.

## Question and scope

Can a bivariate Gaussian random-intercept imputer share cluster uncertainty coherently when the same cluster contains unaudited rows both with and without a recorded value? The existing separate-fit implementation refuses that pattern. The candidate models reference truth `T` and recorded value `W` jointly, conditional on complete `x`, with a bivariate site effect and bivariate row error. It uses the established `pan` engine, flat fixed-effect prior and declared proper inverse-Wishart covariance priors on standardized responses. This is established multilevel multiple imputation, not a new statistical method.

The primary estimand is the population coefficient of `x` in the linear outcome analysis `T ~ x`. We do not infer validity when the repaired quantity is an analysis predictor, or for counts, nonignorable missingness, imperfect references, contextual effects, or a misspecified mean. A coherent imputation posterior does not itself guarantee calibration of Rubin pooling with a cluster jackknife.

## Fixed data-generating process

Each cohort has 24 sites. Balanced cells have 10 rows per site. The unbalanced cell repeats site sizes 5, 8, 12 and 15 six times, also totaling 240 rows. Generate independent standard-normal `x`, then

```
(T, W) = (0.5 + beta*x, 0.8 + 1.05*x) + site_effect + row_error
Omega  = [[0.50, 0.45], [0.45, 0.60]]
Sigma  = [[0.50, 0.55], [0.55, 0.80]]
```

Effects are independent between sites and rows. They are bivariate Gaussian except in the explicit stress cell. The primary slope is `beta = 0.7`; the null cell uses zero. Each replicate randomly chooses six wholly unaudited sites. Gold measurements are exact, of the same quantity and time as `T`. Only selected `T` values enter the public reference column; all others are `NA`. Hidden truth remains outside imputer input and enters only the oracle and scoring.

In the ordinary cells, randomly choose exactly one quarter of all rows for reference measurement from the remaining 18 sites. Within every site, randomly hide `ceiling(0.3 * site_size)` recorded values. Thus every site has both recorded-value patterns, and every wholly unaudited site requires joint stochastic imputations for both. Missingness and audit choices are independent of latent values.

Six cells are fixed:

| Cell | Change from the ordinary balanced Gaussian cohort | Role |
| --- | --- | --- |
| `gaussian_balanced` | None | Primary matched model |
| `gaussian_unbalanced` | Site sizes specified above | Primary matched model |
| `near_zero_icc` | Multiply `Omega` by 0.02, leave `Sigma` fixed | Weak variance-component stress within the model |
| `observed_mar` | Missingness and reference selection described below | Primary ignorable selection |
| `null_slope` | Set `beta = 0` | Primary null calibration |
| `student_t5` | Independent multivariate t with 5 degrees of freedom for site effects and row errors, each scaled to preserve the corresponding covariance above | Distributional misspecification stress |

For `observed_mar`, recorded-value missingness is Bernoulli with log odds `a + 0.8*x`; solve the cohort-specific intercept so the mean probability is 0.30. No pattern is forced afterwards. Among the 18 eligible sites, reference selection is Bernoulli with log odds `b + 0.6*x + 0.6*w_observed`, where `w_observed` is the centered, scaled recorded value when available and zero otherwise. Solve a separate cohort-specific intercept so the expected total audit is one quarter of the cohort. The observed-data missingness mechanism is ignorable for the joint model: unobserved `W` and unobserved `T` do not enter those selection probabilities. Randomly choosing the unaudited sites also does not depend on their latent effects. Bernoulli draws are retained as drawn, even when a fit refuses insufficient support.

The t cell changes the conditional distribution as well as marginal tails. It is not a claim to isolate a single causal reason for any coverage difference.

## Methods and information matching

All imputed analyses use 40 completions, the same cohort and reference mask, and `catmech_pool_lm(..., cluster = "site", cluster_variance = "cv3")`. Target the `x` coefficient at nominal 95% confidence.

- `joint`: `catmech_reference` with normal method, recorded auxiliary, site and explicit `cluster_model = "joint"`. Default prior scale is one. The response covariance priors are inverse-Wishart with 4 degrees of freedom and scale `0.5 * I * prior_scale`, after the engine's recorded standardization. The fixed-effect prior is flat; support checks are part of the tested procedure.
- `plain_recorded`: the existing recorded-value regression without an imputation cluster effect, followed by the same cluster-robust pooled analysis. This is an existing supported baseline that accepts both recorded patterns and uses the same gold, `W` and `x`; it is not claimed to be a correct joint model.
- `oracle`: truth in every row, followed by the same completed-data cluster jackknife. It measures the analysis procedure's own finite-sample limitations and uses extra information unavailable to the imputers. Use the no-imputation branch, with two identical completions only to access the existing analysis helper.
- `joint_prior_quarter` and `joint_prior_four`: the same joint engine with prior scales 0.25 and 4. Run these on the first 100 replicates of `gaussian_balanced`, `near_zero_icc` and `observed_mar`, selected in advance. Report paired differences on exactly that subset. This is a prior-sensitivity check, not a search for a winning prior or a replacement default.

The legacy clustered separate-fit method is not an accuracy comparator because the primary cohorts deliberately trigger its supported-pattern refusal. Add a small explicit contract check that it refuses a mixed stochastic site; do not count that deliberate refusal as a joint-model failure.

## Monte Carlo size, seeds and execution

Run 200 replicates of every cell for the three main methods: 3,600 requests. The prespecified prior checks add 600 requests, for 4,200 scored requests total. At 95% coverage, 200 replicates have approximate Monte Carlo standard error 0.0154; 100 have 0.0218. This is a bounded initial validation, not evidence of one-percentage-point calibration.

The held-out data seed is `1530000000 + cell_index*100000 + replicate`. The held-out imputation seed is `1540000000 + cell_index*100000 + replicate`, shared across the imputation methods. Cell indices follow the table order, starting at one. Seeds are distinct from historical studies. Pilot data and imputation seeds use bases `1510000000` and `1520000000`, respectively. The pilot uses replicates 1 through 3 per cell and is read only for timing, refusals, diagnostic behavior, and engineering defects; do not score or use its coefficient accuracy to set the model or select settings.

The fixed package controls are four independent chains, 1,000 burn-in iterations and initially 2,000 retained iterations per chain. The package may automatically extend the retained diagnostic sample according to the same fixed schedule, 2,000 to 4,000 to 8,000, stopping when its declared R-hat/bulk-ESS/tail-ESS screen passes. At the maximum, a failed screen produces a failed request. Record each request's attempts and retained iterations. Completions are spread over retained iterations and chains; they are not claimed independent. No manual retry, discretionary extension, change of priors, or method-specific exclusion is permitted in the scored run. The engineering pilot may justify a documented global control amendment before held-out scoring; it cannot use coefficient accuracy to select settings. Any follow-up study must use a new plan and seeds.

Hash the plan, runner, installed `pan` version, all production R sources, `DESCRIPTION`, and `NAMESPACE` before and after. Each worker records its source hashes; require exact agreement. Save session information, controls, raw per-request records and errors. Parallel processes may divide whole replicates; output order is restored by cell, replicate and method. No unselected truth is supplied as an auxiliary variable.

## Outputs and interpretation

For every cell and method report requests, returned intervals, failures by error class/message, coverage per request (a failed interval cannot cover), conditional coverage among returned intervals, Wilson Monte Carlo intervals for coverage and failure rates, and their binomial Monte Carlo SEs. Zero observed failures do not imply a zero underlying failure probability.

For returned intervals report bias and its Monte Carlo SE, RMSE and its delta-method Monte Carlo SE, empirical estimate SD, mean model SE, mean interval width and its Monte Carlo SE, plus median/max R-hat and minimum bulk/tail ESS where available. These accuracy summaries are conditional on return and must be labeled if failures occur. Include paired mean error, squared-error, coverage-indicator and width differences against `plain_recorded` on complete pairs, with paired Monte Carlo SEs and pair counts. Report per-request paired coverage differences over all shared requests as well. Prior-sensitivity comparisons use the fixed shared 100-replicate subset, without borrowing precision from other cells.

No post-hoc pass/fail boundary or method selection is introduced. A material numerical error or incoherent shared effects blocks release regardless of coverage. Correct implementation may still be released as an opt-in working model with explicit limitations when stress cells show failure. Undercoverage in the matched cells must be investigated and reported before any claim of calibrated inference; a pilot or a 200-replicate study cannot establish general validity.

## Separate implementation validation required before release

These checks complement the coverage study and do not select the estimator:

1. A fixed-parameter conditional-normal oracle on small mixed-pattern data must check imputed means, variances and cross-pattern within-site covariances against direct multivariate-normal conditioning. Two independently sampled site effects must fail that oracle.
2. Check the actual `pan` adapter against a direct `pan` call with identical standardized data, prior, seed and iteration schedule, including missing `T` only, missing `W` only and both missing. This tests wrapping and indexing; it is not independent validation of `pan`'s mathematics.
3. Preserve gold anchors exactly, original recorded data, caller RNG, row/site alignment, reproducibility and the existing supported default behavior. Include shuffled rows, uneven sizes, all-gold bypass, no-gold/support refusals, and nonconvergence handling.
4. Check the prior parameterization and transformations independently. Sensitivity to covariance prior scale is reported, especially near zero ICC; a proper covariance prior does not establish likelihood identification or remove small-site prior influence.
5. Check that existing complete-W/default requests remain numerically unchanged and that requesting the new model is explicit. The existing warning about using the repaired value as a clustered analysis predictor remains applicable until separately studied.

## Amendments before held-out scoring

Before any pilot or held-out run, review clarified that the baseline is an existing supported comparator, without a claim that it is the strongest possible method. The fixed automatic extension schedule above is part of the tested package procedure; manual retries remain prohibited. Any later engineering changes will be recorded here before the held-out run.

An engineering reproducibility check before the planned pilot found that `pan` retains a hidden random-number cache across calls: an unrelated prior `pan` call could change a nominally same-seed package result, despite restoration of R's RNG state. The adapter will isolate each chain in a fresh R worker process so the recorded chain seed determines its run independently of caller history. This changes execution isolation, not the statistical model, priors, masks, seeds or scored design. The pilot will assess the process-startup cost. No planned-simulation accuracy had been inspected when this amendment was written.
