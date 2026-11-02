# Brief

You are a statistician setting up one family of methods for a planned comparison. Do your honest best for your method: the goal is the setup an expert would actually use, not a weak or a flattering one. Decide everything yourself.

## The problem

Each record has a recorded value, which may be wrong. A random audit of whole clusters (such as facilities) finds the true value for audited rows only. The analyst wants a linear regression coefficient on the true values, with a 95% confidence interval, using the recorded values on all rows to gain precision.

## The design

- Whole clusters are audited, sampled without replacement, at the fraction in `p$fraction`, within strata when `stratum` is present, with at least one per stratum. Sampling probabilities are known from the counts.
- A typical cohort has about 70 clusters with 3 rows each, so about 17 or 35 audited clusters, with 7 or 10 strata or none. Some strata have only one audited cluster.
- Errors may be small and row-level, may shift whole clusters, may under-report, or may include outliers. Some targets are rare counts. Often the recorded value agrees exactly with the truth on many or all audited rows, and it can correlate above 0.99 with the truth. In some real data the values are counts with many zeros.
- The target is the coefficient on the true values in the population of clusters that the cohort represents. Real-data versions of these cohorts, which you do not have, have a similar shape.

## What you write

One file, `setups/<your method label>.R`, defining `fit_method(p, seed)`. It may use base R and packages already installed. Do not install anything. Check `installed.packages()`. Fix any randomness from `seed`.

`p` is a list with:
- `input`: a data frame with `lab_value` (recorded), `gold_value` (true on audited rows, `NA` elsewhere), `cluster_id`, `covariate`, sometimes `stratum` (factor), and sometimes `outcome`.
- `formula`: the analysis model. It names the audited quantity as `lab_value`. It is either `lab_value ~ covariate + stratum` (audited value as outcome), `lab_value ~ covariate` or `outcome ~ lab_value + covariate + stratum` (audited value as predictor).
- `term`: the coefficient wanted, `"covariate"` or `"lab_value"`.
- `audited`: logical by row. `fraction`: the audit fraction. `predictors`: columns you may use in an imputation or working model.

It returns `list(estimate = , se = , df = )`. The interval used will be `estimate +/- qt(.975, df) * se`. Call `stop()` if you cannot give a sound answer for a cohort. A refusal counts against availability. Keep typical run time under about one second per call, since it runs tens of thousands of times, and each call has a hard 60-second limit.

## Known pitfalls

These came up in earlier work. You should avoid them, if they apply.
- `mice` by default drops a predictor correlated at 0.99 or more with the imputed variable, which can silently remove the recorded value. `eps = 0` turns this off, but it passes through to `lmer` under `2l.lmer` and breaks it, so a small wrapper is needed.
- `survey::twophase(method = "full")` with clusters at phase one and no phase-one probabilities or population size gives a phase-one variance of exactly zero when clusters have several rows.
- Cluster-robust variances on few audited clusters can be far too small through leverage, particularly with strata that have one or two audited clusters.

## What you have

`generator.R` simulates cohorts (`simulate_cell`, `simulate_covariate_cell`) with six error types. You may run your own simulations with it to check and improve your setup. `fixtures.rds` holds 12 cohorts, including truth, in the shape you will receive, and `check_interface.R` checks your file. The study hands you only the fields listed above.

## Rules

Work only in this folder and its `setups` folder. Do not open, search or read anything under `/Users/patrickiben/Documents`, and do not look for earlier results, other setups or other methods' performance. Do not search the web for this package. You may read installed package documentation.

Report honestly what you chose and why, what your own simulations showed, what you doubt, and whether you kept to these rules.
