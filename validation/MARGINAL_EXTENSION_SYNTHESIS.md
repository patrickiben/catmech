# Numerical convergence extension: v0.10.1

Both original convergence refusals pass with one additional doubling to 64,000
retained iterations per chain. The complete numerical replay returned all 12,000
MCMC fits and all 27,600 analysis intervals from 3,600 cohorts. All four original
screens pass, including availability. This is a **dependent replay of the same
confirmation cohorts**, not an independent scientific confirmation.

## Change and selected-case evidence

Only the marginal engine's default maximum changes from 32k to 64k. Four chains,
1,000 burn-in iterations, 4,000 initial retained draws, exact seeded restarts, every
monitor, R-hat <= 1.05 and bulk/tail ESS >= 400 remain fixed. Earlier passing attempts
stop as before; explicit 32k user controls and the conditional joint engine keep
their meaning. The model, priors, support restrictions and intervals are unchanged.

| Original case | Old minimum bulk ESS at 32k | New minimum bulk ESS at 64k | New minimum tail ESS | New maximum R-hat |
|---|---:|---:|---:|---:|
| Unequal-site Gaussian, 12 |396.73|835.62|1768.02|1.00648|
| Contextual null, 143 |387.00|805.47|1697.83|1.00774|

The selected [challenge](marginal-extension-challenge/results.json) matches the
recorded first four attempt histories at their stored precision and the 32k
extrema at full precision. It assesses numerical performance, not coverage.

## Full replay and unchanged-result audit

The [plan](MARGINAL_EXTENSION_PLAN.md) was frozen before replay. The final run
uses independently started PSOCK R workers and the same complete seed grid,
prior alternatives, targets and bootstrap seed as v0.10.0. No interrupted-run
records are reused and no cohort is replaced.

All 5,200 baseline and 6800 candidate MCMC fits returned; all 3600 oracle workflows
also returned. Every 27,600 requested interval is present. The contextual,
ordinary-calibration, slope-sensitivity and availability screens pass. The
70.1% measured slope-sensitivity reduction over the declared prior ranges is
retained. Calibration and scientific limitations remain those of the original
[synthesis](JOINT_MARGINAL_SYNTHESIS.md); no broader robustness claim follows.

The [comparison audit](marginal-extension-revalidation-isolated/comparison-with-v0.10.0.json)
checks every field of all 27,596 formerly successful analyses, excluding timing
and source identifiers, at relative tolerance 1e-10 / absolute 1e-11 for numbers.
Their estimates, intervals, uncertainties and diagnostics are unchanged. The
only recovered analyses are direct and CV3 summaries for the two original cases.
[Independent grid/diagnostic/arithmetic and tamper audit](marginal-extension-revalidation-isolated/verification.json).
[Full numerical replay](marginal-extension-revalidation-isolated/summary.csv).

## Interrupted execution and process isolation

The initial forked replay suffered a worker segfault at observed-MAR cohort 34.
It was stopped and preserved as [INCOMPLETE](marginal-extension-revalidation/INCOMPLETE.json),
with its log, checkpoints and unchanged-source hashes. Its missing cohort is
not silently repaired, and it is not a completed study. The crash cause is not
established by the log.

The execution amendment starts independent PSOCK workers and immediately rejects
missing or malformed worker results. The affected cohort passed the selected
[isolated check](marginal-extension-isolated-challenge-final/CHALLENGE.json), with
all 17 estimates/intervals matching the original. The final entire grid was then
run from scratch. The exact 35-file source tree is preserved in
[marginal-extension-isolated-source](marginal-extension-isolated-source/SOURCE.json).

## Package and CI checks

All 1,526 package assertions pass on Ubuntu, Windows and macOS CI. The macOS job
uses the explicit macos-15 image, retaining the same full CRAN-style checks;
macos-latest had remained queued. [Successful three-platform run](https://github.com/patrickiben/catmech/actions/runs/41872451153).
Local source archive and installed-package checks are recorded in
[v0.10.1-checks](v0.10.1-checks/).

Zero observed failures on this reused grid do not guarantee convergence for new
inputs. Unsupported data and failed diagnostic screens must still refuse.

```sh
Rscript validation/run-marginal-extension.R confirmation validation/my-extension
python3 -I tools/verify-marginal-extension.py validation/my-extension --tamper-tests
python3 -I tools/compare-marginal-extension.py validation/joint-marginal-confirmation validation/my-extension
```
