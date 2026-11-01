# Go/no-go plan: audit allocation under interference

Written 2026-11-01 before any evaluation run. This tests one idea that would combine three
manuscripts: in a trial with spillover, a recording error at one site enters the exposure of other
sites, so audits should follow the direction of spread. It decides whether a new paper is worth
pursuing. It is not that paper. The two manuscripts under review are not changed.

Terminology: "spillover" means treatment contamination between arms, as in the CRT manuscript.
"Recording error" means a wrong recorded value, as in catmech.

## Sources of each component

- **CRT manuscript** (demo commit 04ee03b): the 57 public Malaria Atlas Project sites, gravity kernel,
  exposure `s = 1 - exp(-inflow / kappa)`, data-generating model, effect ln(0.5), contamination strength
  c = 0.35, stratified random allocation, and the spillover-adjusted estimator
  (arm, stratum, baseline logit, s, arm by s, dspill).
- **Directed flow:** the manuscript's own directed form, production `p * mass` times attraction
  `(1 - p) * mass` over squared distance.
- **Contagion manuscript:** transmitter against reactive targeting, the symmetrization blend
  `(1 - d) * (K + K') / 2 + d * K`, and fixing predictions in advance. Its VAR is not used, because
  these sites have no long time series.
- **catmech:** the recording-error model, per-site expected error and pooled inference.

Trial assignments and outcomes are simulated. The recorded-mass errors are simulated too. Only the
geography and survey counts are real, and modeled flows are not measured movement.

## Why mass

Recorded site mass is the input studied. In this kernel a site's own mass cancels from its own exposure,
because its inflow and normalizer both scale with it. A mass error therefore acts only through other sites.
The analysis kernel is assumed known in form and direction, as in the manuscript's main tables, and
prevalence, coordinates and binomial denominators are assumed correct.

## Design

Recorded log mass is the true log mass plus an error. Regimes: `moderate` (20% of sites, error Normal(0, 1))
and `gross` (10% of sites, plus or minus log 10). Gross errors are easy to spot, so this control favors
reactive auditing. Cells: moderate at direction 0, 0.5 and 1, and gross at direction 1. Each has 400 replicates.
Allocation, errors and the first audit wave are shared across policies and directions.

Each replicate audits 8 random sites, fits the recording model on those references, audits 8 more by policy,
and analyzes. All scores except random use the same expected recorded error from the first-wave fit:

| Policy | Second-wave score |
|---|---|
| `random` | Uniform |
| `reactive` | Expected recorded error |
| `undirected` | Expected error times exposure leverage on the symmetrized kernel |
| `transmitter` | Expected error times exposure leverage on the directed kernel |
| `influence` | Expected error times the change in the treatment estimate per unit of log mass |
| `oracle` | True absolute error times directed leverage, which is infeasible and bounds the gain |

Baselines: `recorded_only`, `first_wave_only` and `full_audit`. `influence` represents influence-function
sampling (Shepherd et al., 2023) and is the strongest feasible competitor. It looks at outcomes, which the
recording model does not include, so its selection is only approximately ignorable there.

Analyses: `pooled` runs the adjusted estimator on 20 completions from the recording model and uses Rubin's rules.
`plugin` uses audited masses where known and recorded masses elsewhere. The verdict uses `pooled`.
The recording model starts with `iterations = 1200`, `max_iterations = 4800`, because an engineering run found 17%
first-wave convergence refusals with shorter chains. Its predictors are elevation, latitude and baseline logit.

Seeds: replicate base 400,100,000 with offsets of 1,000,000, 2,000,000 and 3,000,000 for outcomes, imputation
and policy ties. Engineering runs used base 500,100,000 and only their statuses, error counts and run time were
examined. A 6-replicate run printed a verdict, which was ignored because the sample is far too small.

## Port check

Before any evaluation, the port reran the manuscript's symmetric kernel for 1,000 replicates under stratified random
allocation. Deposited values are for the manuscript's own run and allocation pool:

| Estimator | Port bias | Deposited bias | Port RMSE | Deposited RMSE | Port coverage | Deposited coverage |
|---|---:|---:|---:|---:|---:|---:|
| Naive | 0.289 | 0.289 | 0.376 | 0.372 | 0.156 | 0.168 |
| Adjusted | -0.004 | 0.005 | 0.158 | 0.149 | 0.946 | 0.963 |

The port is close but not exact. Its adjusted RMSE is about 6% higher, and its coverage about 1.7 points lower.

## Verdict rule, applied by the script to the moderate regime

1. **Audits matter:** at direction 1, `full_audit` MSE at most 0.9 times `recorded_only`, paired z at most -2.5.
2. **Network helps:** at direction 1, `transmitter` MSE at most 0.9 times `reactive` with z at most -2.5, and below `random`
   with z at most -2.5.
3. **Direction helps:** at direction 1, `transmitter` below `undirected` with z at most -2, and that gap larger than at
   direction 0.5.

**GO** needs 1, 2 and 3. **PARTIAL** needs 1 and 2 only: leverage matters but direction does not, so the contagion
component would be dropped. **NO-GO** otherwise.

## Predictions

- At direction 0 `transmitter` and `undirected` should be identical by construction. Any difference is a code error.
- `influence` may beat `transmitter`. If it does clearly, the transmitter framing adds little to existing
  influence-function methods, and that must be reported even under a GO.
- In the gross control, `reactive` should come close to `oracle`, shrinking any network advantage.
- My prior before running is skeptical, roughly one in three for GO, because errors in a randomized trial's covariates
  mostly add noise, and the adjusted estimator may absorb much of the exposure error.

No retuning follows the evaluation run. A change requires new seeds and a new plan. A GO would justify a literature
review and a fuller study, not a claim of a validated method.
