# Prospective test: comparators set up by someone else

Written 2026-11-02 and committed before the setups exist, so version history records the order. Amendments, if any, go in the section at the end.

## Why

In the [comparison with outside methods](OUTSIDE_METHODS_SYNTHESIS.md), the option's author chose every comparator's recipe, and two recipes were visibly hurt by a setting. That limits how much KEEP means. Here, three fresh AI agents each write their best setup of one family of methods, seeing no results and none of the earlier recipes. This reduces that problem and does not remove it: they are still AI, they may share blind spots with the author, and the author writes their brief and this rule. Review by a statistician who knows these packages remains the real fix.

## How the setups are written

Each agent works from a sealed folder outside the repository, holding a brief, a data generator with the study's six simulated error types, 12 fixture cohorts with their truth, and a harness that checks the required interface. Each is told not to open the repository. That is an instruction and is not enforced by the system, and each agent reports whether it kept to it. They may run their own simulations from the generator. The real item pools are not in the folder.

- `ind_mice`: `mice` multilevel imputation, or any other installed imputation package, with its own choice of method and settings.
- `ind_survey`: `survey`, two-phase calibration or raking, with its own choice of auxiliaries and variance.
- `ind_design`: a design-based correction of the recorded analysis, such as a difference or regression estimator, with its own variance.

**What they are told.** The audit design and analysis, that a coefficient and a 95% interval are required, and three known pitfalls, so that the comparators get the benefit of what was learned: `mice` can silently drop the recorded value through its collinearity check; `twophase(method = "full")` with clusters at phase one and no phase-one information gives a phase-one variance of zero; and cluster-robust variances on few audited clusters can be far too small through leverage. They are not told the option's or any method's results, or the rule below.

Each setup is a function `fit_method(p, seed)` returning an estimate, standard error and degrees of freedom. Its code enters the repository as written and is committed before the study. The author may add only a harness around it. If a setup fails the pilot outright, its agent gets one chance to repair it, using statuses only, and that goes in the amendments.

## Data and rule

The same 49,500 cohorts as the [comparison](OUTSIDE_METHODS_PLAN.md), with its saved `oracle`, `linear_mi` and `random_intercept` results under `cv3`. Each request has a 60-second limit, and one over it counts as no interval. Imputation seeds are 1,260,000,000 for `ind_mice`, 1,270,000,000 for `ind_survey` and 1,280,000,000 for `ind_design`, plus cell number x 100,000 plus replicate, numbered as in that comparison.

**The rule is that comparison's rule, unchanged:** conditions N1 to N6 against `random_intercept`, with the same thresholds, scored by the same scorer. A setup is no worse if it meets all six.

- **RETIRE** if any setup is no worse. `catmech_reference(cluster = )` then gets a deprecation warning in the next version, naming the passing method with the smallest N1 ratio.
- **KEEP** otherwise. The option stays experimental.

**Description only:** the six predictor cells, the 11 real items by method, seconds per request, and whole-cluster audits on the [second health system](SECOND_SYSTEM_PLAN.md), `district_33` and `district_50`, with the setups that accept them, on that study's cohorts and seeds.

## Predictions fixed in advance

1. The outcome is KEEP.
2. `ind_mice` has an N1 ratio at most 1.05, and fails N3.
3. `ind_survey` and `ind_design` fail N1.
4. `ind_survey` reaches 92% on more real items at 25% than the earlier `raking` recipe under `approx`, which reached none.

## Limits known in advance

Three agents are not a survey of practice. They share a model with the author. They can read the same generator that this study's simulated cells use, so they could tune to it, as a diligent analyst would. The real items are held back from them, but those items have been seen by the author. No human has reviewed this.

## Amendments

Made after the study ran and before any of its results were read. The original plan and scorer are kept in `independent-setups-source/`.

**A fault in the scorer's control.** The control compares this scorer's values for `mice_2l_lmer` with the earlier study's saved ones. It searched the first mention of that method in the saved file, which is in the status list and not in the decision section, so it stopped although the values matched. It now searches the decision section. No condition or threshold changed.

