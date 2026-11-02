# Comparison with outside methods: results and interpretation

Under the [rule fixed in advance](OUTSIDE_METHODS_PLAN.md) the outcome is **KEEP**. No outside method was no worse than `catmech_reference(cluster = )`, so the option stays, and stays experimental. No package code changed.

**The margin is narrower than the rule's table first suggests.** Much of the gap to `mice` came from one `mice` default, which drops the recorded value when it agrees closely with the audited value. With that check off, in a [rerun made after seeing results](OUTSIDE_METHODS_SENSITIVITY_PLAN.md), `mice_2l_lmer` matched the option on accuracy and failed only one of six conditions. The same model is available in maintained software, and the option's lead over it is small. Its intervals are 1% to 11% wider, and they cover closer to 95% on the real items. It is also about five times faster, and it needs no such setting. That is the case for keeping it, and no more than that.

The first run had 297,000 requests, all returning intervals but 1,016, and took about 49 minutes. The rerun replaced 58,356 of those. [Scored output](outside-methods-decisions.txt), [the rerun's score](outside-methods-sensitivity-decisions.txt) and [post-hoc numbers](outside-methods-posthoc.txt) come from separate scripts. Results from the first run decide the rule. Everything from the rerun is post hoc and changes no outcome.

## Decision rule

Ratios are RMSE against the option's, as a geometric mean across items or cells. Coverage counts use real items at 25% unless stated. The option reaches 92% on 9 of 11 items, with a mean distance from 95% of 1.47 points. At 50% it reaches 10 of 11, at 1.25 points, and 30 of 36 simulated cells are within limits.

| Method | Accuracy, 25% (N1, at most 1.05) | Accuracy, 50% (N2, at most 1.05) | At 92% or more, distance (N3) | At 50% (N4) | Simulations, ratio and cells (N5) | Available (N6) | No worse |
|---|---|---|---|---|---|---|---|
| `audited_only` | **1.618** | **1.349** | 10, 1.65 | Pass | **1.428**, 34 | Pass | No |
| `difference` | **1.199** | **1.136** | **0, 15.86** | **Fail** | **1.083, 2** | Pass | No |
| `raking` | **1.386** | **1.125** | **0, 25.06** | **Fail** | **1.130, 0** | Pass | No |
| `mice_2l_lmer` | **1.348** | **1.121** | **4, 7.16** | **Fail** | 1.015, 30 | Pass | No |
| `mice_2l_pan` | **1.353** | **1.127** | **6, 4.74** | **Fail** | 1.016, 33 | Pass | No |
| *Rerun:* `raking`, `approx` | **1.386** | **1.125** | **0, 13.70** | **Fail** | **1.130, 4** | Pass | No |
| *Rerun:* `mice_2l_lmer`, `eps = 0` | 0.998 | 1.003 | **8, 2.55** | Pass | 0.998, 30 | Pass | No |
| *Rerun:* `mice_2l_pan`, `eps = 0` | 1.036 | 1.016 | **8, 2.60** | **Fail** | 1.000, 33 | Pass | No |

Bold marks a failed condition. N3 also needs at most 98% on at least 10 of 11 items, which `mice_2l_pan` misses in both runs at 9, with `stillbirths` at 99.9%.

## The `mice` default

By default, `mice` removes a predictor whose correlation with the imputed variable is at least 0.99 among observed rows. In an audit, that predictor is the recorded value, and it is the one that matters most. It happened in 4,424 of 49,500 requests, nearly all on items where the records are good. On these cohorts at 25%:

| Item | Requests of 2,000 | RMSE, default | RMSE, `eps = 0` | Coverage, default | Coverage, `eps = 0` |
|---|---|---|---|---|---|
| `deliveries` | 1,571 | 0.2545 | 0.0657 | 58.2 | 94.0 |
| `livebirths` | 951 | 0.2479 | 0.0638 | 57.5 | 95.0 |
| `ancnewreg` | 868 | 0.3640 | 0.1158 | 79.5 | 91.5 |
| `pnc` | 158 | 0.2807 | 0.0600 | 74.7 | 78.5 |

It is not always harmful. On `iron` and `stillbirths`, on the smaller groups of cohorts it touched, the default did better on coverage (86.0 against 75.4, and 96.7 against 87.9). Those are unusual cohorts, and both methods do poorly on them with the check off.

This is the most directly useful finding here for anyone imputing audited values with `mice`: the better the records agree, the likelier the default is to remove them. The only signals are a brief warning that events were logged and an entry in `loggedEvents`. The setting that turns it off, `eps = 0`, cannot be passed straight through with `2l.lmer`, because `mice` hands it on to `lmer`, which rejects it. Then `mice` warns and returns unmodeled draws. The rerun used a small wrapper.

## The option against `mice_2l_lmer` with the check off

These use the same model, and accuracy is the same: the RMSE ratio per item at 25% runs from 0.984 to 1.015. The difference is width. The option's intervals are wider on every item at 25%, by about 1% to 11%, and its coverage is higher on every item, by 0.25 to 2.46 points. The mean distance from 95% is 1.08 points larger for `mice`, with a cohort-resampling range of 0.76 to 1.27 that treats the 11 items as fixed. The condition allowed 0.25.

So wider intervals are all the option has. On these real items, where almost everything undercovers, that helps. In the simulations the two tie at 30 of 36 cells, and both run high in some cells. Wider is not the same as better calibrated in general. Why the option is wider has not been tested, though it draws both variance components from chi-square distributions with few audited clusters, which would be consistent.

Speed: about 0.102 seconds per request for the option, against 0.544 for `mice_2l_lmer` and 0.236 for `mice_2l_pan`, on this machine under load. `lmer` also reported a singular fit in some imputation for 21% of requests in the first run.

`mice_2l_pan` is less stable. It reached 99.9% on `stillbirths`, its width on real items at 25% runs up to +102%, and six of its requests on `iron` hung until the 60-second limit stopped them. These are its defaults, untuned, including its prior.

## Design-based methods and the audited rows alone

- **`audited_only` covers well and is imprecise.** It reaches 92% on more items than the option, at 10, with widths +31% to +205% and an RMSE 1.618 times the option's. The option's real gain over simply analyzing the audit is in accuracy.
- **`raking` and `difference` fail on accuracy, whatever their variance.** Their point estimates do not depend on the variance setting, so these results stand under both runs. With one influence-function auxiliary and about 17 audited clusters, they likely use the recorded value less efficiently than a model does. That reading has not been tested.
- **Their intervals are too narrow, and this part reflects my recipes.** With `raking` under `full`, standard errors ran 0.50 to 0.63 of the spread of the estimate in the shift cells, and with `approx` 0.66 to 0.93. `difference` ran as low as 0.36 in the design with two audited clusters per stratum. That looks like unadjusted leverage: with few audited clusters per stratum, stratum effects shrink the residuals, and these variances have no correction like the jackknife the option gets. It has not been tested. So this comparison does not show those coverage numbers are the best these methods can do. A jackknife might repair them, and it would not change N1.
- **`survey::twophase(method = "full")` did not give consistent results with several rows per cluster.** On one fixture, one total had a standard error of 296 at row level and 589 after aggregating to clusters, against 417 both ways under `approx`. Whether this falls outside the method's intended use has not been settled, and the `survey` authors have not been asked.

## Recorded values alone

At a 25% audit, `recorded` has a lower RMSE than every audit-using method on **6 of 11 items**: `stillbirths`, `deliveries`, `pnc`, `syphilis`, `ancnewreg` and `livebirths`. At 50% that falls to 2 of 11. So a partial audit can add more noise than it removes.

That is no reason to skip the audit, since no one can tell in advance which case applies. On other items the recorded value's RMSE is up to 2.107 times the option's at 25% and 3.722 times at 50%, and its coverage falls as low as 40.6%. It is a reason not to assume that imputation from a small audit improves a good record.

## When the audited value is a predictor

No rule covers this. The earlier finding now has outside support: `mice_2l_lmer` and `mice_2l_pan` attenuate the slope too, at -12.0% and -14.7% with `shift` errors at 25%, against -16.1% for the option and -3.0% for linear MI. So this comes from the random-intercept model, not from catmech's code. It is the same pattern in independent software. The cause has still not been tested directly.

Linear MI without `cluster` stays the advice for this use. `raking` is the notable alternative, at -0.8% to +0.2% bias across the six cells and, under `approx`, 92.4% to 95.0% coverage, with RMSE 1.029 to 1.148 times linear MI's. `audited_only` also has little bias and costs more in accuracy.

## Rare counts

Nothing rescued them. On the rare-count cells, `mice_2l_lmer` with the check off stayed at 83.8% to 91.8%, about where the option is, at 83.2% to 92.0%. `mice_2l_pan` runs higher, from 85.6% to 97.2%, which is the only sign of a difference, and it comes with the instability above. These are normal models applied to counts.

## A correction to an earlier write-up

The [real-audit synthesis](REAL_AUDIT_SYNTHESIS.md) reports a "design-weighted difference estimator" that beat linear MI on recovery RMSE in 15 of 20 cells. That estimator's correction term is zero by the normal equations, as the [check here](outside-methods-check.txt) confirms numerically, so it is least squares on fitted values. It is no difference estimator, and that sentence overstated what was tested. Here, its exploratory form has an RMSE about equal to linear MI's on real items, and higher than the option's by about 5% to 8% ([scored output](outside-methods-decisions.txt)).

## Predictions, scored

| # | Prediction | Result |
|---|---|---|
| 1 | The outcome is RETIRE, through `mice_2l_lmer` | **Failed** |
| 2 | `raking` fails N1 | **Held**, 1.386 |
| 3 | `difference` reaches 92% on at least the option's count at 25%, and fails N1 | **Failed**, it reached none |
| 4 | `audited_only` has an N1 ratio above 1.20 | **Held**, 1.618 |
| 5 | `mice_2l_pan`'s distance at 25% is within 0.50 points of `mice_2l_lmer`'s | **Failed**, 4.74 against 7.16 |
| 6 | `recorded` has lower RMSE than every audit-using method on at least three items at 25% | **Held**, six |
| Rerun 1 | `mice_2l_lmer` with `eps = 0` has an N1 ratio at most 1.05 | **Held**, 0.998 |
| Rerun 2 | It still fails N3, under nine items | **Held**, eight |
| Rerun 3 | `raking` with `approx` still fails N1, and stays under 92% on at least nine items | **Held**, on all 11 |

## What this changes

- **No package code.** `catmech_reference(cluster = )` stays, and stays experimental. Earlier verdicts stand.
- **Documentation** now says the same model is available in `mice`, warns about that default, and states this comparison.
- **The option's case is narrow:** the same accuracy as `mice_2l_lmer`, somewhat wider intervals, greater speed, and no setting to get right, all on one health system.

## Limits

This is not a neutral comparison, and that matters most here. The same assistant built the option, wrote the plan and chose each comparator's recipe, and no reviewer, AI or human, has read this round. Each recipe is a reasonable default use, and none is the best possible use of its method. Two of them were visibly hurt by a setting, and only those two were rerun, after seeing results that looked wrong, while the option has had several rounds of this attention. So KEEP should be read with that in mind. A statistician who knows `mice` or `survey` might well set these up better.

The option's lead over `mice_2l_lmer` rests on one condition, on real items from one health system, where the same items shaped earlier choices about this option, including the variance it uses. The rule's margins favor the comparators, but the option's author chose them. The analysis variance is the jackknife for both imputation methods and was chosen in an earlier round on the option's results, not `mice`'s. The Zambia file has not been used, so none of this has been tested on a second dataset.

On process: the plan was committed before the runner existed, and the runner and scorer were committed before the study. An engineering pilot was read for statuses and timing only. The control reproduced the option's saved results on 660 requests. The fixed-dataset check matched `survey`'s point estimate and second-phase variance, and matched its first-phase variance apart from the factor G / (G - 1). The rerun took three attempts. One stalled for about 30 minutes because `mice_2l_pan` hung, and one failed its control because of the `eps` problem above. Both are recorded in its plan and neither is used. Commit order in a local repository is not third-party proof of timing.

## Sources

See the [plan](OUTSIDE_METHODS_PLAN.md) for references.
