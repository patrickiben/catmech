# Test of separate within-cluster and between-cluster slopes: results and interpretation

Under the [rule fixed in advance](WITHIN_BETWEEN_PLAN.md) the outcome is **NOT ADOPTED**. `catmech_reference(cluster = )` stays as it is and stays experimental, and no package code changes. The thresholds did not move after the results were read.

**The plain reading is that the slopes pay off when the audited value is a predictor and cost accuracy when it is the outcome and the audit is small.** With the audited value as a predictor, coverage moved closer to 95% in every source, and bias fell in both deciding sources. With it as the outcome, coverage was mostly as good or better, but point error rose, by the most on small real audits: 1.258 times the option's on the Rwanda items at 25%, and 1.372 and 1.590 on two of the four Zambia audits. The candidate adds parameters and a small audit has few clusters to estimate them, which fits this pattern, though the test did not isolate that as the cause.

**Separate slopes do not by themselves reproduce `ind_mice`.** That setup, which prompted this test, also has separate slopes, and it had lower error than the candidate on all 11 real items at 25% (on `iron`, 1.458 against 2.092 times the option's). It differs in other ways, including a shrinkage of the imputation variance, and this test does not say which difference matters.

There were 529,500 requests, and the study took about 2 hours. [Scored output](within-between-decisions.txt), [description](within-between-posthoc.txt) and [a later check](within-between-shared-cells.txt) come from separate scripts. The last was written after the results were read, and neither it nor the description enters the rule.

## Decision rule

Each class against the option, on the two deciding sources:

| Source | Class | Cells | Ratio | Distance, candidate and option | Width | Bias %, candidate and option | No worse | Better |
|---|---|---|---|---|---|---|---|---|
| SNAP | Outcome | 8 | 1.003 | 2.60, 2.49 | 1.004 | 0.97, 0.94 | Yes | No |
| SNAP | Predictor | 4 | 0.993 | 9.95, 14.05 | 1.093 | 1.06, 2.52 | Yes | Yes |
| New simulation | Outcome | 66 | **1.026** | 1.78, 2.23 | 1.027 | 0.65, 1.10 | **No** | No |
| New simulation | Predictor | 58 | 0.968 | 1.85, 4.08 | 1.064 | 5.06, 7.62 | Yes | Yes |

Ratio is the candidate's RMSE over the option's, and distance is the mean distance of coverage from 95%. No worse needs a ratio of at most 1.02. SNAP passes. The new simulation fails on the outcome class's ratio, and its cohort-bootstrap interval, 1.021 to 1.030, stays above the limit. Split by audit size, that ratio is 1.048 at 25% and 1.004 at 50%.

The guard, data already seen, could only block adoption, and it did:

| Check | Ratio | Other measures, candidate and option | Result |
|---|---|---|---|
| G1, real items at 25% | **1.258** | At 92% or more on 10 and 9 of 11; distance **1.85** and 1.47 | Fail |
| G2, real items at 50% | **1.059** | At 92% or more on 11 and 10 of 11; distance 1.49 and 1.25 | Fail |
| G3, existing outcome cells | **1.040** | Within limits on 34 and 30 of 36 | Fail |
| G4, existing predictor cells | 0.958 | Bias 5.10% and 6.94% | Pass |
| G5, Zambia `clinic_25` | **1.372** | Distance 7.29 and 7.98 | Fail |
| G5, Zambia `clinic_50` | 1.009 | Distance 10.31 and 10.60 | Pass |
| G5, Zambia `district_33` | **1.590** | Distance 10.00 and 10.12 | Fail |
| G5, Zambia `district_50` | **1.069** | Distance 10.21 and 10.59 | Fail |

Every failed guard check failed on accuracy, and G1 also on distance. The counts of items at 92% or more and of cells within limits improved.

## With the audited value as a predictor

- **Coverage improved.** In the new simulation the candidate reached 92% on 54 of 58 cells, against 44 for the option, 43 for linear MI and 48 for `ind_mice`, which gave no interval on 4 of those cells. In one cell, `p_two_tier_k150` at 25%, the option covered 43.6% and the candidate 95.0%.
- **Bias remained.** The candidate's class bias of 5.06% is well above linear MI's 1.95%. It reached 20.0% on `p_site_dominant` at 25%, where it covered 94.6% and the option's bias was 34.5%. On `covariate_row_only` at 25%, among the existing cells, bias rose from 2.6% to 8.4%.
- **SNAP improved and stayed short.** No method reached 92% on any of the 4 SNAP predictor cells, and neither did the oracle, at 85.5% to 89.3%. The candidate came within about half a point of the oracle in one cell and about 4 points in another.

Linear MI had lower bias than the candidate here, and the candidate had coverage nearer 95%. So the earlier advice stands: do not use the option when the audited value is a predictor. Cluster means added by hand met the predictor checks here but have no rule of their own behind them, so this is not advice to add them.

## With the audited value as the outcome

- **SNAP barely changed,** as predicted: a ratio of 1.003, with most reviews having no error and with errors that hardly cluster (see limits).
- **Real items paid the most.** At 25%, `iron` reached 2.092 times the option's error and `anc1_early` 1.567. At 50%, the ratios fell to 1.126 and 1.190. Coverage moved by less than 3 points on each item.
- **In a few simulated cells it helped outright.** On `o_size_informative` at 25%, the option's bias was 12.4% and it covered 90.6%, where the candidate's bias was 0.4% and it covered 97.6%, with a lower error. In that cell the generator links cluster size to the true value, the covariate and the error.

## Predictions, scored

| # | Prediction | Scored result |
|---|---|---|
| 1 | The outcome is NOT ADOPTED | **Held** |
| 2 | On SNAP, the outcome class has a ratio between 0.98 and 1.02 | **Held**, 1.003 |
| 3 | On the new simulation, the candidate's predictor class has lower bias than the option's | **Held**, 5.06% against 7.62% |
| 4 | On the new simulation, `ind_mice` has a lower predictor-class bias than the candidate | **Printed as failed,** which is a scoring artifact |

Prediction 4 could not be judged as written. `ind_mice` refused on 4 of the 58 predictor cells, so its class mean is undefined, and the scorer counts a comparison with an undefined number as false. The scorer stays frozen, and its line stands in the scored output. On the 54 cells where `ind_mice` gave intervals, in the later check, its bias was 2.10% against the candidate's 3.88%, which is the direction predicted. This is post hoc and is not scored as held.

## Requests with no interval

In the deciding sources 3,003 requests gave no interval. `ind_mice` refused 3,000 of them, on six cells at 25% with too few audited clusters, as the plan's amendment expected. The other 3 are one cohort of `o_offset_x_site`, where linear MI, the option and the candidate all stopped on a design without full rank. In the guard, the candidate's 7 match the option's saved counts on the same items. In Zambia it had 6, which the G5 availability check allowed.

On those same six cells the candidate fell back to the option's model, 3,000 requests in all, which pulls the class ratios toward 1. The guard's failure makes the outcome the same either way. It used different seeds from the option, so those cells' small differences are Monte Carlo noise. In the guard it fell back on both `control_one_row` cells, and in Zambia on 2 requests.

## Limits

SNAP is one paired item with three analyses, in dollars and not in health, and it leaves out cases found ineligible and overissuances that reach the full benefit, which removes some of the largest errors. Its covariates were recorded during the same review. It holds a state reviewer's finding and not a gold standard. The error's intraclass correlation by state was 0.007, so **SNAP tests little of what a clustering method is for**, though 35.5% of its reviews have an error. Its predictor gain may have another source, and this test does not establish one.

The new generator's independence is weaker than planned. Its agent's session loaded an index of the author's notes that named catmech results and methods, though not the candidate or the rule, and it shares a model with the author, who also wrote its brief. Only 5 of its 62 cells have strata. The cells have 500 cohorts each, so single-cell differences of a few points are noise.

The predictor result came from data the rule used, and a version that adds means only for predictor use would be a new hypothesis. SNAP and the new simulation have now been seen, so a fair test of it would need new data. The same assistant wrote the plan, runner, scorer and this write-up. No human has reviewed any of it.

## Sources

See the [plan](WITHIN_BETWEEN_PLAN.md), [the data note](within-between/README.md) and [the generator as written](within-between/generator/generator.R).
