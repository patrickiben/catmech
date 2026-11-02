# Test of what goes wrong with rare counts: results and interpretation

**Rarity is not the problem. Under the [rules fixed in advance](RARE_COUNTS_PLAN.md), the shortfall comes from heavy-tailed cluster effects.** With every row recounted and no cluster effect, rare counts cover 93.67% on average under `cv3`. With a gamma cluster effect of shape 0.5 they cover 88.30%. Counts ten times as common lose just as much, going from 93.77% to 88.07%. So "rare counts are unsolved" was the wrong description. A better one is that cluster-robust intervals fall short when a few clusters are far from the rest, whatever the count level. Real rare items here do not show the shortfall, which fits: `lbw`, `resus` and `stillbirths` cover 94.0%, 96.7% and 95.1% with the option.

**Nothing I tried fixed it.** No other reference distribution, no count-aware imputation and no warning met its rule. No package code changes. Separately, one [independently written setup](INDEPENDENT_SETUPS_SYNTHESIS.md) did reach 94.17% on the six original cells, with wider intervals and after tuning on those same cells, which is described there.

There were 18,000 recount requests, 21,000 imputation requests with 5 giving no interval, and 20,500 diagnostic requests. The study took about 5 minutes. [Scored output](rare-counts-decisions.txt) and [post-hoc numbers](independent-setups-posthoc.txt) come from separate scripts.

## P1: what the shortfall is

Mean coverage under `cv3` with every row recounted:

| Cluster effect | Rare level | Ten times as common |
|---|---|---|
| None | 93.67 | 93.77 |
| Gamma, shape 2 | 92.10 | 91.97 |
| Gamma, shape 0.5 | 88.30 | 88.07 |

The rule's condition for **HEAVY TAILS** was at least 93 with no cluster effect and at most 90 at shape 0.5, and both hold. How heavy real cluster effects are remains unknown. These simulated cells may be more extreme than most audits, though the [second health system](SECOND_SYSTEM_SYNTHESIS.md) shows that single districts can dominate a slope in real data.

## P2: another analysis

| Reference distribution | Rare, shape 0.5 (needed 92) | Other 15 cells (`cv3` 91.91) | Eligible |
|---|---|---|---|
| `cr2` | 88.63 | 91.92 | No |
| `wcr` | 90.50 | 93.14 | No |
| `boot_t` | 86.50 | 89.19 | No |
| `poisson` | 86.00 | 80.80 | No |

The restricted wild bootstrap comes closest, and it does best overall, but it is short of the rule and has no route through imputation here. The bootstrap-t did worse than `cv3`, against my expectation. The Poisson model answers a different question, and under these linear means it is misspecified, which contributes to its low coverage at the common level.

## P3: count-aware imputation

| Candidate | Distance from 95% (linear MI 4.38, option 3.73, recount 4.81) | RMSE over linear MI | At most 98%, of 12 | Eligible |
|---|---|---|---|---|
| `nb_mi` | 4.07 | Not meaningful | 4 | No |
| `pmm_recorded` | 9.25 | 0.855 | 12 | No |

- **`pmm_recorded` is more accurate and undercovers.** It has about 15% lower RMSE but much worse coverage, with intervals about 17% narrower. That matches the earlier finding for this kind of matching.
- **`nb_mi` was not a fair test of count imputation.** My recipe blew up: errors reached the order of 10^11 in some requests, it fell back to Poisson on 71.0% of requests, and even median error was above linear MI's in every cell. The log link and normal coefficient draws are unstable when strata have few events. A count model with partial pooling of strata might behave. This one says only that a plain recipe is not safe.

## P4: a warning

**Not adopted.** The diagnostic, one cluster supplying over half of the variance, pointed the wrong way: flagged requests covered 97.0% and unflagged 88.8%. It also flagged 12.16% of requests where nothing was wrong, against a limit of 5%.

My explanation, which is post hoc and untested, is that an influential cluster in the sample widens the interval, so these flagged intervals are the safe ones. The failures come when such a cluster is missing from the sample, and nothing inside the sample can show that. The [second health system](SECOND_SYSTEM_SYNTHESIS.md) showed the same pattern directly, with coverage of 100.0% with one district and 0.0% without it. So I know of no simple warning for this.

## Predictions, scored

| # | Prediction | Result |
|---|---|---|
| 1 | P1 is HEAVY TAILS | **Held** |
| 2 | No method in P2 is eligible, with `boot_t` the closest | **Failed**: none eligible, but `wcr` was closest and `boot_t` last among the linear ones |
| 3 | Neither candidate in P3 is eligible | **Held** |
| 4 | The warning in P4 is adopted | **Failed** |

## What this changes

- **The documents now say** that rare counts themselves are not the problem, and that heavy-tailed cluster effects are, with no fix found in the package.
- **No package code changes** under these rules.

## Limits

These cells use one covariate pattern and three rows per cluster, and the gamma cluster effect is an assumption. Only three real rare items were available, all from one recount and all seen before. The `D` threshold was a judgment call. The same assistant wrote the plan, runner, scoring and this write-up, and no reviewer has read it.
