# Go/no-go result: audit allocation under interference

**Verdict: NO-GO, by the rule [fixed in advance](PLAN.md).** The idea was that a recording error at one site spreads into other sites'
exposure, so audits should follow the direction of spread. In this test, audits helped, but targeting them by network position did not beat
targeting the most suspicious records, and direction added nothing. I recommend not writing a merged paper on this thesis.

Everything here is simulated except the geography and survey counts. The two manuscripts under review were not changed. The verdict is computed by
the script and [recomputed independently](../../tools/verify-artifacts.py) from replicate records.

| Mean squared error, full direction | Pooled analysis | Plug-in analysis | Sites with errors found (of 11.3) |
|---|---:|---:|---:|
| Recorded values only, no audit | 0.0243 | 0.0243 | 0 |
| First wave only (8 random sites) | 0.0296 | 0.0238 | 1.5 |
| Random second wave | 0.0270 | 0.0231 | 3.1 |
| Reactive: most suspicious | 0.0216 | 0.0223 | 3.9 |
| Undirected network leverage | 0.0212 | 0.0220 | 3.7 |
| Transmitter: directed leverage | 0.0220 | 0.0224 | 3.8 |
| Influence on the estimate | 0.0217 | 0.0222 | 3.7 |
| Oracle, using hidden truth | 0.0238 | 0.0214 | 9.1 |
| Full audit of 57 sites | 0.0212 | 0.0212 | 11.3 |

| Rule, moderate regime at full direction | Result | Met |
|---|---|---|
| Audits matter: full audit against no audit | MSE ratio 0.873, z = -2.53 | Yes |
| Network helps: transmitter against reactive | MSE ratio 1.018, z = +0.81 | **No** |
| Network helps: transmitter against random | MSE ratio 0.817, z = -2.63 | Yes |
| Direction helps: transmitter against undirected | MSE ratio 1.037, z = +2.62, in the wrong direction | **No** |

## Reading the result

- **Little was left to win.** Reactive auditing already reached 0.0216 against 0.0212 for auditing every site. No targeting rule can gain much in that gap.
- **The adjusted estimator tolerated these errors.** With 20% of site masses wrong by a typical factor of 2.7, using recorded values unchanged still covered 95.0%, with MSE 15% above a full audit. My prior was skeptical for this reason.
- **Direction did not help, and was slightly worse.** Undirected leverage beat directed leverage at full direction. At direction 0 the two were identical, as they must be, which passed the plan's code check.
- **Influence-function targeting did no better either** (MSE ratio 1.004 against reactive), so this is not specific to the transmitter rule.
- **Gross errors matter more, and reactive auditing handles them.** With 10% of sites off by a factor of 10, no audit covered 91.5% with MSE 47% above a full audit. Reactive auditing restored 94.6% and matched the transmitter rule (ratio 1.000).
- **The oracle is not a bound for the pooled analysis**, contrary to the plan's wording. It selects on hidden truth, which breaks the recording model's selection assumption and inflates the apparent error rate. In the plug-in analysis it nearly matches a full audit.

## Limits of this verdict

This tested one recorded input (site mass), one estimator, one contamination strength and a kernel whose form and direction the analyst knows. It does not show that errors never travel through a network, or that audits are unimportant in cluster trials.
Another input, such as recorded intervention coverage, or a wrong kernel might behave differently. That would need a new plan and new seeds, and I do not recommend trying variations until one works.

Other caveats: 127 of 1,600 replicates (7.9%) lost every second-wave policy to first-wave convergence refusals, with only 8 references at 57 sites, so policy comparisons are conditional on the rest. Pooled intervals overcovered (96% to 98%) in several cells. The
[port check](reproduce/port-check.csv) was close but not exact: adjusted RMSE 0.158 against the deposited 0.149. The run took 2,227 seconds, against my estimate of about 15 minutes.

Full tables: [summary](study/summary.csv), [paired contrasts](study/paired-contrasts.csv), [verdict](study/verdict.csv).
