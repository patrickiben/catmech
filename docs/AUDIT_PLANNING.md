# Planning an audit: how many records?

**Status: not validated. A [prospective test](../validation/AUDIT_PLAN_SYNTHESIS.md) returned NO-GO under a rule fixed in advance.** Median error in predicting interval width was within 10% at target audits of 10%, 25% and 50%
(7.8% at worst), but 11.3% to 16.8% at a 5% audit, in cohorts of 600. Read the output below as a documented calculation, and keep that limit in mind.

`catmech_audit_plan()` takes a dataset with recorded values and a small random audit, and estimates how much narrower the interval for one regression
coefficient would become with a larger audit. It plans for linear multiple imputation with the recorded value, `catmech_reference(use_recorded = TRUE)`, which is
the [recommended starting point](CONTAMINATION_METHOD.md). It addresses precision only. It says nothing about whether the audit will find a problem, or about bias or coverage.

## The calculation

With audit fraction `f` and a random audit, classical two-phase theory gives

    width(f) = width_full * sqrt(1 + k * (1 - f) / f)

`k` is the share of the analysis variance that the auxiliaries, meaning the recorded value and whatever else is known on every record, cannot recover. This is not a new result. When `k` is small, recorded values
are nearly as good as references and a small audit loses little. When recorded values are discarded, `k` is close to 1 and every unaudited record costs almost its full information.

Solving for `f` answers the planning question directly. To come within a ratio `r` of the full-audit width:

    f = k / (r^2 - 1 + k)

## Running it

```r
plan <- catmech_audit_plan(
  data, value = "recorded", reference = "reference",
  predictors = c("age", "site", "outcome"),
  formula = outcome ~ recorded + age + site, term = "recorded",
  fractions = c(.10, .25, .50, .75), target_ratio = 1.25, seed = 1
)
plan$k        # k with its interval, keeping and discarding recorded values
plan$table    # records against expected width, for each fraction
plan$answer   # records needed to reach target_ratio
plan$flags    # every warning raised, by short name
```

`data` holds every record in the cohort. `reference` is observed on the pilot only. `predictors` should include the analysis outcome and every complete auxiliary you would use in the imputation. A runnable version is in `inst/examples/audit-plan.R`.
The formula must be linear in the audited value, so transform the column first rather than writing `log(value)`. The audited value may be the outcome, a covariate, or part of an interaction, and `term` may be any coefficient.

## Reading the output

| Field | Meaning |
|---|---|
| `k$estimate` | The planning value of `k`. The `recorded` row keeps the recorded value and the `discard` row is the `catmech_reference()` default. |
| `k$lower`, `k$upper` | Bootstrap interval for `k`, seeded, with the caller's random state restored. |
| `k$classical` | The empirical version of the ratio. See below. It is not used for planning width. |
| `k$rubin_check` | The value implied by the pilot fit's own between-imputation variance. A large gap from `estimate` deserves a look, though this value is noisy with few imputations. |
| `table$relative_width` | Expected width over the full-audit width. |
| `table$relative_lower`, `relative_upper` | The same at the limits for `k`. These describe the expected width, not what one audit will realize. |
| `table$predictive_lower`, `predictive_upper` | A range for the width of one audit of that size, which adds that audit's own noise. |
| `table$width` | Width in the coefficient's units, anchored on the pilot's own fit. |
| `table$discard_relative_width`, `discard_width` | The same calculation when recorded values are discarded. |
| `table$small_audit` | That size is below the minimum the function accepts as a pilot. Treat its row as a rough guide at most. |
| `answer$records`, `records_lower`, `records_upper` | Records needed to reach `target_ratio`, from the estimate and from the limits for `k`. |

## Refusals and warnings

| Situation | What happens |
|---|---|
| Fewer than 30 rows with both values, or fewer than 20 residual degrees of freedom in the imputation model | Refused |
| Pilot design rank deficient, for example a category absent from the pilot | Refused |
| Formula not linear in the audited value, or analysis variables with missing values | Refused |
| Fewer than 100 pilot rows | `small_pilot` warning |
| Fewer than 10 disagreeing pilot rows | `few_disagreements` warning |
| No disagreeing pilot row | `no_disagreements` warning, and `k` is estimated as zero. The warning states how high the disagreement rate could still be. |
| Pilot looks selective on the recorded value | `pilot_looks_selective` warning |
| Outcome missing from `predictors` | `outcome_not_in_predictors` warning |
| Pilot discrepancy above 5 robust SDs of the reference, or recorded value above 10 robust SDs from the median | `gross_errors` warning. Do not rely on the plan. Review `catmech_recording_audit()` and consider `catmech_contamination()`. |
| Fewer than 80% of pilot resamples could be fitted | `bootstrap_lost_support` warning |

## Assumptions

- The pilot is a simple random sample of these records, and the later audit is random too. A targeted pilot breaks the formula. The check for this looks only at the recorded value, and it cannot prove randomness.
- The analysis is linear MI with Rubin's rules through `catmech_pool_lm()`, with the same formula and predictors.
- The reference is exact and measures the same quantity at the same time in the same units.
- The relationship between reference and recorded value is close enough to linear. Linear MI [failed under synthetic gross errors](../validation/CONTAMINATION_SYNTHESIS.md), which is why that case gets a warning.
- A pilot says nothing about errors it did not see. If errors are rare and large, a small pilot can miss them and report a `k` that is far too low.

## How k is calculated

Let `psi` be the influence function of the coefficient using true values, and `aux` what is known on every record. Classical theory gives `k = Var(psi - E[psi | aux]) / Var(psi)`.

The planner evaluates that ratio under the working models, because the interval from Rubin's rules comes from them. Write the analysis as linear in the audited value `y`. Under the normal imputation model `y ~ N(mu, tau^2)`, fitted on the pilot's reference rows, `psi` is a quadratic in `y`, so its conditional variance has the closed form

    Var(psi | aux) = tau^2 * (u * E1 + v * A1)^2 + 2 * (A1 * E1)^2 * tau^4

where `u` and `v` are the analysis's design weight and residual at `mu`, and `A1` and `E1` are their slopes in `y`. `k` is the average of this over all records, divided by the model-based `Var(psi)`. Records without a recorded value use an imputation model without it, as `catmech_reference()` does.
It needs no simulation, and with one normal covariate it reduces to the known `(tau^2 / S) * (1 + b^2 * S / s^2 - 2 * b^2 * tau^2 / s^2)`, which the unit tests check.
The exact steps are written out for porting in `validation/audit-plan-study/planner-data.json`.

The `classical` column computes the same ratio with empirical variances on the pilot rows. When residual variance is not constant the two differ, and each has its own use: the working-model value describes the interval Rubin's rules report, and the empirical value describes the true sampling variance.

On the evaluation pools, working-model `k` was 0.042 on ClinicalTrials.gov and 0.070 on NHANES, and it matched the frozen study's interval widths. The empirical value was 0.028 and 0.088, closer to how much RMSE grew in six of seven cells. Use `estimate` to plan interval width. A large gap from `classical` means reported intervals and true sampling error are growing at different rates, and that the analysis interval deserves its own review.

`width` needs a full-audit width to scale. It takes that from the pilot fit's within-imputation variance, which is nearly free of imputation noise, rather than dividing the pilot's total width by its own predicted inflation.

## How accurate it was

From the [synthesis](../validation/AUDIT_PLAN_SYNTHESIS.md), on cohorts of 600 from two public datasets. Median and 90th-percentile absolute error of predicted width, for a pilot of 60 or 150 records:

| Dataset, pilot | Target 30 records | Target 60 | Target 150 | Target 300 |
|---|---:|---:|---:|---:|
| ClinicalTrials.gov, 60 | **11.3%**, 38.8% | | 2.4%, 8.5% | 1.5%, 4.8% |
| ClinicalTrials.gov, 150 | **16.8%**, 43.8% | 6.6%, 20.6% | | 0.9%, 2.6% |
| NHANES, 60 | **13.9%**, 35.3% | | 3.4%, 9.6% | 2.8%, 7.2% |
| NHANES, 150 | **15.5%**, 37.5% | 7.8%, 19.6% | | 1.6%, 3.8% |

What that means in practice:

- **Very small audits are unpredictable, whatever the planner does.** A predictor that knew the median width at 30 records still missed by a median of 16% to 17%, and two fits of the same audit with different imputation seeds differed by 4% on ClinicalTrials.gov and 9% on NHANES. In 7% to 8% of cohorts a 30-record fit was refused because a rare category had no usable row. Do not plan on an audit that small in a design like these.
- **`k` itself is hard to estimate from a small pilot.** On ClinicalTrials.gov, where errors are exact-or-different, 60-record pilots gave `k` from 0.0012 to 0.131 (5th to 95th percentile) around a pool value of 0.042. The answer in records swung from 2 to 114. Its interval contained the pool value in 54% to 74% of pilots across pilot sizes and datasets, not the nominal 95%. Plan on `records_upper`, and treat it as a guide and not a guarantee.
- **Widths at moderate fractions are forgiving.** The formula is not sensitive to `k` at 25% or more, which is why those predictions held up even when `k` was poorly estimated.
- **The `k` interval is not a range for your audit's width.** It covered the actual ratio only 35% to 68% of the time. The predictive range did better, at 60% to 92%, and was still below nominal.
- **The discard columns understate small audits.** Against the frozen study's widths, by 19% and 13% at 5% and 10% on ClinicalTrials.gov and 13% at 10% on NHANES. Their direction is clear even so: on these pools, reaching within 25% of the full-audit width took about 380 of 600 records when discarding recorded values, against 41 and 66 when keeping them.
- **Gross errors break it.** In the synthetic gross-error cell, median error reached 742% and linear MI's own coverage fell to 32% at a 5% audit. The warning fired on all 600 pilots there, on none in the other synthetic cell and on NHANES, and on 28 of 2,991 ClinicalTrials.gov pilots, a dataset where linear MI held up in the frozen study.

## Limits

The evidence comes from two public datasets, in cohorts of 600. Neither is a trial's source-data verification, and the ClinicalTrials.gov reference is sponsor-entered. Resampling a pool says nothing about transport to other data. The planner addresses one coefficient at a time. It does not handle clustered records or targeted audits, and it does not protect against a pilot that happened to miss the errors that matter.
