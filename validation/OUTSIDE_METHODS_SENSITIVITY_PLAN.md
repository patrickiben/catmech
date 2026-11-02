# Follow-up to the comparison with outside methods: two settings

Written 2026-11-02, **after the results of the [comparison](OUTSIDE_METHODS_PLAN.md) were read**, and committed before this rerun. So this is post hoc. It cannot change the outcome of the rule, which was KEEP. It asks whether that outcome depends on two settings in the comparators' recipes.

## What prompted it

- **`mice` removed the recorded value.** By default, `mice` drops a predictor whose correlation with the variable being imputed is at least 0.99 among observed rows. That happened in 4,424 of 49,500 requests, including 1,571 of 2,000 for `deliveries` at 25%. The poor results for `mice_2l_lmer` and `mice_2l_pan` are concentrated on those items. The plan anticipated logged events, and recorded them, but did not anticipate how many.
- **`survey::twophase(method = "full")` understated standard errors with several rows per cluster.** `raking` covered 55.6% to 79.7% at 25%, and worse at 50%, against 92.6% and 93.8% on the one-row control. On fresh seeds, the standard error ran 0.46 to 0.73 of the spread of the estimate without calibration when clusters had three rows, and 0.97 with one row. On one dataset the same total had a standard error of 296 at row level and 589 after aggregating to clusters under `full`, and 417 both ways under `approx` and `simple`. Whether this is outside the method's intended use has not been settled. Point estimates do not change with this setting.

## What runs

- **`mice_2l_lmer` and `mice_2l_pan` with `eps = 0`,** which turns off that check, so the recorded value stays in every model. Same seeds. Only requests with a logged event or no interval in the comparison are rerun, apart from `control_one_row`. As a control, the first three other cohorts in every cell must reproduce the comparison exactly.
- **`raking` with `method = "approx"`,** on all 49,500 requests. Nothing else changes.

`difference` and `audited_only` are not rerun. The other methods keep their comparison results.

## How it is read

The rule's conditions N1 to N6 are applied again by the same scorer, with these results in place of the originals.

Fixed before the rerun:

- The rule's outcome stays KEEP as written, because the rule named those defaults.
- If a rerun method is no worse under N1 to N6, the write-up says plainly that KEEP depended on a software default, and the package documentation names that route as an equal alternative. Whether to retire the option then goes to Patrick as a decision, with a recommendation.
- If none is no worse, the write-up says KEEP holds under these settings too.

## Predictions fixed before the rerun

1. `mice_2l_lmer` with `eps = 0` has an N1 ratio at most 1.05.
2. `mice_2l_lmer` with `eps = 0` still fails N3, reaching 92% on fewer than nine items at 25%.
3. `raking` with `approx` still fails N1, which is certain since its estimates do not change, and stays below 92% on at least nine items at 25%.

## Limits

These two settings were chosen after seeing results that looked wrong, and the option's own recipe has had many more rounds of this attention. `difference` also undercovers, with standard errors as low as 0.36 of the spread in the design with two audited clusters per stratum. That looks like unadjusted leverage from stratum effects, which is not tested here, and it is not corrected. It fails N1 regardless. No reviewer has seen this.
