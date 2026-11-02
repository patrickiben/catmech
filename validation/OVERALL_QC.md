# Overall check before a possible public release

Done 2026-11-02 at commit `423a272`, by five independent AI reviewers working from a fresh clone, covering code, statistics, evidence, documents and release readiness. They saw the code and not the assistant's conclusions. No human took part, so this is not peer review. The assistant that wrote most of this repository then checked each finding against the code and data, and made the fixes. No reviewer has read the fixes.

## What held up

The closed forms and draws matched independent calculations: the jackknife against brute-force refits, Rubin's rules and Barnard-Rubin degrees of freedom, the variance components against explicit projections, the coefficient and intercept draws against exact formulas, and the contamination sampler against a separate one. The numbers in the latest write-ups recomputed from saved results using separate code, and saved rows reproduced through the runners. A search of the full history found no credentials.

## Fixed in version 0.7.1

**Results that were badly wrong with no signal.** All three were reproduced before fixing. Each now warns or refuses, and no number changes in a case that still returns one.

- **Every audited row agrees.** Recorded values were kept with no added uncertainty and no warning. With a 3% audit and large errors in 5% of records, this happened in about a fifth of audits, and those intervals never covered the truth.
- **One cluster determines a coefficient.** The jackknife's pseudoinverse leaves out that cluster, so a site type found at one site covered about a third of the time.
- **One audited cluster per stratum.** With strata in `predictors`, as the documents advised, between-cluster variance could not be estimated, and `spec$cluster_uncertainty` still said it was drawn.

**Smaller code issues.** The five-cluster refusal skipped the all-agree case. Blank cluster ids gave an internal error or merged. An exact offset was refused without `cluster` and accepted with it. `catmech_recording_audit()` told about 30% of truly random small audits that they were not random. Two edge cases in `catmech_holdout_check()`. `catmech_sensitivity()` compared variances of different kinds without saying so. Three steps had no test that would catch their removal, and the two in the clustered draw now have them.

**Claims corrected.**

- The first-run `raking` recipe omitted the phase-one variance. The write-up had blamed inconsistent results in `survey`.
- The documents said the slope bias comes from the model and not this code. The option attenuates more than `mice_2l_lmer` in all six cells.
- The corrected-variance write-up said the plan came before its runner existed. Git shows a 39-second gap with the runner already drafted, so the supported order is plan, then frozen runner, then study.
- The README's table listed one failure for linear MI when the same study recorded others, recommended raking without qualification, implied that the jackknife repairs any method, and quoted coverage from the earlier variance as current.
- The case study's "Nothing was retuned after viewing" overgeneralized, since the clustered option was redesigned on items already seen.
- The AI disclosure said less than the record: an assistant wrote most later plans, rules, scoring and write-ups.
- Several README examples did not run. Stale statements in other documents were corrected.

**Data.** NHANES and ClinicalTrials.gov derived files had no credit, terms or data note, and the only license a reader met was MIT. There is now [a data note](real-data/README.md). The ClinicalTrials.gov terms ask distributors to keep data current, which a frozen snapshot does not do, and the note says so.

## Left open

- **An email address.** `validation/sponsor-clustering/sponsor-map.csv` held one person's institutional email address, in a registry field. It is now a placeholder, as `docs/GITHUB.md` records.
- **Local paths in logs.** Frozen logs show the account name and scratch paths. Two scratch scripts have the path hard-coded and cannot run elsewhere. These are accepted, since editing frozen evidence has its own cost.
- **Smaller items not changed.** `assumption_hash` for `catmech_model` does not depend on the predictors. Seeded draws can differ across collation locales with mixed-case character labels. The joint scale move in the contamination sampler has no test that would catch removing its Jacobian, though a separate sampler agreed with it. The verifier's check on a write-up's decimals accepts a number found anywhere in its sources, and does not check integers.
- **Not checked.** `clubSandwich` and `sandwich` are not installed, so `cv3` was compared with brute-force refits only. No cited source was opened by the statistics reviewer. No secret scanner was run, only pattern searches.
- **Still true of the science.** Linear imputation undercovers on real clustered items. Rare counts are unsolved. All clustered tests use one health system, the comparison was set up by the option's author, and no independent statistician has reviewed any of it.
