# Sponsor clustering screen on ClinicalTrials.gov

Written 2026-11-01 before any evaluation-split analysis and frozen by hash in `tools/sponsor-clustering.py`. It is not changed after this point. A change requires a new versioned run, and the first run stays on record.

## The question

Do real ClinicalTrials.gov discrepancies cluster within sponsors beyond chance, and does finding one predict more in that sponsor's other trials?

This is a cheap screen on real data. It tests clustering only. It cannot show errors spreading, because shared definitions, reporting habits, CROs, SOPs or training could all produce clustering with no spread. Registry entry is also not source-data verification, so a null result closes the registry route only, not clinical-site errors in general.

## What has already been seen

- Overall rates by phase and by sponsor class were already seen in the earlier audit ([synthesis](REAL_DATA_SYNTHESIS.md)): 14.2% disagreement on the evaluation rows, with rates differing across phase (p = 0.0001) and sponsor class (p < 0.0001).
- Before writing this plan I downloaded sponsor names and counted trials per sponsor by split, with no discrepancy data involved. I also saw the covariate levels, the start-year counts of the development split, and the first 11 rows of the pool file, which carry no sponsor names.
- No discrepancy rate by sponsor has been computed on either split. Neither split has been run through any test below.

## Data

- Pool: `validation/real-data/ctgov-pool.csv`, SHA-256 `c2e02ea2542b8291dff54799a9f0dfecfdbcadb71ebf28e90aa798764e0a91ca`. Its `recorded_n`, `reference_n` and covariates are kept exactly.
- Sponsor map: `validation/sponsor-clustering/sponsor-map.csv`, SHA-256 `31a223e87de063d0def6d6befc5f9a4626ddfba9a6fe55291dca56144f46ec81`. It was fetched with the same API filter as `tools/real-data-fetch.py`, requesting only `NCTId` and `LeadSponsorName`, and joined by `nct_id`. Query, retrieval time, match rate and page hashes are in `validation/sponsor-clustering/fetch-manifest.json`. All 34,084 pool rows matched.
- Name normalization is trimming and case-folding only, with no manual merging. That leaves 5,488 distinct names, the same count as before normalization, so it merged nothing. Variants of one organization stay separate, and one name can cover many independent groups.
- Every sponsor has one sponsor class in each split, so sponsors nest within class.
- Discrepancy: D = `recorded_n` differs from `reference_n`. Sensitivity outcome: D2 = |log(`recorded_n` / `reference_n`)| > 0.25, natural log.
- The pool's development split (the `split` column) is used for engineering only. Evaluation rows are the confirmatory data.

Counts on evaluation rows, from names alone: 1,254 sponsors with at least 3 trials (22,912 trials); 931 sponsors with at least 4 trials (21,943 trials); under the split rule below, 465 test sponsors with 11,512 trials and an audit budget of 1,151.

## Design, as fixed by the reviewer

**Test 1, excess clustering.** Among sponsors with at least 3 evaluation trials, the statistic is the sum over sponsors of (observed D - expected D)^2 / expected D, with expected D from the rates within strata of phase x sponsor class x start-year band x trial-size quartile. Permute D across trials within those strata, 2,000 times with a fixed seed, keeping sponsor sizes fixed. Report the one-sided permutation p-value and the observed-to-null ratio.

**Test 2, prediction within held-out sponsors.** Among sponsors with at least 4 evaluation trials, assign each sponsor to train or test by the SHA-256 of its normalized name, half each. In every sponsor, reveal a random half of its trials (fixed seed), as an audit would. Fit two logistic models on training sponsors to predict the unrevealed trials: baseline with phase, sponsor class, start-year band and log trial size, and history with those plus a shrunken revealed discrepancy rate, (revealed discrepancies + 0.5) / (revealed trials + 1), on the logit scale. Score unrevealed trials in test sponsors. Report the AUC difference and Brier difference (history minus baseline) with 95% intervals from a bootstrap that resamples test sponsors, 2,000 times.

**Test 3, audit yield.** In test sponsors, with a budget of 10% of trials, compare discrepancies found per 100 audited under: random, risk (the baseline model's top-ranked trials), and adaptive (spend half the budget at random, then spend the rest on unaudited trials of sponsors where a discrepancy was found, ranked by the baseline model, falling back to the baseline ranking if those run out). Average over 200 fixed seeds for the random elements. Report the adaptive-to-risk yield ratio with a 95% interval from the same sponsor bootstrap, and report both adaptive and risk against random.

**Verdict, computed by the script on D.**

- CLUSTERING-SUPPORTED: Test 1 p < 0.01, AND the AUC difference is at least 0.02 with a 95% lower limit above 0, AND the adaptive-to-risk yield ratio is at least 1.2 with a 95% lower limit above 1.
- WEAK: Test 1 p < 0.01, but the other conditions fail.
- NO-CLUSTERING: Test 1 p >= 0.01.

D2 and results by sponsor class are descriptive and do not enter the verdict. Comparisons use unrounded values.

## Operational details left open above

These fill gaps in the design and do not alter it.

**Variables.**

- Trial size is `recorded_n` throughout. `reference_n` is what an audit would reveal, so using it to choose audits or predict D would leak the outcome.
- Start-year bands: 2013 or earlier, 2014-2015, 2016-2017, 2018 or later. I set them from development-split year counts alone, to make the bands close to equal in size.
- Phase has the pool's four levels. Sponsor class has its two.

**Test 1.**

- The analysis set is the evaluation trials of sponsors with at least 3 evaluation trials. Stratum rates, quartile cut points and permutations all use this set only.
- Quartile cut points are type 7 quantiles of `recorded_n` in the analysis set. A trial is in the lowest group whose upper cut point is at or above its size.
- A sponsor's expected D is the sum of its trials' stratum rates. A sponsor with expected D of 0 contributes 0, and its observed D is then 0 under every permutation.
- Each permutation shuffles D among trials in the same stratum. Every trial keeps its sponsor and stratum.
- p = (1 + number of null statistics at or above the observed statistic) / 2,001, so its floor is about 0.0005. The ratio is the observed statistic over the mean null statistic.

**Test 2.**

- Split: sort the eligible sponsors by the hexadecimal SHA-256 of the UTF-8 normalized name. The first half, rounded up, is train and the rest is test.
- Reveal: draw one `random.Random(seed).random()` value per eligible trial, in `nct_id` order. In a sponsor with m trials, the floor(m/2) trials with the lowest values are revealed. There is one fixed reveal, so the intervals do not include reveal-to-reveal variation.
- Models are fit on the unrevealed trials of training sponsors by Newton's method, with an intercept, dummies for phase, class and start-year band, and natural-log `recorded_n`. History adds logit((revealed D + 0.5) / (revealed trials + 1)).
- AUC uses midranks for ties. Brier is the mean squared difference between prediction and outcome.
- Bootstrap: draw as many test sponsors as there are, with replacement, by `int(random.Random(seed).random() * K)`. A sponsor drawn twice enters twice. Models and the reveal stay fixed. Intervals are type 7 percentile limits at 2.5% and 97.5%.

**Test 3.**

- The population is all evaluation trials of test sponsors, revealed or not, because this is a separate simulated audit. The baseline model is the one fit in Test 2. Revealed status plays no part.
- Budget B = floor((N + 5) / 10) for N trials. The random stage of adaptive audits B1 = floor(B/2) trials and the second stage audits B - B1.
- Ranking is by baseline prediction, highest first, with ties broken by ascending SHA-256 of `nct_id`. Risk audits the top B.
- For each seed, one random order of the trials is drawn. Random audits its first B trials and adaptive's random stage audits its first B1.
- A sponsor is flagged if the random stage found a discrepancy in it. Flagging is settled once, after the random stage, and is not updated during the second stage. The second stage audits the top-ranked unaudited trials of flagged sponsors, then if those run out the top-ranked remaining trials.
- Yield is 100 x found / B. Random and adaptive yields are means over seeds 0 to 199. Ratios are ratios of those means. For the point estimates the random order comes from `random.Random(base + seed)`, one value per trial in `nct_id` order, so a plain-Python check can repeat it exactly.
- Bootstrap: the draws are the same as in Test 2. Each replicate recomputes N, B, the ranking of its trials and all yields. Each copy of a sponsor is a separate sponsor, so a finding in one copy does not flag the other. Copies of a trial tie and stay in copy order. Each replicate's 200 random orders come from NumPy PCG64 seeded with (base, replicate, seed).
- A bootstrap replicate with an undefined value, such as a risk yield of 0, is left out of that interval and the count of such replicates is reported.

**Seeds.** Test 1 permutations 20261101; reveal 20261102; bootstrap draws 20261103; point-estimate audits 20261104 + seed; bootstrap audits (20261105, replicate, seed). D2 uses the same seeds, the same reveal and the same bootstrap draws as D.

**Descriptive results, outside the verdict.**

- D2 goes through Tests 1 to 3 in full, with its own stratum rates and its own refit models.
- By sponsor class: Test 1 is decomposed by class from the same permutations, using the same cut points. Test 2 differences use the same models and bootstrap draws, restricted to each class's unrevealed test trials. Test 3 runs within each class's test-sponsor trials with its own budget, with point estimates over the 200 seeds and no bootstrap interval.

**Engineering.** The same code runs on the development split first, with the same thresholds applied to development trials and with development output marked as engineering only. Code may be fixed after a development run, but this plan may not. The evaluation run happens once into `validation/sponsor-clustering/evaluation/`. A failed fit (no convergence within 100 Newton steps) stops the run and is reported as a failure. No setting is changed to obtain a fit, and no results are retried or retuned.

**Checks.** The script checks the hashes of this plan, the pool and the sponsor map before and after the run, refuses to overwrite, and uses the fixed seeds above. `tools/verify-artifacts.py` recomputes the main results in plain Python from the raw files.

## Prediction, written before running anything

My honest expectation is **WEAK**, at about 75%, with CLUSTERING-SUPPORTED about 20% and NO-CLUSTERING about 5%.

- Test 1: strong excess clustering on D, p at or near the floor and an observed-to-null ratio of about 2 to 4. Sponsors use their own habits when filling in these fields.
- Test 2: baseline AUC about 0.60 to 0.68. History adds about 0.03 to 0.08, with a lower limit above 0 and a small negative Brier difference.
- Test 3: yields per 100 of about 14 for random, 20 to 28 for risk and 18 to 25 for adaptive, for an adaptive-to-risk ratio of about 0.8 to 1.1. That fails the 1.2 condition, because adaptive spends half its budget at random and large sponsors get flagged easily. The risk ranking can already target the trial types that discrepancies concentrate in.
- D2: clustering of the same kind, weaker in absolute count. Class: the industry ratio is higher than the other-class ratio, with low confidence.

## Interpretation limits

- CLUSTERING-SUPPORTED only justifies looking for record-level audit data with shared-resource metadata. It says nothing about direction or spreading.
- WEAK means clustering beyond chance without enough prediction or audit gain to justify adaptive targeting in this registry setting.
- NO-CLUSTERING closes the registry route only, not clinical-site errors in general.
- The result does not cover clinical source-data errors. D is a difference between two sponsor-entered registry fields, and enrolled and started can legitimately differ, so D is not an error rate.
- Clustering here can come from shared definitions, reporting habits, CROs, SOPs or training with no spread at all.
- Lead sponsor is the only shared factor tested. It is not a monitor, clerk, lab or system, and names are not merged, so a large or fragmented name can misstate the real unit.
- The pool holds only completed interventional trials with posted results, so it does not represent all trials.
