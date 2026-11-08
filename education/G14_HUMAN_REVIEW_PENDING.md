# G14 Human Review Pending: catmech 1.0.1 Companion

- **Opened:** `2026-11-07`
- **Status:** PENDING. The automated checks in `verify_educational_companion.py` cover the locks, coverage, evidence classes, readability, links and R labs. The human rows below have not been done for these files.

## Files awaiting review

- Manuscript `paper/manuscript.Rmd`: SHA-256 `e25e6b44561b4c5b8529d87867eda249665ab2269aa406e8217a90e3639c0758` (350 lines)
- Manuscript PDF `paper/manuscript.pdf`: SHA-256 `79c865d1e4f4ef216599aeb7b4ce931c237bb89b2c525c5b97d50ea6f8b56df1` (16 pages)
- Companion `education/COMPANION.md`: SHA-256 `d363ba037e026d34b877fb1ad571e92b4a5c1a84ddec16c5418a2ef2fc5c2390` (30 units)
- `education/COMPANION.pdf`: SHA-256 `d7a996bc23cce5c9752357be85db561c9ed6a7284e66865090b73b3e690099a9` (25 pages)
- `education/STUDY_GUIDE.pdf`: SHA-256 `dc6945a8dc5310adc342fc43d41561fb444bcb74184bde879a44a5200162525c` (9 pages)
- `education/ONE_PAGER.pdf`: SHA-256 `cdad85f32396fc55abecf231f542b5e7eaafc7291279d92ee12efcf08ec436b0` (1 page)
- `education/companion_lab.html`: SHA-256 `f48289c53e6b17a3da18f80515235168254bb0887b40b2664dd31f67ba2966ed`
- Transfer assessment: SHA-256 `6ed39436c9055653b7b7cafbc6cba212ea40a61d8dc34d8b79600f878b10b041` (12 items, unchanged)

[G14_HUMAN_SIGNOFF.md](G14_HUMAN_SIGNOFF.md) records the review of the earlier versions (manuscript SHA-256 `8152288b1c2af2a49103543db3c9324582177b1e842365b5c5ab83cdac8c1fd0`, companion SHA-256 `c47cb3847813a7c1cf88ec930026c53c8c5b0ee50fdb756cc234a03b596d5985`). It does not cover the files above, and it was left unchanged. The companion build no longer writes that file: it reads it, and marks the human rows GREEN only when the hashes it records match the files just built.

## Waiting on the reviewer

1. **Claim review:** no companion statement goes beyond the manuscript (`unrecorded_claims`).
2. **Plain-language and jargon review** of the `In ordinary language` fields.
3. **Semantic equivalence** of `COMPANION.pdf`, `ONE_PAGER.pdf` and `companion_lab.html` with `COMPANION.md`, and the practice feedback in `companion_lab.html` (attempts, correctness, explanations, retry).
4. **Text, render and accessibility QA** of `COMPANION.pdf`, `STUDY_GUIDE.pdf`, `ONE_PAGER.pdf` and `companion_lab.html`.
5. **Transfer assessment and the one-sentence rubric.** The twelve items are unchanged (same assessment SHA-256), so the question is whether a reader of the new text still passes them.
6. **G12:** the author's read of the revised manuscript, recorded in [G12_RESIGN_2026-11-07.md](../paper/G12_RESIGN_2026-11-07.md).

## What changed in the companion

- **C17 and C18** (Section 4 workflow): every analysis outcome in the imputation model, the Nelson-Aalen hazard, the rough audit-size step, the new printed Cox numbers (`0.1877`, hazard ratio `1.2065` with interval `[1.0982, 1.3255]`), and the hash wording in stage 8.
- **C09** (Section 2.3): `catmech_audit_power()` as a rough, unvalidated aid whose `power` column is not the power of a test, with the linear-imputation comparison at the manuscript's settings (`n = 100` against about `n = 257`). It also corrects "28,800 plan rows" to the 7,200 plans recorded in `validation/AUDIT_PLAN_SYNTHESIS.md`.
- **C28** (Section 9): what the hashes are of, and that they are not tamper-evidence.
- **C29** (Section 10): 199 tests with 1,812 expectations.
- **C01, C15 and C16:** matching audit-size and hash wording. **C02 and C08:** manuscript lines changed, companion text unchanged.
- **Study guide lessons S4 and S12, the one-pager, the console tutor (S4 and S12) and the R labs (C09, C17-C18 and C28):** the same changes. The one-pager now says every analysis outcome goes into `predictors`, and the labs check that the GLM and Cox pooling draw no warning.
- **Page locators** for all 30 units and 12 lessons were recomputed from the 16-page PDF. Most of the old locators were about one page off even for the 15-page PDF; for example, the stage 6 to 8 code and output were on page 7, not page 8.

## To sign

Read the files above, then update [G14_HUMAN_SIGNOFF.md](G14_HUMAN_SIGNOFF.md) yourself with the new manuscript and companion SHA-256 values, the date and your verdicts, and rebuild with `education/build_educational_companion.py`. The build then marks the human rows GREEN, records the transfer scores only if the sign-off states 12/12 and 4/4, and removes this file.
