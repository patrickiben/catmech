<!-- gate-basis: README.md@025b4c2030e41f168eaa1b05c58c24750d396eaea9a1647e7477c054ca2eba21 -->
<!-- gate-basis: paper/verify_paper_numbers.py@4dc3caf4abd3e584190e033e6afcffdd3486e61bbcbf290b54bc5318a14e40f4 -->
<!-- gate-basis: paper/paper.bib@e4efb6523c93e3ea4ab40847542c68a6cc6b0e8dce0d33ac370f3d2e1b577cad -->
<!-- gate-basis: paper/paper.md@05b86f68fa8192e5dc33cac32954031ddad3ef99e970cc581b44961303a76deb -->
<!-- gate-basis: paper/manuscript.pdf@9e55292a53694b361fae5bba088659a67df358e5e88c2c8a89f66b64cb09a289 -->
<!-- gate-basis: paper/manuscript.tex@4349cd68b70d01a8182f86cfbeb231790f9eabd4cf1f65e535a107fc926c2ca6 -->
<!-- gate-basis: paper/manuscript.Rmd@e7130db27be8203337a9acd984e5a6744a7dcc743e162d756bb9c61b63a3b1f6 -->
# Pre-Submission Validation Harness Ledger (`catmech` v1.0.0)

- **Date**: 2026-11-06
- **Target**: `catmech` v1.0.0 CRAN package & JSS / JOSS manuscripts (`paper/manuscript.Rmd`, `paper/manuscript.tex`, `paper/manuscript.pdf`, `paper/paper.md`, `paper/paper.bib`)
- **Harness Specification**: `/Users/patrickiben/Documents/Documents - Patrick’s Mac Studio/Research_Validation_Harness/PRE_SUBMISSION_VALIDATION_HARNESS.md`
- **specification_coverage**: `lines 1-590 (100%)`
- **Total Gates in Specification**: `24` (`G0`–`G23`, verified via `grep -cE '^## G[0-9]+:' PRE_SUBMISSION_VALIDATION_HARNESS.md`)
- **Execution Tier**: Inline Deterministic Tier + Self-Audit Pass (multi-agent subagent waves `G5`, `G12`, `G13`, `G18` require separate user cost authorization per harness protocol)

---

## Gate Summary Table (All 24 Gates)

| Gate | Title | Status | Tool / Method | Verdict & Evidence |
| :--- | :--- | :---: | :--- | :--- |
| **G0** | Format-conversion pre-gate | **RUN** | Pandoc / RMarkdown inspection | **PASS**: Repaired 4 list blocks lacking preceding blank lines and institutional BibTeX braces (`fda2003part11`, `ich1998e9`). |
| **G1** | Path, hygiene & repository-surface guard | **RUN** | `check_paths.sh` + `R CMD check` + `verify-release-archive.py` | **PASS**: 0 machine paths in shippable package/paper (`R/`, `src/`, `man/`, `vignettes/`, `inst/`, `tests/`, `tools/`, `paper/`, `README.md`, `DESCRIPTION`, `NAMESPACE`). 2 historical paths in `.Rbuildignore`-excluded `validation/clustered-fix-scratch/` left untouched per Rule 3 (frozen archive). |
| **G2** | Numeric-claim <-> artifact check | **RUN** | `tools/verify-artifacts.py` + `paper/verify_paper_numbers.py` | **PASS**: Replaced ungrounded synthetic table in `manuscript.Rmd` and `paper.md` with exact frozen numbers from `validation/contamination-study/summary.csv`. 343 claims + negative control verified in `paper/verify_paper_numbers.py`. |
| **G3** | Citation & bibliography verification | **RUN** | `check_citations.py` on `manuscript.bbl` and `\nocite{*}` (all 25 `.bib` entries) | **PASS (`0 FLAG`)**: Fixed 4 broken/hijacked DOIs/entries (`huang2013simple`, `breslow2003semiparametric`, `neyman1938contribution`, `little2019statistical`, `nab2021mecor`). Final result: `20 OK | 5 REVIEW (books/manuals verified by hand) | 0 FLAG`. |
| **G4** | Mathematical proof & derivation audit | **RUN** | `tools/check-contamination-oracle.R`, `tools/check-reference-oracle.R`, `tools/check-cluster-algebra.R` + inline audit | **PASS**: Posterior moments, row posteriors, Barnard-Rubin df, CV3 jackknife, Lin-Wei, and Liang-Zeger sandwich formulas match frozen oracle checks (`pass == TRUE`). |
| **G5** | Adversarial referee panel (3 personas) | **PARTIAL** | Inline Qwen peer-review + deterministic self-audit (3-subagent wave deferred pending user cost approval) | **PASS (Inline)**: Addressed Qwen referee findings (Section 7 simulation-to-practice divergence, Section 3 `catmech_pool_coxph`/`catmech_pool_glm` prominence, collapsible README history, exact S3 print output). |
| **G6** | Figure, table & layout inspection | **RUN** | `pdflatex` log + visual page-by-page inspection (`read_file` on `manuscript.pdf` pages 1–10) | **PASS**: `errors=0 undefined=0 overfull=0 underfull=2` across all 10 pages. Author header stacked across 3 centered lines; tables and ASCII decision tree fit within margins. |
| **G7** | AI-artifact & placeholder scrubber | **RUN** | `scrub_ai_artifacts.py` | **PASS**: `0 source hit(s)` across `manuscript.Rmd`, `manuscript.tex`, `paper.md`, `README.md`. |
| **G8** | Anonymization & blind-review check | **NOT RUN** | Target venues (CRAN / JSS / JOSS) are non-blinded | **N/A (Non-blinded venue)**: CRAN, JSS, and JOSS require full author attribution, affiliation, and ORCID. |
| **G9** | Reproducibility & environment audit | **RUN** | `tools/verify-artifacts.py` + `tools/verify-release-archive.py` + `R CMD check` | **PASS**: All 10 frozen study archives, seed streams, SHA-256 source manifests, and `R CMD check` (`0 errors, 0 warnings, 0 notes`) verified. |
| **G10** | Prose style, humanness & US spelling | **RUN** | `check_prose_style_v2.py` + `check_spelling_us.py` | **PASS (`ALL GREEN`)**: `0` British spellings; `0` em-dashes; `0` curly quotes; `0` Title-Case headings; all 12 lexical/artifact categories `[ok]`. Replaced "proves" with "demonstrates" in Section 7. |
| **G11** | Cross-reference, symbol & notation consistency | **RUN** | LaTeX compile (`undefined=0`) + symbol check | **PASS**: Zero undefined references or citations; notation ($T_i, W_i, R_i, \delta_i, \pi_i, \nu_{com}, \nu_{BR}$) consistent across `manuscript.Rmd` and `paper.md`. |
| **G12** | Span-attributed prose verification | **PARTIAL** | Deterministic claim verification via `paper/verify_paper_numbers.py` (2-subagent wave deferred pending user cost approval) | **PASS (Deterministic)**: All numeric, tabular, and code-surface spans verified against frozen CSVs and `NAMESPACE`. |
| **G13** | Cold-reader comprehension check | **NOT RUN** | Requires isolated subagent (deferred pending user cost approval) | **DEFERRED**: Available on request via Lean multi-agent wave. |
| **G14** | Ethical, IRB & data-use compliance | **RUN** | Public dataset manifest & provenance check (`validation/real-data/fetch-manifest.json`, `validation/real-audit/fetch-manifest.json`) | **PASS**: All empirical benchmarks use public, de-identified datasets (ClinicalTrials.gov, CDC NHANES, Dryad/PLOS ONE Rwanda HMIS recount) with frozen SHA-256 manifests. |
| **G15** | Related-work completeness & fairness | **RUN** | Software comparison table audit (`mice`, `simex`, `mecor`, `Amelia`) | **PASS**: Fairly credits existing CRAN packages and explicitly recommends simple linear MI (`use_recorded = TRUE`) when discrepancies are subtle/diffuse. |
| **G16** | Limitations & scope honesty | **RUN** | Section 5, Section 7 & `README.md` (`Limits` section) audit | **PASS**: Explicitly discloses where contamination model loses to linear MI (ClinicalTrials.gov, NHANES), where per-record error probabilities have low AUC (`0.543` / `0.500`), and where `catmech_reference(cluster = )` attenuates predictor slopes. |
| **G17** | Pre-committed analysis-plan honesty | **RUN** | Cross-check against `CONTAMINATION_PLAN.md`, `REAL_DATA_PLAN.md`, `VARIANCE_MI_PLAN.md` | **PASS**: Explicitly discloses in `manuscript.Rmd` and `paper.md` that `catmech_contamination()` passed 3 of 4 primary screen cells and narrowly missed `site_shift_declared` (`93.7%`, Wilson `[0.9202, 0.9505]`). |
| **G18** | Venue fit & AE-surface check | **RUN** | CRAN policy check + JSS/JOSS structure check | **PASS**: `R CMD check` passes with 0 errors/warnings/notes; `paper/manuscript.pdf` (JSS format) and `paper/paper.md` (JOSS format) include all required sections. |
| **G19** | Code-surface <-> paper-surface consistency | **RUN** | `verify_code_surface()` in `paper/verify_paper_numbers.py` | **PASS**: Caught and fixed 2 surface mismatches (`catmech_forest_plot()` -> `plot()`, and non-existent `fact_check.py` -> `verify-artifacts.py`, `verify-release-archive.py`, `verify_paper_numbers.py`). All 30 code/script surface references verified. |
| **G20** | Git & release-tag hygiene | **RUN** | `git status` + `.Rbuildignore` + `.gitignore` check | **PASS**: Clean working tree; build detritus excluded; all frozen archives intact. |
| **G21** | Supplement <-> main-text parity | **RUN** | Parity check across `manuscript.Rmd`, `paper.md`, `README.md`, and `vignettes/` | **PASS**: Synthetic, real-data, and clustered benchmark numbers match across `manuscript.Rmd`, `paper.md`, and `README.md`. |
| **G22** | Literature & empirical-calibration sourcing | **RUN** | `paper/verify_paper_numbers.py` + `tools/verify-artifacts.py` | **PASS**: All simulation and empirical numbers trace to frozen CSVs in `validation/`. |
| **G23** | Final pre-flight sign-off & basis stamping | **RUN** | `check_gate_basis.py --stamp` | **PASS**: Basis hashes stamped below. |

---

## Detailed Gate Audit Notes & `unchecked_directions` (Rule 12b)

### G0: Format-conversion pre-gate
- **Findings fixed**: Added blank lines before 4 markdown lists in `paper/manuscript.Rmd` so Pandoc emits `\begin{enumerate}` / `\begin{itemize}`; double-braced institutional authors in `paper/paper.bib`.
- **unchecked_directions**: HTML-to-Word conversion (not targeted; `.pdf`, `.tex`, and `.html` verified).

### G1: Path, hygiene & repository-surface guard
- **Findings**: `check_paths.sh` passed with 0 hits on `R/ src/ man/ vignettes/ inst/ tests/ tools/ paper/ README.md DESCRIPTION NAMESPACE NEWS.md`. Two historical scratch scripts in `validation/clustered-fix-scratch/` contain `/Users/patrickiben/Documents/catmech` and are excluded from the CRAN tarball via `.Rbuildignore` (`^validation$`), preserved untouched per Rule 3.
- **unchecked_directions**: Windows path separator runtime test (covered via GitHub Actions CI matrix).

### G2 & G22: Numeric-claim <-> artifact check & empirical sourcing
- **Findings fixed**: Replaced the illustrative synthetic table ($\beta = 0.500, N = 1,000$) in `paper/manuscript.Rmd` (Section 4.1) and `paper/paper.md` (Section "1. Synthetic Monte Carlo studies") with the exact frozen benchmark results from `validation/contamination-study/summary.csv` ($\beta = 0.800, N = 600$, 25% reference audit, $m = 20$, 52,950 planned requests), reporting `missed_flags` ($N_{\text{sim}} = 1,000$) and `error_gross_and_small` ($N_{\text{sim}} = 150$) with exact MCSEs. Created `paper/verify_paper_numbers.py` (343 automated assertions + negative control).
- **unchecked_directions**: Secondary stress cells in `validation/contamination-study/summary.csv` not tabulated in main text (all 27 scenarios are verified in `tools/verify-artifacts.py`).

### G3: Citation & bibliography verification
- **Findings fixed**: Running `check_citations.py` on both `paper/manuscript.bbl` and a full `\nocite{*}` build of `paper/paper.bib` caught 5 defective entries:
  1. `huang2013simple`: Fixed journal/title/DOI to *Bayesian Analysis* 8(2):439–452 (`10.1214/13-BA815`).
  2. `breslow2003semiparametric`: Fixed authors/title/journal/DOI to Breslow, McNeney & Wellner (2003), *The Annals of Statistics* 31(4):1110–1139 (`10.1214/aos/1059655907`).
  3. `neyman1938contribution`: Fixed JASA DOI to `10.1080/01621459.1938.10503378`.
  4. `little2019statistical`: Fixed Wiley book DOI to `10.1002/9781119482260`.
  5. `nab2021mecor`: Fixed *Computer Methods and Programs in Biomedicine* page/DOI to `106238` (`10.1016/j.cmpb.2021.106238`).
- **Result**: `20 OK | 5 REVIEW | 0 FLAG` across all 25 entries in `paper/paper.bib`.
- **unchecked_directions**: Page-level quotation check inside cited textbooks (`cochran1977sampling`, `rubin1987multiple`).

### G6: Figure, table & layout inspection
- **Findings fixed**: Stacked the author block across 3 centered lines; adjusted proportional pipe-table column widths in Sections 3, 4.1, and 5; narrowed the ASCII decision tree to 68 columns; achieved `errors=0 undefined=0 overfull=0 underfull=2` in `manuscript.log` and visually verified all 10 pages of `paper/manuscript.pdf`.
- **unchecked_directions**: Print grayscale contrast of HTML interactive output.

### G7 & G10: AI-artifact scrubber, prose style & US spelling
- **Findings fixed**: Converted all headings in `paper/manuscript.Rmd` and `paper/paper.md` to sentence case; removed all em-dashes (`—`) from prose; replaced "proves" with "demonstrates" in Section 7.
- **Result**: `check_prose_style_v2.py paper/manuscript.tex` exits `0` (`within thresholds on every implemented category`, `0` em-dashes, `0` curly quotes, `0` Title-Case headings); `scrub_ai_artifacts.py` exits `0` (`0` hits); `check_spelling_us.py` exits `0` (`0` British spellings).
- **unchecked_directions**: Wikipedia-specific categories listed as unimplemented by `check_prose_style_v2.py`.

### G17: Pre-committed analysis-plan honesty
- **Findings fixed**: Added explicit disclosure in Section 4.1 of `paper/manuscript.Rmd` and `paper/paper.md` that across the four primary 1,000-cohort cells (`clean`, `flagged_errors`, `missed_flags`, `site_shift_declared`), `catmech_contamination()` passed the pre-specified `[0.925, 0.975]` Wilson coverage screen in three cells (`96.5%`, `96.3%`, `94.9%`) and narrowly missed in `site_shift_declared` at `93.7%` (`[0.9202, 0.9505]`, bias `-0.0085`, MCSE `0.0016`).
- **unchecked_directions**: None.

### G19: Code-surface <-> paper-surface consistency
- **Findings fixed**: Automated check in `paper/verify_paper_numbers.py` caught `catmech_forest_plot()` (not an exported function; forest plots use S3 `plot()`) and `fact_check.py` (temporary script name replaced by committed `verify-artifacts.py`, `verify-release-archive.py`, and `verify_paper_numbers.py`).
- **unchecked_directions**: Internal non-exported C++ helper names (not referenced in manuscript).
FAIL: paper/manuscript.tex/paper/manuscript.pdf does not exist
