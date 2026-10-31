# Validation record

Date: 2026-10-31. Package: catmech 0.1.0. Environment: R 4.6.0, mice 3.19.0,
macOS 27.0.1 on Apple silicon. This is a scoped software release and exploratory
statistical assessment. It is not a completed manuscript submission certification.

## Verified release behavior

| Check | Result | Evidence |
|---|---|---|
| Original supplied prototype | Reproduced in a temporary copy; original files unchanged | `original-sha256.csv`, `tools/verify-artifacts.py --original PATH` |
| Original reported smoke result | Reproduced coefficient 0.57390237 and its 24 assertions; does not rescue the identified counterexamples | `docs/ASSESSMENT.md` |
| Package source build | Successful | `R CMD build .` |
| Installed package checks | 0 errors, 0 warnings, 0 notes; 113 test assertions passed | `R-CMD-check.log`, `tests.log` |
| Fresh temporary installation | Namespace-only calls, completion adapter, pooling, export and supplied example succeeded | `installed-smoke.log`, `tools/smoke-installed.R` |
| Numeric core | Exact endpoint membership, singleton donor selection, posterior oracle, type-1 matching, ties, rank/support checks, sparse support, q endpoints and auditing covered | `tests/testthat/test-imputer.R`, `tests/testthat/test-spec.R` |
| Workflow | Originals unchanged, seeded repeatability, RNG restoration, YAML normalization, distinct infinity hashes, named-vector hashing, installed mice integration and output refusal covered | `tests/testthat/test-api.R`, `tests/testthat/test-spec.R` |
| Unit transformation | Positive affine target/threshold transformation preserves donor choices; wrong unchanged thresholds rejected | `tests/testthat/test-api.R` |
| Repeated simulation | 5 scenarios x 100 datasets x 6 methods; all 3,000 method fits succeeded | `benchmark/replicates.csv`, `benchmark/RESULTS.md` |
| Exact local reproduction | Six scientific/report artifacts byte-identical in a second clean R process | `benchmark/reproducibility.txt` |
| Independent numeric correspondence | Python recomputes saved bias, RMSE, coverage, MCSE and widths from replicate rows; detects a deliberately corrupted bias | `tools/verify-artifacts.py`, `artifact-verification.log` |
| Baseline audit | No clear general advantage over flag-aware mice PMM; stress cases show serious failure | `benchmark/paired_comparisons.csv`, `docs/ASSESSMENT.md` |
| Dependency record | Exact installed dependency versions captured | `renv.lock`; restoration has not been tested in an empty external environment |
| GitHub CI | Three-platform workflow supplied; remote execution not yet performed | `.github/workflows/R-CMD-check.yaml` |

The 100 simulation warning records all concern ordinary PMM removing the constant
all-No flag in the missing-only scenario. They are retained in the raw results.
The package checks and test suite had no warnings.

The normal fallback has boundary tests but no comparative coverage validation.
Neither package checks nor reproducibility establish proper MI or scientific
superiority. The uncertain-flag and category-error results are substantive
limitations, not software-check failures to be hidden.

## Scoped use of the research-validation harness

The validation skill and a partially read personal master protocol informed these
checks. The located master protocol is more extensive than the skill summary.
It was not read end-to-end; no full named gate is represented as passed. The table
records adapted evidence and explicit gaps. Human author sign-off is not inferred.

| Gate or topic | Scope/status here | Evidence or remaining work |
|---|---|---|
| G0 estimand | Scoped definition and pre-run plan | `BENCHMARK_PLAN.md`; no author novelty attestation |
| G1 reproducibility | Partial | Built-package install and clean-process repeats; no fresh container or restored lockfile run |
| G2 numeric grounding | Adapted checks completed | `tools/verify-artifacts.py`, MCSEs and raw replicates |
| G3 citations | Targeted primary-source review | `docs/REFERENCES.bib`; no comprehensive retraction-service audit |
| G4 robustness | Exploratory only | Five scenarios; failed assumptions reported; no full calibration study |
| G5 adversarial review | Multiple agent lenses plus source/probe verification | Novelty, statistics and implementation reviews; not independent human refereeing |
| G6 reporting | ADEMP-style plan provided | `BENCHMARK_PLAN.md`; no journal checklist certification |
| G7 AI/authorship | Disclosure supplied | `CONTRIBUTING.md`; author approval remains a human action |
| G8 strong alternatives | Matched PMM baseline run | No clear gain; hmi not executed, as documented in assessment |
| G9 claim/artifact correspondence | Adapted numeric checks completed | Source hashes, raw rows and independent recomputation; not a universal grounding proof |
| G10 prose | Detector within thresholds and editorial read | No human writing-style certification |
| G11 dimensions | Applicable unit-change invariant checked | Target/threshold affine transformation and a known-wrong control |
| G12 read-through | Agent review completed | README, method and assessment inspected; no signed author manuscript freeze |
| G13 external contact | Prior-art research completed within stated scope | Primary articles, package docs and archived source; no human cold-reader study |
| G14 teaching transfer | Not applicable | No claimed educational companion or learning study |
| G15 artifact integrity | Adapted source/build checks completed | Source hashes, extracted package checks and release checksums |
| G16 text/artifact correspondence | Partial | Generated benchmark report and recomputed numbers; no frozen manuscript/PDF |
| G17 derivative manuscript | Not applicable | Software revision, no manuscript split |
| G18 venue fit | Not applicable | GitHub preparation, no journal submission |
| G19 pretrained provenance | Not applicable to runtime method | No pretrained model used for imputation |
| G20 external review | Original QC report reassessed | Confirmed errors and corrections in `docs/ASSESSMENT.md` |
| G21 harness effectiveness | Not evaluated | No comparative study of the harness itself |
| G22 ability to fail | Adapted negative controls | Bad donor sampling, wrong matching, unsupported predictor level and corrupted result tests |
| G23 premise/scope | Explicit assumptions | `docs/METHOD.md`; no chain of formal computed bounds claimed |

## Reproduction commands

```sh
Rscript -e 'testthat::test_local()'
R CMD build .
R CMD check --no-manual catmech_0.1.0.tar.gz
Rscript tools/benchmark.R 100 10 validation/new-run
python3 tools/verify-artifacts.py
```

`validation/new-run` must be new or empty. To provision the recorded package
versions with renv, install renv and use
`renv::restore(lockfile = "validation/renv.lock", prompt = FALSE)` in a disposable
checkout, with a compatible R version. This restoration route is provided from the
saved lockfile; it was not exercised in a fresh container during this assessment.
