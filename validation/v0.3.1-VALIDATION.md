# Validation record: catmech 0.3.1

Date: 2026-10-31. This is an experimental software release with independent
implementation checks and prespecified synthetic studies. It is not clinical
validation, a general coverage guarantee, or a claim to new imputation methodology.

## Main findings

The [extension synthesis](EXTENSION_SYNTHESIS.md) reports the full new study.
The [four-panel evidence figure](extension-figures/extension-evidence.pdf) shows
the main benefits, equal-information comparisons and failure conditions;
[figure provenance and interpretation](extension-figures/README.md) accompany it.
With a random, accurate 25% reference subset, default reference-normal MI covered
95.2% of nominal 95% intervals (MCSE 0.68 percentage points). The old flag-masking
workflow covered 58.1% under missed errors and 12.7% under the compound site/error
stress. The reference method receives additional measured information; it does
not identify missed errors from flags alone. Reference-only analysis covered
95.3%, with mean interval width .3304 versus .2782 for reference-normal MI.
Independent standard normal MI performed similarly.

The study also exposed failures that remain visible:

- Partial-reference PMM covered only 85.1% with random 25% references and exactly
  matched its independent mice comparator. The release requires explicit sensitivity
  opt-in for its pooled analysis. It remains an experimental completion option.
- PMM improved heavy-tail coverage at equal information (.940 versus normal .917),
  but worsened strong-MAR coverage (.911 versus normal .949). It is not a general
  replacement for the normal default.
- Reference-normal coverage fell to .055 under nonignorable reference selection
  and .678 when the supplied reference had measurement error. Gold-standard status
  and selection are assumptions, not something the API can establish.
- Five of six Gaussian no-category normal settings passed the prespecified
  Wilson-containment screen. Missing-only narrowly failed at .940, Wilson
  [.9235,.9531]; the independent normal comparator gave the same result.

All five primary clinical recording scenarios reuse the same cohorts, reference
subsets and imputation seeds. Their reference results are identical by design.
This is one paired 1,000-cohort reference experiment, not five independent
replications of its coverage. Smaller reference-fraction stresses have 250 cohorts
and explicitly wider Monte Carlo uncertainty.

## Software and independent checks

| Check | Result | Evidence |
|---|---|---|
| Unit suite | 318 assertions passed, no warnings or skips | `tests.log` |
| Source build and installed package check | 0 errors, 0 warnings, 0 notes | `build.log`, `R-CMD-check.log` |
| Fresh archive installation | Namespace-only APIs, all three examples, reference/PMM exports, mids alignment and pooling guards | `installed-smoke.log`, `tools/smoke-installed.R` |
| Gaussian category posterior oracle | Four moments within 1.19 MCSE of independent quadrature; closed-form and wrong-likelihood controls passed | `model-oracle-final/README.md`, frozen v0.2 source |
| Reference joint predictive oracle | 14 mean/covariance checks from 12,000 draws, largest discrepancy 1.91 MCSE; 12 contract controls passed | `reference-oracle/README.md` |
| Reference negative controls | Incorrect zero cross-record covariance and plug-in variance detected | `reference-oracle/negative-controls.csv` |
| Pooling patch equivalence | 54 supported scalar old/new comparisons bit-identical; 15 guard controls passed | `pooling-guards/README.md` |
| Full/no-reference endpoints | Ordinary-lm targets/estimate/SE/df/CI agree to 2.66e-15; 40 expected refusals when no references exist | `extension-controls/replicates.csv` |
| Engineering pilot | 805 requests; serial/parallel scientific files byte-identical; only expected no-reference refusals | `extension-pilot-check.txt` |
| Exact PMM comparator | 23,500 new/independent PMM numerical pairs agree exactly | `extension-evidence-check.txt` |
| Independent numeric correspondence | Python recomputes summaries, Wilson intervals, paired MCSEs and decomposition identities; requires every planned cell and seed | `artifact-verification.log`, `tools/verify-artifacts.py` |
| Artifact negative controls | Corrupted bias and deletion of an entire planned scenario are detected | `artifact-verification.log` |
| Scientific figure | All 21 plotted rows, Monte Carlo intervals and source hashes verified; PNG and rendered PDF inspected | `extension-figures/visual-QA.txt` |
| Original supplied prototype | Original hashes unchanged | `original-sha256.csv`, `artifact-verification.log` |
| Dependency record | Installed dependency versions retained; no new dependency added in v0.3 | `renv.lock` |
| GitHub Actions | Three-platform workflow supplied; not run remotely | `.github/workflows/R-CMD-check.yaml` |

The local environment is R 4.6.0, mice 3.19.0 and posterior 1.7.0 on Apple silicon
macOS. Exact session records accompany each study. Minimum-version R, other
platforms and a fresh-container dependency restore have not been tested locally.
The independent oracles check known examples and data contracts; they do not
establish model adequacy or repeated-sampling coverage in an arbitrary application.

## Frozen studies and the v0.3.1 guards

The [new protocol](EXTENSION_PLAN.md) preceded all new runs, with streams disjoint
from prior studies and engineering pilots. The completed study has **132,580
requests: 132,540 intervals, 40 expected no-reference refusals, zero other errors**:

| Bundle | Requests | Scope |
|---|---:|---|
| `extension-distribution/` | 80,000 | Normal/PMM by category-information factorial, Gaussian, heavy-tail, null and strong-MAR settings |
| `extension-reference/` | 35,000 | Random 25% accurate reference subsets, paired recording-error scenarios |
| `extension-stress/` | 17,500 | 10%/50% reference fractions, nonignorable selection and noisy references |
| `extension-controls/` | 80 | All-reference equality and no-reference refusal |

The numerical study used frozen v0.3.0 sources from commit `ea5ec6a`, preserved
under `v0.3.0-source/`. Every run checked source hashes before and after execution.
Version 0.3.1 subsequently added two scoped safeguards: explicit rejection of
unsupported multivariate/empty-coefficient analyses, and sensitivity-only pooling
for partial-reference PMM after its demonstrated failure. Neither change alters
imputation values, scalar estimates or interval calculations when explicitly
requested. The 54 paired scalar comparisons verify that claim, and the reference
oracle was rerun against the released sources. The adverse study results were
not retuned or rerun to obtain favorable estimates.

Historical v0.1 and v0.2 results retain their original source snapshots and logs.
The [v0.2 validation record](v0.2.0-VALIDATION.md) and
[earlier coverage synthesis](v0.2.0-COVERAGE_SYNTHESIS.md) remain available. Source
hash verification uses the appropriate frozen version rather than substituting
current code into historical fingerprints. Large held-out studies were run once;
serial/parallel exact reproducibility was checked on independent engineering pilots.

Coverage is conditional on a returned interval. Refusals have no confidence
interval and do not count as successful 95% inference. Report request yields and
failures separately. Retrospective bias-centered coverage uses hidden truth and
is diagnostic, not a real-data correction. The primary Wilson-containment criterion
applies only to prespecified 1,000-replicate primary cells; descriptive stresses
and smaller studies are not forced to pass it.

## Evidence limits

Reference values must measure the same quantity, time and units and be sufficiently
accurate for the stated model. Targeted verification, unsupported strata, nonlinear
or heteroskedastic conditional models and unknown reference error remain important
limits. Normal MI is not a general solution for survival models, interactions or
cluster-correlated outcomes. Category-error matrices remain fixed assumptions;
uncertainty in them is not automatically propagated by pooling separate scenarios.
[Reference method](../docs/REFERENCE_METHOD.md) and [Gaussian method](../docs/METHOD.md)
state the assumptions and established alternatives.

The author's research-validation skill and partially read master protocol informed
these adapted checks. The complete named gate system was not executed or certified.
Adversarial agent review is not independent human refereeing, author sign-off,
clinical adjudication, a comprehensive retraction audit or a validation of the
harness itself. AI-assisted work is disclosed in `CONTRIBUTING.md`.

## Reproduce

```sh
Rscript -e 'testthat::test_local()'
R CMD build .
R CMD check --no-manual catmech_0.3.1.tar.gz
Rscript tools/check-reference-oracle.R validation/my-reference-oracle
Rscript tools/check-pooling-guards.R validation/my-pooling-guards
Rscript tools/extension-benchmark.R distribution validation/my-distribution
Rscript tools/extension-benchmark.R reference validation/my-reference-study
Rscript tools/extension-benchmark.R stress validation/my-stress
Rscript tools/extension-benchmark.R controls validation/my-controls
Rscript tools/plot-extension.R validation/my-extension-figures
python3 tools/verify-artifacts.py
```

Use new output directories. To reproduce the exact historical source fingerprint,
run the benchmark script from its frozen source directory. The verifier checks
published bundle locations against their recorded source versions. Current scripts
produce fresh runs whose environment/version fingerprints may differ.
`CATMECH_BENCHMARK_WORKERS=1` uses the verified sequential path; Windows uses that
path automatically. Exact cross-platform numerical equality is not promised.

Dependency restoration is available through
`renv::restore(lockfile="validation/renv.lock", prompt=FALSE)` in a disposable
checkout with compatible R. Empty-container restoration was not independently
exercised. Installable source and repository archives have separate release checksums.
