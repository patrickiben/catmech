# Validation record: catmech 0.8.0

Date: 2026-11-03. This is an experimental research-software release. It adds an explicit joint clustered model and validates its implementation separately from its inferential performance. It is not a new statistical methodology, clinical validation, or a general coverage guarantee. The complete [v0.7.2 record](v0.7.2-VALIDATION.md) preserves the earlier studies and release checks, with links to their preceding versions.

## New joint model

`catmech_reference(..., use_recorded = TRUE, cluster = "site", cluster_model = "joint")` models reference truth and recorded measurements jointly with a bivariate site effect. It supports recorded-present and recorded-missing unaudited rows in the same site, and predicts sites with neither measurement observed. Gold reference values remain exact and raw input is preserved. The default `cluster_model = "moment"` retains its supported behavior and mixed-pattern guard.

The joint model uses the established `pan` sampler, declared inverse-Wishart covariance priors on empirically standardized responses and a flat fixed-effect prior. Conditional prediction uses the original observed components to draw one shared site vector per completion. See the [model and limitations](../docs/JOINT_CLUSTER_METHOD.md).

The [prospective study plan](JOINT_CLUSTER_PLAN.md) was committed before execution (`33ef227`). A reproducibility amendment (`2a4880a`) preceded the first pilot; the initial implementation was frozen at `554b0c9`. An execution-only worker-limit correction and its amendment were frozen at `4467ca1` before the final pilot and held-out run. The full final study consists of 4,200 requests: six scenarios with 200 replications for three main methods, plus 600 prespecified prior-sensitivity requests. Its [synthesis](JOINT_CLUSTER_SYNTHESIS.md) reports coverage including failed requests, conditional coverage, Monte Carlo uncertainty, bias, error, width and paired comparisons. Earlier Rwanda, rare-count and predictor-attenuation results concern different models and do not validate this option.

## Scientific findings

The full grid contains 4,200 requests: 4,178 returned intervals and 22 refused convergence. The default covariance prior with `m = 40` returned 1,185 of 1,200 fits. Returned-interval coverage across the six scenarios ranged from 96.5% to 99.5%; balanced Gaussian coverage was conservative. Coverage per request, counting refusals as noncoverage, fell to 93.0% in the t5 stress case with eight refusals. Average slope biases were small in these cohorts. Width and coverage changed materially under the prespecified prior sensitivity checks, so the results do not establish exact 95% calibration or prior-insensitive inference.

Against the same-information plain recorded-value comparator, the joint option had lower observed error and narrower intervals in most scenarios, but not near zero cluster variance. Some summaries condition on successful fits; paired comparisons and their Monte Carlo errors are supplied in the synthesis. The oracle uses additional information and is a check on the analysis procedure, not a fair accuracy competitor. The study supports the implementation as an opt-in working model, not uniform superiority or repair of the historical rare-count failures. It used 40 completions and does not directly validate the public API's default of 20.

## Numerical and reproducibility checks

The new test file has 140 passing assertions. A separate full multivariate-normal conditioning calculation checks conditional means and every covariance entry against 10,000 simulated completions. It includes both recording patterns in one site, gold-only anchors, W-only information, a completely unobserved site, shuffled rows and a negative control that detects losing the cross-pattern site covariance. An independent direct `pan` replay checks adapter inputs and parameter draws; it does not independently validate `pan`'s mathematics.

Other checks cover exact gold preservation, unchanged raw input, `mids` alignment, same-seed reproducibility after unrelated `pan` use, caller RNG restoration, all-reference and all-agree endpoints, support refusals and failed convergence with retained diagnostic evidence. The test file passes under `_R_CHECK_LIMIT_CORES_=TRUE`.

Two execution defects were caught before release. `pan` 2.0 retains a cached normal variate that its public seed does not reset; each chain now runs in a fresh R worker. The first CRAN-style package check then rejected four simultaneous workers. Chains now run in batches of at most two, with fresh processes for every chain and retry. The complete development fit was bit-identical before and after batching. The initial pilot and 350 saved requests from the interrupted first study are retained in `joint-cluster-pilot/` and `joint-cluster-study/`; the latter has an `ABORTED.json` marker and no completed-study claim. Accuracy was not inspected to choose the repair or settings. The final run uses the original prespecified seeds, with an exact comparison to those 350 saved requests.

## Release verification

All 886 package assertions passed without warnings or skips. The final `R CMD check --as-cran --no-manual` finished with zero errors, zero warnings and one routine CRAN incoming “New submission” note. Current logs are in `v0.8.0-checks/`. The final study verifier independently recomputes the request grid, seeds, scalar summaries, Wilson intervals, Monte Carlo standard errors, paired differences, diagnostic thresholds and source hashes. Negative controls must reject a missing cell or altered summary. Passing software checks and convergence screens does not establish calibration outside the tested scenarios.

| Check | Evidence |
|---|---|
| Full package tests and CRAN-style source check | `v0.8.0-checks/R-CMD-check.log`, `v0.8.0-checks/tests.log` |
| Fresh archive, namespace API smoke and six bundled examples | `v0.8.0-checks/install.log`, `v0.8.0-checks/installed-smoke.log` |
| Exact compatibility: 18 v0.7.2, 18 v0.3.1 and 12 v0.4.0 workflows | `v0.8.0-checks/v0.7.2-unchanged-check.txt`, `v0.8.0-checks/v0.3-unchanged-check.txt`, `v0.8.0-checks/v0.4.0-postfreeze-check.txt` |
| Historical saved studies and learning-page checks | `v0.8.0-checks/artifact-verification.log` |
| README claims and dated historical case study | `v0.8.0-checks/fact-check.log`, `../docs/case-study/FACT_CHECK.md` |
| Archive contents and checksum | `v0.8.0-checks/build-integrity.txt`, `v0.8.0-checks/source-package-sha256.txt` |

The frozen v0.7.2 package sources are in `v0.7.2-source/`, copied from public commit `b52e31b` for compatibility testing. Other historical source manifests and all earlier numerical results remain unchanged. The dated case-study HTML/PDF retains its November 2, v0.7.1 status. Current GitHub Actions checks run on Ubuntu, macOS and Windows; see the [workflow results](https://github.com/patrickiben/catmech/actions).

## Reproduction and limits

Run `Rscript validation/run-joint-cluster.R heldout validation/my-joint-study` from this release's source, then `python3 -I tools/verify-joint-cluster.py validation/my-joint-study`. Use a fresh output directory. The final saved run can additionally be checked with `--compare-aborted validation/joint-cluster-study`. The initial pilot and aborted run belong to source commit `554b0c9`; their original hashes are intentionally not rewritten to match the repaired source.

The new study targets the audited value as a regression outcome, with exact references and broad paired support. It does not validate nonignorable missingness, imperfect gold standards, concentrated paired references, count models, random slopes, influential-site inference or use of the audited value as an analysis predictor. AI agents performed implementation and review; no independent human statistician has reviewed this release.
