# Integration map from the supplied v0.9 draft

The supplied `catmech-local` directory was inspected and reproduced in a temporary
copy. Its files were not edited. This repository is the proposed replacement
software package. It contains no patient data; examples and benchmarks are synthetic.
`validation/original-sha256.csv` records hashes of the original files reviewed.

| Draft component | Replacement | Reason |
|---|---|---|
| Root runner with fixed lab names | `catmech_impute()` and `inst/examples/lab.R` | Reusable one-target API; explicit predictors; preserve originals |
| Permissive nested YAML | `catmech_spec()`, `read_catmech_spec()`, `inst/specs/lab.yaml` | Reject unsupported/inert assumptions; no implied joint or contamination model |
| Conflicting interval conventions | `R/spec.R` | One right-closed convention, including exact thresholds |
| Hash canonicalizer | `catmech_hash()` | Distinguish infinity signs and missing values; preserve numerical precision |
| Donor selection with `sample(index,1)` | `R/imputer.R` | Fix singleton-index bug; random ties; enforce support |
| Broadening to suspect fit rows | Strict original trusted sample | Refuse unsupported reconstruction |
| Generic pseudoinverse draw | Rank-aware QR and recipient estimability checks | Avoid invented uncertainty in unidentified directions |
| Bernoulli corruption narrative | Explicit sensitivity replacement rule | No fitted corruption likelihood or posterior is claimed |
| Flag-only conflict checks | All-row validation and explicit conflict policy | Detect contradictions in unflagged observations too |
| Ambiguous changed-rate summary | Per-cell decisions plus separate changes | A reconstruction need not change a value |
| Saved `violation_rate` actually measuring compliance | Count of verified threshold violations | Correct diagnostic meaning |
| Single-seed smoke comparison | Repeated paired benchmark with stronger baselines | Measure bias, RMSE, coverage and Monte Carlo uncertainty |
| Source-only tests | Standard package tests, source build, R CMD check and CI | Test installed behavior and reproducible entry points |
| Broad speculative roadmap | Narrow experimental scope | No unimplemented topology, causal intervention or MNAR features |

The initial refined package version was `0.1.0`; the old `0.9` label described a prototype and
is not a claim of package maturity. No manuscript source was supplied, so there
are no proposed paper-text replacements. The new assessment and method document
are the proposed publication-facing explanation.

## Version 0.2 strengthening

| Addition | Reason |
|---|---|
| `catmech_model()` | Replace arbitrary retention with a declared Gaussian model; use optional exact/noisy category likelihoods |
| Trusted values remain unchanged | Category contradictions cannot silently override numeric truth assumptions |
| R-hat and bulk/tail ESS refusal | Unresolved MCMC does not become a successful analysis |
| `catmech_pool_lm()` | Explicit pooling checks and separate legacy-sensitivity interpretation |
| `catmech_sensitivity()` | Same-data comparisons, planned-scenario completeness and retained failures |
| Independent posterior quadrature | Verify one full sampler posterior without reusing Gibbs code |
| Held-out coverage and clinical fixtures | Strong baselines, residual-bias diagnosis and explicit remaining undercoverage |

The public package version is now `0.2.0`. Legacy API behavior remains available;
new model assumptions require the new entry point. See [NEWS](../NEWS.md).

## Version 0.3 extension

| Addition | Reason |
|---|---|
| `catmech_reference()` | Treat independently measured reference values as anchors; impute every unvalidated true target |
| `catmech_pmm()` | Expose ordinary flag-masking PMM without legacy retention or forced categories |
| Honest donor identity audit | Distinguish an identified row from duplicate-value ambiguity; separate available donors from nearest-pool size |
| Fully referenced analysis endpoint | Return ordinary regression intervals when explicit metadata confirms no imputation |
| Predictive-t oracle | Check shared posterior parameter uncertainty through cross-record covariance |
| Equal-information study | Separate PMM versus normal from the effect of additional category predictors |
| Reference selection/noise stresses | Expose assumptions that additional measurements still cannot bypass |

Version `0.3.0` retains the old engines and evidence. It adds no automatic
corruption detector and no estimated category-confusion likelihood.
