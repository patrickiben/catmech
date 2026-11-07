# Complete-Companion Coverage Map

Basis: Patrick Iben, *Statistical Inference and Two-Phase Audit Design for Contaminated Measurements in Clinical Research: The catmech R Package*, canonical manuscript [paper/manuscript.Rmd](../paper/manuscript.Rmd) (`350` lines, SHA-256 `8152288b1c2af2a49103543db3c9324582177b1e842365b5c5ab83cdac8c1fd0`) and 15-page compiled PDF [paper/manuscript.pdf](../paper/manuscript.pdf) (SHA-256 `1d17480886d4f38c38e7cf5f900ea3c5f1ef48f3c4484fa801a72a83de05e0db`). Unit IDs refer to [COMPANION.md](COMPANION.md) ([PDF](COMPANION.pdf)). The [concise study guide](STUDY_GUIDE.md) ([PDF](STUDY_GUIDE.pdf)) provides a twelve-lesson route (`S1`-`S12`).

Every single line (`1-350`) of `paper/manuscript.Rmd` is mapped to a unique unit below with zero gaps and zero overlaps, alongside its paired executable R check in [run_companion_labs.R](run_companion_labs.R) and interactive web surface.

| Unit | Manuscript lines | PDF pages | Evidence class | Source passage / R function / benchmark covered |
| :--- | ---: | ---: | :--- | :--- |
| **C01** | 1-29 | 1 | `SCOPE` | Clinical audit designs, two-phase inference, and the simulation-to-practice synthesis (Title, YAML header, structured abstract, and keywords) |
| **C02** | 30-41 | 1-2 | `DEFINITION` | The three-way clinical audit dilemma: naive attenuation, complete-case variance inflation, and discard-recorded loss (Section 1 Introduction) |
| **C03** | 42-56 | 2 | `DEFINITION` | Data generating and contamination framework: latent truth, corruption indicator, recorded value, and Phase II reference (Section 2.1, Equations (1)-(3)) |
| **C04** | 57-60 | 2-3 | `ANALYTIC` | Priors and collapsed slice-within-Gibbs MCMC architecture in catmech_contamination (Section 2.2 opening paragraph) |
| **C05** | 61-64 | 3 | `ANALYTIC` | Marginalized corruption probability for unaudited records via Gaussian convolution (Section 2.2, Equation (4)) |
| **C06** | 65-70 | 3 | `ANALYTIC` | Conjugate Gaussian posterior for contaminated unaudited records and reference anchoring (Section 2.2, Equations (5)-(6)) |
| **C07** | 71-72 | 3-4 | `DEFINITION` | Gaussian reference imputation in C++17: use_recorded = TRUE (Simple Linear MI) versus FALSE (Section 2.2 closing paragraph) |
| **C08** | 73-80 | 4 | `ANALYTIC` | Two-phase influence function variance decomposition and the interval-width precision ratio (Section 2.3, Equations (7)-(8)) |
| **C09** | 81-84 | 4 | `MEASURED_DATED` | Monte Carlo audit power planning and the pre-specified NO-GO boundary on catmech_audit_plan (Section 2.3, Equation (9) and validation disclosure) |
| **C10** | 85-90 | 4-5 | `ANALYTIC` | Rubin multiple imputation combining rules and the fraction of missing information (Section 2.4, Equations (10)-(11)) |
| **C11** | 91-94 | 5 | `ANALYTIC` | Barnard and Rubin small-sample degrees of freedom adjustment for finite samples, clusters, and events (Section 2.4, Equation (12)) |
| **C12** | 95-98 | 5 | `ANALYTIC` | Downstream pooling for continuous (CV3 jackknife) and generalized linear (CR1S sandwich) endpoints (Section 2.4, Endpoint items 1 and 2, Equation (13)) |
| **C13** | 99-100 | 5-6 | `ANALYTIC` | Time-to-event Cox proportional hazards pooling with event-based degrees of freedom and Lin-Wei sandwich (Section 2.4, Endpoint item 3) |
| **C14** | 101-106 | 6 | `ANALYTIC` | Out-of-sample cross-validated predictive calibration checking via rank PIT and its slope-bias limit (Section 2.5, Equation (14)) |
| **C15** | 107-118 | 6 | `DEFINITION` | Modular S3 software architecture across planning, diagnostics, imputation, pooling, and provenance (Section 3 opening and module table) |
| **C16** | 119-133 | 6-7 | `DEFINITION` | Comparison with mice, simex, mecor, and Amelia, plus broom generics and strategy comparison (Section 3 capability comparison table and closing paragraph) |
| **C17** | 134-173 | 7-8 | `DEFINITION` | Reproducible 8-stage clinical trial workflow (Stages 1-5): simulation, power, screening, mixture fit, and calibration (Section 4 code chunk (Stages 1-5)) |
| **C18** | 174-210 | 8 | `MODELED_LOCKED` | Reproducible 8-stage clinical trial workflow (Stages 6-8): multi-endpoint pooling, strategy comparison, and provenance (Section 4 code chunk (Stages 6-8) and console output) |
| **C19** | 211-218 | 8-9 | `MODELED_LOCKED` | Design of the 52,950-request synthetic contamination study across 27 scenarios and 7 strategies (Section 5 opening and Section 5.1 introduction) |
| **C20** | 219-240 | 9-10 | `MODELED_LOCKED` | Synthetic benchmark results in missed_flags and error_gross_and_small, plus the four primary cells (Section 5.1 benchmark table and primary cell analysis) |
| **C21** | 241-249 | 10 | `MODELED_LOCKED` | Pre-specified stress-cell failure boundaries and holdout calibration alarm power (Section 5.1 stress-cell failure inventory) |
| **C22** | 250-258 | 10-11 | `MEASURED_DATED` | Clustered variance study (148,500 requests): why CV3 beats CR1S/CR2 and the predictor-attenuation boundary (Section 5.2) |
| **C23** | 259-277 | 11-12 | `MEASURED_DATED` | Public paired-discrepancy benchmarks on ClinicalTrials.gov and NHANES (33,800 requests) (Section 6.1 and empirical benchmark table) |
| **C24** | 278-279 | 12 | `MEASURED_DATED` | Why per-record corruption probabilities fail on diffuse real-world discrepancies (AUC 0.543 and 0.500) (Section 6.1 closing paragraph) |
| **C25** | 280-289 | 12-13 | `MEASURED_DATED` | Real-audit slope recovery benchmark (45,600 requests) on Rwanda HMIS and psychology recounts (Section 6.2) |
| **C26** | 290-314 | 13 | `INTERPRETIVE` | Methodological decision framework: choosing between Simple Linear MI and Contamination Mixture (Section 7 decision tree and rules) |
| **C27** | 315-328 | 13-14 | `INTERPRETIVE` | Discussion: simulation-to-practice divergence and why catmech centers on a diagnostic workflow (Section 8 Discussion) |
| **C28** | 329-337 | 14 | `SCOPE` | Regulatory documentation and data governance: what catmech_audit_trail records and where its scope stops (Section 9 Regulatory documentation and data governance) |
| **C29** | 338-345 | 14 | `MEASURED_DATED` | Quality assurance and verification: 1,754 unit tests, R CMD check, and 500,000+ replicate checks (Section 10 Quality assurance and verification) |
| **C30** | 346-350 | 14-15 | `DEFINITION` | Software availability, MIT license, R/C++17 dependencies, and bibliographic references (Section 11 Availability and Section 12 References) |

## Multimodal & Executable R Cross-Reference

- **Document Modality (`DOCUMENT`):** [COMPANION.md](COMPANION.md), [COMPANION.pdf](COMPANION.pdf), [STUDY_GUIDE.md](STUDY_GUIDE.md), [STUDY_GUIDE.pdf](STUDY_GUIDE.pdf)
- **One-Pager Modality (`ONE_PAGER`):** [ONE_PAGER.pdf](ONE_PAGER.pdf) and [catmech-case-study.pdf](../docs/case-study/catmech-case-study.pdf)
- **Interactive Modality (`INTERACTIVE`):** [companion_lab.html](companion_lab.html), [catmech Guide & Method Picker](../docs/learn/index.html), [Live Linear-MI Simulator](../docs/learn/simulator.html), and [Record Audit Planner](../docs/audit-planner/index.html)
- **Executable R Verification & Console Tutor:** [run_companion_labs.R](run_companion_labs.R) (`16/16` live checks) and [interactive_tutor.R](interactive_tutor.R) (`S1`-`S12` console tutor)
