# Generate an auditable compliance and integrity trail for a catmech analysis

Extracts complete provenance, input and completion data hashes, model
specifications, convergence monitors, imputation statistics, and
downstream pooled analyses into a tamper-evident audit document.
Supports Markdown, JSON, and standalone HTML formats.

## Usage

``` r
catmech_audit_trail(
  object,
  pooled = NULL,
  sensitivity = NULL,
  format = c("markdown", "json", "html"),
  path = NULL,
  title = "Catmech Audit Trail & Compliance Report"
)
```

## Arguments

- object:

  A \`catmech_result\` object, or a list containing a \`result\` and
  optional downstream objects (\`pooled\`, \`sensitivity\`,
  \`audit_plan\`, \`holdout_check\`).

- pooled:

  Optional pooled analysis object (\`catmech_pooled_lm\`,
  \`catmech_pooled_glm\`, or \`catmech_pooled_coxph\`).

- sensitivity:

  Optional \`catmech_sensitivity\` object.

- format:

  Output format: \`"markdown"\` (default), \`"json"\`, or \`"html"\`.

- path:

  Optional file path to write the report to. If \`NULL\`, the character
  string is returned invisibly.

- title:

  Document title (default: \`"Catmech Audit Trail & Compliance
  Report"\`).

## Value

A character vector of class \`c("catmech_audit_trail", "character")\`
containing the generated audit document.

## Examples

``` r
set.seed(123)
df <- data.frame(
  lab = stats::rnorm(30, 10, 2),
  age = stats::rnorm(30, 50, 10),
  y = stats::rnorm(30, 5, 1)
)
df$ref <- df$lab
df$ref[11:30] <- NA
res <- catmech_reference(df, value = "lab", reference = "ref",
                         predictors = "age", m = 2, seed = 1)
pooled <- catmech_pool_lm(res, y ~ lab + age)
trail <- catmech_audit_trail(res, pooled = pooled)
print(trail)
#> # Catmech Audit Trail & Compliance Report
#> 
#> **Generated:** 2026-11-07T03:04:44Z UTC | **Catmech Version:** 1.0.0
#> **R Version:** 4.6.2 | **Platform:** x86_64-pc-linux-gnu
#> 
#> ---
#> 
#> ## 1. Cryptographic Data & Specification Integrity (SHA-256)
#> 
#> - **Original Input Data Hash:** `b0dbfaef433a8845a59cd2740963aae87c551da520912d26b69d94ebe24f9fe9`
#> - **Assumption Specification Hash:** `7ba7a2acef9b95325ceacfeaf8413bd993aaeb6afe1e7c9de47e1e0c33f7a75d`
#> - **Completed Dataset Checksums:**
#> 
#>   - Completion `m=1`: `fd9e0211aa706a46d0e63e085181c831144480ec10b12580688cf67fbaa0dfa1`
#>   - Completion `m=2`: `a833b76a6fae53ebab5588a3691804afd77f7583d7b151fc50167c287e726f42`
#> 
#> ## 2. Cohort Specification & Imputation Profile
#> 
#> | Parameter | Value |
#> | :--- | :--- |
#> | Engine | `reference_normal` |
#> | Target Value Column | `lab` |
#> | Reference Column | `ref` |
#> | Predictors | age |
#> | Cluster Variable | `none` |
#> | Total Observations | 30 |
#> | Completions (m) | 2 |
#> | Imputed Records | 20 (66.67%) |
#> | Convergence Status | Not Applicable |
#> 
#> ## 3. Downstream Pooled Inference
#> 
#> **Model Class:** `catmech_pooled_lm`
#> 
#> | term | estimate | std.error | df | lower | upper | interpretation |
#> | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
#> | (Intercept) |  7.6465 | 1.2881 |  9.1610 |  4.7404 | 10.5526 | conditional_reference_model |
#> | lab |  0.2246 | 0.1256 |  7.4031 | -0.0692 |  0.5184 | conditional_reference_model |
#> | age | -0.0915 | 0.0355 | 24.9164 | -0.1647 | -0.0184 | conditional_reference_model |
#> 
#> ---
#> 
#> ### Regulatory Audit & Compliance Certification
#> > This document provides an auditable, verifiable record of data imputation and downstream analysis.
#> > All completions and statistical intervals are cryptographically tied to the recorded input hash and assumption specification.
#>  
```
