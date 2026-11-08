# Generate a provenance record for a catmech analysis

Collects the R and package versions, SHA-256 hashes of the input data,
the assumption specification and each completed dataset, the imputation
settings, the convergence status and any pooled or sensitivity results
into one document, in Markdown, JSON or standalone HTML.

## Usage

``` r
catmech_audit_trail(
  object,
  pooled = NULL,
  sensitivity = NULL,
  format = c("markdown", "json", "html"),
  path = NULL,
  title = "Catmech Audit Trail & Provenance Report"
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

  Document title (default: \`"Catmech Audit Trail & Provenance
  Report"\`).

## Value

A character vector of class \`c("catmech_audit_trail", "character")\`
containing the generated audit document.

## Details

The input-data and completion hashes are taken by \[digest::digest()\]
over the serialized R data frames, so they can be recomputed only in R,
not with a file checksum tool. The specification hash is of its
canonical JSON text, which \[write_catmech()\] saves as
\`assumptions.canonical.json\`. The hashes let a rerun be compared with
the record; they do not make the document tamper-evident, because anyone
who edits the data can also regenerate the record. The record can
support documentation, for example under FDA 21 CFR Part 11 or ICH E9,
but compliance is a property of the validated system, access controls
and procedures around it, not of this function.

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
                         predictors = c("age", "y"), m = 2, seed = 1)
pooled <- catmech_pool_lm(res, y ~ lab + age)
trail <- catmech_audit_trail(res, pooled = pooled)
print(trail)
#> # Catmech Audit Trail & Provenance Report
#> 
#> **Generated:** 2026-11-08T05:40:18Z UTC | **Catmech Version:** 1.0.1
#> **R Version:** 4.6.2 | **Platform:** x86_64-pc-linux-gnu
#> 
#> ---
#> 
#> ## 1. Data and Specification Hashes (SHA-256)
#> 
#> - **Original Input Data Hash:** `b0dbfaef433a8845a59cd2740963aae87c551da520912d26b69d94ebe24f9fe9`
#> - **Assumption Specification Hash:** `0d01fc948329eb1374ac3e7a605d2d654b5029a26059bf334c5523074d21396d`
#> - **Completed Dataset Checksums:**
#> 
#>   - Completion `m=1`: `66e25c1df2d231868f37c76cb6c2dc8ff51e6c6e4c2d15e86faedbaa713958f0`
#>   - Completion `m=2`: `502dfc01798a93b727b4356ea205b8e8e44040933c474441e2fa40dd5fc5032f`
#> 
#> ## 2. Cohort Specification & Imputation Profile
#> 
#> | Parameter | Value |
#> | :--- | :--- |
#> | Engine | `reference_normal` |
#> | Target Value Column | `lab` |
#> | Reference Column | `ref` |
#> | Predictors | age, y |
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
#> | (Intercept) |  8.8500 | 4.0865 | 0.7515 | -107.8258 | 125.5258 | conditional_reference_model |
#> | lab |  0.0865 | 0.3030 | 0.7883 |   -7.3156 |   7.4885 | conditional_reference_model |
#> | age | -0.0873 | 0.1282 | 0.6015 |   -8.6900 |   8.5154 | conditional_reference_model |
#> 
#> ---
#> 
#> ### Audit Trail & Provenance Record
#> > This document records the data, settings and results of one imputation and analysis run.
#> > The data and completion hashes are of the serialized R data frames, and the specification hash is of its canonical JSON text. A rerun can be compared with them. They do not cover the pooled results, and they do not make this document tamper-evident.
#>  
```
