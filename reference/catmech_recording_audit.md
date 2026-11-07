# Summarize what reference measurements reveal about recording errors

A descriptive first step for any data with both recorded and reference
values. It fits no model, so it can be read before a modeling assumption
is chosen.

## Usage

``` r
catmech_recording_audit(data, value, reference, by = NULL, agreement_tol = 0)
```

## Arguments

- data:

  A data frame.

- value:

  Name of the numeric recorded measurement.

- reference:

  Name of the numeric reference measurement.

- by:

  Optional categorical columns to compare, such as flag or site.

- agreement_tol:

  Recorded values within this distance of the reference count as
  correct. See \[catmech_contamination()\].

## Value

A list with \`status\`, \`overall\`, \`by\`, \`heterogeneity\`,
\`selection\` and \`notes\`. Rates use Jeffreys intervals.
\`heterogeneity\` holds exact tests of whether the disagreement rate
differs across each \`by\` column; it is a screen for unmodeled
structure, not a confirmatory test. \`selection\` compares recorded
values on referenced and unreferenced rows; a large standardized
difference shows the reference sample is selective, and agreement cannot
show the reverse.

## Examples

``` r
set.seed(123)
df <- data.frame(
  recorded = c(10, 12, 14, 100, 15, 16),
  ref = c(10, 12, 14, 10, NA, NA)
)
audit <- catmech_recording_audit(df, value = "recorded", reference = "ref", agreement_tol = 1)
print(audit)
#> -- Catmech Recording Audit --
#> Status:             few_errors_observed
#> Disagreements:      1 / 4 calibration rows (25.0%)
#> Reference subset:   4 / 6 total rows (66.7%)
#> Selection SMD:      0.5941
#> 
#> Notes:
#>   * Fewer than five disagreements were observed, so the error distribution is weakly determined.
```
