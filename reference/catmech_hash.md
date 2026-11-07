# Compute a canonical assumption hash

Named maps are sorted, sequence order is retained, and numeric values
use full round-trip precision. Missing values, signed infinities and
strings have distinct type tags. No numerical rounding or null-dropping
is applied.

## Usage

``` r
catmech_hash(spec)
```

## Arguments

- spec:

  A list, usually a validated catmech_spec.

## Value

A SHA-256 hex string. This hashes assumptions, not the input dataset.

## Examples

``` r
sp <- catmech_spec(
  value = "lab", category = "stage",
  levels = c("low", "med", "high"),
  thresholds = c(-Inf, 5, 10, Inf)
)
catmech_hash(sp)
#> [1] "e90472a0424908caef0771842c035cea95b798b7345d7bb34038f12c0d9db9f1"
```
