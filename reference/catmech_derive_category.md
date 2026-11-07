# Derive exact interval categories

Derive exact interval categories

## Usage

``` r
catmech_derive_category(x, thresholds, levels)
```

## Arguments

- x:

  Numeric values. Missing values remain missing.

- thresholds, levels:

  Endpoints and category labels as in \[catmech_spec()\].

## Value

A factor with the specified levels. Out-of-domain values give NA.

## Examples

``` r
vals <- c(3.2, 5.5, 7.8, 12.0)
catmech_derive_category(
  vals,
  thresholds = c(-Inf, 5, 10, Inf),
  levels = c("low", "med", "high")
)
#> [1] low  med  med  high
#> Levels: low med high
```
