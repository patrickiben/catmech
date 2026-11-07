# Declare assumptions for interval-constrained sensitivity completion

The interval convention is \[a,b\] for the first category and (a,b\] for
subsequent categories. Finite endpoints may define a restricted domain.
The probability is a sensitivity weight for replacing flagged, in-bound
observed values, not an estimated posterior corruption probability.

## Usage

``` r
catmech_spec(
  value,
  category,
  levels,
  thresholds,
  flag = NULL,
  flag_yes = "Yes",
  flag_no = "No",
  corruption_probability = NULL,
  conflict_policy = c("error", "reconstruct"),
  min_donors = 1L,
  n_donors = 5L,
  fallback = c("error", "normal")
)
```

## Arguments

- value, category:

  Names of the continuous measurement and interval category.

- levels:

  Ordered, unique category labels.

- thresholds:

  Strictly increasing numeric endpoints, length levels + 1.

- flag:

  Optional quality-flag column name; NULL means no flags.

- flag_yes, flag_no:

  Labels for suspect and non-suspect flags. Missing flags do not mark an
  observed value as suspect; this is an explicit assumption.

- corruption_probability:

  Fixed scalar replacement probability in \[0,1\]. Must be supplied
  explicitly when a flag column is configured. With no flag, NULL is
  normalized to zero.

- conflict_policy:

  Either "error" or "reconstruct". The latter declares the category
  authoritative for all observed category/value conflicts.

- min_donors:

  Minimum eligible donors before using the fallback policy.

- n_donors:

  Number of nearest eligible donors from which to sample.

- fallback:

  "error" (default) or "normal" for an explicit, experimental
  truncated-normal fallback when too few category donors exist.

## Value

A validated catmech_spec list.

## Examples

``` r
spec <- catmech_spec("value", "band", c("low", "high"), c(-Inf, 0, Inf))
spec$conflict_policy
#> [1] "error"
```
