# Plot holdout check calibration diagnostics

Displays a 2-panel diagnostic plot of predictive calibration: 1.
Uniformity histogram of randomized Probability Integral Transform (PIT)
values. 2. Scatter plot of predictive means versus observed reference
values.

## Usage

``` r
# S3 method for class 'catmech_holdout_check'
plot(x, y, ...)
```

## Arguments

- x:

  A \`catmech_holdout_check\` object.

- y:

  Unused.

- ...:

  Additional arguments passed to \[graphics::par()\].

## Value

Invisible \`NULL\`.
