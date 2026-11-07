# Plot a marginal prior sensitivity grid

Displays a 2-panel diagnostic plot showing posterior coefficient
estimates and credible intervals across covariance prior scales and
regression prior scales, compared against baseline estimates and Monte
Carlo error bands.

## Usage

``` r
# S3 method for class 'catmech_marginal_sensitivity_grid'
plot(x, y, col_points = "navy", col_base = "firebrick", ...)
```

## Arguments

- x:

  A \`catmech_marginal_sensitivity_grid\` object.

- y:

  Unused.

- col_points:

  Color for scenario points (default \`"navy"\`).

- col_base:

  Color for baseline reference line (default \`"firebrick"\`).

- ...:

  Additional arguments passed to \[graphics::par()\].

## Value

Invisible \`NULL\`.
