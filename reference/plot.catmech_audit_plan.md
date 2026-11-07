# Plot a catmech audit plan

Displays the planned precision curve showing relative interval width
across audit sampling fractions, bootstrap confidence bands, and the
target precision ratio.

## Usage

``` r
# S3 method for class 'catmech_audit_plan'
plot(x, y, col_curve = "royalblue", col_target = "firebrick", ...)
```

## Arguments

- x:

  A \`catmech_audit_plan\` object.

- y:

  Unused.

- col_curve:

  Line color for relative width curve (default \`"royalblue"\`).

- col_target:

  Line color for target ratio threshold (default \`"firebrick"\`).

- ...:

  Additional arguments passed to \[graphics::par()\].

## Value

Invisible \`NULL\`.

## Examples

``` r
set.seed(123)
n <- 120
df <- data.frame(
  outcome = stats::rnorm(n, 10, 2),
  lab = stats::rnorm(n, 5, 1),
  age = stats::rnorm(n, 50, 10)
)
df$ref <- df$lab + stats::rnorm(n, 0, 0.2)
df$ref[41:n] <- NA
plan <- catmech_audit_plan(df, value = "lab", reference = "ref",
                           predictors = c("age", "outcome"), formula = outcome ~ lab + age,
                           term = "lab", bootstrap = 100L, m = 2L, seed = 1L)
#> Warning: Small pilot: k rests on 40 rows, so its interval is wide. Below 100 rows, plan on the upper limit.
plot(plan)
```
