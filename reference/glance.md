# Summarize a pooled catmech model with glance statistics

Summarize a pooled catmech model with glance statistics

## Usage

``` r
glance(x, ...)
```

## Arguments

- x:

  An object to be summarized.

- ...:

  Additional arguments.

## Value

A one-row data.frame.

## Examples

``` r
set.seed(123)
n <- 40
df <- data.frame(
  lab = stats::rnorm(n, 10, 2),
  age = stats::rnorm(n, 50, 10),
  y = stats::rnorm(n, 5, 1)
)
df$lab_ref <- df$lab
df$lab_ref[11:n] <- NA
res <- catmech_reference(df, value = "lab", reference = "lab_ref",
                         predictors = "age", m = 2, seed = 1)
pooled <- catmech_pool_lm(res, y ~ lab + age)
glance(pooled)
#>   nobs m
#> 1   40 2
```
