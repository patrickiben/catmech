# Turn a pooled catmech model into a tidy data frame

Turn a pooled catmech model into a tidy data frame

## Usage

``` r
tidy(x, ...)
```

## Arguments

- x:

  An object to be tidied.

- ...:

  Additional arguments.

## Value

A data.frame with standardized model columns.

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
                         predictors = c("age", "y"), m = 2, seed = 1)
pooled <- catmech_pool_lm(res, y ~ lab + age)
tidy(pooled)
#>          term     estimate  std.error  statistic       df     p.value
#> 1 (Intercept)  4.955481037 1.37993451  3.5910987 10.37974 0.004638164
#> 2         lab  0.049906412 0.08013362  0.6227900 15.47335 0.542493120
#> 3         age -0.008163162 0.01731801 -0.4713683 15.73525 0.643846948
#>      conf.low  conf.high
#> 1  1.89597482 8.01498725
#> 2 -0.12044039 0.22025321
#> 3 -0.04492597 0.02859965
```
