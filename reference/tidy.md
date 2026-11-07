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
                         predictors = "age", m = 2, seed = 1)
pooled <- catmech_pool_lm(res, y ~ lab + age)
tidy(pooled)
#>          term    estimate  std.error  statistic       df      p.value
#> 1 (Intercept)  5.48654328 0.97808156  5.6094947 25.86205 6.906928e-06
#> 2         lab  0.01301963 0.07946479  0.1638415 15.12639 8.720251e-01
#> 3         age -0.01208767 0.01447075 -0.8353171 35.09559 4.091900e-01
#>      conf.low  conf.high
#> 1  3.47554580 7.49754077
#> 2 -0.15623237 0.18227163
#> 3 -0.04146199 0.01728666
```
