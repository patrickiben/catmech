# Simulate measurements, quality flags, and interval categories

Synthetic data only. The coefficient of lab_true in the outcome model is
0.8. Flags occur among observed values; a configurable fraction indicate
actual corruption. category_error_rate deliberately violates the
package's exact-category assumption. Truth fields must never enter
imputation predictors.

## Usage

``` r
simulate_catmech_lab(
  n = 300L,
  seed = 1L,
  corruption_probability = 1,
  missing_rate = 0.2,
  flag_rate = 0.1,
  category_error_rate = 0,
  error_distribution = c("normal", "heavy_tailed"),
  thresholds = c(-Inf, 5.5, 7, Inf)
)
```

## Arguments

- n:

  Number of rows, at least 20.

- seed:

  Nonnegative integer seed, with caller RNG restored.

- corruption_probability:

  Probability of corruption conditional on a flag.

- missing_rate:

  Probability the continuous value is missing (MCAR).

- flag_rate:

  Flag probability among observed measurements.

- category_error_rate:

  Probability of deliberately miscoding a category.

- error_distribution:

  "normal" or variance-standardized "heavy_tailed" Student t residuals
  with 3 degrees of freedom for the latent measurement.

- thresholds:

  Four interval endpoints for Low, Medium, High.

## Value

A data.frame with lab_true and corrupted for validation only.

## Examples

``` r
dat <- simulate_catmech_lab(n = 50, seed = 123)
head(dat)
#>   patient_id      age sex inflammation lab_value lab_category outlier_flag
#> 1          1 44.39524   F    3.0255714  6.080598       Medium           No
#> 2          2 47.69823   M    1.7152270        NA       Medium           No
#> 3          3 65.58708   M    0.7792823  6.746956       Medium           No
#> 4          4 50.70508   F    2.1813035 14.770935       Medium          Yes
#> 5          5 51.29288   M    1.8611086  6.423890       Medium           No
#> 6          6 67.15065   M    2.0057642  8.804965         High           No
#>    outcome lab_true corrupted
#> 1 3.197003 6.080598     FALSE
#> 2 5.918825 6.477911     FALSE
#> 3 5.783591 6.746956     FALSE
#> 4 4.265588 5.664049      TRUE
#> 5 5.257257 6.423890     FALSE
#> 6 7.178011 8.804965     FALSE
```
