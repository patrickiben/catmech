# Write a self-contained result bundle

Write a self-contained result bundle

## Usage

``` r
write_catmech(result, directory)
```

## Arguments

- result:

  A catmech_result.

- directory:

  New or empty output directory. Existing files are never overwritten.

## Value

The normalized output path, invisibly. Files include a lossless RDS,
long CSV, per-cell audit, diagnostics, JSON manifest and session
information.

## Examples

``` r
set.seed(123)
df <- data.frame(y = stats::rnorm(20), x = stats::rnorm(20))
df$ref <- df$y; df$ref[6:20] <- NA
res <- catmech_reference(df, value = "y", reference = "ref", predictors = "x", m = 2, seed = 1)
tmp_dir <- file.path(tempdir(), "test_catmech_bundle")
write_catmech(res, tmp_dir)
unlink(tmp_dir, recursive = TRUE)
```
