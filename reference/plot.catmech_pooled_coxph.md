# Plot pooled Cox PH hazard ratios (Forest Plot)

Plot pooled Cox PH hazard ratios (Forest Plot)

## Usage

``` r
# S3 method for class 'catmech_pooled_coxph'
plot(
  x,
  terms = NULL,
  intercept = FALSE,
  exponentiate = NULL,
  col_point = "#1f77b4",
  col_ci = "#1f77b4",
  col_ref = "darkred",
  lwd = 2,
  pch = 19,
  cex = 1.2,
  main = NULL,
  xlab = NULL,
  ...
)
```

## Arguments

- x:

  A \`catmech_pooled_coxph\` object.

- terms:

  Character vector of terms to plot (default \`NULL\` for all).

- intercept:

  Logical: include intercept (default \`FALSE\`).

- exponentiate:

  Logical: exponentiate estimates (default \`TRUE\` if available).

- col_point:

  Point color (default \`"#1f77b4"\`).

- col_ci:

  CI color (default \`"#1f77b4"\`).

- col_ref:

  Null reference line color (default \`"darkred"\`).

- lwd:

  Line width (default 2).

- pch:

  Point symbol (default 19).

- cex:

  Point size (default 1.2).

- main:

  Plot title.

- xlab:

  X-axis label.

- ...:

  Additional arguments passed to \[graphics::plot()\].

## Value

\`x\`, invisibly.

## Examples

``` r
set.seed(123)
n <- 40
df <- data.frame(
  time = stats::rexp(n, 0.1), status = stats::rbinom(n, 1, 0.5),
  lab = stats::rnorm(n, 10, 2), age = stats::rnorm(n, 50, 10)
)
df$lab_ref <- df$lab
df$lab_ref[11:n] <- NA
res <- catmech_reference(df, value = "lab", reference = "lab_ref",
                         predictors = "age", m = 2, seed = 1)
pooled <- catmech_pool_coxph(res, survival::Surv(time, status) ~ lab + age)
plot(pooled)
```
