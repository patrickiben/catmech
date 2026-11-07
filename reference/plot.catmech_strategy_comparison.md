# Plot audit strategy comparison (Forest Plot)

Draws a comparative forest plot showing point estimates and confidence
intervals across all evaluated audit strategies for a chosen coefficient
term.

## Usage

``` r
# S3 method for class 'catmech_strategy_comparison'
plot(
  x,
  term = NULL,
  colors = c(validation_subsample = "#898781", naive_trust_recorded = "#d03b3b",
    linear_mi_recorded = "#2a78d6", contamination_model = "#eb6834"),
  lwd = 2,
  pch = 19,
  cex = 1.3,
  main = NULL,
  xlab = NULL,
  ...
)
```

## Arguments

- x:

  A \`catmech_strategy_comparison\` object.

- term:

  Character: which coefficient term to plot. If \`NULL\`, the first
  non-intercept term is plotted.

- colors:

  Named vector of colors for the strategies.

- lwd:

  Line width (default 2).

- pch:

  Point symbol (default 19).

- cex:

  Point size (default 1.3).

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
n <- 50
df <- data.frame(
  outcome = stats::rnorm(n, 10, 2),
  lab = stats::rnorm(n, 5, 1),
  age = stats::rnorm(n, 50, 10)
)
df$ref <- df$lab + stats::rnorm(n, 0, 0.2)
df$ref[16:n] <- NA
comp <- catmech_compare_strategies(
  df, value = "lab", reference = "ref",
  formula = outcome ~ lab + age, m = 2, seed = 1
)
plot(comp, term = "lab")
```
