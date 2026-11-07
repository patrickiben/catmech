# Plot an audit power curve

Displays statistical power versus audit sample size with reference
threshold lines.

## Usage

``` r
# S3 method for class 'catmech_audit_power'
plot(x, col = "#1f77b4", lwd = 2, ...)
```

## Arguments

- x:

  A \`catmech_audit_power\` object.

- col:

  Line color for the power curve (default \`"#1f77b4"\`).

- lwd:

  Line width (default 2).

- ...:

  Additional arguments passed to \[graphics::plot()\].

## Value

\`x\`, invisibly.

## Examples

``` r
pwr <- catmech_audit_power(N = 500, target_ratio = 1.3, n_sim = 200, seed = 1)
plot(pwr)
```
