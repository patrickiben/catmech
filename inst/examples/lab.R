library(catmech)
dat <- simulate_catmech_lab(n = 300, seed = 1)
dat <- dat[setdiff(names(dat), c("lab_true", "corrupted"))]
spec <- read_catmech_spec(system.file("specs", "lab.yaml", package = "catmech"))
predictors <- c("age", "sex", "inflammation", "outcome")
fit <- catmech_impute(dat, spec, predictors, m = 20, seed = 123)
print(fit$diagnostics)
print(summary(mice::pool(with(fit$mids,
  lm(outcome ~ lab_value + age + sex + inflammation)))))

# Compare explicitly chosen sensitivity policies on the same synthetic input.
grid <- do.call(rbind, lapply(c(0, 0.5, 0.8, 1), function(q) {
  s <- spec
  s$corruption_probability <- q
  f <- catmech_impute(dat, s, predictors, m = 20, seed = 123)
  pooled <- summary(mice::pool(with(f$mids,
    lm(outcome ~ lab_value + age + sex + inflammation))))
  row <- pooled[pooled$term == "lab_value", c("estimate", "std.error")]
  data.frame(q = q, row, mean_reconstructed = mean(f$diagnostics$reconstructed))
}))
print(grid)

# To export, use a new or empty directory:
# write_catmech(fit, "output/example")
