#!/usr/bin/env Rscript
# Design calibration only: within-cluster correlation of the imputation model's residuals on very large draws. No method is run.
a <- list(b0 = 2, b1 = .6, size_log_mean = 3, size_log_sd = .5, x_noise_sd = 3, cluster_variance = 20, row_variance = 80, stratum_sd = 2)
DESIGNS <- list(rwanda = c(14L, 9L, 8L, 13L, 5L, 9L, 11L), pairs = rep(8L, 10L), none = 70L)
icc_anova <- function(e, id) {
  id <- factor(id); k <- tapply(e, id, length); g <- nlevels(id); n <- length(e)
  means <- tapply(e, id, mean)
  between <- sum(k * (means - mean(e))^2) / (g - 1); within <- sum((e - means[id])^2) / (n - g)
  k0 <- (n - sum(k^2) / n) / (g - 1)
  (between - within) / (between + (k0 - 1) * within)
}
draw <- function(errors, design, scale) {
  sizes <- DESIGNS[[design]] * scale; rows <- 3L; g <- sum(sizes); n <- g * rows
  id <- rep(seq_len(g), each = rows); stratum_of <- rep(seq_along(sizes), sizes); row_stratum <- stratum_of[id]
  size <- rep(rlnorm(g, a$size_log_mean, a$size_log_sd), each = rows)
  x <- size + rnorm(n, 0, a$x_noise_sd)
  if (errors == "rare") {
    level <- if (length(sizes) > 1L) runif(length(sizes), 0, .2)[row_stratum] else 0
    frailty <- rep(rgamma(g, .5, .5), each = rows)
    truth <- rpois(n, pmax(.01, .1 + .02 * x + level) * frailty)
    slip <- runif(n) < .2
    recorded <- pmax(0, truth + slip * sample(c(-1, 1), n, replace = TRUE))
  } else {
    level <- if (length(sizes) > 1L) rnorm(length(sizes), 0, a$stratum_sd)[row_stratum] else 0
    if (errors == "row_only") {
      truth <- a$b0 + a$b1 * x + level + rnorm(n, 0, sqrt(a$cluster_variance + a$row_variance))
    } else {
      truth <- a$b0 + a$b1 * x + level + rep(rnorm(g, 0, sqrt(a$cluster_variance)), each = rows) + rnorm(n, 0, sqrt(a$row_variance))
    }
    recorded <- switch(errors,
      shift_low = truth + rep(rnorm(g, 0, 3), each = rows) + rnorm(n, 0, 6),
      shift = truth + rep(rnorm(g, 0, 5), each = rows) + rnorm(n, 0, 3),
      under = (1 - rep(rbeta(g, 2, 6), each = rows)) * truth + rnorm(n, 0, 2),
      outlier = truth + rep((runif(g) < .1) * rnorm(g, 0, 20), each = rows) + rnorm(n, 0, 3),
      row_only = truth + rnorm(n, 0, 6))
  }
  d <- data.frame(truth = truth, recorded = recorded, x = x, s = factor(row_stratum), id = id)
  fit <- if (length(sizes) > 1L) lm(truth ~ recorded + x + s, data = d) else lm(truth ~ recorded + x, data = d)
  c(icc = icc_anova(residuals(fit), id), mean_truth = mean(truth), wrong = mean(abs(truth - recorded) > 1e-8), slope_ref = unname(coef(lm(truth ~ x + s, data = d))["x"]), slope_rec = unname(coef(lm(recorded ~ x + s, data = d))["x"]))
}
set.seed(20261102)
for (errors in c("shift_low", "shift", "under", "outlier", "rare", "row_only")) for (design in c("rwanda", "pairs")) {
  scale <- if (design == "rwanda") 1500L else 1300L
  v <- draw(errors, design, scale)
  cat(sprintf("%-10s %-7s residual ICC %6.3f; mean target %6.2f; wrong rows %.2f; slope, reference %.4f, recorded %.4f\n", errors, design, v[["icc"]], v[["mean_truth"]], v[["wrong"]], v[["slope_ref"]], v[["slope_rec"]]))
}
