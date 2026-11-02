#!/usr/bin/env Rscript
# Independent numerical checks of catmech_reference(cluster = ) and catmech_pool_lm(cluster = ) on fixed datasets:
# lme4 and mitools where they apply, and hand-written matrix code that shares nothing with the package otherwise.
# Exact checks use a tolerance of 1e-6 (fitted models) or 1e-9 (algebra). Sampling checks use five Monte Carlo errors.
# Usage: Rscript tools/check-cluster-algebra.R [OUTPUT_FILE]
args <- commandArgs(trailingOnly = TRUE)
file_arg <- grep("^--file=", commandArgs(), value = TRUE)
root <- normalizePath(file.path(dirname(sub("^--file=", "", file_arg[1])), ".."))
output <- if (length(args)) args[1] else file.path(root, "validation", "cluster-algebra-check.txt")
pkgload::load_all(root, quiet = TRUE)
ns <- asNamespace("catmech")
moments <- get("catmech_reference_cluster_moments", ns)
gls_fit <- get("catmech_reference_cluster_gls", ns)
fill <- get("catmech_reference_cluster_fill", ns)
cr1s <- get("catmech_cluster_vcov", ns)
lines <- character(); failures <- 0L
report <- function(name, value, limit) {
  ok <- is.finite(value) && value <= limit
  if (!ok) failures <<- failures + 1L
  lines <<- c(lines, sprintf("%-4s %-82s largest %.3g, limit %.3g", if (ok) "PASS" else "FAIL", name, value, limit))
}

make_data <- function(sizes, seed, tau = 2, sigma = 1) {
  set.seed(seed)
  G <- length(sizes); id <- rep(seq_len(G), sizes); n <- length(id)
  x <- stats::rlnorm(G, 1, .5)[id] + stats::rnorm(n, 0, .5)
  w <- x + stats::rnorm(n)
  y <- 1 + .5 * x + .4 * w + stats::rnorm(G, 0, tau)[id] + stats::rnorm(n, 0, sigma)
  list(X = cbind(`(Intercept)` = 1, x = x, w = w), y = y, id = id, cl = factor(as.character(id)), n = n, G = G)
}
DATA <- list(balanced = make_data(rep(4L, 15L), 101L), unequal = make_data(rep(c(1L, 2L, 3L, 5L, 8L), 4L), 102L),
             three_rows = make_data(rep(3L, 17L), 103L), low_icc = make_data(rep(3L, 25L), 104L, tau = .2))

# Variance components in a balanced one-way layout, against the ANOVA formula and REML.
set.seed(105L)
each <- 4L; G <- 15L; id <- rep(seq_len(G), each = each); y0 <- stats::rnorm(G, 0, 2)[id] + stats::rnorm(G * each)
cl0 <- factor(as.character(id)); got <- moments(matrix(1, G * each, 1L), y0, cl0)
means <- tapply(y0, id, mean)
ms_b <- each * sum((means - mean(y0))^2) / (G - 1L); ms_w <- sum((y0 - means[id])^2) / (G * each - G)
tau_hand <- (ms_b - ms_w) / each
tau_moment <- got$df_between * (got$ss_between / got$df_between - got$sse_within / got$df_within) / got$trace
report("Moments, balanced one-way: within mean square against ANOVA", abs(got$sse_within / got$df_within - ms_w), 1e-9)
report("Moments, balanced one-way: intercept variance against ANOVA", abs(tau_moment - tau_hand), 1e-9)
reml <- lme4::lmer(y0 ~ 1 + (1 | id), REML = TRUE)
report("Moments, balanced one-way: intercept variance against lme4 REML", abs(tau_moment - as.data.frame(lme4::VarCorr(reml))$vcov[1L]), 1e-4)
report("Moments, balanced one-way: residual variance against lme4 REML", abs(ms_w - stats::sigma(reml)^2), 1e-4)

for (name in names(DATA)) {
  d <- DATA[[name]]
  m <- moments(d$X, d$y, d$cl)
  D <- stats::model.matrix(~ cl - 1, data = list(cl = d$cl))
  M <- diag(d$n) - d$X %*% solve(crossprod(d$X), t(d$X))
  report(sprintf("%s: moment trace against tr(Z'MZ) by explicit matrices", name), abs(m$trace - sum(diag(t(D) %*% M %*% D))), 1e-8)

  fit <- lme4::lmer(d$y ~ d$X[, "x"] + d$X[, "w"] + (1 | d$id), REML = TRUE)
  comps <- as.data.frame(lme4::VarCorr(fit))
  tau2 <- comps$vcov[1L]; sigma2 <- stats::sigma(fit)^2
  if (tau2 < 1e-6) { lines <- c(lines, sprintf("NOTE %s: lme4 estimates a zero intercept variance, so its GLS checks reduce to OLS", name)) }
  gls <- gls_fit(d$X, d$y, d$cl, sigma2, tau2)
  report(sprintf("%s: GLS coefficients at REML components against lme4 fixed effects", name), max(abs(unname(gls$coefficients) - unname(lme4::fixef(fit)))), 1e-5)
  covariance <- sigma2 * tcrossprod(backsolve(qr.R(gls$qr), diag(3L)))
  report(sprintf("%s: GLS covariance at REML components against lme4 vcov", name), max(abs(covariance - as.matrix(stats::vcov(fit)))), 1e-5)
  V <- sigma2 * diag(d$n) + tau2 * tcrossprod(D); Vi <- solve(V)
  report(sprintf("%s: GLS coefficients against explicit inverse of V", name), max(abs(unname(gls$coefficients) - as.vector(solve(t(d$X) %*% Vi %*% d$X, t(d$X) %*% Vi %*% d$y)))), 1e-8)

  # Conditional intercepts, as the fill computes them, against lme4's conditional modes and variances.
  sizes <- tabulate(d$cl)
  residual_mean <- as.vector(rowsum(d$y - as.vector(d$X %*% gls$coefficients), d$cl, reorder = TRUE)) / sizes
  blup <- sizes * tau2 / (sigma2 + sizes * tau2) * residual_mean
  ran <- lme4::ranef(fit, condVar = TRUE)$`d$id`
  order_ids <- match(levels(d$cl), rownames(ran))
  report(sprintf("%s: conditional intercept means against lme4 ranef", name), max(abs(blup - ran[order_ids, 1L])), 1e-5)
  conditional <- tau2 * sigma2 / (sigma2 + sizes * tau2)
  report(sprintf("%s: conditional intercept variances against lme4 condVar", name), max(abs(conditional - attr(ran, "postVar")[1L, 1L, order_ids])), 1e-5)

  # The package's cluster-robust variance, against a hand-written double sum.
  frame <- data.frame(y = d$y, x = d$X[, "x"], w = d$X[, "w"], id = d$id)
  model <- stats::lm(y ~ x + w, data = frame)
  e <- stats::residuals(model); B <- solve(crossprod(d$X)); total <- matrix(0, 3L, 3L)
  for (g in unique(d$id)) for (i in which(d$id == g)) for (j in which(d$id == g)) total <- total + e[i] * e[j] * tcrossprod(d$X[i, ], d$X[j, ])
  by_hand <- B %*% total %*% B * (d$G / (d$G - 1)) * ((d$n - 1) / (d$n - 3))
  report(sprintf("%s: CR1S variance against a hand-written double sum", name), max(abs(cr1s(model, d$id) - by_hand)), 1e-8)
}

# The fill's draws against a separate implementation with explicit V matrices, sharing no package code.
hand_draws <- function(X, y, cl, x_new, cl_new, reps, seed) {
  set.seed(seed)
  n <- nrow(X); p <- ncol(X); Z <- stats::model.matrix(~ cl - 1, data = list(cl = cl)); levels_all <- levels(cl)
  P_x <- X %*% solve(crossprod(X), t(X)); XZ <- cbind(X, Z); q <- qr(XZ); rank_all <- q$rank
  P_all <- qr.Q(q)[, seq_len(rank_all), drop = FALSE]; P_all <- tcrossprod(P_all)
  sse_w <- sum(((diag(n) - P_all) %*% y)^2); sse_x <- sum(((diag(n) - P_x) %*% y)^2)
  df_w <- n - rank_all; df_b <- rank_all - p; ss_b <- sse_x - sse_w
  tr <- sum(diag(t(Z) %*% (diag(n) - P_x) %*% Z))
  out <- matrix(NA_real_, nrow(x_new), reps)
  units <- unique(cl_new)
  for (r in seq_len(reps)) {
    s2 <- sse_w / stats::rchisq(1L, df_w)
    t2 <- max(0, df_b * (ss_b / stats::rchisq(1L, df_b) - s2) / tr)
    Vi <- solve(s2 * diag(n) + t2 * tcrossprod(Z))
    cov_b <- solve(t(X) %*% Vi %*% X)
    b <- as.vector(cov_b %*% t(X) %*% Vi %*% y) + as.vector(t(chol(cov_b)) %*% stats::rnorm(p))
    effect <- stats::rnorm(length(units), 0, sqrt(t2)); names(effect) <- units
    for (u in units[units %in% levels_all]) {
      rows <- cl == u; k <- sum(rows)
      effect[u] <- k * t2 / (s2 + k * t2) * mean(y[rows] - X[rows, , drop = FALSE] %*% b) + stats::rnorm(1L, 0, sqrt(t2 * s2 / (s2 + k * t2)))
    }
    out[, r] <- as.vector(x_new %*% b) + unname(effect[cl_new]) + stats::rnorm(nrow(x_new), 0, sqrt(s2))
  }
  out
}
compare_draws <- function(name, data, audited_clusters, keep_partly) {
  audited <- data$id %in% audited_clusters
  fit_rows <- which(audited & !(data$id == keep_partly & duplicated(data$id, fromLast = TRUE) == FALSE))
  fill_rows <- setdiff(seq_len(data$n), fit_rows)
  reps <- 6000L
  set.seed(301L)
  package <- fill(data$X, data$y, fit_rows, fill_rows, as.character(data$id), reps)
  hand <- hand_draws(data$X[fit_rows, , drop = FALSE], data$y[fit_rows], factor(as.character(data$id[fit_rows])),
                     data$X[fill_rows, , drop = FALSE], as.character(data$id[fill_rows]), reps, 302L)
  mean_gap <- (rowMeans(package) - rowMeans(hand)) / sqrt((apply(package, 1L, stats::var) + apply(hand, 1L, stats::var)) / reps)
  report(sprintf("%s: mean draw per row against the separate implementation, in Monte Carlo errors", name), max(abs(mean_gap)), 5)
  ratio <- apply(package, 1L, stats::var) / apply(hand, 1L, stats::var)
  report(sprintf("%s: variance of draws per row, largest departure of the ratio from 1", name), max(abs(ratio - 1)), .15)
  within <- split(seq_along(fill_rows), data$id[fill_rows])
  pairs <- within[lengths(within) >= 2L][seq_len(min(8L, sum(lengths(within) >= 2L)))]
  gaps <- vapply(pairs, function(u) abs(stats::cov(package[u[1L], ], package[u[2L], ]) - stats::cov(hand[u[1L], ], hand[u[2L], ])) / sqrt(stats::var(package[u[1L], ]) * stats::var(package[u[2L], ])), numeric(1))
  report(sprintf("%s: covariance of two rows in one cluster, largest gap as a correlation", name), max(gaps), .08)
}
DRAWS <- list(three_rows_17_of_34 = make_data(rep(3L, 34L), 201L), unequal_low_icc = make_data(rep(c(2L, 3L, 5L), 12L), 202L, tau = .2))
compare_draws("three rows, 17 of 34 clusters audited", DRAWS[[1L]], 1:17, 3L)
compare_draws("unequal clusters, low correlation", DRAWS[[2L]], 1:18, 2L)

# Pooling: catmech_pool_lm against a hand-written Barnard-Rubin calculation and mitools.
set.seed(401L)
ids <- rep(seq_len(30L), each = 3L); x <- stats::rnorm(90L); truth <- x + stats::rnorm(30L)[ids] + stats::rnorm(90L)
recorded <- truth + stats::rnorm(90L, 0, .5)
d <- data.frame(id = ids, x = x, recorded = recorded, gold = ifelse(ids <= 12L, truth, NA_real_))
result <- catmech_reference(d, "recorded", "gold", "x", use_recorded = TRUE, cluster = "id", m = 10L, seed = 7L)
pooled <- catmech_pool_lm(result, recorded ~ x)
pooled <- pooled[pooled$term == "x", ]
fits <- lapply(result$completed, function(z) stats::lm(recorded ~ x, data = z))
q <- vapply(fits, function(f) unname(stats::coef(f)["x"]), numeric(1))
u_mats <- lapply(seq_along(fits), function(j) cr1s(fits[[j]], result$completed[[j]]$id))
u <- vapply(u_mats, function(v) v["x", "x"], numeric(1))
m_imp <- length(q); ubar <- mean(u); b_var <- stats::var(q); total <- ubar + (1 + 1 / m_imp) * b_var
gamma <- (1 + 1 / m_imp) * b_var / total; nu_old <- (m_imp - 1) / gamma^2; nu_com <- 29
nu_obs <- (nu_com + 1) / (nu_com + 3) * nu_com * (1 - gamma); nu <- 1 / (1 / nu_old + 1 / nu_obs)
report("Pooling: estimate against hand-written Rubin's rules", abs(pooled$estimate - mean(q)), 1e-10)
report("Pooling: standard error against hand-written Rubin's rules", abs(pooled$std.error - sqrt(total)), 1e-10)
report("Pooling: degrees of freedom against hand-written Barnard-Rubin", abs(pooled$df - nu), 1e-8)
mit <- mitools::MIcombine(lapply(fits, stats::coef), u_mats, df.complete = nu_com)
report("Pooling: estimate against mitools", abs(pooled$estimate - unname(mit$coefficients["x"])), 1e-10)
report("Pooling: standard error against mitools", abs(pooled$std.error - sqrt(mit$variance["x", "x"])), 1e-10)
mit_df <- unname(mit$df["x"])
lines <- c(lines, sprintf("NOTE Pooling: degrees of freedom %.4f in catmech and %.4f in mitools, which may use a different formula, so this is reported and not scored", pooled$df, mit_df))

header <- c("Independent numerical checks of the random-intercept imputation and cluster-robust pooling.",
            "Tools: lme4, mitools and hand-written matrix code. clubSandwich, sandwich, pan, jomo and brms were not used.",
            sprintf("Failures: %d", failures), "")
writeLines(c(header, lines), output)
cat(c(header, lines), sep = "\n")
if (failures) quit(status = 1L)
