#!/usr/bin/env Rscript
# Scores the complete-data analysis under corrected cluster-robust variances, on the replacement study's cohorts.
# Settings are locked in validation/VARIANCE_CHECK_PLAN.md. No imputation is involved.
# Usage: Rscript tools/oracle-variance-check.R [NEW_DIRECTORY] [REPS]
args <- commandArgs(trailingOnly = TRUE)
file_arg <- grep("^--file=", commandArgs(), value = TRUE)
root <- normalizePath(file.path(dirname(sub("^--file=", "", file_arg[1])), ".."))
outdir <- if (length(args)) args[1] else file.path(root, "validation", "variance-check")
reps <- if (length(args) >= 2L) as.integer(args[2]) else 500L
if (dir.exists(outdir) && length(list.files(outdir, all.files = TRUE, no.. = TRUE)))
  stop("Use a new or empty output directory; existing evidence is never overwritten.")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
study <- file.path(root, "validation", "clustered-fix-study")
config <- jsonlite::read_json(file.path(study, "settings.json"), simplifyVector = TRUE)
expressions <- parse(file.path(root, "tools/clustered-fix-benchmark.R"))
assigned <- vapply(expressions, function(e) if (is.call(e) && identical(e[[1]], as.name("<-")) && is.name(e[[2]])) as.character(e[[2]]) else "", character(1))
first <- which(assigned == "rwanda")[1]; last <- which(assigned == "real_case")[1]
for (e in expressions[first:last]) eval(e, envir = globalenv())
grid <- utils::read.csv(file.path(study, "simulation-grid.csv"), stringsAsFactors = FALSE)
BOOT <- 9999L

sym_power <- function(A, power) {
  e <- eigen((A + t(A)) / 2, symmetric = TRUE)
  keep <- e$values > 1e-10
  e$vectors[, keep, drop = FALSE] %*% (t(e$vectors[, keep, drop = FALSE]) * e$values[keep]^power)
}

# Every variance for coefficient `term`, from the design matrix, response and cluster labels.
score_variances <- function(X, y, ids, term, target, seed, boot = BOOT) {
  n <- nrow(X); p <- ncol(X); k <- match(term, colnames(X))
  units <- split(seq_len(n), ids); G <- length(units)
  XtXi <- chol2inv(chol(crossprod(X)))
  beta <- as.vector(XtXi %*% crossprod(X, y)); e <- as.vector(y - X %*% beta)
  b <- XtXi[, k]
  scores <- vapply(units, function(u) sum(b * crossprod(X[u, , drop = FALSE], e[u])), numeric(1))
  se_cr1s <- sqrt((G / (G - 1)) * ((n - 1) / (n - p)) * sum(scores^2))
  w <- matrix(0, G, n); cr2 <- 0; cv3 <- 0; leverage <- numeric(G)
  for (g in seq_len(G)) {
    u <- units[[g]]
    Xg <- X[u, , drop = FALSE]
    I_H <- diag(length(u)) - Xg %*% XtXi %*% t(Xg)
    A <- sym_power(I_H, -.5)
    AXb <- as.vector(A %*% (Xg %*% b))
    cr2 <- cr2 + sum(AXb * e[u])^2
    wg <- -as.vector(X %*% (XtXi %*% crossprod(Xg, AXb)))
    wg[u] <- wg[u] + AXb
    w[g, ] <- wg
    cv3 <- cv3 + sum(as.vector(Xg %*% b) * as.vector(sym_power(I_H, -1) %*% e[u]))^2
    leverage[g] <- sum(as.vector(Xg %*% b)^2) / sum(b * crossprod(X, X %*% b))
  }
  cr2_df <- { WW <- tcrossprod(w); sum(diag(WW))^2 / sum(WW^2) }
  se_cr2 <- sqrt(cr2); se_cv3 <- sqrt((G - 1) / G * cv3)
  t1 <- stats::qt(.975, G - 1); t2 <- stats::qt(.975, cr2_df)
  covers <- function(se, q) abs(beta[k] - target) <= q * se
  list(estimate = beta[k], se_cr1s = se_cr1s, se_cr2 = se_cr2, se_cv3 = se_cv3, cr2_df = cr2_df, clusters = G,
       leverage = max(leverage), width_cr1s = 2 * t1 * se_cr1s, width_cr2 = 2 * t2 * se_cr2, width_cv3 = 2 * t1 * se_cv3,
       cover_cr1s = covers(se_cr1s, t1), cover_cr2 = covers(se_cr2, t2), cover_cv3 = covers(se_cv3, t1),
       wcr_p = wild_p(X, y, ids, k, target, seed, boot))
}

# Restricted wild cluster bootstrap p-value for H0: coefficient k = target, by cluster-level sums.
wild_p <- function(X, y, ids, k, target, seed, boot = BOOT) {
  Z <- X[, -k, drop = FALSE]; x <- X[, k]
  x_tilde <- stats::lm.fit(Z, x)$residuals
  u_tilde <- stats::lm.fit(Z, y - x * target)$residuals
  groups <- as.integer(factor(ids)); G <- max(groups)
  XtXi <- chol2inv(chol(crossprod(X)))
  s <- as.vector(rowsum(x_tilde * u_tilde, groups, reorder = TRUE))
  L <- rowsum(x_tilde * X, groups, reorder = TRUE)
  R <- rowsum(u_tilde * X, groups, reorder = TRUE)
  A <- L %*% XtXi %*% t(R)
  stat <- function(xi) {
    numerator <- as.vector(xi %*% s)
    cluster <- sweep(xi, 2L, s, `*`) - xi %*% t(A)
    numerator / sqrt(rowSums(cluster^2))
  }
  observed <- stat(matrix(1, 1L, G))
  old <- if (exists(".Random.seed", envir = globalenv())) get(".Random.seed", envir = globalenv()) else NULL
  set.seed(seed)
  xi <- matrix(sample(c(-1, 1), boot * G, replace = TRUE), boot, G)
  if (!is.null(old)) assign(".Random.seed", old, envir = globalenv())
  (1 + sum(abs(stat(xi)) >= abs(observed))) / (boot + 1)
}

# Controls on small datasets, against brute force and explicit matrices.
set.seed(20261102)
control <- function() {
  g <- 12L; id <- rep(seq_len(g), times = rep(c(2L, 3L, 4L), 4L)); n <- length(id)
  X <- cbind(`(Intercept)` = 1, covariate = stats::rexp(g)[id] + stats::rnorm(n, 0, .3), z = stats::rnorm(n))
  y <- as.vector(X %*% c(1, .5, -.2)) + stats::rnorm(g)[id] + stats::rnorm(n)
  got <- score_variances(X, y, id, "covariate", .5, 1L, boot = 19999L)
  full <- stats::lm.fit(X, y)$coefficients[["covariate"]]
  brute <- vapply(unique(id), function(i) stats::lm.fit(X[id != i, ], y[id != i])$coefficients[["covariate"]], numeric(1))
  stopifnot(abs(got$se_cv3 - sqrt((g - 1) / g * sum((brute - full)^2))) < 1e-8)
  # Explicit CR2 and Bell-McCaffrey with full matrices.
  P <- X %*% solve(crossprod(X), t(X)); M <- diag(n) - P; b <- solve(crossprod(X))[, 2L]
  W <- t(vapply(unique(id), function(i) {
    u <- which(id == i); ev <- eigen(diag(length(u)) - P[u, u], symmetric = TRUE)
    A <- ev$vectors %*% diag(ev$values^-.5, length(u)) %*% t(ev$vectors)
    as.vector(M[, u, drop = FALSE] %*% A %*% X[u, , drop = FALSE] %*% b)
  }, numeric(n)))
  stopifnot(abs(got$se_cr2 - sqrt(sum((W %*% y)^2))) < 1e-8)
  stopifnot(abs(got$cr2_df - sum(diag(W %*% t(W)))^2 / sum((W %*% t(W))^2)) < 1e-8)
  # One row per cluster gives HC2.
  one <- seq_len(40L); X1 <- cbind(`(Intercept)` = 1, covariate = stats::rexp(40L)); y1 <- X1[, 2L] + stats::rnorm(40L)
  got1 <- score_variances(X1, y1, one, "covariate", 1, 2L, boot = 99L)
  fit1 <- stats::lm(y1 ~ X1[, 2L]); h <- stats::hatvalues(fit1); xt <- X1[, 2L] - mean(X1[, 2L])
  stopifnot(abs(got1$se_cr2 - sqrt(sum(xt^2 * stats::residuals(fit1)^2 / (1 - h)) / sum(xt^2)^2)) < 1e-8)
  # Wild bootstrap against brute-force refits with the same weights.
  Z <- X[, -2L]; restricted <- stats::lm.fit(Z, y - X[, 2L] * .5); u <- restricted$residuals; fitted0 <- y - u
  t_of <- function(yy) {
    f <- stats::lm.fit(X, yy); XtXi <- solve(crossprod(X)); sc <- rowsum(X * f$residuals, id, reorder = TRUE)
    (f$coefficients[[2L]] - .5) / sqrt((XtXi %*% crossprod(sc) %*% XtXi)[2L, 2L])
  }
  set.seed(5L); xi <- matrix(sample(c(-1, 1), 400L * g, replace = TRUE), 400L, g)
  brute_t <- apply(xi, 1L, function(w) t_of(fitted0 + w[as.integer(factor(id))] * u))
  fast <- wild_p(X, y, id, 2L, .5, 5L, boot = 400L)
  stopifnot(abs(fast - (1 + sum(abs(brute_t) >= abs(t_of(y)))) / 401) < 1e-12)
  "cv3 against refits, cr2 and its degrees of freedom against explicit matrices, cr2 against HC2 with one row per cluster, and the wild bootstrap against brute-force refits: PASS"
}
control_line <- control(); cat(control_line, "\n")

cells <- list()
for (z in real) cells[[length(cells) + 1L]] <- list(kind = "real", z = z, analysis = z$id, level = z$level, fraction = .25, cell = 100L + length(cells) + 1L)
for (i in seq_len(nrow(grid))) cells[[length(cells) + 1L]] <- list(kind = "simulation", spec = as.list(grid[i, ]), index = i, analysis = grid$label[i], level = "simulation", fraction = grid$fraction[i], cell = i)
one_request <- function(s, r) {
  p <- if (s$kind == "real") real_case(s$z, 1L, r) else simulate_cell(s$spec, as.integer(config$simulation$cohort + s$index * 1e5 + r), s$fraction)
  d <- p$input; d$lab_value <- p$truth
  X <- stats::model.matrix(p$formula, data = d)
  got <- score_variances(X, d$lab_value, d$cluster_id, "covariate", p$target, as.integer(1600000000 + s$cell * 1e5 + r))
  data.frame(level = s$level, analysis = s$analysis, fraction = s$fraction, replicate = r, target = p$target, as.data.frame(got))
}
started <- Sys.time()
records <- list()
for (s in cells) {
  batch <- parallel::mclapply(seq_len(reps), function(r) one_request(s, r), mc.cores = 14L, mc.preschedule = TRUE, mc.set.seed = FALSE)
  if (any(vapply(batch, inherits, logical(1), what = "try-error"))) stop("A worker failed on ", s$analysis, ": ", conditionMessage(attr(batch[[which(vapply(batch, inherits, logical(1), what = "try-error"))[1]]], "condition")))
  records <- c(records, batch)
}
results <- do.call(rbind, records)
out <- gzfile(file.path(outdir, "replicates.csv.gz"), "wt")
utils::write.csv(results, out, row.names = FALSE)
close(out)
writeLines(c(control_line, paste("Requests:", nrow(results)), paste("Elapsed seconds:", round(as.numeric(difftime(Sys.time(), started, units = "secs")), 1)), capture.output(sessionInfo())), file.path(outdir, "session-info.txt"))
cat("Completed", nrow(results), "requests in", round(as.numeric(difftime(Sys.time(), started, units = "secs"))), "seconds\n")
