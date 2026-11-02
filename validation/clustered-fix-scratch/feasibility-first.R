#!/usr/bin/env Rscript
# Scratch feasibility check, simulated data only, separate seeds. Not evidence and no real item is touched.
root <- "/Users/patrickiben/Documents/catmech"
out <- commandArgs(trailingOnly = TRUE)[1]
pkgload::load_all(root, quiet = TRUE)
workers <- 14L; reps <- 500L; m <- 20L
namespace <- asNamespace("catmech")
package_fill <- get("catmech_reference_cluster_fill", envir = namespace)

rescaled_fill <- function(x, y, fit_rows, fill_rows, ids, strata, m) {
  by_cluster <- split(fit_rows, as.character(ids[fit_rows]))
  groups <- split(names(by_cluster), as.character(strata[vapply(by_cluster, `[`, integer(1), 1L)]))
  out <- matrix(NA_real_, length(fill_rows), m)
  for (j in seq_len(m)) {
    draw <- NULL
    for (attempt in seq_len(200L)) {
      rows <- integer(); w <- numeric()
      for (u in groups) {
        nh <- length(u)
        pick <- if (nh == 1L) u else u[sample.int(nh, nh - 1L, replace = TRUE)]
        got <- by_cluster[pick]
        rows <- c(rows, unlist(got, use.names = FALSE)); w <- c(w, rep(if (nh == 1L) 1 else nh / (nh - 1), sum(lengths(got))))
      }
      fit <- stats::lm.wfit(x[rows, , drop = FALSE], y[rows], w)
      if (fit$rank != ncol(x) || sum(w) <= ncol(x)) next
      variance <- sum(w * fit$residuals^2) / (sum(w) - ncol(x))
      if (!is.finite(variance) || variance <= 0) next
      draw <- as.vector(x[fill_rows, , drop = FALSE] %*% fit$coefficients) + stats::rnorm(length(fill_rows), 0, sqrt(variance))
      break
    }
    if (is.null(draw)) stop("lost design support")
    out[, j] <- draw
  }
  out
}

# Normal draw around the full audit's fit, with the stratified cluster-robust covariance.
robust_fill <- function(x, y, fit_rows, fill_rows, ids, strata, m) {
  X <- x[fit_rows, , drop = FALSE]; n <- nrow(X); p <- ncol(X)
  fit <- stats::lm.fit(X, y[fit_rows])
  if (fit$rank != p || n <= p) stop("no full rank")
  bread <- chol2inv(qr.R(fit$qr))
  clusters <- as.character(ids[fit_rows])
  scores <- rowsum(X * fit$residuals, clusters, reorder = FALSE)
  stratum <- if (is.null(strata)) rep("all", nrow(scores)) else as.character(strata[fit_rows])[match(rownames(scores), clusters)]
  meat <- matrix(0, p, p)
  for (h in unique(stratum)) {
    z <- scores[stratum == h, , drop = FALSE]; nh <- nrow(z)
    meat <- meat + if (nh == 1L) crossprod(z) else nh / (nh - 1) * crossprod(sweep(z, 2L, colMeans(z)))
  }
  V <- (n - 1) / (n - p) * bread %*% meat %*% bread
  e <- eigen((V + t(V)) / 2, symmetric = TRUE)
  root_v <- e$vectors %*% (t(e$vectors) * sqrt(pmax(e$values, 0)))
  s2 <- sum(fit$residuals^2) / (n - p)
  out <- matrix(NA_real_, length(fill_rows), m)
  for (j in seq_len(m)) {
    beta <- fit$coefficients + as.vector(root_v %*% stats::rnorm(p))
    sigma2 <- s2 * (n - p) / stats::rchisq(1L, n - p)
    out[, j] <- as.vector(x[fill_rows, , drop = FALSE] %*% beta) + stats::rnorm(length(fill_rows), 0, sqrt(sigma2))
  }
  out
}

set_fill <- function(f) {
  unlockBinding("catmech_reference_cluster_fill", namespace)
  assign("catmech_reference_cluster_fill", f, envir = namespace)
  lockBinding("catmech_reference_cluster_fill", namespace)
}

SIZES <- c(14L, 9L, 8L, 13L, 5L, 9L, 11L)
simulate <- function(dgp, seed, fraction) {
  set.seed(seed)
  rows <- 3L; g <- sum(SIZES)
  id <- rep(seq_len(g), each = rows)
  stratum_of <- rep(seq_along(SIZES), SIZES)
  stratum <- factor(paste0("s", stratum_of[id]))
  size <- rep(stats::rlnorm(g, 3, .5), each = rows)
  x <- size + stats::rnorm(g * rows, 0, 3)
  if (dgp %in% c("shift", "under")) {
    effect <- stats::rnorm(length(SIZES), 0, 2)[stratum_of[id]]
    u <- rep(stats::rnorm(g, 0, sqrt(.2 * 100)), each = rows)
    truth <- 2 + .6 * x + effect + u + stats::rnorm(g * rows, 0, sqrt(.8 * 100))
    target <- .6
    recorded <- if (dgp == "shift") truth + rep(stats::rnorm(g, 0, 5), each = rows) + stats::rnorm(g * rows, 0, 3) else {
      d <- rep(stats::rbeta(g, 2, 6), each = rows)
      (1 - d) * truth + stats::rnorm(g * rows, 0, 2)
    }
  } else {
    frailty <- rep(stats::rgamma(g, .5, .5), each = rows)
    effect <- stats::runif(length(SIZES), 0, .2)[stratum_of[id]]
    truth <- stats::rpois(g * rows, pmax(.01, (.1 + .02 * x + effect)) * frailty)
    target <- .02
    slip <- stats::runif(g * rows) < .2
    recorded <- pmax(0, truth + slip * sample(c(-1, 1), g * rows, replace = TRUE))
  }
  chosen <- logical(g); start <- 0L
  for (nh in SIZES) {
    chosen[start + sample.int(nh, max(1L, as.integer(round(fraction * nh))))] <- TRUE
    start <- start + nh
  }
  audited <- chosen[id]
  gold <- truth; gold[!audited] <- NA_real_
  list(input = data.frame(covariate = x, stratum = stratum, cluster_id = id, lab_value = recorded, gold_value = gold), truth = truth, target = target,
       formula = lab_value ~ covariate + stratum, predictors = c("covariate", "stratum"))
}

VARIANTS <- list(oracle = list(), existing = list(),
                 within_strata = list(cluster = "cluster_id", strata = "stratum", fill = package_fill),
                 no_strata = list(cluster = "cluster_id", strata = NULL, fill = package_fill),
                 rao_wu = list(cluster = "cluster_id", strata = "stratum", fill = rescaled_fill),
                 robust_draw = list(cluster = "cluster_id", strata = "stratum", fill = robust_fill))
run_variant <- function(name, p, seed) {
  v <- VARIANTS[[name]]; input <- p$input
  if (name == "oracle") input$gold_value <- p$truth
  fit <- if (is.null(v$cluster)) {
    catmech_reference(input, "lab_value", "gold_value", p$predictors, method = "normal", m = m, seed = seed, use_recorded = TRUE)
  } else {
    set_fill(v$fill); on.exit(set_fill(package_fill), add = TRUE)
    catmech_reference(input, "lab_value", "gold_value", p$predictors, method = "normal", m = m, seed = seed, use_recorded = TRUE, cluster = v$cluster, cluster_strata = v$strata)
  }
  pooled <- catmech_pool_lm(fit, p$formula, cluster = "cluster_id")
  pooled[pooled$term == "covariate", ]
}
icc_anova <- function(e, id) {
  id <- factor(id); k <- tapply(e, id, length); g <- nlevels(id); n <- length(e)
  means <- tapply(e, id, mean)
  between <- sum(k * (means - mean(e))^2) / (g - 1); within <- sum((e - means[id])^2) / (n - g)
  k0 <- (n - sum(k^2) / n) / (g - 1)
  (between - within) / (between + (k0 - 1) * within)
}
one <- function(dgp, index, fraction, r) {
  p <- simulate(dgp, as.integer(1100000000 + index * 1e5 + r), fraction)
  seed <- as.integer(1200000000 + index * 1e5 + r)
  full <- stats::residuals(stats::lm(p$truth ~ p$input$lab_value + p$input$covariate + p$input$stratum))
  rows <- lapply(names(VARIANTS), function(name) {
    z <- tryCatch(suppressWarnings(run_variant(name, p, seed)), error = identity)
    bad <- inherits(z, "error")
    data.frame(dgp = dgp, fraction = fraction, r = r, variant = name, ok = !bad,
               estimate = if (bad) NA else z$estimate, lower = if (bad) NA else z$lower, upper = if (bad) NA else z$upper,
               df = if (bad) NA else z$df, b = if (bad) NA else z$between_variance, target = p$target, icc = icc_anova(full, p$input$cluster_id))
  })
  do.call(rbind, rows)
}
started <- Sys.time()
cells <- expand.grid(dgp = c("shift", "under", "rare"), fraction = c(.25, .5), stringsAsFactors = FALSE)
got <- list()
for (i in seq_len(nrow(cells))) got[[i]] <- do.call(rbind, parallel::mclapply(seq_len(reps), function(r) one(cells$dgp[i], i, cells$fraction[i], r), mc.cores = workers, mc.preschedule = TRUE, mc.set.seed = FALSE))
d <- do.call(rbind, got)
utils::write.csv(d, out, row.names = FALSE)
oracle <- d[d$variant == "oracle", c("dgp", "fraction", "r", "estimate")]; names(oracle)[4] <- "oracle_estimate"
d <- merge(d, oracle)
s <- do.call(rbind, lapply(split(d, interaction(d$dgp, d$fraction, d$variant, drop = TRUE)), function(g) {
  k <- g$ok
  data.frame(dgp = g$dgp[1], fraction = g$fraction[1], variant = g$variant[1], returned = sum(k), coverage = round(100 * mean(g$lower[k] <= g$target[k] & g$upper[k] >= g$target[k]), 1),
             median_width = signif(stats::median(g$upper[k] - g$lower[k]), 3), wide = sum((g$upper[k] - g$lower[k]) > 5 * stats::median(g$upper[k] - g$lower[k])),
             rmse = signif(sqrt(mean((g$estimate[k] - g$target[k])^2)), 3), recovery = signif(sqrt(mean((g$estimate[k] - g$oracle_estimate[k])^2)), 3),
             median_b = signif(stats::median(g$b[k]), 3), icc = round(stats::median(g$icc), 2))
}))
order_v <- match(s$variant, names(VARIANTS))
s <- s[order(s$dgp, s$fraction, order_v), ]; rownames(s) <- NULL
print(s, row.names = FALSE)
cat("Seconds:", round(as.numeric(difftime(Sys.time(), started, units = "secs"))), "\n")
