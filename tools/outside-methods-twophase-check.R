#!/usr/bin/env Rscript
# Why raking's first-run standard errors were too small: what survey::twophase does with phase-one clusters and no phase-one
# probabilities. Simulated cohorts with seeds used nowhere else; no study result is read. See validation/OUTSIDE_METHODS_SYNTHESIS.md.
# Usage: Rscript tools/outside-methods-twophase-check.R [OUTPUT_FILE]
args <- commandArgs(trailingOnly = TRUE)
file_arg <- grep("^--file=", commandArgs(), value = TRUE)
root <- normalizePath(file.path(dirname(sub("^--file=", "", file_arg[1])), ".."))
output <- if (length(args)) args[1] else file.path(root, "validation", "outside-methods-twophase-check.txt")
suppressMessages({pkgload::load_all(root, quiet = TRUE); library(survey)})
options(survey.lonely.psu = "adjust")
config <- jsonlite::read_json(file.path(root, "validation/clustered-fix-study/settings.json"), simplifyVector = TRUE)
expressions <- parse(file.path(root, "tools/clustered-fix-benchmark.R"))
assigned <- vapply(expressions, function(e) if (is.call(e) && identical(e[[1]], as.name("<-")) && is.name(e[[2]])) as.character(e[[2]]) else "", character(1))
for (e in expressions[which(assigned == "rwanda")[1]:which(assigned == "real_case")[1]]) eval(e, envir = globalenv())
prep <- function(cfg, r) {
  p <- simulate_cell(list(design = cfg$design, errors = cfg$errors, rows = cfg$rows), as.integer(2000000 + r + cfg$off), cfg$fraction)
  d <- p$input; d$audited <- p$audited; d$value_r <- d$gold_value; d$bigN <- 1e6
  if (!"stratum" %in% names(d)) d$stratum <- factor("all")
  fit <- stats::lm(p$formula, data = p$input); x <- stats::model.matrix(fit)
  d$if_k <- as.vector(x %*% chol2inv(qr.R(fit$qr))[, 2]) * stats::resid(fit)
  list(p = p, d = d)
}
lines <- c("Simulated cohorts, seeds used nowhere else.", "")
z <- prep(list(design = "none", rows = 3L, errors = "shift", fraction = .5, off = 1), 1); d <- z$d; d$y <- d$value_r
a <- stats::aggregate(y ~ cluster_id, d[d$audited, ], sum)
g <- d[!duplicated(d$cluster_id), c("cluster_id", "audited", "stratum", "bigN")]; g$y <- a$y[match(g$cluster_id, a$cluster_id)]
lines <- c(lines, "ONE TOTAL, ONE DATASET, three rows per cluster: standard error, then its phase one and phase two variance")
parts <- function(design) {
  v <- stats::vcov(survey::svytotal(~y, design)); ph <- attr(v, "phases")
  sprintf("SE %.0f, phase one %s, phase two %s", sqrt(as.numeric(v)), if (is.null(ph)) "not reported" else format(as.numeric(ph$phase1), digits = 4), if (is.null(ph)) "not reported" else format(as.numeric(ph$phase2), digits = 4))
}
for (m in c("full", "approx", "simple")) {
  row <- survey::twophase(id = list(~cluster_id, ~cluster_id), subset = ~audited, data = d, method = m)
  cl <- survey::twophase(id = list(~cluster_id, ~cluster_id), subset = ~audited, data = g, method = m)
  lines <- c(lines, sprintf("  %s, row level: %s", m, parts(row)), sprintf("  %s, one row per cluster: %s", m, parts(cl)))
}
row <- survey::twophase(id = list(~cluster_id, ~cluster_id), fpc = list(~bigN, NULL), subset = ~audited, data = d, method = "full")
lines <- c(lines, sprintf("  full with a phase-one population of one million clusters, row level: %s", parts(row)), "")
one <- function(cfg, r, how) {
  z <- prep(cfg, r); d <- z$d; st <- if (nlevels(d$stratum) > 1L) ~stratum else NULL
  des <- if (how == "full_fpc") survey::twophase(id = list(~cluster_id, ~cluster_id), strata = list(NULL, st), fpc = list(~bigN, NULL), subset = ~audited, data = d) else
    survey::twophase(id = list(~cluster_id, ~cluster_id), strata = list(NULL, st), subset = ~audited, data = d, method = how)
  f <- stats::update(z$p$formula, value_r ~ .)
  g0 <- survey::svyglm(f, design = des); cal <- survey::calibrate(des, phase = 2, formula = ~if_k, calfun = "raking"); g1 <- survey::svyglm(f, design = cal)
  c(unname(stats::coef(g0)[2]), unname(survey::SE(g0)[2]), unname(stats::coef(g1)[2]), unname(survey::SE(g1)[2]))
}
cfgs <- list(list(design = "none", rows = 3L, errors = "shift", fraction = .5, off = 0), list(design = "none", rows = 3L, errors = "shift", fraction = .25, off = 1000),
             list(design = "rwanda", rows = 3L, errors = "shift", fraction = .5, off = 2000), list(design = "rwanda", rows = 3L, errors = "shift", fraction = .25, off = 4000),
             list(design = "none_210", rows = 1L, errors = "row_only", fraction = .5, off = 3000))
lines <- c(lines, "STANDARD ERROR OVER THE SPREAD OF THE ESTIMATE, 300 cohorts: without calibration, then raked")
for (how in c("full", "approx", "full_fpc")) for (cfg in cfgs) {
  res <- do.call(rbind, parallel::mclapply(1:300, function(r) tryCatch(one(cfg, r, how), error = function(e) rep(NA_real_, 4)), mc.cores = 14L, mc.set.seed = FALSE))
  lines <- c(lines, sprintf("  %s, %s, %d rows, %d%%: %.2f, %.2f", how, cfg$design, cfg$rows, 100 * cfg$fraction,
    mean(res[, 2], na.rm = TRUE) / stats::sd(res[, 1], na.rm = TRUE), mean(res[, 4], na.rm = TRUE) / stats::sd(res[, 3], na.rm = TRUE)))
}
writeLines(lines, output)
cat(lines, sep = "\n"); cat("\n")
