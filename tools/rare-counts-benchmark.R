#!/usr/bin/env Rscript
# Prospective test of what goes wrong with rare counts; settings locked in validation/RARE_COUNTS_PLAN.md.
# Usage: Rscript tools/rare-counts-benchmark.R pilot|study [NEW_DIRECTORY] [PILOT_REPS]
args <- commandArgs(trailingOnly = TRUE)
mode <- if (length(args)) match.arg(args[1], c("pilot", "study")) else "pilot"
file_arg <- grep("^--file=", commandArgs(), value = TRUE)
root <- normalizePath(file.path(dirname(sub("^--file=", "", file_arg[1])), ".."))
outdir <- if (length(args) >= 2L) args[2] else file.path(root, "validation", paste0("rare-counts-", mode))
if (dir.exists(outdir) && length(list.files(outdir, all.files = TRUE, no.. = TRUE)))
  stop("Use a new or empty output directory; existing evidence is never overwritten.")
hash_sources <- function() {
  paths <- c(sort(list.files(file.path(root, "R"), "\\.R$", full.names = TRUE)),
    file.path(root, c("DESCRIPTION", "NAMESPACE", "tools/rare-counts-benchmark.R", "tools/rare-counts-score.py", "tools/variance-mi-benchmark.R", "tools/clustered-fix-benchmark.R",
                      "validation/RARE_COUNTS_PLAN.md", "validation/real-audit/rwanda-pool.csv", "validation/variance-mi-study/replicates.csv.gz")))
  data.frame(file = substring(paths, nchar(root) + 2L),
    sha256 = vapply(paths, function(p) digest::digest(file = p, algo = "sha256"), character(1)), row.names = NULL)
}
hashes_before <- hash_sources()
pkgload::load_all(root, quiet = TRUE)
suppressPackageStartupMessages({library(mice); library(MASS)})
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
pilot <- mode == "pilot"
pilot_reps <- if (length(args) >= 3L) as.integer(args[3]) else 3L
config <- jsonlite::read_json(file.path(root, "validation/clustered-fix-study/settings.json"), simplifyVector = TRUE)
if (pilot) {
  config$held_out <- list(cohort = 1060100000, selection = 1070100000, imputation = 1080000000)
  config$development <- list(cohort = 1090100000, selection = 1100100000, imputation = 1105000000)
  config$simulation <- list(cohort = 1106100000, imputation = 1107000000)
}
SEEDS <- if (pilot) list(cohort = 1010000000, boot = 1020000000, wild = 1030000000, target = 1035000000, nb = 1040000000, pmm = 1050000000) else
  list(cohort = 1210000000, boot = 1220000000, wild = 1230000000, target = 1235000000, nb = 1240000000, pmm = 1250000000)
settings <- list(mode = mode, m = 20L, part1_reps = if (pilot) pilot_reps else 1000L, part2_real_25 = if (pilot) pilot_reps else 2000L, part2_other = if (pilot) pilot_reps else 500L,
                 part3_reps = if (pilot) pilot_reps else 500L, boot = 399L, wild = 999L, seeds = SEEDS)
jsonlite::write_json(settings, file.path(outdir, "settings.json"), pretty = TRUE, auto_unbox = TRUE)
fix <- parse(file.path(root, "tools/clustered-fix-benchmark.R"))
fix_names <- vapply(fix, function(e) if (is.call(e) && identical(e[[1]], as.name("<-")) && is.name(e[[2]])) as.character(e[[2]]) else "", character(1))
for (e in fix[which(fix_names == "rwanda")[1]:which(fix_names == "real_case")[1]]) eval(e, envir = globalenv())
variance <- parse(file.path(root, "tools/variance-mi-benchmark.R"))
variance_names <- vapply(variance, function(e) if (is.call(e) && identical(e[[1]], as.name("<-")) && is.name(e[[2]])) as.character(e[[2]]) else "", character(1))
for (n in c("sym_power", "analysis_parts", "one_analysis", "wild_pooled_p")) eval(variance[[which(variance_names == n)[1]]], envir = globalenv())
BOOT <- settings$wild
grid <- utils::read.csv(file.path(root, "validation/clustered-fix-study/simulation-grid.csv"), stringsAsFactors = FALSE)
FRAILTY <- c(none = Inf, shape_2 = 2, shape_05 = .5)
part1 <- expand.grid(design = c("rwanda", "pairs", "none"), frailty = names(FRAILTY), scale = c(1, 10), stringsAsFactors = FALSE)
part1$cell <- seq_len(nrow(part1))
part1$label <- paste(part1$design, part1$frailty, ifelse(part1$scale == 1, "rare", "common"), sep = "_")

generate <- function(cell, seed, clusters = NULL) {
  set.seed(seed)
  a <- config$simulation_parameters
  sizes <- DESIGNS[[cell$design]]; if (!is.null(clusters)) sizes <- rep(sizes, each = 1L) * clusters
  rows <- 3L; g <- sum(sizes); n <- g * rows
  stratified <- length(sizes) > 1L
  id <- rep(seq_len(g), each = rows)
  row_stratum <- rep(seq_along(sizes), sizes)[id]
  x <- rep(stats::rlnorm(g, a$size_log_mean, a$size_log_sd), each = rows) + stats::rnorm(n, 0, a$x_noise_sd)
  level <- if (stratified) stats::runif(length(sizes), 0, .2)[row_stratum] else 0
  shape <- FRAILTY[[cell$frailty]]
  frailty <- if (is.finite(shape)) rep(stats::rgamma(g, shape, shape), each = rows) else 1
  y <- as.numeric(stats::rpois(n, pmax(.01, .1 + .02 * x + level) * frailty * cell$scale))
  d <- data.frame(y = y, x = x, id = id)
  if (stratified) d$stratum <- factor(paste0("s", row_stratum))
  d
}
design_of <- function(d) if ("stratum" %in% names(d)) stats::model.matrix(~ x + stratum, data = d) else stats::model.matrix(~ x, data = d)
log_targets <- vapply(seq_len(nrow(part1)), function(i) {
  big <- generate(part1[i, ], as.integer(SEEDS$target + i), clusters = if (pilot) 200L else ceiling(200000 / sum(DESIGNS[[part1$design[i]]])))
  unname(stats::coef(stats::glm(y ~ . - id, family = stats::quasipoisson(), data = big))["x"])
}, numeric(1))
part1$log_target <- log_targets
utils::write.csv(part1, file.path(outdir, "part1-cells.csv"), row.names = FALSE)

boot_t_cover <- function(X, y, ids, strata, k, estimate, se, target, seed) {
  set.seed(seed)
  units <- split(seq_len(nrow(X)), ids)
  by_stratum <- split(seq_along(units), strata[!duplicated(ids)][order(unique(ids))])
  stats <- rep(NA_real_, settings$boot)
  for (b in seq_len(settings$boot)) {
    chosen <- unlist(lapply(by_stratum, function(u) u[sample.int(length(u), length(u), replace = TRUE)]), use.names = FALSE)
    rows <- unlist(units[chosen], use.names = FALSE)
    Xb <- X[rows, , drop = FALSE]; yb <- y[rows]; idb <- rep(seq_along(chosen), lengths(units[chosen]))
    fit <- stats::lm.fit(Xb, yb)
    if (fit$rank < ncol(X)) next
    bread <- chol2inv(qr.R(fit$qr))
    scores <- rowsum(Xb * fit$residuals, idb, reorder = FALSE)
    G <- length(chosen); n <- nrow(Xb)
    v <- (bread %*% crossprod(scores) %*% bread)[k, k] * (G / (G - 1)) * ((n - 1) / (n - ncol(X)))
    if (v > 0) stats[b] <- (fit$coefficients[k] - estimate) / sqrt(v)
  }
  q <- stats::quantile(stats, c(.025, .975), na.rm = TRUE, names = FALSE)
  list(cover = as.numeric(estimate - q[2] * se <= target && target <= estimate - q[1] * se), width = (q[2] - q[1]) * se)
}
diagnostic <- function(parts, y) {
  delta <- as.vector(parts$C %*% y)
  max(delta^2) / sum(delta^2)
}
poisson_interval <- function(d, target) {
  fit <- suppressWarnings(stats::glm(y ~ . - id, family = stats::poisson(), data = d))
  if (!fit$converged) stop("Poisson fit did not converge.")
  X <- stats::model.matrix(fit); mu <- stats::fitted(fit)
  scores <- rowsum(X * (d$y - mu), d$id, reorder = FALSE)
  bread <- chol2inv(chol(crossprod(X * sqrt(mu))))
  G <- nrow(scores)
  k <- match("x", colnames(X))
  v <- (bread %*% crossprod(scores) %*% bread)[k, k] * G / (G - 1)
  half <- stats::qt(.975, G - 1) * sqrt(v)
  c(cover = as.numeric(abs(unname(stats::coef(fit)["x"]) - target) <= half), width = 2 * half)
}
part1_replicate <- function(i, r) {
  cell <- part1[i, ]
  d <- generate(cell, as.integer(SEEDS$cohort + i * 1e5 + r))
  X <- design_of(d); k <- match("x", colnames(X)); target <- .02 * cell$scale
  row <- list(label = cell$label, replicate = r, status = "ok", estimate = NA_real_, cover_cv3 = NA_real_, cover_cr2 = NA_real_, cover_wcr = NA_real_, cover_boot = NA_real_, cover_poisson = NA_real_,
              width_cv3 = NA_real_, width_cr2 = NA_real_, width_boot = NA_real_, width_poisson = NA_real_, D = NA_real_)
  if (qr(X)$rank < ncol(X)) { row$status <- "failed"; return(as.data.frame(row)) }
  parts <- analysis_parts(X, d$id, k)
  a <- one_analysis(parts, matrix(d$y, ncol = 1L))
  G <- parts$G
  half_cv3 <- stats::qt(.975, G - 1) * sqrt(a$cv3); half_cr2 <- stats::qt(.975, parts$df) * sqrt(a$cr2)
  row$estimate <- a$q
  row$cover_cv3 <- as.numeric(abs(a$q - target) <= half_cv3); row$width_cv3 <- 2 * half_cv3
  row$cover_cr2 <- as.numeric(abs(a$q - target) <= half_cr2); row$width_cr2 <- 2 * half_cr2
  row$D <- diagnostic(parts, d$y)
  se_cr1 <- sqrt(a$cr1s)
  p <- tryCatch(wild_pooled_p(X, matrix(d$y, ncol = 1L), d$id, k, target, as.integer(SEEDS$wild + i * 1e5 + r), (a$q - target) / se_cr1, parts$f), error = function(e) NA_real_)
  row$cover_wcr <- as.numeric(p > .05)
  strata <- if ("stratum" %in% names(d)) d$stratum else factor(rep("all", nrow(d)))
  b <- tryCatch(boot_t_cover(X, d$y, d$id, strata, k, a$q, se_cr1, target, as.integer(SEEDS$boot + i * 1e5 + r)), error = function(e) list(cover = NA_real_, width = NA_real_))
  row$cover_boot <- b$cover; row$width_boot <- b$width
  pois <- tryCatch(poisson_interval(d, cell$log_target), error = function(e) c(cover = NA_real_, width = NA_real_))
  row$cover_poisson <- unname(pois["cover"]); row$width_poisson <- unname(pois["width"])
  as.data.frame(row, stringsAsFactors = FALSE)
}

# Part 2: count-aware imputation on the cohorts of the corrected-variance study.
nb_draws <- function(input, predictors, m) {
  audited <- !is.na(input$gold_value)
  frame <- input[, c(predictors), drop = FALSE]; frame$log_recorded <- log1p(pmax(0, input$lab_value))
  X <- stats::model.matrix(stats::reformulate(names(frame)), data = frame)
  fit_data <- data.frame(y = input$gold_value[audited], X[audited, -1L, drop = FALSE], check.names = TRUE)
  if (any(fit_data$y < 0) || any(fit_data$y != round(fit_data$y))) stop("The audited values are not counts.")
  fit <- tryCatch(suppressWarnings(MASS::glm.nb(y ~ ., data = fit_data)), error = function(e) NULL)
  use_nb <- !is.null(fit) && isTRUE(fit$converged) && is.finite(fit$theta) && fit$theta < 1e4 && all(is.finite(stats::coef(fit)))
  if (!use_nb) fit <- suppressWarnings(stats::glm(y ~ ., family = stats::poisson(), data = fit_data))
  beta <- stats::coef(fit); V <- stats::vcov(fit)
  if (anyNA(beta) || any(!is.finite(V))) stop("The audited rows do not have full column rank.")
  chol_V <- chol(V + diag(1e-10, ncol(V)))
  fill <- which(!audited)
  out <- matrix(rep(input$gold_value, m), nrow(input), m)
  for (j in seq_len(m)) {
    b <- beta + as.vector(crossprod(chol_V, stats::rnorm(length(beta))))
    mu <- exp(pmin(30, as.vector(X[fill, , drop = FALSE] %*% b)))
    out[fill, j] <- if (use_nb) {
      theta <- exp(stats::rnorm(1L, log(fit$theta), fit$SE.theta / fit$theta))
      stats::rnbinom(length(fill), size = theta, mu = mu)
    } else stats::rpois(length(fill), mu)
  }
  list(draws = out, fallback = !use_nb)
}
pmm_draws <- function(input, predictors, m, seed) {
  use <- unique(c(predictors, "lab_value", "gold_value"))
  d <- input[, use, drop = FALSE]
  matrix_p <- matrix(0L, length(use), length(use), dimnames = list(use, use)); matrix_p["gold_value", setdiff(use, "gold_value")] <- 1L
  chosen <- stats::setNames(rep("", length(use)), use); chosen["gold_value"] <- "pmm"
  fit <- suppressWarnings(mice::mice(d, m = m, maxit = 1, method = chosen, predictorMatrix = matrix_p, printFlag = FALSE, seed = seed, eps = 0))
  list(draws = vapply(seq_len(m), function(j) mice::complete(fit, j)$gold_value, numeric(nrow(d))), fallback = FALSE)
}
pool_draws <- function(p, draws, ids) {
  fits <- lapply(seq_len(ncol(draws)), function(j) {
    d <- p$input; d$lab_value <- draws[, j]
    stats::lm(p$formula, data = d, na.action = stats::na.fail)
  })
  q <- vapply(fits, function(f) unname(stats::coef(f)["covariate"]), numeric(1))
  u <- vapply(fits, function(f) catmech_cluster_jackknife(f, ids)["covariate", "covariate"], numeric(1))
  mice::pool.scalar(q, u, n = length(unique(ids)), k = 1)
}
limited <- function(expr) {
  job <- parallel::mcparallel(tryCatch(expr, error = identity), silent = TRUE)
  got <- parallel::mccollect(job, wait = FALSE, timeout = 60)
  if (is.null(got)) {
    tools::pskill(job$pid, tools::SIGKILL); parallel::mccollect(job, wait = TRUE)
    return(simpleError("exceeded the 60 second limit"))
  }
  parallel::mccollect(job, wait = TRUE)
  out <- got[[1L]]
  if (is.null(out)) simpleError("the child returned nothing") else out
}
part2_cells <- list()
for (i in which(grid$errors == "rare")) part2_cells[[length(part2_cells) + 1L]] <- list(kind = "simulation", spec = as.list(grid[i, ]), index = i, analysis = grid$label[i], fraction = grid$fraction[i], reps = settings$part2_other, number = i)
for (z in real) if (z$id %in% c("lbw", "resus", "stillbirths")) for (j in seq_along(config$fractions))
  part2_cells[[length(part2_cells) + 1L]] <- list(kind = "real", z = z, fraction_index = j, analysis = z$id, fraction = config$fractions[j], reps = if (j == 1L) settings$part2_real_25 else settings$part2_other, number = 100L + z$index * 2L + j + if (z$level == "development") 20L else 0L)
part2_cohort <- function(s, r) {
  if (s$kind == "real") real_case(s$z, s$fraction_index, r) else simulate_cell(s$spec, as.integer(config$simulation$cohort + s$index * 1e5 + r), s$fraction)
}
part2_replicate <- function(s, r) {
  p <- part2_cohort(s, r)
  ids <- p$input$cluster_id
  do.call(rbind, lapply(c("nb_mi", "pmm_recorded"), function(method) {
    seed <- as.integer(SEEDS[[if (method == "nb_mi") "nb" else "pmm"]] + s$number * 1e5 + r)
    started <- Sys.time()
    got <- limited({
      set.seed(seed)
      drawn <- if (method == "nb_mi") nb_draws(p$input, p$predictors, settings$m) else pmm_draws(p$input, p$predictors, settings$m, seed)
      c(pool_draws(p, drawn$draws, ids), list(fallback = drawn$fallback))
    })
    failed <- inherits(got, "error")
    half <- if (failed) NA_real_ else stats::qt(.975, got$df) * sqrt(got$t)
    data.frame(analysis = s$analysis, fraction = s$fraction, replicate = r, method = method, status = if (failed) "failed" else "ok", message = if (failed) substr(conditionMessage(got), 1L, 90L) else "",
               estimate = if (failed) NA_real_ else got$qbar, target = p$target, width = 2 * half, cover = if (failed) NA_real_ else as.numeric(abs(got$qbar - p$target) <= half), fallback = if (failed) NA else got$fallback,
               seconds = round(as.numeric(difftime(Sys.time(), started, units = "secs")), 4), stringsAsFactors = FALSE)
  }))
}
# Part 3: the diagnostic on the recount, for cells that are not rare and for real items at 25%.
part3_replicate <- function(s, r) {
  p <- part2_cohort(s, r)
  d <- p$input; d$lab_value <- p$truth
  X <- stats::model.matrix(p$formula, data = d)
  k <- match("covariate", colnames(X))
  parts <- analysis_parts(X, d$cluster_id, k)
  a <- one_analysis(parts, matrix(d$lab_value, ncol = 1L))
  half <- stats::qt(.975, parts$G - 1) * sqrt(a$cv3)
  data.frame(analysis = s$analysis, fraction = s$fraction, kind = s$kind, replicate = r, D = diagnostic(parts, d$lab_value), cover_cv3 = as.numeric(abs(a$q - p$target) <= half))
}
part3_cells <- list()
for (i in which(grid$errors != "rare" & grid$rows == 3L)) part3_cells[[length(part3_cells) + 1L]] <- list(kind = "simulation", spec = as.list(grid[i, ]), index = i, analysis = grid$label[i], fraction = grid$fraction[i])
for (z in real) part3_cells[[length(part3_cells) + 1L]] <- list(kind = "real", z = z, fraction_index = 1L, analysis = z$id, fraction = config$fractions[1])

start <- Sys.time()
run_all <- function(n, f) {
  batch <- parallel::mclapply(seq_len(n), f, mc.cores = 14L, mc.preschedule = TRUE, mc.set.seed = FALSE)
  bad <- vapply(batch, inherits, logical(1), what = "try-error")
  if (any(bad)) stop("A worker failed: ", conditionMessage(attr(batch[[which(bad)[1]]], "condition")))
  do.call(rbind, batch)
}
write_gz <- function(x, name) { out <- gzfile(file.path(outdir, name), "wt"); utils::write.csv(x, out, row.names = FALSE); close(out) }
results1 <- do.call(rbind, lapply(seq_len(nrow(part1)), function(i) {
  got <- run_all(settings$part1_reps, function(r) part1_replicate(i, r))
  cat("part 1", part1$label[i], "; elapsed", round(as.numeric(difftime(Sys.time(), start, units = "secs"))), "seconds\n"); flush.console()
  got
}))
write_gz(results1, "part1.csv.gz")
results2 <- do.call(rbind, lapply(part2_cells, function(s) {
  got <- run_all(s$reps, function(r) part2_replicate(s, r))
  cat("part 2", s$analysis, s$fraction, "; elapsed", round(as.numeric(difftime(Sys.time(), start, units = "secs"))), "seconds\n"); flush.console()
  got
}))
write_gz(results2, "part2.csv.gz")
results3 <- do.call(rbind, lapply(part3_cells, function(s) run_all(settings$part3_reps, function(r) part3_replicate(s, r))))
write_gz(results3, "part3.csv.gz")
if (!identical(hashes_before, hash_sources())) stop("Source, plan, pool or reused results changed during the run; evidence is not verified.")
utils::write.csv(hashes_before, file.path(outdir, "source-hashes.csv"), row.names = FALSE)
writeLines(c(paste("Completed UTC:", format(Sys.time(), tz = "UTC", usetz = TRUE)), paste("Elapsed seconds:", round(as.numeric(difftime(Sys.time(), start, units = "secs")), 1)), capture.output(sessionInfo())), file.path(outdir, "session-info.txt"))
cat("Completed", mode, ":", nrow(results1), "+", nrow(results2), "+", nrow(results3), "requests;", paste(names(table(results2$status)), table(results2$status), collapse = ", "), "\n")
