#!/usr/bin/env Rscript
# Prospective comparison of catmech_reference(cluster = ) with outside methods; settings locked in validation/OUTSIDE_METHODS_PLAN.md.
# Cohorts come from the earlier studies' own functions and seeds. Saved oracle, linear MI and option results are reused, not rerun.
# Usage: Rscript tools/outside-methods-benchmark.R check|pilot|study [NEW_DIRECTORY] [PILOT_REPS]
args <- commandArgs(trailingOnly = TRUE)
mode <- if (length(args)) match.arg(args[1], c("check", "pilot", "study")) else "pilot"
file_arg <- grep("^--file=", commandArgs(), value = TRUE)
root <- normalizePath(file.path(dirname(sub("^--file=", "", file_arg[1])), ".."))
outdir <- if (length(args) >= 2L) args[2] else file.path(root, "validation", paste0("outside-methods-", mode))
if (dir.exists(outdir) && length(list.files(outdir, all.files = TRUE, no.. = TRUE)))
  stop("Use a new or empty output directory; existing evidence is never overwritten.")
hash_sources <- function() {
  paths <- c(sort(list.files(file.path(root, "R"), "\\.R$", full.names = TRUE)),
    file.path(root, c("DESCRIPTION", "NAMESPACE", "tools/outside-methods-benchmark.R", "tools/outside-methods-score.py", "tools/variance-mi-benchmark.R", "tools/clustered-fix-benchmark.R",
                      "validation/OUTSIDE_METHODS_PLAN.md", "validation/real-audit/rwanda-pool.csv", "validation/variance-mi-study/replicates.csv.gz")))
  data.frame(file = substring(paths, nchar(root) + 2L),
    sha256 = vapply(paths, function(p) digest::digest(file = p, algo = "sha256"), character(1)), row.names = NULL)
}
hashes_before <- hash_sources()
pkgload::load_all(root, quiet = TRUE)
suppressPackageStartupMessages({library(survey); library(mice)})
options(survey.lonely.psu = "adjust")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
pilot <- mode == "pilot"
pilot_reps <- if (length(args) >= 3L) as.integer(args[3]) else 3L
config <- jsonlite::read_json(file.path(root, "validation/clustered-fix-study/settings.json"), simplifyVector = TRUE)
earlier <- jsonlite::read_json(file.path(root, "validation/variance-mi-study/settings.json"), simplifyVector = TRUE)
if (pilot) {
  config$held_out <- list(cohort = 1610100000, selection = 1620100000, imputation = 1630000000)
  config$development <- list(cohort = 1640100000, selection = 1650100000, imputation = 1660000000)
  config$simulation <- list(cohort = 1600100000, imputation = 1670000000)
}
expressions <- parse(file.path(root, "tools/clustered-fix-benchmark.R"))
assigned <- vapply(expressions, function(e) if (is.call(e) && identical(e[[1]], as.name("<-")) && is.name(e[[2]])) as.character(e[[2]]) else "", character(1))
first <- which(assigned == "rwanda")[1]; last <- which(assigned == "real_case")[1]
for (e in expressions[first:last]) eval(e, envir = globalenv())
SEEDS <- if (pilot) list(covariate_cohort = 1680100000, covariate_imputation = 1690000000, lmer = 1550000000, pan = 1570000000) else
  list(covariate_cohort = earlier$seeds$covariate_cohort, covariate_imputation = earlier$seeds$covariate_imputation, lmer = 1950000000, pan = 1970000000)
runner <- parse(file.path(root, "tools/variance-mi-benchmark.R"))
named <- vapply(runner, function(e) if (is.call(e) && identical(e[[1]], as.name("<-")) && is.name(e[[2]])) as.character(e[[2]]) else "", character(1))
eval(runner[[which(named == "simulate_covariate_cell")[1]]], envir = globalenv())
grid <- utils::read.csv(file.path(root, "validation/clustered-fix-study/simulation-grid.csv"), stringsAsFactors = FALSE)
covariate_grid <- expand.grid(design = "rwanda", errors = c("shift", "under", "row_only"), fraction = config$fractions, stringsAsFactors = FALSE)
covariate_grid$rows <- 3L
settings <- list(mode = mode, m = 20L, reps_real_25 = if (pilot) pilot_reps else 2000L, reps_other = if (pilot) pilot_reps else 500L, control_reps = if (mode == "study") 10L else 0L, seeds = SEEDS)
jsonlite::write_json(settings, file.path(outdir, "settings.json"), pretty = TRUE, auto_unbox = TRUE)

# The analysis design with `value` standing in for the audited quantity; the analysis formula always names it lab_value.
design_for <- function(p, value, rows = seq_len(nrow(p$input))) {
  d <- p$input[rows, , drop = FALSE]; d$lab_value <- value[rows]
  frame <- stats::model.frame(p$formula, data = d, na.action = stats::na.fail)
  list(x = stats::model.matrix(p$formula, data = frame), y = stats::model.response(frame), ids = d$cluster_id)
}
weighted_fit <- function(a, k, w = rep(1, length(a$y))) {
  f <- stats::lm.wfit(a$x, a$y, w)
  if (f$rank < ncol(a$x)) stop("The audited rows do not have full column rank.")
  bread <- chol2inv(qr.R(f$qr))
  bread[f$qr$pivot, f$qr$pivot] <- bread
  e <- a$y - as.vector(a$x %*% f$coefficients)
  list(estimate = unname(f$coefficients[k]), row_score = as.vector(a$x %*% bread[, k]) * e, residual = e)
}
term_of <- function(p) if (is.null(p$term)) "covariate" else p$term
strata_of <- function(p) if ("stratum" %in% names(p$input)) as.character(p$input$stratum) else rep("all", nrow(p$input))
# Clusters by stratum, and d = N_h / m_h for each cluster.
cluster_design <- function(p) {
  audited <- p$audited; ids <- p$input$cluster_id; strata <- strata_of(p)
  first_row <- !duplicated(ids)
  units <- data.frame(id = ids[first_row], stratum = strata[first_row], audited = audited[first_row], stringsAsFactors = FALSE)
  N <- table(units$stratum); m <- table(units$stratum[units$audited])
  units$N_h <- as.numeric(N[units$stratum]); units$m_h <- as.numeric(m[units$stratum]); units$d <- units$N_h / units$m_h
  units
}
interval <- function(estimate, se, df, target) {
  half <- stats::qt(.975, df) * se
  list(estimate = estimate, se = se, df = df, width = 2 * half, cover = as.numeric(abs(estimate - target) <= half))
}

fit_recorded <- function(p) {
  a <- design_for(p, p$input$lab_value)
  lmfit <- stats::lm(a$y ~ a$x - 1); names(lmfit$coefficients) <- colnames(a$x)
  if (anyNA(lmfit$coefficients)) stop("The recorded rows do not have full column rank.")
  V <- catmech_cluster_jackknife(lmfit, a$ids); dimnames(V) <- list(colnames(a$x), colnames(a$x))
  G <- length(unique(a$ids)); term <- term_of(p)
  interval(unname(lmfit$coefficients[term]), sqrt(V[term, term]), G - 1, p$target)
}
fit_audited_only <- function(p) {
  a <- design_for(p, p$input$gold_value, which(p$audited))
  lmfit <- stats::lm(a$y ~ a$x - 1); names(lmfit$coefficients) <- colnames(a$x)
  if (anyNA(lmfit$coefficients)) stop("The audited rows do not have full column rank.")
  V <- catmech_cluster_jackknife(lmfit, a$ids); dimnames(V) <- list(colnames(a$x), colnames(a$x))
  G <- length(unique(a$ids)); term <- term_of(p)
  if (!is.finite(V[term, term]) || V[term, term] <= 0) stop("The audited rows have no residual variance.")
  interval(unname(lmfit$coefficients[term]), sqrt(V[term, term]), G - 1, p$target)
}
difference_parts <- function(p) {
  units <- cluster_design(p); term <- term_of(p)
  rows <- which(p$audited); weights <- units$d[match(p$input$cluster_id[rows], units$id)]
  all_w <- design_for(p, p$input$lab_value)
  a_w <- design_for(p, p$input$lab_value, rows); a_r <- design_for(p, p$input$gold_value, rows)
  k <- match(term, colnames(all_w$x))
  f_all <- weighted_fit(all_w, k); f_w <- weighted_fit(a_w, k, weights); f_r <- weighted_fit(a_r, k, weights)
  ids <- as.character(p$input$cluster_id[rows])
  psi_r <- rowsum(f_r$row_score, ids, reorder = FALSE)[, 1]
  delta <- psi_r - rowsum(f_w$row_score, ids, reorder = FALSE)[, 1]
  a_units <- units[units$audited, ][match(names(psi_r), as.character(units$id[units$audited])), ]
  G <- nrow(units); H <- length(unique(units$stratum))
  phase_two <- function(score) {
    scaled <- a_units$d * score
    out <- 0
    for (h in unique(a_units$stratum)) {
      i <- a_units$stratum == h; f_h <- a_units$m_h[i][1] / a_units$N_h[i][1]; m_h <- sum(i)
      out <- out + (1 - f_h) * if (m_h > 1L) (m_h / (m_h - 1)) * sum((scaled[i] - mean(scaled[i]))^2) else scaled[i]^2
    }
    out
  }
  V1 <- (G / (G - 1)) * sum(a_units$d * psi_r^2)
  V2 <- phase_two(delta)
  list(estimate = f_all$estimate + f_r$estimate - f_w$estimate, V1 = V1, V2 = V2, G = G, df = max(1, sum(units$audited) - H), psi_r = psi_r, delta = delta, ids = names(psi_r))
}
fit_difference <- function(p) {
  z <- difference_parts(p)
  if (!is.finite(z$V1 + z$V2) || z$V1 + z$V2 <= 0) stop("The audited rows have no residual variance.")
  interval(z$estimate, sqrt(z$V1 + z$V2), z$df, p$target)
}
# Least squares on fitted values; exploratory, accuracy only.
fitted_value_parts <- function(p) {
  units <- cluster_design(p); rows <- which(p$audited); term <- term_of(p)
  weights <- units$d[match(p$input$cluster_id[rows], units$id)]
  predictors <- stats::model.matrix(stats::reformulate(c("lab_value", p$predictors)), data = p$input)
  f <- stats::lm.wfit(predictors[rows, , drop = FALSE], p$input$gold_value[rows], weights)
  if (f$rank < ncol(predictors)) stop("The audited rows do not have full column rank.")
  proxy <- as.vector(predictors %*% f$coefficients)
  all_p <- design_for(p, proxy); k <- match(term, colnames(all_p$x))
  list(estimate = weighted_fit(all_p, k)$estimate, correction = weighted_fit(design_for(p, p$input$gold_value - proxy, rows), k, weights)$estimate)
}
raking_design <- function(p) {
  d <- p$input; d$audited <- p$audited
  if (!"stratum" %in% names(d)) d$stratum <- factor("all")
  all_w <- design_for(p, p$input$lab_value); k <- match(term_of(p), colnames(all_w$x))
  d$if_k <- weighted_fit(all_w, k)$row_score
  d$value_r <- d$gold_value
  list(d = d, design = survey::twophase(id = list(~cluster_id, ~cluster_id), strata = list(NULL, if (nlevels(d$stratum) > 1L) ~stratum else NULL), subset = ~audited, data = d))
}
fit_raking <- function(p) {
  r <- raking_design(p)
  calibrated <- survey::calibrate(r$design, phase = 2, formula = ~if_k, calfun = "raking")
  formula <- stats::update(p$formula, if (identical(all.vars(p$formula)[1], "lab_value")) value_r ~ . else . ~ . - lab_value + value_r)
  f <- survey::svyglm(formula, design = calibrated)
  term <- if (identical(term_of(p), "lab_value")) "value_r" else term_of(p)
  interval(unname(stats::coef(f)[term]), unname(survey::SE(f)[term]), survey::degf(calibrated), p$target)
}
fit_mice <- function(p, method, seed) {
  use <- unique(c(p$predictors, "lab_value", "gold_value", "cluster_id"))
  d <- p$input[, use, drop = FALSE]
  matrix_p <- matrix(0L, length(use), length(use), dimnames = list(use, use))
  matrix_p["gold_value", setdiff(use, c("gold_value", "cluster_id"))] <- 1L; matrix_p["gold_value", "cluster_id"] <- -2L
  chosen <- stats::setNames(rep("", length(use)), use); chosen["gold_value"] <- method
  lmer_failed <- FALSE; singular <- 0L
  fit <- withCallingHandlers(mice::mice(d, m = settings$m, maxit = 1, method = chosen, predictorMatrix = matrix_p, printFlag = FALSE, seed = seed),
    warning = function(w) { if (grepl("lmer does not run", conditionMessage(w))) lmer_failed <<- TRUE; invokeRestart("muffleWarning") },
    message = function(m) { if (grepl("singular", conditionMessage(m))) singular <<- singular + 1L; invokeRestart("muffleMessage") })
  if (lmer_failed) stop("lmer does not run.")
  term <- term_of(p); ids <- p$input$cluster_id; G <- length(unique(ids))
  pieces <- lapply(seq_len(settings$m), function(j) {
    a <- design_for(p, mice::complete(fit, j)$gold_value)
    lmfit <- stats::lm(a$y ~ a$x - 1); names(lmfit$coefficients) <- colnames(a$x)
    V <- catmech_cluster_jackknife(lmfit, ids); dimnames(V) <- list(colnames(a$x), colnames(a$x))
    c(q = unname(lmfit$coefficients[term]), u = V[term, term])
  })
  pooled <- mice::pool.scalar(vapply(pieces, `[[`, numeric(1), "q"), vapply(pieces, `[[`, numeric(1), "u"), n = G, k = 1)
  out <- interval(pooled$qbar, sqrt(pooled$t), pooled$df, p$target)
  out$logged <- if (is.null(fit$loggedEvents)) 0L else nrow(fit$loggedEvents)
  out$singular <- singular
  out
}
METHODS <- c("recorded", "audited_only", "difference", "raking", "mice_2l_lmer", "mice_2l_pan")
run_method <- function(method, p, s, r) {
  switch(method,
    recorded = fit_recorded(p), audited_only = fit_audited_only(p), difference = fit_difference(p), raking = fit_raking(p),
    mice_2l_lmer = fit_mice(p, "2l.lmer", as.integer(SEEDS$lmer + s$number * 1e5 + r)),
    mice_2l_pan = fit_mice(p, "2l.pan", as.integer(SEEDS$pan + s$number * 1e5 + r)))
}
cohort <- function(s, r) {
  if (s$kind == "real") real_case(s$z, s$fraction_index, r) else if (s$kind == "simulation")
    simulate_cell(s$spec, as.integer(config$simulation$cohort + s$index * 1e5 + r), s$fraction) else
    simulate_covariate_cell(s$spec, as.integer(SEEDS$covariate_cohort + s$index * 1e5 + r), s$fraction)
}
catmech_seed <- function(s, r) {
  if (s$kind == "real") as.integer(s$z$seeds$imputation + s$z$index * 1e6 + s$fraction_index * 1e5 + r) else if (s$kind == "simulation")
    as.integer(config$simulation$imputation + s$index * 1e5 + r) else as.integer(SEEDS$covariate_imputation + s$index * 1e5 + r)
}
one_replicate <- function(s, r) {
  p <- cohort(s, r)
  fitted_estimate <- tryCatch(fitted_value_parts(p)$estimate, error = function(e) NA_real_)
  do.call(rbind, lapply(METHODS, function(method) {
    warned <- 0L
    started <- Sys.time()
    got <- tryCatch(withCallingHandlers(run_method(method, p, s, r), warning = function(w) { warned <<- warned + 1L; invokeRestart("muffleWarning") }), error = identity)
    seconds <- as.numeric(difftime(Sys.time(), started, units = "secs"))
    failed <- inherits(got, "error")
    status <- if (!failed) "ok" else if (grepl("full column rank|no residual variance", conditionMessage(got))) "refused" else "failed"
    data.frame(level = s$level, analysis = s$analysis, fraction = s$fraction, replicate = r, method = method, status = status, message = if (failed) substr(conditionMessage(got), 1L, 120L) else "",
               warned = warned, estimate = if (failed) NA_real_ else got$estimate, se = if (failed) NA_real_ else got$se, df = if (failed) NA_real_ else got$df, width = if (failed) NA_real_ else got$width,
               cover = if (failed) NA_real_ else got$cover, target = p$target, logged = if (failed || is.null(got$logged)) 0L else got$logged, singular = if (failed || is.null(got$singular)) 0L else got$singular, fitted_estimate = if (method == "difference") fitted_estimate else NA_real_,
               seconds = round(seconds, 4), stringsAsFactors = FALSE)
  }))
}

cells <- list()
for (z in real) for (j in seq_along(config$fractions)) cells[[length(cells) + 1L]] <- list(kind = "real", z = z, fraction_index = j, level = z$level, analysis = z$id, fraction = config$fractions[j],
  reps = if (j == 1L) settings$reps_real_25 else settings$reps_other, number = 100L + length(cells) + 1L)
for (i in seq_len(nrow(grid))) cells[[length(cells) + 1L]] <- list(kind = "simulation", spec = as.list(grid[i, ]), index = i, level = "simulation", analysis = grid$label[i], fraction = grid$fraction[i], reps = settings$reps_other, number = i)
for (i in seq_len(nrow(covariate_grid))) cells[[length(cells) + 1L]] <- list(kind = "covariate", spec = as.list(covariate_grid[i, ]), index = i, level = "imputed_covariate", analysis = paste0("covariate_", covariate_grid$errors[i]), fraction = covariate_grid$fraction[i], reps = settings$reps_other, number = 200L + i)

# The plan's checks on a fixed dataset, before any study run.
if (mode == "check") {
  lines <- character()
  note <- function(...) { lines <<- c(lines, paste0(...)); cat(paste0(..., "\n")) }
  tests <- list(list(i = 1L, r = 900001L), list(i = 7L, r = 900002L), list(i = 13L, r = 900003L), list(i = 2L, r = 900004L))
  worst <- c(weights = 0, estimate = 0, phase1 = 0, phase2 = 0, correction = 0)
  for (t in tests) {
    s <- list(kind = "simulation", spec = as.list(grid[t$i, ]), index = t$i, fraction = grid$fraction[t$i])
    p <- cohort(s, t$r)
    r <- raking_design(p); z <- difference_parts(p); units <- cluster_design(p)
    expected <- units$d[match(r$d$cluster_id[r$d$audited], units$id)]
    worst["weights"] <- max(worst["weights"], max(abs(stats::weights(r$design) - expected)))
    # Point estimate: recorded fit on all rows plus the design-weighted fit of (reference - recorded) from survey.
    dd <- r$d; dd$delta_value <- dd$gold_value - dd$lab_value
    twop <- survey::twophase(id = list(~cluster_id, ~cluster_id), strata = list(NULL, if (nlevels(dd$stratum) > 1L) ~stratum else NULL), subset = ~audited, data = dd)
    via_survey <- unname(stats::coef(stats::lm(p$formula, data = p$input))["covariate"]) + unname(stats::coef(survey::svyglm(stats::update(p$formula, delta_value ~ .), design = twop))["covariate"])
    worst["estimate"] <- max(worst["estimate"], abs(via_survey - z$estimate))
    # Variance components: survey's phase one on psi_R and phase two on delta, at cluster level.
    cl <- data.frame(id = units$id, stratum = units$stratum, audited = units$audited, stringsAsFactors = FALSE)
    cl$psi_r <- z$psi_r[match(as.character(cl$id), z$ids)]; cl$delta <- z$delta[match(as.character(cl$id), z$ids)]
    cl$stratum <- factor(cl$stratum)
    twoc <- survey::twophase(id = list(~id, ~id), strata = list(NULL, if (nlevels(cl$stratum) > 1L) ~stratum else NULL), subset = ~audited, data = cl)
    v1 <- attr(stats::vcov(survey::svytotal(~psi_r, twoc)), "phases")$phase1
    v2 <- attr(stats::vcov(survey::svytotal(~delta, twoc)), "phases")$phase2
    # The plan's phase one has the factor G / (G - 1), which survey's does not.
    worst["phase1"] <- max(worst["phase1"], abs(as.numeric(v1) * (z$G / (z$G - 1)) / z$V1 - 1))
    worst["phase2"] <- max(worst["phase2"], abs(as.numeric(v2) / z$V2 - 1))
    worst["correction"] <- max(worst["correction"], abs(fitted_value_parts(p)$correction))
  }
  note("Check of the plan's fixed-dataset conditions, on four simulated cohorts with seeds used nowhere else.")
  note("Largest gap: twophase weights against N_h / m_h ", format(worst["weights"], digits = 3), "; difference estimate against survey ", format(worst["estimate"], digits = 3),
       "; relative gap in phase one variance, apart from the factor G / (G - 1) in the plan's form, ", format(worst["phase1"], digits = 3), "; relative gap in phase two variance ", format(worst["phase2"], digits = 3), "; fitted_value correction term ", format(worst["correction"], digits = 3), ".")
  writeLines(lines, file.path(outdir, "check.txt"))
  quit(save = "no")
}

start <- Sys.time(); records <- list()
for (s in cells) {
  batch <- parallel::mclapply(seq_len(s$reps), function(r) one_replicate(s, r), mc.cores = 14L, mc.preschedule = TRUE, mc.set.seed = FALSE)
  bad <- vapply(batch, inherits, logical(1), what = "try-error")
  if (any(bad)) stop("A worker failed on ", s$analysis, ": ", conditionMessage(attr(batch[[which(bad)[1]]], "condition")))
  records <- c(records, batch)
  cat(s$level, s$analysis, s$fraction, "; elapsed", round(as.numeric(difftime(Sys.time(), start, units = "secs"))), "seconds\n"); flush.console()
}
results <- do.call(rbind, records)
out <- gzfile(file.path(outdir, "replicates.csv.gz"), "wt")
utils::write.csv(results, out, row.names = FALSE)
close(out)

# Control on the first cohorts of every cell: the option run afresh must reproduce the saved study.
if (settings$control_reps) {
  saved <- utils::read.csv(gzfile(file.path(root, "validation/variance-mi-study/replicates.csv.gz")), stringsAsFactors = FALSE)
  saved <- saved[saved$imputer == "random_intercept" & saved$replicate <= settings$control_reps, ]
  control <- do.call(rbind, lapply(cells, function(s) do.call(rbind, parallel::mclapply(seq_len(settings$control_reps), function(r) {
    p <- cohort(s, r)
    started <- Sys.time()
    got <- tryCatch(suppressWarnings(catmech_pool_lm(catmech_reference(p$input, "lab_value", "gold_value", p$predictors, method = "normal", m = settings$m, seed = catmech_seed(s, r), use_recorded = TRUE, cluster = "cluster_id"), p$formula, cluster = "cluster_id")), error = identity)
    seconds <- as.numeric(difftime(Sys.time(), started, units = "secs"))
    row <- if (inherits(got, "error")) NULL else got[got$term == term_of(p), ]
    data.frame(analysis = s$analysis, fraction = s$fraction, replicate = r, status = if (is.null(row)) "failed" else "ok", estimate = if (is.null(row)) NA_real_ else row$estimate, se = if (is.null(row)) NA_real_ else row$std.error, df = if (is.null(row)) NA_real_ else row$df, target = p$target, seconds = round(seconds, 4))
  }, mc.cores = 14L, mc.preschedule = TRUE, mc.set.seed = FALSE))))
  utils::write.csv(control, file.path(outdir, "control.csv"), row.names = FALSE)
  m <- saved[match(paste(control$analysis, control$fraction, control$replicate), paste(saved$analysis, saved$fraction, saved$replicate)), ]
  ok <- control$status == "ok"
  gaps <- c(estimate = max(abs(control$estimate - m$estimate)[ok]), se = max(abs(control$se - m$se_cv3)[ok]), target = max(abs(control$target - m$target)), statuses = sum(ok != (m$status == "ok")))
  floored <- ok & 1.05 * m$between / m$se_cv3^2 < 1e-4
  floored[is.na(floored)] <- FALSE
  df_gap <- max(abs(control$df - m$df_cv3)[ok & !floored], 0)
  writeLines(c(paste0("Control: the option run afresh through the package on the first ", settings$control_reps, " cohorts of every cell, ", nrow(control), " requests."),
    paste0("Largest gap against the saved study: estimate ", format(gaps["estimate"], digits = 3), ", standard error ", format(gaps["se"], digits = 3), ", degrees of freedom ", format(df_gap, digits = 3), ", target ", format(gaps["target"], digits = 3), "; status mismatches ", gaps["statuses"], ".")),
    file.path(outdir, "control.txt"))
  cat(readLines(file.path(outdir, "control.txt")), sep = "\n"); cat("\n")
  if (gaps["estimate"] > 1e-8 || gaps["se"] > 1e-8 || df_gap > 1e-6 || gaps["target"] > 1e-10 || gaps["statuses"]) stop("The control failed; results are kept but are not verified.")
}
if (!identical(hashes_before, hash_sources())) stop("Source, plan, pool or reused results changed during the run; evidence is not verified.")
utils::write.csv(hashes_before, file.path(outdir, "source-hashes.csv"), row.names = FALSE)
writeLines(c(paste("Completed UTC:", format(Sys.time(), tz = "UTC", usetz = TRUE)), paste("Elapsed seconds:", round(as.numeric(difftime(Sys.time(), start, units = "secs")), 1)), capture.output(sessionInfo())), file.path(outdir, "session-info.txt"))
cat("Completed", mode, ":", nrow(results), "requests;", paste(names(table(results$status)), table(results$status), collapse = ", "), "\n")
