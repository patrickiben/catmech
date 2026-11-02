#!/usr/bin/env Rscript
# Prospective test on a second health system; settings locked in validation/SECOND_SYSTEM_PLAN.md.
# Usage: Rscript tools/second-system-benchmark.R pilot|study [NEW_DIRECTORY] [PILOT_REPS]
args <- commandArgs(trailingOnly = TRUE)
mode <- if (length(args)) match.arg(args[1], c("pilot", "study")) else "pilot"
file_arg <- grep("^--file=", commandArgs(), value = TRUE)
root <- normalizePath(file.path(dirname(sub("^--file=", "", file_arg[1])), ".."))
outdir <- if (length(args) >= 2L) args[2] else file.path(root, "validation", paste0("second-system-", mode))
if (dir.exists(outdir) && length(list.files(outdir, all.files = TRUE, no.. = TRUE)))
  stop("Use a new or empty output directory; existing evidence is never overwritten.")
hash_sources <- function() {
  paths <- c(sort(list.files(file.path(root, "R"), "\\.R$", full.names = TRUE)),
    file.path(root, c("DESCRIPTION", "NAMESPACE", "tools/second-system-benchmark.R", "tools/second-system-score.py", "tools/second-system-prepare.R",
                      "validation/SECOND_SYSTEM_PLAN.md", "validation/second-system/zambia-pool.csv", "validation/second-system/raw/Zambia_data.xlsx")))
  data.frame(file = substring(paths, nchar(root) + 2L),
    sha256 = vapply(paths, function(p) digest::digest(file = p, algo = "sha256"), character(1)), row.names = NULL)
}
hashes_before <- hash_sources()
pkgload::load_all(root, quiet = TRUE)
suppressPackageStartupMessages(library(mice))
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
pilot <- mode == "pilot"
pilot_reps <- if (length(args) >= 3L) as.integer(args[3]) else 5L
SEEDS <- if (pilot) list(cohort = 1110000000, audit = 1120000000, imputation = 1130000000, mice = 1140000000) else
  list(cohort = 2010000000, audit = 2020000000, imputation = 2030000000, mice = 2040000000)
settings <- list(mode = mode, m = 20L, reps = if (pilot) pilot_reps else 1000L, sensitivity_reps = if (pilot) pilot_reps else 500L, seeds = SEEDS)
jsonlite::write_json(settings, file.path(outdir, "settings.json"), pretty = TRUE, auto_unbox = TRUE)
pool_all <- utils::read.csv(file.path(root, "validation/second-system/zambia-pool.csv"), stringsAsFactors = FALSE)
ITEMS <- c("opd", "deliveries", "anc", "postnatal", "immunization", "ipt", "fpmethods", "counselling", "exposedbabies")
AUDITS <- list(clinic_25 = list(unit = "clinic", fraction = .25), clinic_50 = list(unit = "clinic", fraction = .5),
               district_33 = list(unit = "district", clusters = 6L), district_50 = list(unit = "district", clusters = 9L))
mice.impute.2l.lmer.all <- function(y, ry, x, type, wy = NULL, intercept = TRUE, eps = NULL, ...) mice::mice.impute.2l.lmer(y, ry, x, type, wy, intercept)

draw_cohort <- function(pool, item_index, r) {
  set.seed(as.integer(SEEDS$cohort + item_index * 1e6 + r))
  by_district <- split(seq_len(nrow(pool)), pool$district)
  names_all <- sort(names(by_district))
  picked <- names_all[sample.int(length(names_all), length(names_all), replace = TRUE)]
  rows <- unlist(by_district[picked], use.names = FALSE)
  ids <- rep(seq_along(picked), lengths(by_district[picked]))
  list(rows = rows, ids = ids)
}
draw_audit <- function(cohort, audit, audit_index, item_index, r) {
  set.seed(as.integer(SEEDS$audit + item_index * 1e6 + audit_index * 1e5 + r))
  n <- length(cohort$rows)
  if (audit$unit == "clinic") {
    seq_len(n) %in% sample.int(n, max(1L, as.integer(round(audit$fraction * n))))
  } else {
    G <- max(cohort$ids)
    cohort$ids %in% sample.int(G, min(G, audit$clusters))
  }
}
pool_fits <- function(values, x, ids, df_only = FALSE) {
  fits <- lapply(values, function(v) { d <- data.frame(y = v, covariate = x); stats::lm(y ~ covariate, data = d, na.action = stats::na.fail) })
  list(q = vapply(fits, function(f) unname(stats::coef(f)[2]), numeric(1)),
       u = vapply(fits, function(f) catmech_cluster_jackknife(f, ids)[2, 2], numeric(1)))
}
as_interval <- function(estimate, se, df, target) {
  half <- stats::qt(.975, df) * se
  list(estimate = estimate, se = se, df = df, width = 2 * half, cover = as.numeric(abs(estimate - target) <= half))
}
one_value <- function(y, x, ids, target) {
  if (length(unique(ids)) < 2L) stop("Fewer than two clusters.")
  f <- pool_fits(list(y), x, ids)
  if (!is.finite(f$u) || f$u <= 0) stop("no residual variance")
  as_interval(f$q, sqrt(f$u), length(unique(ids)) - 1L, target)
}
fit_method <- function(method, d, target, seed_imputation, seed_mice) {
  switch(method,
    oracle = one_value(d$truth, d$covariate, d$cluster_id, target),
    recorded = one_value(d$lab_value, d$covariate, d$cluster_id, target),
    audited_only = one_value(d$truth[d$audited], d$covariate[d$audited], d$cluster_id[d$audited], target),
    linear_mi = , random_intercept = {
      input <- data.frame(covariate = d$covariate, cluster_id = d$cluster_id, lab_value = d$lab_value, gold_value = ifelse(d$audited, d$truth, NA_real_))
      options <- list(input, "lab_value", "gold_value", "covariate", method = "normal", m = settings$m, seed = seed_imputation, use_recorded = TRUE)
      if (method == "random_intercept") options$cluster <- "cluster_id"
      fit <- do.call(catmech_reference, options)
      row <- catmech_pool_lm(fit, lab_value ~ covariate, cluster = "cluster_id")
      row <- row[row$term == "covariate", ]
      out <- as_interval(row$estimate, row$std.error, row$df, target)
      out$between <- row$between_variance
      out$kept <- as.numeric(grepl("_recorded_kept$", row$interpretation))
      out
    },
    mice_2l_lmer = {
      input <- data.frame(covariate = d$covariate, lab_value = d$lab_value, gold_value = ifelse(d$audited, d$truth, NA_real_), cluster_id = d$cluster_id)
      use <- names(input)
      matrix_p <- matrix(0L, 4L, 4L, dimnames = list(use, use))
      matrix_p["gold_value", c("covariate", "lab_value")] <- 1L; matrix_p["gold_value", "cluster_id"] <- -2L
      chosen <- stats::setNames(rep("", 4L), use); chosen["gold_value"] <- "2l.lmer.all"
      lmer_failed <- FALSE
      fit <- withCallingHandlers(mice::mice(input, m = settings$m, maxit = 1, method = chosen, predictorMatrix = matrix_p, printFlag = FALSE, seed = seed_mice, eps = 0),
        warning = function(w) { if (grepl("lmer does not run", conditionMessage(w))) lmer_failed <<- TRUE; invokeRestart("muffleWarning") },
        message = function(m) invokeRestart("muffleMessage"))
      if (lmer_failed) stop("lmer does not run.")
      f <- pool_fits(lapply(seq_len(settings$m), function(j) mice::complete(fit, j)$gold_value), d$covariate, d$cluster_id)
      pooled <- mice::pool.scalar(f$q, f$u, n = length(unique(d$cluster_id)), k = 1)
      as_interval(pooled$qbar, sqrt(pooled$t), pooled$df, target)
    })
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
METHODS <- c("oracle", "recorded", "audited_only", "linear_mi", "random_intercept", "mice_2l_lmer")
SENSITIVITY <- c("oracle", "linear_mi", "random_intercept")
one_replicate <- function(item_index, audit_index, r, pool, target, methods, run) {
  cohort <- draw_cohort(pool, item_index, r)
  audited <- draw_audit(cohort, AUDITS[[audit_index]], audit_index, item_index, r)
  d <- data.frame(covariate = pool$covariate[cohort$rows], lab_value = pool$reported[cohort$rows], truth = pool$verified[cohort$rows], cluster_id = cohort$ids, audited = audited)
  seed_i <- as.integer(SEEDS$imputation + item_index * 1e6 + audit_index * 1e5 + r)
  seed_m <- as.integer(SEEDS$mice + item_index * 1e6 + audit_index * 1e5 + r)
  do.call(rbind, lapply(methods, function(method) {
    warned <- character()
    started <- Sys.time()
    got <- limited(withCallingHandlers(fit_method(method, d, target, seed_i, seed_m), warning = function(w) { warned <<- c(warned, conditionMessage(w)); invokeRestart("muffleWarning") }))
    seconds <- as.numeric(difftime(Sys.time(), started, units = "secs"))
    failed <- inherits(got, "error")
    data.frame(run = run, item = ITEMS[item_index], audit = names(AUDITS)[audit_index], replicate = r, method = method, status = if (failed) "failed" else "ok", message = if (failed) substr(conditionMessage(got), 1L, 100L) else "",
               estimate = if (failed) NA_real_ else got$estimate, se = if (failed) NA_real_ else got$se, df = if (failed) NA_real_ else got$df, width = if (failed) NA_real_ else got$width,
               cover = if (failed) NA_real_ else got$cover, target = target, kept = if (failed || is.null(got$kept)) 0 else got$kept, audited_clusters = length(unique(d$cluster_id[d$audited])), audited_rows = sum(d$audited), seconds = round(seconds, 4), stringsAsFactors = FALSE)
  }))
}
target_of <- function(pool) unname(stats::coef(stats::lm(verified ~ covariate, data = pool))[2])
start <- Sys.time(); records <- list(); targets <- list()
for (run in c("primary", "sensitivity")) {
  methods <- if (run == "primary") METHODS else SENSITIVITY
  reps <- if (run == "primary") settings$reps else settings$sensitivity_reps
  for (i in seq_along(ITEMS)) {
    pool <- pool_all[pool_all$item == ITEMS[i], ]
    if (run == "primary") pool <- pool[pool$both_zero == 0L, ]
    target <- target_of(pool)
    targets[[length(targets) + 1L]] <- data.frame(run = run, item = ITEMS[i], clinics = nrow(pool), districts = length(unique(pool$district)), target = target, recorded_slope = unname(stats::coef(stats::lm(reported ~ covariate, data = pool))[2]))
    for (a in seq_along(AUDITS)) {
      batch <- parallel::mclapply(seq_len(reps), function(r) one_replicate(i, a, r, pool, target, methods, run), mc.cores = 14L, mc.preschedule = TRUE, mc.set.seed = FALSE)
      bad <- vapply(batch, inherits, logical(1), what = "try-error")
      if (any(bad)) stop("A worker failed on ", ITEMS[i], ": ", conditionMessage(attr(batch[[which(bad)[1]]], "condition")))
      records <- c(records, batch)
    }
    cat(run, ITEMS[i], "; elapsed", round(as.numeric(difftime(Sys.time(), start, units = "secs"))), "seconds\n"); flush.console()
  }
}
if (!identical(hashes_before, hash_sources())) stop("Source, plan or pool changed during the run; evidence is not verified.")
results <- do.call(rbind, records)
out <- gzfile(file.path(outdir, "replicates.csv.gz"), "wt")
utils::write.csv(results, out, row.names = FALSE)
close(out)
utils::write.csv(do.call(rbind, targets), file.path(outdir, "targets.csv"), row.names = FALSE)
# The error's intraclass correlation by district, by moments, on each item's pool. Description only.
icc <- do.call(rbind, lapply(ITEMS, function(item) {
  p <- pool_all[pool_all$item == item & pool_all$both_zero == 0L, ]
  e <- p$reported - p$verified; g <- factor(p$district)
  a <- stats::anova(stats::lm(e ~ g)); k <- (length(e) - sum(table(g)^2) / length(e)) / (nlevels(g) - 1)
  msb <- a[["Mean Sq"]][1]; msw <- a[["Mean Sq"]][2]
  data.frame(item = item, clinics = length(e), share_disagree = mean(e != 0), mean_error = mean(e), sd_error = stats::sd(e), error_icc = (msb - msw) / (msb + (k - 1) * msw), outcome_mean = mean(p$verified))
}))
utils::write.csv(icc, file.path(outdir, "error-icc.csv"), row.names = FALSE)
utils::write.csv(hashes_before, file.path(outdir, "source-hashes.csv"), row.names = FALSE)
writeLines(c(paste("Completed UTC:", format(Sys.time(), tz = "UTC", usetz = TRUE)), paste("Elapsed seconds:", round(as.numeric(difftime(Sys.time(), start, units = "secs")), 1)), capture.output(sessionInfo())), file.path(outdir, "session-info.txt"))
cat("Completed", mode, ":", nrow(results), "requests;", paste(names(table(results$status)), table(results$status), collapse = ", "), "\n")
