#!/usr/bin/env Rscript
# Post hoc rerun of two comparator settings; see validation/OUTSIDE_METHODS_SENSITIVITY_PLAN.md. Same cohorts and seeds as the comparison.
# Usage: Rscript tools/outside-methods-sensitivity.R [NEW_DIRECTORY]
args <- commandArgs(trailingOnly = TRUE)
file_arg <- grep("^--file=", commandArgs(), value = TRUE)
root <- normalizePath(file.path(dirname(sub("^--file=", "", file_arg[1])), ".."))
outdir <- if (length(args)) args[1] else file.path(root, "validation", "outside-methods-sensitivity")
if (dir.exists(outdir) && length(list.files(outdir, all.files = TRUE, no.. = TRUE)))
  stop("Use a new or empty output directory; existing evidence is never overwritten.")
hash_sources <- function() {
  paths <- c(sort(list.files(file.path(root, "R"), "\\.R$", full.names = TRUE)),
    file.path(root, c("DESCRIPTION", "NAMESPACE", "tools/outside-methods-sensitivity.R", "tools/outside-methods-benchmark.R", "tools/outside-methods-score.py",
                      "validation/OUTSIDE_METHODS_SENSITIVITY_PLAN.md", "validation/real-audit/rwanda-pool.csv", "validation/outside-methods-study/replicates.csv.gz")))
  data.frame(file = substring(paths, nchar(root) + 2L),
    sha256 = vapply(paths, function(p) digest::digest(file = p, algo = "sha256"), character(1)), row.names = NULL)
}
hashes_before <- hash_sources()
pkgload::load_all(root, quiet = TRUE)
suppressPackageStartupMessages({library(survey); library(mice)})
options(survey.lonely.psu = "adjust")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
config <- jsonlite::read_json(file.path(root, "validation/clustered-fix-study/settings.json"), simplifyVector = TRUE)
earlier <- jsonlite::read_json(file.path(root, "validation/variance-mi-study/settings.json"), simplifyVector = TRUE)
study <- jsonlite::read_json(file.path(root, "validation/outside-methods-study/settings.json"), simplifyVector = TRUE)
settings <- list(m = study$m); SEEDS <- study$seeds
named_from <- function(path) {
  e <- parse(path)
  list(e = e, names = vapply(e, function(x) if (is.call(x) && identical(x[[1]], as.name("<-")) && is.name(x[[2]])) as.character(x[[2]]) else "", character(1)))
}
fix <- named_from(file.path(root, "tools/clustered-fix-benchmark.R"))
for (e in fix$e[which(fix$names == "rwanda")[1]:which(fix$names == "real_case")[1]]) eval(e, envir = globalenv())
variance <- named_from(file.path(root, "tools/variance-mi-benchmark.R"))
eval(variance$e[[which(variance$names == "simulate_covariate_cell")[1]]], envir = globalenv())
main <- named_from(file.path(root, "tools/outside-methods-benchmark.R"))
for (n in c("design_for", "weighted_fit", "term_of", "interval", "cohort")) eval(main$e[[which(main$names == n)[1]]], envir = globalenv())
grid <- utils::read.csv(file.path(root, "validation/clustered-fix-study/simulation-grid.csv"), stringsAsFactors = FALSE)
covariate_grid <- expand.grid(design = "rwanda", errors = c("shift", "under", "row_only"), fraction = config$fractions, stringsAsFactors = FALSE)
covariate_grid$rows <- 3L
cells <- list()
for (z in real) for (j in seq_along(config$fractions)) cells[[length(cells) + 1L]] <- list(kind = "real", z = z, fraction_index = j, level = z$level, analysis = z$id, fraction = config$fractions[j],
  reps = if (j == 1L) study$reps_real_25 else study$reps_other, number = 100L + length(cells) + 1L)
for (i in seq_len(nrow(grid))) cells[[length(cells) + 1L]] <- list(kind = "simulation", spec = as.list(grid[i, ]), index = i, level = "simulation", analysis = grid$label[i], fraction = grid$fraction[i], reps = study$reps_other, number = i)
for (i in seq_len(nrow(covariate_grid))) cells[[length(cells) + 1L]] <- list(kind = "covariate", spec = as.list(covariate_grid[i, ]), index = i, level = "imputed_covariate", analysis = paste0("covariate_", covariate_grid$errors[i]), fraction = covariate_grid$fraction[i], reps = study$reps_other, number = 200L + i)

fit_raking_approx <- function(p) {
  d <- p$input; d$audited <- p$audited
  if (!"stratum" %in% names(d)) d$stratum <- factor("all")
  all_w <- design_for(p, p$input$lab_value); k <- match(term_of(p), colnames(all_w$x))
  d$if_k <- weighted_fit(all_w, k)$row_score
  d$value_r <- d$gold_value
  design <- survey::twophase(id = list(~cluster_id, ~cluster_id), strata = list(NULL, if (nlevels(d$stratum) > 1L) ~stratum else NULL), subset = ~audited, data = d, method = "approx")
  calibrated <- survey::calibrate(design, phase = 2, formula = ~if_k, calfun = "raking")
  formula <- stats::update(p$formula, if (identical(all.vars(p$formula)[1], "lab_value")) value_r ~ . else . ~ . - lab_value + value_r)
  f <- survey::svyglm(formula, design = calibrated)
  term <- if (identical(term_of(p), "lab_value")) "value_r" else term_of(p)
  interval(unname(stats::coef(f)[term]), unname(survey::SE(f)[term]), survey::degf(calibrated), p$target)
}
fit_mice_all <- function(p, method, seed) {
  use <- unique(c(p$predictors, "lab_value", "gold_value", "cluster_id"))
  d <- p$input[, use, drop = FALSE]
  matrix_p <- matrix(0L, length(use), length(use), dimnames = list(use, use))
  matrix_p["gold_value", setdiff(use, c("gold_value", "cluster_id"))] <- 1L; matrix_p["gold_value", "cluster_id"] <- -2L
  chosen <- stats::setNames(rep("", length(use)), use); chosen["gold_value"] <- method
  lmer_failed <- FALSE; singular <- 0L
  fit <- withCallingHandlers(mice::mice(d, m = settings$m, maxit = 1, method = chosen, predictorMatrix = matrix_p, printFlag = FALSE, seed = seed, eps = 0),
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
run_one <- function(s, r, method, kept) {
  p <- cohort(s, r)
  started <- Sys.time()
  got <- tryCatch(switch(method,
    raking = fit_raking_approx(p),
    mice_2l_lmer = fit_mice_all(p, "2l.lmer", as.integer(SEEDS$lmer + s$number * 1e5 + r)),
    mice_2l_pan = fit_mice_all(p, "2l.pan", as.integer(SEEDS$pan + s$number * 1e5 + r))), error = identity)
  seconds <- as.numeric(difftime(Sys.time(), started, units = "secs"))
  failed <- inherits(got, "error")
  data.frame(level = s$level, analysis = s$analysis, fraction = s$fraction, replicate = r, method = method, status = if (!failed) "ok" else "failed", message = if (failed) substr(conditionMessage(got), 1L, 120L) else "",
             warned = 0L, estimate = if (failed) NA_real_ else got$estimate, se = if (failed) NA_real_ else got$se, df = if (failed) NA_real_ else got$df, width = if (failed) NA_real_ else got$width,
             cover = if (failed) NA_real_ else got$cover, target = p$target, logged = if (failed || is.null(got$logged)) 0L else got$logged, singular = if (failed || is.null(got$singular)) 0L else got$singular, fitted_estimate = NA_real_,
             seconds = round(seconds, 4), stringsAsFactors = FALSE)
}
saved <- utils::read.csv(gzfile(file.path(root, "validation/outside-methods-study/replicates.csv.gz")), stringsAsFactors = FALSE)
key <- function(a, f, r, m) paste(a, f, r, m)
saved_keys <- key(saved$analysis, saved$fraction, saved$replicate, saved$method)
start <- Sys.time(); records <- list(); control_records <- list()
for (s in cells) {
  here <- saved[saved$analysis == s$analysis & saved$fraction == s$fraction, ]
  jobs <- data.frame(r = seq_len(s$reps), method = "raking", stringsAsFactors = FALSE)
  control <- data.frame(r = integer(), method = character(), stringsAsFactors = FALSE)
  for (m in c("mice_2l_lmer", "mice_2l_pan")) {
    h <- here[here$method == m, ]
    again <- h$logged > 0 | h$status != "ok"
    if (s$analysis != "control_one_row") jobs <- rbind(jobs, data.frame(r = h$replicate[again], method = m, stringsAsFactors = FALSE))
    usable <- h$replicate[!again]
    if (length(usable)) control <- rbind(control, data.frame(r = utils::head(sort(usable), 3L), method = m, stringsAsFactors = FALSE))
  }
  all_jobs <- rbind(jobs, control)
  batch <- parallel::mclapply(seq_len(nrow(all_jobs)), function(i) run_one(s, all_jobs$r[i], all_jobs$method[i]), mc.cores = 14L, mc.preschedule = TRUE, mc.set.seed = FALSE)
  if (any(vapply(batch, inherits, logical(1), what = "try-error"))) stop("A worker failed on ", s$analysis)
  out <- do.call(rbind, batch)
  n <- nrow(jobs)
  records[[length(records) + 1L]] <- out[seq_len(n), , drop = FALSE]
  if (nrow(control)) control_records[[length(control_records) + 1L]] <- out[(n + 1L):nrow(out), , drop = FALSE]
  cat(s$level, s$analysis, s$fraction, nrow(jobs), "requests; elapsed", round(as.numeric(difftime(Sys.time(), start, units = "secs"))), "seconds\n"); flush.console()
}
results <- do.call(rbind, records)
file <- gzfile(file.path(outdir, "replicates.csv.gz"), "wt")
utils::write.csv(results, file, row.names = FALSE)
close(file)
checks <- do.call(rbind, control_records)
utils::write.csv(checks, file.path(outdir, "control.csv"), row.names = FALSE)
m <- saved[match(key(checks$analysis, checks$fraction, checks$replicate, checks$method), saved_keys), ]
gaps <- c(estimate = max(abs(checks$estimate - m$estimate)), se = max(abs(checks$se - m$se)), df = max(abs(checks$df - m$df)))
raking_new <- results[results$method == "raking" & results$status == "ok", ]
raking_old <- saved[match(key(raking_new$analysis, raking_new$fraction, raking_new$replicate, "raking"), saved_keys), ]
ok_old <- raking_old$status == "ok"
estimate_gap <- max(abs(raking_new$estimate - raking_old$estimate)[ok_old])
writeLines(c(paste0("Control: ", nrow(checks), " requests with no logged event in the comparison, run again with eps = 0, reproduce it. Largest gap: estimate ", format(gaps["estimate"], digits = 3), ", standard error ", format(gaps["se"], digits = 3), ", degrees of freedom ", format(gaps["df"], digits = 3), "."),
  paste0("Raking under approx against full, ", sum(ok_old), " requests: largest gap in the point estimate ", format(estimate_gap, digits = 3), ".")), file.path(outdir, "control.txt"))
cat(readLines(file.path(outdir, "control.txt")), sep = "\n"); cat("\n")
if (any(gaps > 1e-8) || estimate_gap > 1e-8) stop("The control failed; results are kept but are not verified.")
if (!identical(hashes_before, hash_sources())) stop("Source, plan, pool or reused results changed during the run; evidence is not verified.")
utils::write.csv(hashes_before, file.path(outdir, "source-hashes.csv"), row.names = FALSE)
writeLines(c(paste("Completed UTC:", format(Sys.time(), tz = "UTC", usetz = TRUE)), paste("Elapsed seconds:", round(as.numeric(difftime(Sys.time(), start, units = "secs")), 1)), capture.output(sessionInfo())), file.path(outdir, "session-info.txt"))
cat("Completed:", nrow(results), "requests;", paste(names(table(results$status)), table(results$status), collapse = ", "), "\n")
