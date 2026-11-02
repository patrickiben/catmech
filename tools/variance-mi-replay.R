#!/usr/bin/env Rscript
# Replays the study's intervals under the adopted variance through catmech_pool_lm(), after that variance moved into the
# package, using the runner's own functions read from its source without running its main loop.
# Usage: Rscript tools/variance-mi-replay.R [NEW_DIRECTORY]
args <- commandArgs(trailingOnly = TRUE)
file_arg <- grep("^--file=", commandArgs(), value = TRUE)
root <- normalizePath(file.path(dirname(sub("^--file=", "", file_arg[1])), ".."))
outdir <- if (length(args)) args[1] else file.path(root, "validation", "variance-mi-replay")
if (dir.exists(outdir) && length(list.files(outdir, all.files = TRUE, no.. = TRUE)))
  stop("Use a new or empty output directory; existing evidence is never overwritten.")
pkgload::load_all(root, quiet = TRUE)
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
study <- file.path(root, "validation", "variance-mi-study")
settings <- jsonlite::read_json(file.path(study, "settings.json"), simplifyVector = TRUE)
config <- jsonlite::read_json(file.path(root, "validation/clustered-fix-study/settings.json"), simplifyVector = TRUE)
expressions <- parse(file.path(root, "tools/clustered-fix-benchmark.R"))
assigned <- vapply(expressions, function(e) if (is.call(e) && identical(e[[1]], as.name("<-")) && is.name(e[[2]])) as.character(e[[2]]) else "", character(1))
first <- which(assigned == "rwanda")[1]; last <- which(assigned == "real_case")[1]
for (e in expressions[first:last]) eval(e, envir = globalenv())
runner <- parse(file.path(root, "tools/variance-mi-benchmark.R"))
named <- vapply(runner, function(e) if (is.call(e) && identical(e[[1]], as.name("<-")) && is.name(e[[2]])) as.character(e[[2]]) else "", character(1))
SEEDS <- settings$seeds
eval(runner[[which(named == "simulate_covariate_cell")[1]]], envir = globalenv())
grid <- utils::read.csv(file.path(root, "validation/clustered-fix-study/simulation-grid.csv"), stringsAsFactors = FALSE)
covariate_grid <- expand.grid(design = "rwanda", errors = c("shift", "under", "row_only"), fraction = config$fractions, stringsAsFactors = FALSE)
covariate_grid$rows <- 3L
cells <- list()
for (z in real) for (j in seq_along(config$fractions)) cells[[length(cells) + 1L]] <- list(kind = "real", z = z, fraction_index = j, analysis = z$id, fraction = config$fractions[j], reps = if (j == 1L) settings$reps_real_25 else settings$reps_other)
for (i in seq_len(nrow(grid))) cells[[length(cells) + 1L]] <- list(kind = "simulation", spec = as.list(grid[i, ]), index = i, analysis = grid$label[i], fraction = grid$fraction[i], reps = settings$reps_other)
for (i in seq_len(nrow(covariate_grid))) cells[[length(cells) + 1L]] <- list(kind = "covariate", spec = as.list(covariate_grid[i, ]), index = i, analysis = paste0("covariate_", covariate_grid$errors[i]), fraction = covariate_grid$fraction[i], reps = settings$reps_other)
one <- function(s, r) {
  if (s$kind == "real") {
    p <- real_case(s$z, s$fraction_index, r)
    seed <- as.integer(s$z$seeds$imputation + s$z$index * 1e6 + s$fraction_index * 1e5 + r)
  } else if (s$kind == "simulation") {
    p <- simulate_cell(s$spec, as.integer(config$simulation$cohort + s$index * 1e5 + r), s$fraction)
    seed <- as.integer(config$simulation$imputation + s$index * 1e5 + r)
  } else {
    p <- simulate_covariate_cell(s$spec, as.integer(SEEDS$covariate_cohort + s$index * 1e5 + r), s$fraction)
    seed <- as.integer(SEEDS$covariate_imputation + s$index * 1e5 + r)
  }
  term <- if (is.null(p$term)) "covariate" else p$term
  do.call(rbind, lapply(c("oracle", "linear_mi", "random_intercept"), function(imputer) {
    input <- p$input
    if (imputer == "oracle") input$gold_value <- p$truth
    options <- list(input, "lab_value", "gold_value", p$predictors, method = "normal", m = settings$m, seed = seed, use_recorded = TRUE)
    if (imputer == "random_intercept") options$cluster <- "cluster_id"
    got <- tryCatch(suppressWarnings(catmech_pool_lm(do.call(catmech_reference, options), p$formula, cluster = "cluster_id")), error = identity)
    failed <- inherits(got, "error")
    row <- if (failed) NULL else got[got$term == term, ]
    data.frame(analysis = s$analysis, fraction = s$fraction, replicate = r, imputer = imputer, status = if (failed) "failed" else "ok",
               estimate = if (failed) NA_real_ else row$estimate, se = if (failed) NA_real_ else row$std.error, df = if (failed) NA_real_ else row$df,
               lower = if (failed) NA_real_ else row$lower, upper = if (failed) NA_real_ else row$upper)
  }))
}
started <- Sys.time(); records <- list()
for (s in cells) {
  batch <- parallel::mclapply(seq_len(s$reps), function(r) one(s, r), mc.cores = 14L, mc.preschedule = TRUE, mc.set.seed = FALSE)
  if (any(vapply(batch, inherits, logical(1), what = "try-error"))) stop("A replay worker failed.")
  records <- c(records, batch)
}
replayed <- do.call(rbind, records)
out <- gzfile(file.path(outdir, "replay.csv.gz"), "wt")
utils::write.csv(replayed, out, row.names = FALSE)
close(out)
saved <- utils::read.csv(gzfile(file.path(study, "replicates.csv.gz")), stringsAsFactors = FALSE)
key <- function(d, method) paste(d$analysis, d$fraction, d$replicate, d[[method]])
m <- saved[match(key(replayed, "imputer"), key(saved, "imputer")), ]
ok <- replayed$status == "ok"
stopifnot(all((m$status == "ok") == ok))
floored <- ok & m$imputer != "oracle" & 1.05 * m$between / m$se_cv3^2 < 1e-4
estimate_gap <- max(abs(replayed$estimate - m$estimate)[ok]); se_gap <- max(abs(replayed$se - m$se_cv3)[ok])
df_gap <- max(abs(replayed$df - m$df_cv3)[ok & !floored], 0); floored_gap <- max(abs(replayed$df - m$df_cv3)[floored], 0)
half <- (replayed$upper - replayed$lower) / 2
flips <- sum(as.numeric(abs(replayed$estimate - m$target) <= half)[ok] != m$cover_cv3[ok])
writeLines(c(paste0("Replayed all ", nrow(replayed), " study requests under `cv3` through the package's own `catmech_reference()` and `catmech_pool_lm()`, using the runner's cohort and audit functions and the study's seeds."),
  paste0(sum(ok), " returned intervals, and the rest match the study's refusals. Largest difference: estimate ", format(estimate_gap, digits = 3), ", standard error ", format(se_gap, digits = 3), ", degrees of freedom ", format(df_gap, digits = 3), ". In ", sum(floored), " requests with almost no between-imputation variance, degrees of freedom differ by up to ", format(floored_gap, digits = 3), ", from the runner fault the plan describes. Coverage flips: ", flips, ".")),
  file.path(outdir, "README.md"))
cat(readLines(file.path(outdir, "README.md")), sep = "\n")
cat("\nSeconds:", round(as.numeric(difftime(Sys.time(), started, units = "secs"))), "\n")
if (estimate_gap > 1e-8 || se_gap > 1e-8 || df_gap > 1e-6 || floored_gap > .01 || flips) stop("Replay does not match the saved study.")
