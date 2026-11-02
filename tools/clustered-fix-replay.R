#!/usr/bin/env Rscript
# Replays every study request for the selected candidate through the package itself, after the candidate moved into
# R/reference.R, using the runner's own functions read from its source without running its main loop. The saved results
# let tools/verify-artifacts.py confirm that the package gives what the study scored.
# Usage: Rscript tools/clustered-fix-replay.R [NEW_DIRECTORY]
args <- commandArgs(trailingOnly = TRUE)
file_arg <- grep("^--file=", commandArgs(), value = TRUE)
root <- normalizePath(file.path(dirname(sub("^--file=", "", file_arg[1])), ".."))
outdir <- if (length(args)) args[1] else file.path(root, "validation", "clustered-fix-replay")
if (dir.exists(outdir) && length(list.files(outdir, all.files = TRUE, no.. = TRUE)))
  stop("Use a new or empty output directory; existing evidence is never overwritten.")
pkgload::load_all(root, quiet = TRUE)
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
study <- file.path(root, "validation", "clustered-fix-study")
config <- jsonlite::read_json(file.path(study, "settings.json"), simplifyVector = TRUE)
decisions <- readLines(file.path(root, "validation", "clustered-fix-decisions.txt"))
selected <- sub("^.*SELECTED: (\\S+).*$", "\\1", grep("SELECTED: ", decisions, value = TRUE)[1])
expressions <- parse(file.path(root, "tools/clustered-fix-benchmark.R"))
assigned <- vapply(expressions, function(e) if (is.call(e) && identical(e[[1]], as.name("<-")) && is.name(e[[2]])) as.character(e[[2]]) else "", character(1))
first <- which(assigned == "rwanda")[1]; last <- which(assigned == "real_case")[1]
for (e in expressions[first:last]) eval(e, envir = globalenv())
grid <- utils::read.csv(file.path(study, "simulation-grid.csv"), stringsAsFactors = FALSE)
scenarios <- list()
for (z in real) for (j in seq_along(config$fractions))
  scenarios[[length(scenarios) + 1L]] <- list(kind = "real", z = z, fraction_index = j, analysis = z$id, fraction = config$fractions[j])
for (i in seq_len(nrow(grid)))
  scenarios[[length(scenarios) + 1L]] <- list(kind = "simulation", cell = as.list(grid[i, ]), index = i, analysis = grid$label[i], fraction = grid$fraction[i])
one <- function(s, r) {
  if (s$kind == "real") {
    p <- real_case(s$z, s$fraction_index, r)
    seed <- as.integer(s$z$seeds$imputation + s$z$index * 1e6 + s$fraction_index * 1e5 + r)
  } else {
    p <- simulate_cell(s$cell, as.integer(config$simulation$cohort + s$index * 1e5 + r), s$fraction)
    seed <- as.integer(config$simulation$imputation + s$index * 1e5 + r)
  }
  options <- list(p$input, "lab_value", "gold_value", p$predictors, method = "normal", m = config$m, seed = seed, use_recorded = TRUE, cluster = "cluster_id")
  if ("cluster_strata" %in% names(formals(catmech_reference))) options$cluster_strata <- p$strata
  result <- tryCatch(suppressWarnings(catmech_pool_lm(do.call(catmech_reference, options), p$formula, cluster = "cluster_id")), error = identity)
  failed <- inherits(result, "error")
  row <- if (failed) NULL else result[result$term == "covariate", ]
  data.frame(analysis = s$analysis, fraction = s$fraction, replicate = r, status = if (failed) "refused" else "ok",
             estimate = if (failed) NA_real_ else row$estimate, se = if (failed) NA_real_ else row$std.error, df = if (failed) NA_real_ else row$df,
             lower = if (failed) NA_real_ else row$lower, upper = if (failed) NA_real_ else row$upper)
}
started <- Sys.time()
records <- list()
for (s in scenarios) {
  batch <- parallel::mclapply(seq_len(config$reps), function(r) one(s, r), mc.cores = 14L, mc.preschedule = TRUE, mc.set.seed = FALSE)
  if (any(vapply(batch, inherits, logical(1), what = "try-error"))) stop("A replay worker failed.")
  records <- c(records, batch)
}
replayed <- do.call(rbind, records)
out <- gzfile(file.path(outdir, "replay.csv.gz"), "wt")
utils::write.csv(replayed, out, row.names = FALSE)
close(out)
saved <- utils::read.csv(gzfile(file.path(study, "replicates.csv.gz")), stringsAsFactors = FALSE)
saved <- saved[saved$method == selected, ]
key <- function(d) paste(d$analysis, d$fraction, d$replicate)
match_saved <- saved[match(key(replayed), key(saved)), ]
near <- function(a, b) (is.na(a) & is.na(b)) | (!is.na(a) & !is.na(b) & abs(a - b) < 1e-10)
same <- replayed$status == match_saved$status & near(replayed$estimate, match_saved$estimate) & near(replayed$se, match_saved$se) & near(replayed$df, match_saved$df) &
  near(replayed$lower, match_saved$lower) & near(replayed$upper, match_saved$upper)
largest <- max(abs(replayed$estimate - match_saved$estimate), abs(replayed$se - match_saved$se), na.rm = TRUE)
writeLines(c(paste0("Replayed all ", nrow(replayed), " study requests for `", selected, "` through the package's own `catmech_reference()` and `catmech_pool_lm()`, using the runner's cohort and audit functions and the study's seeds."),
             paste0(sum(same), " match the saved study in status, estimate, standard error, degrees of freedom and limits to within 1e-10. Largest difference in an estimate or standard error: ", format(largest, digits = 3), ".")),
           file.path(outdir, "README.md"))
cat(readLines(file.path(outdir, "README.md")), sep = "\n")
cat("Seconds:", round(as.numeric(difftime(Sys.time(), started, units = "secs"))), "\n")
if (!all(same)) stop("Replay does not match the saved study.")
