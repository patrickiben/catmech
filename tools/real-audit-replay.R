#!/usr/bin/env Rscript
# Replays a sample of the real-audit study from its seeds, using the runner's own functions read from its source
# without running its main loop, and saves the results so tools/verify-artifacts.py can compare them with the saved
# study. This catches relabeled or rescaled results, which recomputing summaries from saved rows cannot.
# Usage: Rscript tools/real-audit-replay.R [NEW_DIRECTORY] [REPLICATES_PER_CELL]
args <- commandArgs(trailingOnly = TRUE)
file_arg <- grep("^--file=", commandArgs(), value = TRUE)
root <- normalizePath(file.path(dirname(sub("^--file=", "", file_arg[1])), ".."))
outdir <- if (length(args)) args[1] else file.path(root, "validation", "real-audit-replay")
per_cell <- if (length(args) >= 2L) as.integer(args[2]) else 5L
if (dir.exists(outdir) && length(list.files(outdir, all.files = TRUE, no.. = TRUE)))
  stop("Use a new or empty output directory; existing evidence is never overwritten.")
pkgload::load_all(root, quiet = TRUE)
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
study <- file.path(root, "validation", "real-audit-study")
config <- jsonlite::read_json(file.path(study, "settings.json"), simplifyVector = TRUE)
expressions <- parse(file.path(root, "tools/real-audit-benchmark.R"))
assigned <- vapply(expressions, function(e) if (is.call(e) && identical(e[[1]], as.name("<-")) && is.name(e[[2]])) as.character(e[[2]]) else "", character(1))
first <- which(assigned == "rwanda")[1]; last <- which(assigned == "run_replicate")[1]
for (e in expressions[first:last]) eval(e, envir = globalenv())
cells <- list()
for (i in seq_along(analyses)) for (j in seq_along(config$fractions)) {
  z <- analyses[[i]]
  cells[[length(cells) + 1L]] <- list(index = i, fraction_index = j, level = z$level, analysis = z$id, fraction = config$fractions[j],
    scenario = paste0(z$level, ":", z$id, ":audit_", 100 * config$fractions[j], "pct"))
}
jobs <- expand.grid(cell = seq_along(cells), replicate = seq_len(per_cell))
replayed <- parallel::mclapply(seq_len(nrow(jobs)), function(k) run_replicate(cells[[jobs$cell[k]]], jobs$replicate[k])$records,
                               mc.cores = 14L, mc.preschedule = FALSE, mc.set.seed = FALSE)
if (any(vapply(replayed, inherits, logical(1), what = "try-error"))) stop("A replay worker failed.")
replayed <- do.call(rbind, replayed)
saved <- utils::read.csv(file.path(study, "replicates.csv"), stringsAsFactors = FALSE)
key <- function(d) paste(d$scenario, d$replicate, d$method)
match_saved <- saved[match(key(replayed), key(saved)), ]
same <- replayed$status == match_saved$status & (is.na(replayed$estimate) | (abs(replayed$estimate - match_saved$estimate) < 1e-8 &
  abs(replayed$se - match_saved$se) < 1e-8 & abs(replayed$se_cluster - match_saved$se_cluster) < 1e-8))
utils::write.csv(replayed[c("scenario", "replicate", "method", "status", "estimate", "se", "se_cluster")], file.path(outdir, "replay.csv"), row.names = FALSE)
writeLines(paste0("Replayed ", nrow(replayed), " requests from ", per_cell, " cohorts per cell with the runner's own functions; ", sum(same), " match the saved study in status, estimate and both standard errors to within 1e-8."), file.path(outdir, "README.md"))
cat(readLines(file.path(outdir, "README.md")), sep = "\n")
if (!all(same)) stop("Replay does not match the saved study.")
