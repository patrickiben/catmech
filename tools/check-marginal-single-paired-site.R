#!/usr/bin/env Rscript
# Compare preparation only; never fit a model or inspect coefficient accuracy.
# Usage: Rscript tools/check-marginal-single-paired-site.R STAGE OLD_SAMPLER OLD_RUNNER OUTPUT_JSON
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) stop("Expected STAGE OLD_SAMPLER OLD_RUNNER OUTPUT_JSON.")
script <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1L])
root <- normalizePath(file.path(dirname(script), ".."), mustWork = TRUE)
stage <- normalizePath(args[1L], mustWork = TRUE)
old_sampler <- normalizePath(args[2L], mustWork = TRUE)
old_runner <- normalizePath(args[3L], mustWork = TRUE)
output <- args[4L]
if (file.exists(output)) stop("Do not overwrite an existing check record.")
settings <- jsonlite::read_json(file.path(stage, "settings.json"), simplifyVector = TRUE)
source_hashes <- utils::read.csv(file.path(stage, "source-hashes-before.csv"), stringsAsFactors = FALSE)
expected <- function(file) {
  z <- source_hashes$sha256[source_hashes$file == file]
  stopifnot(length(z) == 1L)
  z
}
stopifnot(identical(digest::digest(file = old_sampler, algo = "sha256"), expected("R/marginal-sampler.R")),
  identical(digest::digest(file = old_runner, algo = "sha256"), expected("validation/run-joint-marginal.R")))
pkgload::load_all(root, quiet = TRUE)
old <- new.env(parent = asNamespace("catmech"))
sys.source(old_sampler, envir = old)
generator <- new.env(parent = globalenv())
generator$config <- settings
generator$contextual_cells <- settings$contextual_cells
required <- c("rmv", "logit_intercept", "make_cohort")
found <- character()
for (expression in parse(old_runner)) {
  if (is.call(expression) && identical(expression[[1L]], as.name("<-")) &&
      is.symbol(expression[[2L]]) && as.character(expression[[2L]]) %in% required) {
    eval(expression, generator)
    found <- c(found, as.character(expression[[2L]]))
  }
}
stopifnot(setequal(found, required))
records <- utils::read.csv(file.path(stage, "records.csv"), stringsAsFactors = FALSE)
cohorts <- unique(records[c("cell", "replicate", "input_sha256")])
checks <- list()
for (index in seq_along(settings$cells)) for (replicate in seq_len(settings$primary_reps)) {
  p <- generator$make_cohort(settings$cells[index], index, replicate)
  known <- cohorts[cohorts$cell == settings$cells[index] & cohorts$replicate == replicate, ]
  stopifnot(nrow(known) == 1L, identical(known$input_sha256, p$meta$input_sha256))
  d <- p$input
  before <- old$catmech_marginal_cluster_prepare(d, "W", d$gold, "x", "site")
  after <- catmech:::catmech_marginal_cluster_prepare(d, "W", d$gold, "x", "site")
  paired <- !is.na(d$W) & !is.na(d$gold)
  paired_sites <- length(unique(d$site[paired]))
  stopifnot(paired_sites > 1L, identical(before, after),
    identical(serialize(before, NULL, version = 2), serialize(after, NULL, version = 2)))
  checks[[length(checks) + 1L]] <- data.frame(cell = settings$cells[index], replicate = replicate,
    data_seed = p$meta$data_seed, paired_sites = paired_sites,
    input_sha256 = p$meta$input_sha256,
    preparation_sha256 = digest::digest(serialize(after, NULL, version = 2), algo = "sha256", serialize = FALSE))
}
checks <- do.call(rbind, checks)
report <- list(status = "pass", stage = settings$mode,
  cohorts = nrow(checks), minimum_paired_sites = min(checks$paired_sites),
  all_cohort_input_hashes_match = TRUE, all_preparation_objects_byte_identical = TRUE,
  single_paired_site_branch_entered = FALSE, fitting_or_accuracy_evaluation_performed = FALSE,
  frozen_source_fingerprint = settings$source_fingerprint,
  old_sampler_sha256 = digest::digest(file = old_sampler, algo = "sha256"),
  current_sampler_sha256 = digest::digest(file = file.path(root, "R/marginal-sampler.R"), algo = "sha256"),
  frozen_runner_sha256 = digest::digest(file = old_runner, algo = "sha256"),
  check_script_sha256 = digest::digest(file = script, algo = "sha256"), cohort_checks = checks)
jsonlite::write_json(report, output, auto_unbox = TRUE, pretty = TRUE, digits = NA)
cat(nrow(checks), "cohorts regenerated; minimum", min(checks$paired_sites),
  "paired sites; inputs match; full preparations byte-identical; no fits or accuracy evaluated.\n")
