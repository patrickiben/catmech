#!/usr/bin/env Rscript
# Engineering check only: no model fits, estimates, intervals or accuracy scoring.
# Usage: Rscript tools/check-marginal-release-equivalence.R STAGE FROZEN_SOURCE OUTPUT_JSON
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) stop("Expected STAGE FROZEN_SOURCE OUTPUT_JSON.")
script <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1L])
root <- normalizePath(file.path(dirname(script), ".."), mustWork = TRUE)
stage <- normalizePath(args[1L], mustWork = TRUE)
frozen <- normalizePath(args[2L], mustWork = TRUE)
output <- args[3L]
if (file.exists(output)) stop("Do not overwrite an existing equivalence receipt.")
if (!file.exists(file.path(stage, "COMPLETE.json")))
  stop("A completed study is required; do not inspect or alter a running study.")
hash <- function(file) digest::digest(file = file, algo = "sha256")
bytes <- function(file) readBin(file, "raw", n = file.info(file)$size)
settings <- jsonlite::read_json(file.path(stage, "settings.json"), simplifyVector = TRUE)
complete <- jsonlite::read_json(file.path(stage, "COMPLETE.json"), simplifyVector = TRUE)
read_manifest <- function(file) utils::read.csv(file, stringsAsFactors = FALSE,
  colClasses = c(file = "character", bytes = "numeric", sha256 = "character"))
before <- read_manifest(file.path(stage, "source-hashes-before.csv"))
after <- read_manifest(file.path(stage, "source-hashes-after.csv"))
stopifnot(identical(before, after), !anyDuplicated(before$file),
  identical(complete$status, "complete"), identical(complete$mode, settings$mode),
  identical(settings$mode, "confirmation"), length(settings$cells) == 9L,
  settings$primary_reps == 400L, settings$prior_reps == 200L,
  length(settings$prior_cells) == 4L, complete$analysis_requests == 27600L,
  complete$mcmc_fit_requests == 12000L, complete$fitted_workflow_requests == 15600L,
  identical(complete$source_fingerprint, settings$source_fingerprint),
  identical(digest::digest(before, algo = "sha256"), settings$source_fingerprint))
for (i in seq_len(nrow(before))) {
  path <- file.path(frozen, before$file[i])
  stopifnot(file.exists(path), file.info(path)$size == before$bytes[i],
    identical(hash(path), before$sha256[i]))
}

# Reject any additional scientific source file or any unrelated source change.
source_set <- function(base) sort(c(
  "DESCRIPTION", "NAMESPACE",
  paste0("R/", list.files(file.path(base, "R"), "\\.R$")),
  paste0("src/", list.files(file.path(base, "src"), "\\.(c|cpp|h)$|^Makevars"))))
scientific <- before$file[grepl("^(R/|src/)", before$file) |
  before$file %in% c("DESCRIPTION", "NAMESPACE")]
stopifnot(identical(source_set(root), sort(scientific)))
sampler <- "R/marginal-sampler.R"
old_call <- "ave(x[, j], ids, FUN = mean)"
new_call <- "stats::ave(x[, j], ids, FUN = mean)"
old_text <- rawToChar(bytes(file.path(frozen, sampler)))
positions <- gregexpr(old_call, old_text, fixed = TRUE)[[1L]]
stopifnot(length(positions) == 1L, positions[1L] > 0L,
  !grepl(new_call, old_text, fixed = TRUE))
expected <- sub(old_call, new_call, old_text, fixed = TRUE)
stopifnot(identical(bytes(file.path(root, sampler)), charToRaw(expected)))
release_hashes <- before
for (i in seq_len(nrow(before))) {
  path <- file.path(root, before$file[i])
  stopifnot(file.exists(path))
  release_hashes$bytes[i] <- file.info(path)$size
  release_hashes$sha256[i] <- hash(path)
  if (before$file[i] != sampler)
    stopifnot(identical(bytes(path), bytes(file.path(frozen, before$file[i]))))
}
stopifnot(identical(before$file[before$sha256 != release_hashes$sha256], sampler))

pkgload::load_all(root, quiet = TRUE, compile = FALSE)
namespace <- asNamespace("catmech")
old <- new.env(parent = namespace)
sys.source(file.path(frozen, sampler), envir = old)
current_prepare <- get("catmech_marginal_cluster_prepare", envir = namespace)
stopifnot(identical(environment(current_prepare), namespace),
  identical(get("ave", envir = old, inherits = TRUE), stats::ave))

# Evaluate only the three frozen generator definitions, never the study runner.
generator <- new.env(parent = globalenv())
generator$config <- settings
generator$contextual_cells <- settings$contextual_cells
required <- c("rmv", "logit_intercept", "make_cohort")
found <- character()
for (expression in parse(file.path(frozen, "validation/run-joint-marginal.R"))) {
  if (is.call(expression) && identical(expression[[1L]], as.name("<-")) &&
      is.symbol(expression[[2L]]) && as.character(expression[[2L]]) %in% required) {
    eval(expression, generator)
    found <- c(found, as.character(expression[[2L]]))
  }
}
stopifnot(identical(sort(found), sort(required)))
record_file <- file.path(stage, "records.csv")
columns <- names(utils::read.csv(record_file, nrows = 0L, check.names = FALSE))
classes <- stats::setNames(rep("NULL", length(columns)), columns)
classes[c("cell", "replicate", "input_sha256")] <- "character"
cohorts <- unique(utils::read.csv(record_file, colClasses = classes,
  stringsAsFactors = FALSE))
expected_count <- length(settings$cells) * settings$primary_reps
stopifnot(nrow(cohorts) == expected_count, complete$cohorts == expected_count,
  !anyDuplicated(paste(cohorts$cell, cohorts$replicate, sep = ":")))
checks <- vector("list", expected_count)
counter <- 0L
for (index in seq_along(settings$cells)) for (replicate in seq_len(settings$primary_reps)) {
  p <- generator$make_cohort(settings$cells[index], index, replicate)
  known <- cohorts[cohorts$cell == settings$cells[index] & cohorts$replicate == replicate, ]
  stopifnot(nrow(known) == 1L, identical(known$input_sha256, p$meta$input_sha256))
  d <- p$input
  old_problem <- old$catmech_marginal_cluster_prepare(d, "W", d$gold, "x", "site")
  new_problem <- current_prepare(d, "W", d$gold, "x", "site")
  old_serial <- serialize(old_problem, NULL, version = 2)
  new_serial <- serialize(new_problem, NULL, version = 2)
  stopifnot(identical(old_problem, new_problem), identical(old_serial, new_serial))
  if (counter == 0L) { example <- d; example_problem <- new_problem }
  counter <- counter + 1L
  checks[[counter]] <- data.frame(cell = settings$cells[index], replicate = replicate,
    data_seed = p$meta$data_seed, input_sha256 = p$meta$input_sha256,
    old_preparation_sha256 = digest::digest(old_serial, algo = "sha256", serialize = FALSE),
    new_preparation_sha256 = digest::digest(new_serial, algo = "sha256", serialize = FALSE))
}

# Namespace-only behavior is the actual portability correction, not a new model.
stopifnot("package:stats" %in% search())
detach("package:stats", unload = FALSE, character.only = TRUE)
stopifnot(!"package:stats" %in% search(),
  !exists("ave", envir = namespace, inherits = TRUE))
old_failure <- tryCatch(old$catmech_marginal_cluster_prepare(example, "W",
  example$gold, "x", "site"), error = identity)
stopifnot(inherits(old_failure, "error"), grepl("ave", conditionMessage(old_failure), fixed = TRUE),
  identical(current_prepare(example, "W", example$gold, "x", "site"), example_problem))
for (i in seq_len(nrow(release_hashes)))
  stopifnot(identical(hash(file.path(root, release_hashes$file[i])), release_hashes$sha256[i]))

changed <- which(before$file == sampler)
report <- list(schema_version = 1L, status = "pass",
  correction_id = "qualify_stats_ave", correction_kind = "post_confirmation_namespace_portability",
  stage = settings$mode, study_directory = stage, frozen_source_directory = frozen,
  frozen_source_fingerprint = settings$source_fingerprint,
  study_complete_sha256 = hash(file.path(stage, "COMPLETE.json")),
  study_manifest_sha256 = hash(file.path(stage, "source-hashes-before.csv")),
  study_records_sha256 = hash(record_file), check_script_sha256 = hash(script),
  permitted_source_delta = list(file = sampler, replacement_count = 1L,
    old_text = old_call, new_text = new_call,
    old_bytes = before$bytes[changed], new_bytes = release_hashes$bytes[changed],
    old_sha256 = before$sha256[changed], new_sha256 = release_hashes$sha256[changed]),
  all_frozen_sources_verified = TRUE, all_other_frozen_sources_byte_identical = TRUE,
  release_source_file_set_unchanged = TRUE, cohorts = counter,
  all_cohort_input_hashes_match = TRUE, all_preparation_objects_byte_identical = TRUE,
  legacy_ave_resolved_to_stats_ave = TRUE,
  stats_detached_old_prepare_failed = TRUE, stats_detached_current_prepare_passed = TRUE,
  stats_detached_current_preparation_unchanged = TRUE,
  fitting_or_accuracy_evaluation_performed = FALSE,
  interpretation = paste("Namespace portability correction only; release source differs from confirmation.",
    "The original confirmation evidence is retained, and all its public inputs and model preparations",
    "are byte-identical under the correction. This check does not generate new coverage evidence."),
  frozen_source_checks = before, release_source_checks = release_hashes,
  cohort_checks = do.call(rbind, checks))
jsonlite::write_json(report, output, auto_unbox = TRUE, pretty = TRUE, digits = NA)
cat(counter, "cohorts regenerated; inputs match; preparations byte-identical;",
  "one exact stats::ave qualification; stats-detached namespace check passed; no fits.\n")
