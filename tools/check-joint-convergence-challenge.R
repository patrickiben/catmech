#!/usr/bin/env Rscript
# Selected numerical stress check only. No interval or point-estimate analysis.
args <- commandArgs(trailingOnly = TRUE)
script <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1L])
root <- normalizePath(file.path(dirname(script), ".."), mustWork = TRUE)
out <- if (length(args)) args[1L] else file.path(root, "validation/v0.9.0-development/challenge")
if (dir.exists(out) && length(list.files(out, all.files = TRUE, no.. = TRUE)))
  stop("Use a new or empty output directory; challenge evidence is not overwritten.")
relative_sources <- c("R/joint-cluster.R", "R/reference.R", "R/model.R", "R/spec.R", "R/impute.R",
  "tools/check-joint-convergence-challenge.R",
  "validation/v0.8.0-source/validation/run-joint-cluster.R",
  "validation/joint-cluster-study-final/records.csv",
  "validation/joint-cluster-study-final/settings.json")
hash_sources <- function() data.frame(file = relative_sources,
  sha256 = vapply(file.path(root, relative_sources), function(f)
    digest::digest(file = f, algo = "sha256"), character(1)), row.names = NULL)
before <- hash_sources()
fingerprint <- digest::digest(before, algo = "sha256")
pkgload::load_all(root, quiet = TRUE)
baseline <- utils::read.csv(file.path(root, "validation/joint-cluster-study-final/records.csv"),
  stringsAsFactors = FALSE)
selected <- baseline[baseline$error_class == "catmech_convergence_error", , drop = FALSE]
if (nrow(selected) != 22L || any(!startsWith(selected$method, "joint")) ||
    anyDuplicated(selected[c("cell", "replicate", "method")]))
  stop("The selected challenge must contain exactly the 22 recorded v0.8 convergence refusals.")
old_settings <- jsonlite::read_json(file.path(root, "validation/joint-cluster-study-final/settings.json"),
  simplifyVector = TRUE)
generator <- new.env(parent = globalenv())
generator$settings <- old_settings
# Read only the three generator functions. Evaluating the entire historical
# runner would launch another study and is deliberately avoided.
required <- c("rmv", "logistic_intercept", "make_cohort")
loaded <- character()
for (expression in parse(file.path(root, "validation/v0.8.0-source/validation/run-joint-cluster.R"))) {
  if (is.call(expression) && identical(expression[[1L]], as.name("<-")) &&
      is.symbol(expression[[2L]]) && as.character(expression[[2L]]) %in% required) {
    eval(expression, generator)
    loaded <- c(loaded, as.character(expression[[2L]]))
  }
}
stopifnot(setequal(loaded, required))
control <- catmech:::catmech_joint_cluster_control()
stopifnot(identical(control, list(chains = 4L, burnin = 1000L, iterations = 4000L,
  max_iterations = 32000L, prior_scale = 1)))
dir.create(out, recursive = TRUE, showWarnings = FALSE)
utils::write.csv(before, file.path(out, "source-hashes-before.csv"), row.names = FALSE)
jsonlite::write_json(list(kind = "selected_previous_convergence_failures", requests = nrow(selected),
  controls = control, m = old_settings$m, pan_version = as.character(utils::packageVersion("pan")),
  catmech_version = as.character(utils::packageVersion("catmech")),
  source_fingerprint = fingerprint, accuracy_computed = FALSE,
  limitation = "Selected numerical challenge; it cannot estimate a population refusal rate or interval calibration."),
  file.path(out, "settings.json"), auto_unbox = TRUE, pretty = TRUE)
capture.output(utils::sessionInfo(), file = file.path(out, "session-info.txt"))
records <- summaries <- list()
started <- Sys.time()
for (i in seq_len(nrow(selected))) {
  if (!identical(hash_sources(), before)) stop("Challenge inputs or numerical sources changed before a request.")
  old <- selected[i, , drop = FALSE]
  p <- generator$make_cohort(old$cell, old$cell_index, old$replicate)
  if (!identical(p$meta$input_sha256, old$input_sha256) ||
      p$meta$data_seed != old$data_seed || p$meta$imputation_seed != old$imputation_seed)
    stop("Historical input or seeds were not reproduced exactly.")
  prior_scale <- switch(old$method, joint_prior_quarter = .25, joint_prior_four = 4, joint = 1,
    stop("Unknown historical joint method."))
  warnings <- character()
  start <- proc.time()[["elapsed"]]
  result <- tryCatch(withCallingHandlers(catmech_reference(p$input, "W", "gold", "x",
    m = as.integer(old_settings$m), seed = old$imputation_seed, use_recorded = TRUE,
    cluster = "site", cluster_model = "joint", joint_control = list(prior_scale = prior_scale)),
    warning = function(w) {warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning")}),
    error = identity)
  elapsed <- proc.time()[["elapsed"]] - start
  failed <- inherits(result, "error")
  diagnostics <- result$model_diagnostics
  s <- diagnostics$summary
  aggregate <- function(z, op) if (length(z) && all(is.finite(z))) op(z) else NA_real_
  row <- data.frame(cell = old$cell, replicate = old$replicate, method = old$method,
    data_seed = old$data_seed, imputation_seed = old$imputation_seed,
    input_sha256 = old$input_sha256, prior_scale = prior_scale,
    original_iterations = old$iterations_used, original_max_rhat = old$max_rhat,
    original_min_ess_bulk = old$min_ess_bulk, original_min_ess_tail = old$min_ess_tail,
    status = if (failed) "error" else "ok",
    error_class = if (failed) class(result)[1L] else "",
    message = if (failed) conditionMessage(result) else "",
    warnings = paste(unique(warnings), collapse = " | "),
    convergence_pass = isTRUE(diagnostics$convergence_pass),
    iterations_used = if (length(diagnostics$iterations_used) == 1L) diagnostics$iterations_used else NA_integer_,
    max_rhat = aggregate(s$rhat, max), min_ess_bulk = aggregate(s$ess_bulk, min),
    min_ess_tail = aggregate(s$ess_tail, min),
    attempts = if (is.null(diagnostics$attempts)) "" else as.character(jsonlite::toJSON(diagnostics$attempts, auto_unbox = TRUE)),
    elapsed_seconds = elapsed, source_fingerprint = fingerprint, stringsAsFactors = FALSE)
  records[[i]] <- row
  if (is.data.frame(s) && nrow(s))
    summaries[[i]] <- cbind(cell = old$cell, replicate = old$replicate, method = old$method, s)
  utils::write.csv(do.call(rbind, records), file.path(out, "records.csv"), row.names = FALSE, na = "")
  if (length(summaries)) utils::write.csv(do.call(rbind, summaries),
    file.path(out, "diagnostics.csv"), row.names = FALSE, na = "")
  cat(sprintf("%d/%d %s replicate %d %s: %s; retained %s; %.3fs\n", i, nrow(selected),
    old$cell, old$replicate, old$method, row$status, row$iterations_used, elapsed))
  if (!identical(hash_sources(), before)) stop("Challenge inputs or numerical sources changed during a request.")
}
after <- hash_sources()
stopifnot(identical(before, after))
utils::write.csv(after, file.path(out, "source-hashes-after.csv"), row.names = FALSE)
records <- do.call(rbind, records)
utils::write.csv(records[records$status == "error", , drop = FALSE], file.path(out, "errors.csv"),
  row.names = FALSE, na = "")
jsonlite::write_json(list(status = "complete", requests = nrow(records),
  returned = sum(records$status == "ok"), errors = sum(records$status == "error"),
  accuracy_computed = FALSE, source_fingerprint = fingerprint,
  elapsed_seconds = as.numeric(difftime(Sys.time(), started, units = "secs")),
  limitation = "These previously failed cases were selected. This is not an estimate of a population failure rate or interval accuracy."),
  file.path(out, "COMPLETE.json"), auto_unbox = TRUE, pretty = TRUE)
cat("Selected numerical challenge complete. No point estimates, intervals or coverage were computed.\n")
