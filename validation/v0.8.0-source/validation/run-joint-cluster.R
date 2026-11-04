#!/usr/bin/env Rscript
# Protocol: validation/JOINT_CLUSTER_PLAN.md. Pilot accuracy is not exported.
args <- commandArgs(trailingOnly = TRUE)
mode <- if (length(args)) match.arg(args[1L], c("pilot", "heldout")) else "pilot"
script <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1L])
root <- normalizePath(file.path(dirname(script), ".."), mustWork = TRUE)
out <- if (length(args) >= 2L) args[2L] else file.path(root, "validation", paste0("joint-cluster-", mode))
if (dir.exists(out) && length(list.files(out, all.files = TRUE, no.. = TRUE)))
  stop("Use a new or empty output directory; evidence is not overwritten.")
workers <- suppressWarnings(as.integer(Sys.getenv("CATMECH_BENCHMARK_WORKERS", "4")))
if (length(workers) != 1L || is.na(workers) || workers < 1L || workers > 16L)
  stop("CATMECH_BENCHMARK_WORKERS must be between 1 and 16.")
if (.Platform$OS.type == "windows") workers <- 1L
hash_sources <- function() {
  rel <- c("DESCRIPTION", "NAMESPACE", sort(list.files(file.path(root, "R"), "\\.R$", full.names = FALSE)),
    "validation/JOINT_CLUSTER_PLAN.md", "validation/run-joint-cluster.R")
  rel[grepl("\\.R$", rel) & !grepl("/", rel)] <- paste0("R/", rel[grepl("\\.R$", rel) & !grepl("/", rel)])
  data.frame(file = rel, sha256 = vapply(file.path(root, rel), function(f)
    digest::digest(file = f, algo = "sha256"), character(1)), row.names = NULL)
}
before <- hash_sources()
source_fingerprint <- digest::digest(before, algo = "sha256")
pkgload::load_all(root, quiet = TRUE)
if (!"cluster_model" %in% names(formals(catmech_reference))) stop("The joint model API is not available.")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
control <- list(chains = 4L, burnin = 1000L, iterations = 2000L, max_iterations = 8000L, prior_scale = 1)
cells <- c("gaussian_balanced", "gaussian_unbalanced", "near_zero_icc", "observed_mar", "null_slope", "student_t5")
sensitivity_cells <- c("gaussian_balanced", "near_zero_icc", "observed_mar")
reps <- if (mode == "pilot") 3L else 200L
settings <- list(mode = mode, m = 40L, cells = cells, reps = reps,
  sensitivity_cells = sensitivity_cells, sensitivity_reps = if (mode == "pilot") 3L else 100L,
  control = control, workers = workers, pan_version = as.character(utils::packageVersion("pan")),
  data_seed_base = if (mode == "pilot") 1510000000L else 1530000000L,
  imputation_seed_base = if (mode == "pilot") 1520000000L else 1540000000L,
  source_fingerprint = source_fingerprint,
  accuracy_exported = mode == "heldout")
utils::write.csv(before, file.path(out, "source-hashes-before.csv"), row.names = FALSE)
jsonlite::write_json(settings, file.path(out, "settings.json"), auto_unbox = TRUE, pretty = TRUE)
capture.output(utils::sessionInfo(), file = file.path(out, "session-info.txt"))

rmv <- function(n, covariance, heavy = FALSE) {
  z <- matrix(stats::rnorm(n * 2L), nrow = n, ncol = 2L) %*% chol(covariance)
  if (heavy) z <- z * sqrt(3 / stats::rchisq(n, df = 5))
  z
}
logistic_intercept <- function(eta, target) {
  stats::uniroot(function(a) mean(stats::plogis(a + eta)) - target, c(-50, 50), tol = 1e-12)$root
}
make_cohort <- function(cell, cell_index, replicate) {
  data_seed <- as.integer(settings$data_seed_base + 100000L * cell_index + replicate)
  imputation_seed <- as.integer(settings$imputation_seed_base + 100000L * cell_index + replicate)
  set.seed(data_seed)
  sizes <- if (cell == "gaussian_unbalanced") rep(c(5L, 8L, 12L, 15L), 6L) else rep(10L, 24L)
  ids <- rep(seq_len(24L), times = sizes)
  n <- length(ids)
  x <- stats::rnorm(n)
  Omega <- matrix(c(.5, .45, .45, .6), 2L, 2L)
  if (cell == "near_zero_icc") Omega <- Omega * .02
  Sigma <- matrix(c(.5, .55, .55, .8), 2L, 2L)
  u <- rmv(24L, Omega, cell == "student_t5")
  e <- rmv(n, Sigma, cell == "student_t5")
  truth_beta <- if (cell == "null_slope") 0 else .7
  truth <- .5 + truth_beta * x + u[ids, 1L] + e[, 1L]
  recorded <- .8 + 1.05 * x + u[ids, 2L] + e[, 2L]
  unaudited_sites <- sample(seq_len(24L), size = 6L, replace = FALSE)
  eligible <- !ids %in% unaudited_sites
  if (cell == "observed_mar") {
    miss_eta <- .8 * x
    miss_intercept <- logistic_intercept(miss_eta, .3)
    missing_recorded <- stats::runif(n) < stats::plogis(miss_intercept + miss_eta)
    recorded[missing_recorded] <- NA_real_
    w <- rep(0, n)
    w[!missing_recorded] <- as.numeric(scale(recorded[!missing_recorded]))
    audit_eta <- .6 * x + .6 * w
    audit_intercept <- logistic_intercept(audit_eta[eligible], (.25 * n) / sum(eligible))
    known <- eligible & (stats::runif(n) < stats::plogis(audit_intercept + audit_eta))
  } else {
    known <- rep(FALSE, n)
    known[sample(which(eligible), size = as.integer(n / 4), replace = FALSE)] <- TRUE
    missing_recorded <- rep(FALSE, n)
    for (g in seq_len(24L)) {
      rows <- which(ids == g)
      missing_recorded[sample(rows, size = ceiling(.3 * length(rows)), replace = FALSE)] <- TRUE
    }
    recorded[missing_recorded] <- NA_real_
    miss_intercept <- audit_intercept <- NA_real_
  }
  gold <- truth
  gold[!known] <- NA_real_
  input <- data.frame(W = recorded, gold = gold, x = x, site = sprintf("site%02d", ids))
  stopifnot(all(is.na(input$gold[!known])), identical(input$gold[known], truth[known]))
  # Private truth is never a column of input and is passed only to the oracle.
  mixed_stochastic_sites <- sum(vapply(split(seq_len(n), ids), function(rows)
    any(!known[rows] & !missing_recorded[rows]) && any(!known[rows] & missing_recorded[rows]), logical(1)))
  meta <- data.frame(cell = cell, cell_index = cell_index, replicate = replicate,
    data_seed = data_seed, imputation_seed = imputation_seed, n = n, sites = 24L,
    gold_n = sum(known), paired_n = sum(known & !missing_recorded),
    missing_recorded_n = sum(missing_recorded),
    gold_sites = length(unique(ids[known])), wholly_unaudited_sites = 24L - length(unique(ids[known])),
    mixed_stochastic_sites = mixed_stochastic_sites, truth = truth_beta,
    missing_logit_intercept = miss_intercept, audit_logit_intercept = audit_intercept,
    input_sha256 = digest::digest(input, algo = "sha256"), stringsAsFactors = FALSE)
  list(input = input, truth = truth, meta = meta)
}
methods_for <- function(cell, replicate) {
  c("joint", "plain_recorded", "oracle", if (cell %in% sensitivity_cells && replicate <= settings$sensitivity_reps)
    c("joint_prior_quarter", "joint_prior_four"))
}
extract_diagnostics <- function(z) {
  d <- if (inherits(z, "error")) z$model_diagnostics else z$model_diagnostics
  if (is.null(d)) return(list(max_rhat = NA_real_, min_ess_bulk = NA_real_, min_ess_tail = NA_real_,
    iterations_used = NA_integer_, attempts = "", summary = NULL))
  s <- d$summary
  max_valid <- function(v) if (length(v) && any(is.finite(v))) max(v[is.finite(v)]) else NA_real_
  min_valid <- function(v) if (length(v) && any(is.finite(v))) min(v[is.finite(v)]) else NA_real_
  list(max_rhat = max_valid(s$rhat), min_ess_bulk = min_valid(s$ess_bulk), min_ess_tail = min_valid(s$ess_tail),
    iterations_used = if (length(d$iterations_used) == 1L) as.integer(d$iterations_used) else NA_integer_,
    attempts = if (is.null(d$attempts)) "" else as.character(jsonlite::toJSON(d$attempts, auto_unbox = TRUE)), summary = s)
}
run_method <- function(method, p) {
  d <- p$input
  call <- list(data = d, value = "W", reference = "gold", predictors = "x", method = "normal",
    m = settings$m, seed = p$meta$imputation_seed, use_recorded = TRUE)
  if (method == "oracle") {
    call$data$gold <- p$truth
    call$m <- 2L
  } else if (startsWith(method, "joint")) {
    call$cluster <- "site"
    call$cluster_model <- "joint"
    call$joint_control <- control
    call$joint_control$prior_scale <- switch(method, joint_prior_quarter = .25, joint_prior_four = 4, 1)
  }
  fit <- do.call(catmech_reference, call)
  if (method != "oracle") {
    known <- !is.na(d$gold)
    if (!identical(fit$original, d)) stop("Original input was modified.")
    if (!all(vapply(fit$completed, function(z) identical(z$W[known], d$gold[known]), logical(1))))
      stop("Gold anchors were changed.")
  }
  ans <- catmech_pool_lm(fit, W ~ x, cluster = "site", cluster_variance = "cv3")
  ans <- ans[ans$term == "x", , drop = FALSE]
  if (nrow(ans) != 1L || any(!is.finite(unlist(ans[c("estimate", "std.error", "lower", "upper")]))))
    stop("Missing or nonfinite primary coefficient.")
  if (ans$std.error <= 0 || ans$upper <= ans$lower) stop("Degenerate interval.")
  list(answer = ans, model_diagnostics = fit$model_diagnostics)
}
run_replicate <- function(cell, index, replicate) {
  if (!identical(hash_sources(), before)) stop("Source changed before worker task.")
  p <- make_cohort(cell, index, replicate)
  records <- list()
  diagnostics <- list()
  for (method in methods_for(cell, replicate)) {
    warnings <- character()
    start <- proc.time()[["elapsed"]]
    result <- tryCatch(withCallingHandlers(run_method(method, p), warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning")
    }), error = identity)
    elapsed <- proc.time()[["elapsed"]] - start
    failed <- inherits(result, "error")
    diag <- extract_diagnostics(result)
    row <- cbind(p$meta, data.frame(method = method, status = if (failed) "error" else "ok",
      error_class = if (failed) class(result)[1L] else "", message = if (failed) conditionMessage(result) else "",
      warnings = paste(unique(warnings), collapse = " | "), estimate = NA_real_, se = NA_real_, df = NA_real_,
      lower = NA_real_, upper = NA_real_, width = NA_real_, error = NA_real_, covered = NA_integer_,
      max_rhat = diag$max_rhat, min_ess_bulk = diag$min_ess_bulk, min_ess_tail = diag$min_ess_tail,
      iterations_used = diag$iterations_used, attempts = diag$attempts, elapsed_seconds = elapsed,
      source_fingerprint = source_fingerprint, stringsAsFactors = FALSE))
    if (!failed && mode == "heldout") {
      a <- result$answer
      row$estimate <- a$estimate; row$se <- a$std.error; row$df <- a$df
      row$lower <- a$lower; row$upper <- a$upper; row$width <- a$upper - a$lower
      row$error <- a$estimate - p$meta$truth
      row$covered <- as.integer(a$lower <= p$meta$truth && a$upper >= p$meta$truth)
    }
    if (mode == "heldout" && failed) row$covered <- 0L
    if (!is.null(diag$summary) && nrow(diag$summary)) {
      diagnostics[[method]] <- cbind(cell = cell, replicate = replicate, method = method, diag$summary)
    }
    records[[method]] <- row
  }
  if (!identical(hash_sources(), before)) stop("Source changed during worker task.")
  list(records = do.call(rbind, records), diagnostics = if (length(diagnostics)) do.call(rbind, diagnostics) else NULL)
}
se_mean <- function(x) if (length(x) > 1L) stats::sd(x) / sqrt(length(x)) else NA_real_
wilson <- function(k, n) {
  if (!n) return(c(lower = NA_real_, upper = NA_real_))
  z <- stats::qnorm(.975); p <- k / n; denominator <- 1 + z^2 / n
  center <- (p + z^2 / (2 * n)) / denominator
  half <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / denominator
  c(lower = max(0, center - half), upper = min(1, center + half))
}
summarize_records <- function(records) {
  groups <- split(records, interaction(records$cell, records$method, drop = TRUE))
  do.call(rbind, lapply(groups, function(z) {
    good <- z$status == "ok"; n <- nrow(z); ok <- z[good, , drop = FALSE]; nk <- nrow(ok)
    ncover <- sum(z$covered); pc <- ncover / n; cc <- if (nk) ncover / nk else NA_real_
    f <- n - nk; pw <- wilson(ncover, n); cw <- wilson(ncover, nk); fw <- wilson(f, n)
    rmse <- if (nk) sqrt(mean(ok$error^2)) else NA_real_
    data.frame(cell = z$cell[1L], method = z$method[1L], requests = n, intervals = nk, failures = f,
      coverage_per_request = pc, coverage_per_request_mcse = sqrt(pc * (1 - pc) / n),
      coverage_per_request_lower = pw[1L], coverage_per_request_upper = pw[2L],
      coverage_conditional = cc, coverage_conditional_mcse = if (nk) sqrt(cc * (1 - cc) / nk) else NA_real_,
      coverage_conditional_lower = cw[1L], coverage_conditional_upper = cw[2L],
      failure_rate = f / n, failure_rate_mcse = sqrt((f / n) * (1 - f / n) / n),
      failure_rate_lower = fw[1L], failure_rate_upper = fw[2L],
      bias = if (nk) mean(ok$error) else NA_real_, bias_mcse = se_mean(ok$error),
      rmse = rmse, rmse_mcse = if (is.finite(rmse) && rmse > 0) se_mean(ok$error^2) / (2 * rmse) else NA_real_,
      empirical_sd = if (nk > 1L) stats::sd(ok$estimate) else NA_real_,
      mean_se = if (nk) mean(ok$se) else NA_real_, mean_width = if (nk) mean(ok$width) else NA_real_,
      mean_width_mcse = se_mean(ok$width), total_seconds = sum(z$elapsed_seconds),
      median_seconds = stats::median(z$elapsed_seconds),
      max_rhat = if (all(is.na(z$max_rhat))) NA_real_ else max(z$max_rhat, na.rm = TRUE),
      median_max_rhat = if (all(is.na(z$max_rhat))) NA_real_ else stats::median(z$max_rhat, na.rm = TRUE),
      min_ess_bulk = if (all(is.na(z$min_ess_bulk))) NA_real_ else min(z$min_ess_bulk, na.rm = TRUE),
      min_ess_tail = if (all(is.na(z$min_ess_tail))) NA_real_ else min(z$min_ess_tail, na.rm = TRUE),
      stringsAsFactors = FALSE, row.names = NULL)
  }))
}
paired_records <- function(records) {
  all <- list()
  for (cell in cells) {
    z <- records[records$cell == cell, , drop = FALSE]
    for (method in setdiff(unique(z$method), "plain_recorded")) {
      baseline <- if (method %in% c("joint_prior_quarter", "joint_prior_four")) "joint" else "plain_recorded"
      a <- z[z$method == method, ]; b <- z[z$method == baseline, ]
      pair <- merge(a, b, by = "replicate", suffixes = c("_method", "_baseline"))
      stopifnot(all(pair$input_sha256_method == pair$input_sha256_baseline))
      complete <- pair$status_method == "ok" & pair$status_baseline == "ok"
      for (metric in c("coverage_per_request", "coverage_complete_pairs", "error", "squared_error", "width")) {
        q <- if (metric == "coverage_per_request") pair else pair[complete, , drop = FALSE]
        diff <- switch(metric,
          coverage_per_request = q$covered_method - q$covered_baseline,
          coverage_complete_pairs = q$covered_method - q$covered_baseline,
          error = q$error_method - q$error_baseline,
          squared_error = q$error_method^2 - q$error_baseline^2,
          width = q$width_method - q$width_baseline)
        all[[length(all) + 1L]] <- data.frame(cell = cell, method = method, baseline = baseline, metric = metric,
          pairs = nrow(q), difference = if (length(diff)) mean(diff) else NA_real_, mcse = se_mean(diff))
      }
    }
  }
  do.call(rbind, all)
}

records <- diagnostics <- list()
started <- Sys.time()
for (index in seq_along(cells)) {
  cell <- cells[index]
  for (chunk in split(seq_len(reps), ceiling(seq_len(reps) / 10L))) {
    batch <- parallel::mclapply(chunk, function(r) run_replicate(cell, index, r),
      mc.cores = workers, mc.preschedule = TRUE, mc.set.seed = FALSE)
    if (any(vapply(batch, inherits, logical(1), what = "try-error")))
      stop("Worker infrastructure failure; partial artifacts are not a completed study.")
    records <- c(records, lapply(batch, `[[`, "records"))
    diagnostics <- c(diagnostics, lapply(batch, `[[`, "diagnostics"))
    raw <- do.call(rbind, records)
    raw <- raw[order(raw$cell_index, raw$replicate, raw$method), ]
    rownames(raw) <- NULL
    utils::write.csv(raw, file.path(out, "records.csv"), row.names = FALSE, na = "")
    diagnostic_rows <- do.call(rbind, Filter(Negate(is.null), diagnostics))
    if (!is.null(diagnostic_rows)) utils::write.csv(diagnostic_rows,
      file.path(out, "diagnostics.csv"), row.names = FALSE, na = "")
    cat(sprintf("%s: replicates through %d/%d; %d requests saved; %.1f seconds elapsed\n",
      cell, max(chunk), reps, nrow(raw), as.numeric(difftime(Sys.time(), started, units = "secs"))))
  }
}
expected <- length(cells) * reps * 3L + length(sensitivity_cells) * min(reps, settings$sensitivity_reps) * 2L
stopifnot(nrow(raw) == expected, !anyDuplicated(raw[c("cell", "replicate", "method")]))
if (mode == "heldout") {
  utils::write.csv(summarize_records(raw), file.path(out, "summary.csv"), row.names = FALSE, na = "")
  utils::write.csv(paired_records(raw), file.path(out, "paired.csv"), row.names = FALSE, na = "")
}
errors <- raw[raw$status != "ok", c("cell", "replicate", "method", "error_class", "message", "warnings")]
utils::write.csv(errors, file.path(out, "errors.csv"), row.names = FALSE, na = "")
after <- hash_sources()
utils::write.csv(after, file.path(out, "source-hashes-after.csv"), row.names = FALSE)
if (!identical(before, after)) stop("Sources changed during study; no COMPLETE marker is written.")
jsonlite::write_json(list(status = "complete", mode = mode, requests = nrow(raw),
  intervals = sum(raw$status == "ok"), errors = sum(raw$status != "ok"),
  started_utc = format(started, tz = "UTC", usetz = TRUE),
  finished_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  elapsed_seconds = as.numeric(difftime(Sys.time(), started, units = "secs")),
  source_fingerprint = source_fingerprint), file.path(out, "COMPLETE.json"), auto_unbox = TRUE, pretty = TRUE)
cat(sprintf("Complete: %d requests, %d intervals, %d errors.\n", nrow(raw), sum(raw$status == "ok"), sum(raw$status != "ok")))
