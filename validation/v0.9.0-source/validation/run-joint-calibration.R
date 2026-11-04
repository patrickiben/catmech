#!/usr/bin/env Rscript
# Frozen stages and estimands: validation/JOINT_CALIBRATION_PLAN.md.
args <- commandArgs(trailingOnly = TRUE)
mode <- if (length(args)) match.arg(args[1L], c("pilot", "development", "confirmation")) else "pilot"
script <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1L])
root <- normalizePath(file.path(dirname(script), ".."), mustWork = TRUE)
out <- if (length(args) >= 2L) args[2L] else file.path(root, "validation", paste0("joint-calibration-", mode))
if (dir.exists(out) && length(list.files(out, all.files = TRUE, no.. = TRUE)))
  stop("Use a new or empty output directory; evidence is not overwritten.")
workers <- suppressWarnings(as.integer(Sys.getenv("CATMECH_BENCHMARK_WORKERS", "8")))
if (length(workers) != 1L || is.na(workers) || workers < 1L || workers > 16L) stop("Invalid worker count.")
if (.Platform$OS.type == "windows") workers <- 1L
baseline_root <- file.path(root, "validation/v0.8.0-source")
hash_sources <- function() {
  rel <- c("DESCRIPTION", "NAMESPACE", "validation/JOINT_CALIBRATION_PLAN.md",
    "validation/run-joint-calibration.R", paste0("R/", sort(list.files(file.path(root, "R"), "\\.R$"))),
    "validation/v0.8.0-source/DESCRIPTION", "validation/v0.8.0-source/NAMESPACE",
    paste0("validation/v0.8.0-source/R/", sort(list.files(file.path(baseline_root, "R"), "\\.R$"))))
  data.frame(file = rel, bytes = as.numeric(file.info(file.path(root, rel))$size),
    sha256 = vapply(file.path(root, rel), function(f) digest::digest(file = f, algo = "sha256"), character(1)),
    row.names = NULL)
}
before <- hash_sources()
source_fingerprint <- digest::digest(before, algo = "sha256")
pkgload::load_all(root, quiet = TRUE)
if (!exists("catmech_joint_summary", mode = "function")) stop("Candidate posterior-summary API is unavailable.")
old <- new.env(parent = globalenv())
for (f in sort(list.files(file.path(baseline_root, "R"), "\\.R$", full.names = TRUE))) sys.source(f, envir = old)
dir.create(out, recursive = TRUE, showWarnings = FALSE)
cells <- c("gaussian_balanced", "gaussian_unbalanced", "near_zero_icc", "observed_mar", "null_slope", "student_t5", "contextual_x")
prior_cells <- c("gaussian_balanced", "near_zero_icc", "observed_mar")
primary_cells <- c("gaussian_balanced", "gaussian_unbalanced", "observed_mar", "null_slope")
new_control <- list(chains = 4L, burnin = 1000L, iterations = 4000L, max_iterations = 32000L, prior_scale = 1)
old_control <- list(chains = 4L, burnin = 1000L, iterations = 2000L, max_iterations = 8000L, prior_scale = 1)
config <- list(mode = mode, workers = workers, m = 40L, cells = cells,
  primary_cells = primary_cells, prior_cells = prior_cells,
  primary_reps = switch(mode, pilot = 1L, development = 30L, confirmation = 400L),
  contextual_reps = switch(mode, pilot = 1L, development = 30L, confirmation = 200L),
  prior_reps = switch(mode, pilot = 1L, development = 30L, confirmation = 200L),
  data_seed_base = switch(mode, pilot = 1590000000L, development = 1610000000L, confirmation = 1630000000L),
  imputation_seed_base = switch(mode, pilot = 1600000000L, development = 1620000000L, confirmation = 1640000000L),
  bootstrap_seed = switch(mode, pilot = NA_integer_, development = 1650000000L, confirmation = 1660000000L),
  bootstrap_reps = 2000L, new_control = new_control, old_control = old_control,
  baseline_source_commit = "b236293", pan_version = as.character(utils::packageVersion("pan")),
  source_fingerprint = source_fingerprint, accuracy_exported = mode != "pilot")
utils::write.csv(before, file.path(out, "source-hashes-before.csv"), row.names = FALSE)
jsonlite::write_json(config, file.path(out, "settings.json"), auto_unbox = TRUE, pretty = TRUE, na = "null")
capture.output(utils::sessionInfo(), file = file.path(out, "session-info.txt"))

rmv <- function(n, covariance, heavy = FALSE) {
  z <- matrix(stats::rnorm(n * 2L), n, 2L) %*% chol(covariance)
  if (heavy) z <- z * sqrt(3 / stats::rchisq(n, 5))
  z
}
logit_intercept <- function(eta, target) stats::uniroot(function(a)
  mean(stats::plogis(a + eta)) - target, c(-50, 50), tol = 1e-12)$root
make_cohort <- function(cell, index, replicate) {
  seed <- as.integer(config$data_seed_base + index * 100000L + replicate)
  imp_seed <- as.integer(config$imputation_seed_base + index * 100000L + replicate)
  set.seed(seed)
  sizes <- if (cell == "gaussian_unbalanced") rep(c(5L, 8L, 12L, 15L), 6L) else rep(10L, 24L)
  ids <- rep(seq_len(24L), sizes); n <- length(ids)
  if (cell == "contextual_x") {
    a <- stats::rnorm(24L)
    x <- a[ids] + stats::rnorm(n)
  } else { a <- rep(0, 24L); x <- stats::rnorm(n) }
  Omega <- matrix(c(.5, .45, .45, .6), 2L, 2L)
  if (cell == "near_zero_icc") Omega <- .02 * Omega
  Sigma <- matrix(c(.5, .55, .55, .8), 2L, 2L)
  u <- rmv(24L, Omega, cell == "student_t5")
  e <- rmv(n, Sigma, cell == "student_t5")
  structural <- if (cell == "null_slope") 0 else .7
  target <- if (cell == "contextual_x") .9 else structural
  truth <- .5 + structural * x + .4 * a[ids] + u[ids, 1L] + e[, 1L]
  W <- .8 + 1.05 * x + u[ids, 2L] + e[, 2L]
  eligible <- !ids %in% sample(seq_len(24L), 6L)
  if (cell == "observed_mar") {
    miss_intercept <- logit_intercept(.8 * x, .3)
    missing_W <- stats::runif(n) < stats::plogis(miss_intercept + .8 * x)
    W[missing_W] <- NA_real_
    w <- rep(0, n); w[!missing_W] <- as.numeric(scale(W[!missing_W]))
    eta <- .6 * x + .6 * w
    audit_intercept <- logit_intercept(eta[eligible], .25 * n / sum(eligible))
    known <- eligible & (stats::runif(n) < stats::plogis(audit_intercept + eta))
  } else {
    known <- rep(FALSE, n); known[sample(which(eligible), as.integer(n / 4L))] <- TRUE
    missing_W <- rep(FALSE, n)
    for (g in seq_len(24L)) {
      rows <- which(ids == g)
      missing_W[sample(rows, ceiling(.3 * length(rows)))] <- TRUE
    }
    W[missing_W] <- NA_real_
    miss_intercept <- audit_intercept <- NA_real_
  }
  gold <- truth; gold[!known] <- NA_real_
  input <- data.frame(W = W, gold = gold, x = x, site = sprintf("site%02d", ids))
  stopifnot(identical(input$gold[known], truth[known]), all(is.na(input$gold[!known])))
  meta <- data.frame(cell = cell, cell_index = index, replicate = replicate,
    data_seed = seed, imputation_seed = imp_seed, n = n, sites = 24L,
    gold_n = sum(known), paired_n = sum(known & !missing_W), missing_recorded_n = sum(missing_W),
    wholly_unaudited_sites = 24L - length(unique(ids[known])), target = target,
    structural_slope = structural, missing_logit_intercept = miss_intercept, audit_logit_intercept = audit_intercept,
    input_sha256 = digest::digest(input, algo = "sha256"), stringsAsFactors = FALSE)
  list(input = input, truth = truth, meta = meta)
}
fit_specs <- function(cell, replicate) {
  c("old", "new", "oracle", if (cell %in% prior_cells && replicate <= config$prior_reps)
    c("new_prior_quarter", "new_prior_four"))
}
diag_info <- function(result) {
  d <- result$model_diagnostics
  summary <- d$summary
  extremum <- function(x, op) if (length(x) && any(is.finite(x))) op(x[is.finite(x)]) else NA_real_
  list(summary = summary, max_rhat = extremum(summary$rhat, max),
    min_bulk = extremum(summary$ess_bulk, min), min_tail = extremum(summary$ess_tail, min),
    iterations = if (length(d$iterations_used) == 1L) as.integer(d$iterations_used) else NA_integer_,
    attempts = if (is.null(d$attempts)) "" else as.character(jsonlite::toJSON(d$attempts, auto_unbox = TRUE)))
}
capture <- function(fun) {
  warnings <- character(); start <- proc.time()[["elapsed"]]
  value <- tryCatch(withCallingHandlers(fun(), warning = function(w) {
    warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning")
  }), error = identity)
  list(value = value, warnings = paste(unique(warnings), collapse = " | "),
    elapsed = proc.time()[["elapsed"]] - start)
}
fit_one <- function(id, p) {
  d <- p$input
  call <- list(data = d, value = "W", reference = "gold", predictors = "x", method = "normal",
    m = config$m, seed = p$meta$imputation_seed, use_recorded = TRUE)
  fun <- if (id %in% c("old", "oracle")) old$catmech_reference else catmech_reference
  if (id == "oracle") {
    call$data$gold <- p$truth; call$m <- 2L
  } else {
    call$cluster <- "site"; call$cluster_model <- "joint"
    call$joint_control <- if (id == "old") old_control else new_control
    call$joint_control$prior_scale <- switch(id, new_prior_quarter = .25, new_prior_four = 4, 1)
  }
  fit <- do.call(fun, call)
  known <- !is.na(call$data$gold)
  if (!identical(fit$original, call$data)) stop("Original input was changed.")
  if (!all(vapply(fit$completed, function(d) identical(d$W[known], call$data$gold[known]), logical(1))))
    stop("Reference anchors were changed.")
  fit
}
analyze_one <- function(fit, id, analysis) {
  if (analysis == "direct") {
    z <- catmech_joint_summary(fit, level = .95)
    z <- z[z$term == "x", , drop = FALSE]
    if (nrow(z) != 1L) stop("Posterior coefficient x is absent.")
    ans <- list(estimate = z$estimate, spread = z$posterior_sd, lower = z$lower, upper = z$upper,
      df = NA_real_, within_variance = NA_real_, between_variance = NA_real_,
      posterior_rhat = z$rhat, posterior_ess_bulk = z$ess_bulk, posterior_ess_tail = z$ess_tail,
      mcse_mean = z$mcse_mean, mcse_lower = z$mcse_lower, mcse_upper = z$mcse_upper,
      posterior_draws = z$draws, interpretation = z$interpretation)
  } else {
    fun <- if (id %in% c("old", "oracle")) old$catmech_pool_lm else catmech_pool_lm
    z <- fun(fit, W ~ x, cluster = "site", cluster_variance = "cv3")
    z <- z[z$term == "x", , drop = FALSE]
    if (nrow(z) != 1L) stop("Pooled coefficient x is absent.")
    ans <- list(estimate = z$estimate, spread = z$std.error, lower = z$lower, upper = z$upper,
      df = z$df, within_variance = z$within_variance, between_variance = z$between_variance,
      posterior_rhat = NA_real_, posterior_ess_bulk = NA_real_, posterior_ess_tail = NA_real_,
      mcse_mean = NA_real_, mcse_lower = NA_real_, mcse_upper = NA_real_, posterior_draws = NA_integer_,
      interpretation = z$interpretation)
  }
  if (any(!is.finite(unlist(ans[c("estimate", "spread", "lower", "upper")])))) stop("Nonfinite analysis summary.")
  if (ans$spread <= 0 || ans$upper <= ans$lower) stop("Degenerate interval.")
  ans
}
run_replicate <- function(cell, index, replicate) {
  if (!identical(hash_sources(), before)) stop("Source changed before worker task.")
  p <- make_cohort(cell, index, replicate)
  records <- fits <- diagnostics <- list()
  for (id in fit_specs(cell, replicate)) {
    fit <- capture(function() fit_one(id, p))
    fit_failed <- inherits(fit$value, "error")
    d <- diag_info(fit$value)
    fit_row <- cbind(p$meta, data.frame(fit_id = id,
      prior_scale = switch(id, new_prior_quarter = .25, new_prior_four = 4, 1),
      fit_status = if (fit_failed) "error" else "ok", error_class = if (fit_failed) class(fit$value)[1L] else "",
      message = if (fit_failed) conditionMessage(fit$value) else "", warnings = fit$warnings,
      elapsed_seconds = fit$elapsed, max_rhat = d$max_rhat, min_ess_bulk = d$min_bulk,
      min_ess_tail = d$min_tail, iterations_used = d$iterations, attempts = d$attempts,
      source_fingerprint = source_fingerprint, stringsAsFactors = FALSE))
    fits[[id]] <- fit_row
    if (!is.null(d$summary) && nrow(d$summary)) diagnostics[[id]] <- cbind(
      cell = cell, replicate = replicate, fit_id = id, d$summary)
    analyses <- if (startsWith(id, "new")) c("cv3", "direct") else "cv3"
    for (analysis in analyses) {
      method <- if (analysis == "direct") sub("^new", "direct", id) else
        if (startsWith(id, "new")) sub("^new", "new_cv3", id) else paste0(id, "_cv3")
      analyzed <- if (fit_failed) list(value = fit$value, elapsed = 0, warnings = "") else
        capture(function() analyze_one(fit$value, id, analysis))
      failed <- inherits(analyzed$value, "error")
      row <- cbind(p$meta, data.frame(method = method, fit_id = id, analysis = analysis,
        prior_scale = fit_row$prior_scale, status = if (failed) "error" else "ok",
        failure_stage = if (!failed) "" else if (fit_failed) "fit" else "analysis",
        error_class = if (failed) class(analyzed$value)[1L] else "",
        message = if (failed) conditionMessage(analyzed$value) else "",
        warnings = paste(Filter(nzchar, c(fit$warnings, analyzed$warnings)), collapse = " | "),
        estimate = NA_real_, spread = NA_real_, spread_type = if (analysis == "direct") "posterior_sd" else "pooled_se",
        lower = NA_real_, upper = NA_real_, width = NA_real_, df = NA_real_,
        within_variance = NA_real_, between_variance = NA_real_, error = NA_real_, covered = NA_integer_,
        interval_score = NA_real_, posterior_rhat = NA_real_, posterior_ess_bulk = NA_real_,
        posterior_ess_tail = NA_real_, mcse_mean = NA_real_, mcse_lower = NA_real_, mcse_upper = NA_real_,
        posterior_draws = NA_integer_, interpretation = "", analysis_seconds = analyzed$elapsed,
        source_fingerprint = source_fingerprint, stringsAsFactors = FALSE))
      if (!failed && mode != "pilot") {
        a <- analyzed$value
        for (field in names(a)) row[[field]] <- a[[field]]
        row$error <- a$estimate - p$meta$target; row$width <- a$upper - a$lower
        row$covered <- as.integer(a$lower <= p$meta$target && a$upper >= p$meta$target)
        row$interval_score <- row$width + 40 * max(0, a$lower - p$meta$target) + 40 * max(0, p$meta$target - a$upper)
      }
      if (failed && mode != "pilot") row$covered <- 0L
      records[[method]] <- row
    }
  }
  if (!identical(hash_sources(), before)) stop("Source changed during worker task.")
  list(records = do.call(rbind, records), fits = do.call(rbind, fits),
    diagnostics = if (length(diagnostics)) do.call(rbind, diagnostics) else NULL)
}
se_mean <- function(x) if (length(x) > 1L) stats::sd(x) / sqrt(length(x)) else NA_real_
wilson <- function(k, n) {
  if (!n) return(c(lower = NA_real_, upper = NA_real_))
  z <- stats::qnorm(.975); p <- k / n; d <- 1 + z^2 / n
  center <- (p + z^2 / (2 * n)) / d; half <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / d
  c(lower = max(0, center - half), upper = min(1, center + half))
}
summary_rows <- function(records) do.call(rbind, lapply(split(records,
  interaction(records$cell, records$method, drop = TRUE)), function(z) {
  ok <- z[z$status == "ok", , drop = FALSE]; n <- nrow(z); nk <- nrow(ok)
  covered <- sum(z$covered); pc <- covered / n; cc <- if (nk) covered / nk else NA_real_
  pw <- wilson(covered, n); cw <- wilson(covered, nk); fw <- wilson(n - nk, n)
  rmse <- if (nk) sqrt(mean(ok$error^2)) else NA_real_
  row <- data.frame(cell = z$cell[1L], method = z$method[1L], requests = n, intervals = nk, failures = n - nk,
    coverage_per_request = pc, coverage_per_request_mcse = sqrt(pc * (1 - pc) / n),
    coverage_per_request_lower = pw[1L], coverage_per_request_upper = pw[2L],
    coverage_conditional = cc, coverage_conditional_mcse = if (nk) sqrt(cc * (1 - cc) / nk) else NA_real_,
    coverage_conditional_lower = cw[1L], coverage_conditional_upper = cw[2L],
    failure_rate = (n - nk) / n, failure_rate_mcse = sqrt((n - nk) / n * nk / n / n),
    failure_rate_lower = fw[1L], failure_rate_upper = fw[2L],
    bias = if (nk) mean(ok$error) else NA_real_, bias_mcse = se_mean(ok$error), rmse = rmse,
    rmse_mcse = if (is.finite(rmse) && rmse > 0) se_mean(ok$error^2) / (2 * rmse) else NA_real_,
    empirical_sd = if (nk > 1L) stats::sd(ok$estimate) else NA_real_,
    mean_spread = if (nk) mean(ok$spread) else NA_real_,
    rms_spread = if (nk) sqrt(mean(ok$spread^2)) else NA_real_, stringsAsFactors = FALSE, row.names = NULL)
  for (name in c("width", "interval_score", "within_variance", "between_variance", "df", "mcse_mean", "mcse_lower", "mcse_upper")) {
    v <- ok[[name]]; v <- v[is.finite(v)]
    row[[paste0("mean_", name)]] <- if (length(v)) mean(v) else NA_real_
    row[[paste0("mean_", name, "_mcse")]] <- se_mean(v)
  }
  row
}))
paired_rows <- function(records) {
  result <- list()
  for (cell in cells) {
    z <- records[records$cell == cell, , drop = FALSE]
    comparisons <- list(c("new_cv3", "old_cv3"), c("direct", "new_cv3"), c("direct", "old_cv3"))
    if (cell %in% prior_cells) for (analysis in c("new_cv3", "direct"))
      for (suffix in c("_prior_quarter", "_prior_four")) comparisons[[length(comparisons) + 1L]] <- c(paste0(analysis, suffix), analysis)
    for (comparison in comparisons) {
      a <- z[z$method == comparison[1L], ]; b <- z[z$method == comparison[2L], ]
      pair <- merge(a, b, by = "replicate", suffixes = c("_method", "_baseline"))
      stopifnot(all(pair$input_sha256_method == pair$input_sha256_baseline))
      complete <- pair$status_method == "ok" & pair$status_baseline == "ok"
      for (metric in c("coverage_per_request", "coverage_complete_pairs", "error", "squared_error", "width", "interval_score")) {
        q <- if (metric == "coverage_per_request") pair else pair[complete, , drop = FALSE]
        delta <- switch(metric, coverage_per_request = q$covered_method - q$covered_baseline,
          coverage_complete_pairs = q$covered_method - q$covered_baseline,
          error = q$error_method - q$error_baseline,
          squared_error = q$error_method^2 - q$error_baseline^2,
          width = q$width_method - q$width_baseline, interval_score = q$interval_score_method - q$interval_score_baseline)
        result[[length(result) + 1L]] <- data.frame(cell = cell, method = comparison[1L], baseline = comparison[2L],
          metric = metric, pairs = nrow(q), difference = if (length(delta)) mean(delta) else NA_real_, mcse = se_mean(delta))
      }
    }
  }
  do.call(rbind, result)
}
prior_rows <- function(records) {
  triples <- list(); summary <- list()
  for (cell in prior_cells) {
    z <- records[records$cell == cell & records$replicate <= config$prior_reps, , drop = FALSE]
    for (replicate in sort(unique(z$replicate))) {
      r <- z[z$replicate == replicate, , drop = FALSE]
      row <- data.frame(cell = cell, replicate = replicate, direct_complete = FALSE, cv3_complete = FALSE,
        direct_sensitivity = NA_real_, cv3_sensitivity = NA_real_,
        direct_estimate_shift = NA_real_, cv3_estimate_shift = NA_real_, difference = NA_real_)
      for (analysis in c("direct", "new_cv3")) {
        selected <- r[match(paste0(analysis, c("", "_prior_quarter", "_prior_four")), r$method), ]
        ok <- nrow(selected) == 3L && all(!is.na(selected$status) & selected$status == "ok")
        prefix <- if (analysis == "direct") "direct" else "cv3"
        row[[paste0(prefix, "_complete")]] <- ok
        if (ok) {
          row[[paste0(prefix, "_sensitivity")]] <- mean(abs(log(selected$width[-1L] / selected$width[1L])))
          row[[paste0(prefix, "_estimate_shift")]] <- mean(abs(selected$estimate[-1L] - selected$estimate[1L]))
        }
      }
      if (row$direct_complete && row$cv3_complete) row$difference <- row$direct_sensitivity - row$cv3_sensitivity
      triples[[length(triples) + 1L]] <- row
    }
  }
  triples <- do.call(rbind, triples)
  for (cell in prior_cells) {
    z <- triples[triples$cell == cell, ]; good <- z[is.finite(z$difference), ]
    summary[[cell]] <- data.frame(cell = cell, requested_triples = nrow(z), direct_triples = sum(z$direct_complete),
      cv3_triples = sum(z$cv3_complete), complete_pairs = nrow(good),
      mean_direct_sensitivity = if (nrow(good)) mean(good$direct_sensitivity) else NA_real_,
      mean_cv3_sensitivity = if (nrow(good)) mean(good$cv3_sensitivity) else NA_real_,
      difference = if (nrow(good)) mean(good$difference) else NA_real_, mcse = se_mean(good$difference))
  }
  summary <- do.call(rbind, summary); rownames(summary) <- NULL
  overall <- mean(summary$difference); mcse <- sqrt(sum(summary$mcse^2)) / length(prior_cells)
  list(triples = triples, summary = summary, gate = list(difference = overall, mcse = mcse,
    lower = overall - stats::qnorm(.975) * mcse, upper = overall + stats::qnorm(.975) * mcse,
    no_positive_cell = all(is.finite(summary$difference) & summary$difference <= 0),
    pass = is.finite(overall) && is.finite(mcse) && overall + stats::qnorm(.975) * mcse < 0 &&
      all(is.finite(summary$difference) & summary$difference <= 0)))
}
calibration_screen <- function(records) {
  methods <- c("old_cv3", "new_cv3", "direct")
  matrices <- lapply(primary_cells, function(cell) {
    z <- records[records$cell == cell, ]; reps <- sort(unique(z$replicate))
    coverage <- matrix(NA_real_, length(reps), length(methods), dimnames = list(reps, methods))
    for (method in methods) {
      r <- z[z$method == method, ]; r <- r[match(reps, r$replicate), ]
      coverage[, method] <- ifelse(r$status == "ok", r$covered, NA_real_)
    }
    coverage
  })
  names(matrices) <- primary_cells
  conditional <- vapply(matrices, function(z) colMeans(z, na.rm = TRUE), numeric(length(methods)))
  distances <- rowMeans(abs(conditional - .95))
  direct_rows <- records[records$method == "direct" & records$cell %in% primary_cells, ]
  biases <- vapply(primary_cells, function(cell) mean(direct_rows$error[direct_rows$cell == cell & direct_rows$status == "ok"]), numeric(1))
  RNGkind("Mersenne-Twister", "Inversion", "Rejection"); set.seed(config$bootstrap_seed)
  bootstrap <- matrix(NA_real_, config$bootstrap_reps, 2L,
    dimnames = list(NULL, c("direct_minus_old", "direct_minus_new_cv3")))
  for (b in seq_len(config$bootstrap_reps)) {
    d <- matrix(NA_real_, length(methods), length(primary_cells))
    for (j in seq_along(primary_cells)) {
      z <- matrices[[j]]; chosen <- sample.int(nrow(z), nrow(z), replace = TRUE)
      d[, j] <- abs(colMeans(z[chosen, , drop = FALSE], na.rm = TRUE) - .95)
    }
    score <- rowMeans(d)
    bootstrap[b, ] <- c(score[3L] - score[1L], score[3L] - score[2L])
  }
  bounds <- apply(bootstrap, 2L, function(v) if (all(is.finite(v)))
    stats::quantile(v, probs = c(.025, .975), type = 7, names = FALSE) else c(NA_real_, NA_real_))
  finite <- all(is.finite(c(conditional, biases, bootstrap, bounds)))
  gate <- list(conditional_coverage = setNames(lapply(seq_along(methods), function(k) setNames(as.list(conditional[k, ]), primary_cells)), methods),
    methods = methods, mean_absolute_distance = setNames(as.list(distances), methods),
    direct_bias = setNames(as.list(biases), primary_cells),
    coverage_band_pass = all(is.finite(conditional[3L, ]) & conditional[3L, ] >= .925 & conditional[3L, ] <= .975),
    bias_bound_pass = all(is.finite(biases) & abs(biases) <= .02),
    direct_minus_old = list(difference = distances[3L] - distances[1L], lower = bounds[1L, 1L], upper = bounds[2L, 1L]),
    direct_minus_new_cv3 = list(difference = distances[3L] - distances[2L], lower = bounds[1L, 2L], upper = bounds[2L, 2L]))
  gate$pass <- finite && gate$coverage_band_pass && gate$bias_bound_pass &&
    all(distances[3L] < distances[1:2]) && all(bounds[2L, ] < 0)
  list(bootstrap = data.frame(bootstrap = seq_len(config$bootstrap_reps), bootstrap), gate = gate)
}
availability_screen <- function(fits) {
  rows <- lapply(setdiff(cells, "contextual_x"), function(cell) {
    z <- fits[fits$cell == cell, ]; a <- z[z$fit_id == "new", ]; b <- z[z$fit_id == "old", ]
    pair <- merge(a, b, by = "replicate", suffixes = c("_new", "_old"))
    new_fail <- pair$fit_status_new != "ok"; old_fail <- pair$fit_status_old != "ok"
    delta <- as.numeric(new_fail) - as.numeric(old_fail)
    data.frame(cell = cell, requests = nrow(pair), old_failures = sum(old_fail), new_failures = sum(new_fail),
      old_failure_rate = mean(old_fail), new_failure_rate = mean(new_fail), difference = mean(delta), mcse = se_mean(delta))
  })
  rows <- do.call(rbind, rows)
  delta <- mean(rows$difference); mcse <- sqrt(sum(rows$mcse^2)) / nrow(rows)
  list(cells = rows, gate = list(difference = delta, mcse = mcse,
    lower = delta - stats::qnorm(.975) * mcse, upper = delta + stats::qnorm(.975) * mcse,
    every_return_rate_at_least_99 = all(rows$new_failure_rate <= .01),
    no_observed_failure_increase = all(rows$new_failure_rate <= rows$old_failure_rate),
    pass = all(rows$new_failure_rate <= .01) && all(rows$new_failure_rate <= rows$old_failure_rate) &&
      is.finite(mcse) && delta + stats::qnorm(.975) * mcse < 0))
}

all_records <- all_fits <- all_diagnostics <- list()
started <- Sys.time()
for (index in seq_along(cells)) {
  cell <- cells[index]
  reps <- if (cell == "contextual_x") config$contextual_reps else config$primary_reps
  for (chunk in split(seq_len(reps), ceiling(seq_len(reps) / 16L))) {
    batch <- parallel::mclapply(chunk, function(r) run_replicate(cell, index, r),
      mc.cores = workers, mc.preschedule = TRUE, mc.set.seed = FALSE)
    if (any(vapply(batch, inherits, logical(1), what = "try-error")))
      stop("Worker infrastructure error; saved checkpoints are not a completed study.")
    all_records <- c(all_records, lapply(batch, `[[`, "records"))
    all_fits <- c(all_fits, lapply(batch, `[[`, "fits"))
    all_diagnostics <- c(all_diagnostics, lapply(batch, `[[`, "diagnostics"))
    records <- do.call(rbind, all_records)
    records <- records[order(records$cell_index, records$replicate, records$method), ]; rownames(records) <- NULL
    fits <- do.call(rbind, all_fits)
    fits <- fits[order(fits$cell_index, fits$replicate, fits$fit_id), ]; rownames(fits) <- NULL
    diagnostics <- do.call(rbind, Filter(Negate(is.null), all_diagnostics))
    utils::write.csv(records, file.path(out, "records.csv"), row.names = FALSE, na = "")
    utils::write.csv(fits, file.path(out, "fit-records.csv"), row.names = FALSE, na = "")
    if (!is.null(diagnostics)) utils::write.csv(diagnostics, file.path(out, "diagnostics.csv"), row.names = FALSE, na = "")
    cat(sprintf("%s: %d/%d cohorts; %d analysis requests, %d fitted-workflow requests; %.1fs elapsed\n",
      cell, max(chunk), reps, nrow(records), nrow(fits), as.numeric(difftime(Sys.time(), started, units = "secs"))))
  }
}
cohorts <- 6L * config$primary_reps + config$contextual_reps
expected_analyses <- 4L * cohorts + 12L * config$prior_reps
expected_fits <- 3L * cohorts + 6L * config$prior_reps
stopifnot(nrow(records) == expected_analyses, nrow(fits) == expected_fits,
  !anyDuplicated(records[c("cell", "replicate", "method")]), !anyDuplicated(fits[c("cell", "replicate", "fit_id")]))
if (mode != "pilot") {
  utils::write.csv(summary_rows(records), file.path(out, "summary.csv"), row.names = FALSE, na = "")
  utils::write.csv(paired_rows(records), file.path(out, "paired.csv"), row.names = FALSE, na = "")
  prior <- prior_rows(records)
  utils::write.csv(prior$triples, file.path(out, "prior-triples.csv"), row.names = FALSE, na = "")
  utils::write.csv(prior$summary, file.path(out, "prior-summary.csv"), row.names = FALSE, na = "")
  calibration <- calibration_screen(records)
  utils::write.csv(calibration$bootstrap, file.path(out, "calibration-bootstrap.csv"), row.names = FALSE, na = "")
  availability <- availability_screen(fits)
  utils::write.csv(availability$cells, file.path(out, "availability.csv"), row.names = FALSE, na = "")
  jsonlite::write_json(list(stage = mode, calibration = calibration$gate, prior_sensitivity = prior$gate,
    availability = availability$gate, activate_prior_development_branch = mode == "development" && !isTRUE(prior$gate$pass),
    interpretation = if (mode == "development") "development_only_no_confirmation_claim" else "separate_prespecified_screens_no_universal_validity"),
    file.path(out, "adoption.json"), auto_unbox = TRUE, pretty = TRUE, na = "null", digits = 16)
}
utils::write.csv(records[records$status != "ok", c("cell", "replicate", "method", "fit_id", "failure_stage", "error_class", "message", "warnings")],
  file.path(out, "errors.csv"), row.names = FALSE, na = "")
after <- hash_sources()
utils::write.csv(after, file.path(out, "source-hashes-after.csv"), row.names = FALSE)
if (!identical(before, after)) stop("Sources changed during stage; no COMPLETE marker is written.")
jsonlite::write_json(list(status = "complete", mode = mode, cohorts = cohorts, analysis_requests = nrow(records),
  fitted_workflow_requests = nrow(fits), mcmc_fit_requests = sum(fits$fit_id != "oracle"),
  returned_intervals = sum(records$status == "ok"), analysis_failures = sum(records$status != "ok"),
  fit_failures = sum(fits$fit_status != "ok"), source_fingerprint = source_fingerprint,
  started_utc = format(started, tz = "UTC", usetz = TRUE), finished_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  elapsed_seconds = as.numeric(difftime(Sys.time(), started, units = "secs"))),
  file.path(out, "COMPLETE.json"), auto_unbox = TRUE, pretty = TRUE)
cat(sprintf("Complete: %d analysis requests, %d intervals, %d analysis refusals; %d fit refusals.\n",
  nrow(records), sum(records$status == "ok"), sum(records$status != "ok"), sum(fits$fit_status != "ok")))
