#!/usr/bin/env Rscript
# Frozen stages and estimands: validation/MARGINAL_EXTENSION_PLAN.md.
args <- commandArgs(trailingOnly = TRUE)
mode <- if (length(args)) match.arg(args[1L], c("pilot", "development", "confirmation")) else "pilot"
script <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1L])
root <- normalizePath(file.path(dirname(script), ".."), mustWork = TRUE)
out <- if (length(args) >= 2L) args[2L] else file.path(root, "validation", paste0("joint-marginal-", mode))
if (dir.exists(out) && length(list.files(out, all.files = TRUE, no.. = TRUE)))
  stop("Use a new or empty output directory; evidence is not overwritten.")
workers <- suppressWarnings(as.integer(Sys.getenv("CATMECH_BENCHMARK_WORKERS", "8")))
if (length(workers) != 1L || is.na(workers) || workers < 1L || workers > 16L) stop("Invalid worker count.")
if (.Platform$OS.type == "windows") workers <- 1L
baseline_root <- file.path(root, "validation/v0.9.0-source")
hash_sources <- function() {
  rel <- c("DESCRIPTION", "NAMESPACE", "validation/MARGINAL_EXTENSION_PLAN.md",
    "validation/run-marginal-extension.R",
    paste0("src/", sort(list.files(file.path(root, "src"), "\\.(c|cpp|h)$|^Makevars"))),
    paste0("R/", sort(list.files(file.path(root, "R"), "\\.R$"))),
    "validation/v0.9.0-source/DESCRIPTION", "validation/v0.9.0-source/NAMESPACE",
    paste0("validation/v0.9.0-source/R/", sort(list.files(file.path(baseline_root, "R"), "\\.R$"))))
  data.frame(file = rel, bytes = as.numeric(file.info(file.path(root, rel))$size),
    sha256 = vapply(file.path(root, rel), function(f) digest::digest(file = f, algo = "sha256"), character(1)),
    row.names = NULL)
}
before <- hash_sources()
source_fingerprint <- digest::digest(before, algo = "sha256")
pkgload::load_all(root, quiet = TRUE)
if (!exists("catmech_marginal_summary", mode = "function")) stop("Candidate posterior-summary API is unavailable.")
old <- new.env(parent = globalenv())
for (f in sort(list.files(file.path(baseline_root, "R"), "\\.R$", full.names = TRUE))) sys.source(f, envir = old)
dir.create(out, recursive = TRUE, showWarnings = FALSE)
cells <- c("gaussian_balanced", "gaussian_unbalanced", "near_zero_icc", "observed_mar", "null_slope", "student_t5", "contextual_x", "contextual_unbalanced", "contextual_null")
prior_cells <- c("gaussian_balanced", "near_zero_icc", "observed_mar", "contextual_x")
contextual_cells <- c("contextual_x", "contextual_unbalanced", "contextual_null")
primary_cells <- c("gaussian_balanced", "gaussian_unbalanced", "observed_mar", "null_slope")
new_control <- list(chains = 4L, burnin = 1000L, iterations = 4000L, max_iterations = 64000L, prior_scale = 1, regression_scale = 10)
old_control <- list(chains = 4L, burnin = 1000L, iterations = 4000L, max_iterations = 32000L, prior_scale = 1)
compiler_fields <- c("CC", "CFLAGS", "CXX17", "CXX17STD", "CXX17FLAGS", "CPPFLAGS", "LDFLAGS", "FLIBS")
compiler <- setNames(lapply(compiler_fields, function(field)
  paste(system2(file.path(R.home("bin"), "R"), c("CMD", "config", field),
    stdout = TRUE, stderr = TRUE), collapse = "\n")), compiler_fields)
makevars_user <- Sys.getenv("R_MAKEVARS_USER", unset = "")
compiler$R_MAKEVARS_USER <- makevars_user
compiler$makevars_user_sha256 <- if (nzchar(makevars_user) && file.exists(makevars_user))
  digest::digest(file = makevars_user, algo = "sha256") else NA_character_
config <- list(mode = mode, workers = workers, m = 40L, cells = cells,
  primary_cells = primary_cells, contextual_cells = contextual_cells, prior_cells = prior_cells,
  primary_reps = switch(mode, pilot = 1L, development = 30L, confirmation = 400L),
  prior_reps = switch(mode, pilot = 1L, development = 30L, confirmation = 200L),
  data_seed_base = switch(mode, pilot = 1670000000L, development = 1750000000L, confirmation = 1710000000L),
  imputation_seed_base = switch(mode, pilot = 1680000000L, development = 1760000000L, confirmation = 1720000000L),
  bootstrap_seed = switch(mode, pilot = NA_integer_, development = 1770000000L, confirmation = 1740000000L),
  bootstrap_reps = 2000L, new_control = new_control, old_control = old_control,
  baseline_source_commit = "ae8f94086710ea1acc0df937f45ee847b6a0976b", pan_version = as.character(utils::packageVersion("pan")),
  candidate_sampler = "factorized_hw_matrixnormal_gibbs_logscale_mh_2_per_component_sd1",
  source_fingerprint = source_fingerprint, compiler = compiler, accuracy_exported = mode != "pilot")
config$validation_kind <- "dependent_numerical_replay"
config$replays_source_fingerprint <- "607f0585be72b5e116a6d384f1c1df2f1294241737540e80deaf4aa413d2814e"
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
  sizes <- if (cell %in% c("gaussian_unbalanced", "contextual_unbalanced")) rep(c(5L, 8L, 12L, 15L), 6L) else rep(10L, 24L)
  ids <- rep(seq_len(24L), sizes); n <- length(ids)
  if (cell %in% contextual_cells) {
    a <- stats::rnorm(24L)
    x <- a[ids] + stats::rnorm(n)
  } else { a <- rep(0, 24L); x <- stats::rnorm(n) }
  Omega <- matrix(c(.5, .45, .45, .6), 2L, 2L)
  if (cell == "near_zero_icc") Omega <- .02 * Omega
  Sigma <- matrix(c(.5, .55, .55, .8), 2L, 2L)
  u <- rmv(24L, Omega, cell == "student_t5")
  e <- rmv(n, Sigma, cell == "student_t5")
  structural <- if (cell == "null_slope") 0 else if (cell == "contextual_null") -.2 else .7
  target <- if (cell == "contextual_null") 0 else if (cell %in% contextual_cells) .9 else structural
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
  c("baseline", "marginal", "oracle", if (cell %in% prior_cells && replicate <= config$prior_reps)
    c("baseline_prior_quarter", "baseline_prior_four", "marginal_prior_quarter", "marginal_prior_four",
      "marginal_regression_half", "marginal_regression_double"))
}
prior_scale <- function(id) if (endsWith(id, "_quarter")) .25 else if (endsWith(id, "_four")) 4 else 1
regression_scale <- function(id) if (!startsWith(id, "marginal")) NA_real_ else
  if (endsWith(id, "_half")) 5 else if (endsWith(id, "_double")) 20 else 10
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
  fun <- if (startsWith(id, "baseline") || id == "oracle") old$catmech_reference else catmech_reference
  if (id == "oracle") {
    call$data$gold <- p$truth; call$m <- 2L
  } else {
    call$cluster <- "site"
    call$cluster_model <- if (startsWith(id, "baseline")) "joint" else "joint_marginal"
    call$joint_control <- if (startsWith(id, "baseline")) old_control else new_control
    call$joint_control$prior_scale <- prior_scale(id)
    if (startsWith(id, "marginal")) call$joint_control$regression_scale <- regression_scale(id)
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
    fun <- if (startsWith(id, "baseline")) old$catmech_joint_summary else catmech_marginal_summary
    z <- fun(fit, level = .95)
    z <- z[z$term == "x", , drop = FALSE]
    if (nrow(z) != 1L) stop("Posterior coefficient x is absent.")
    ans <- list(estimate = z$estimate, spread = z$posterior_sd, lower = z$lower, upper = z$upper,
      df = NA_real_, within_variance = NA_real_, between_variance = NA_real_,
      posterior_rhat = z$rhat, posterior_ess_bulk = z$ess_bulk, posterior_ess_tail = z$ess_tail,
      mcse_mean = z$mcse_mean, mcse_lower = z$mcse_lower, mcse_upper = z$mcse_upper,
      posterior_draws = z$draws, interpretation = z$interpretation)
  } else {
    fun <- if (startsWith(id, "baseline") || id == "oracle") old$catmech_pool_lm else catmech_pool_lm
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
      prior_scale = prior_scale(id), regression_scale = regression_scale(id),
      fit_status = if (fit_failed) "error" else "ok", error_class = if (fit_failed) class(fit$value)[1L] else "",
      message = if (fit_failed) conditionMessage(fit$value) else "", warnings = fit$warnings,
      elapsed_seconds = fit$elapsed, max_rhat = d$max_rhat, min_ess_bulk = d$min_bulk,
      min_ess_tail = d$min_tail, iterations_used = d$iterations, attempts = d$attempts,
      source_fingerprint = source_fingerprint, stringsAsFactors = FALSE))
    fits[[id]] <- fit_row
    if (!is.null(d$summary) && nrow(d$summary)) diagnostics[[id]] <- cbind(
      cell = cell, replicate = replicate, fit_id = id, d$summary)
    analyses <- if (id == "oracle") "cv3" else c("cv3", "direct")
    for (analysis in analyses) {
      family <- if (startsWith(id, "baseline")) "baseline" else if (startsWith(id, "marginal")) "marginal" else "oracle"
      method <- sub(family, paste0(family, "_", analysis), id, fixed = TRUE)
      analyzed <- if (fit_failed) list(value = fit$value, elapsed = 0, warnings = "") else
        capture(function() analyze_one(fit$value, id, analysis))
      failed <- inherits(analyzed$value, "error")
      row <- cbind(p$meta, data.frame(method = method, fit_id = id, analysis = analysis,
        prior_scale = fit_row$prior_scale, regression_scale = fit_row$regression_scale, status = if (failed) "error" else "ok",
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
    comparisons <- list(c("marginal_direct", "baseline_direct"),
      c("marginal_cv3", "baseline_cv3"), c("marginal_direct", "marginal_cv3"),
      c("baseline_direct", "baseline_cv3"))
    if (cell %in% prior_cells) for (analysis in c("baseline_direct", "baseline_cv3", "marginal_direct", "marginal_cv3"))
      for (suffix in c("_prior_quarter", "_prior_four"))
        comparisons[[length(comparisons) + 1L]] <- c(paste0(analysis, suffix), analysis)
    if (cell %in% prior_cells) for (analysis in c("marginal_direct", "marginal_cv3"))
      for (suffix in c("_regression_half", "_regression_double"))
        comparisons[[length(comparisons) + 1L]] <- c(paste0(analysis, suffix), analysis)
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
    for (replicate in sort(unique(z$replicate))) for (analysis in c("direct", "cv3")) {
      r <- z[z$replicate == replicate, , drop = FALSE]
      row <- data.frame(cell = cell, replicate = replicate, analysis = analysis,
        baseline_complete = FALSE, marginal_complete = FALSE,
        marginal_hw_complete = FALSE, marginal_regression_complete = FALSE,
        baseline_coefficient_sensitivity = NA_real_, marginal_coefficient_sensitivity = NA_real_,
        baseline_width_sensitivity = NA_real_, marginal_width_sensitivity = NA_real_,
        marginal_hw_coefficient_sensitivity = NA_real_, marginal_regression_coefficient_sensitivity = NA_real_,
        marginal_hw_width_sensitivity = NA_real_, marginal_regression_width_sensitivity = NA_real_,
        baseline_quarter_shift = NA_real_, baseline_four_shift = NA_real_,
        marginal_quarter_shift = NA_real_, marginal_four_shift = NA_real_,
        marginal_regression_half_shift = NA_real_, marginal_regression_double_shift = NA_real_,
        coefficient_difference = NA_real_, width_difference = NA_real_)
      for (family in c("baseline", "marginal")) {
        methods <- paste0(family, "_", analysis, c("", "_prior_quarter", "_prior_four"))
        selected <- r[match(methods, r$method), ]
        ok <- nrow(selected) == 3L && all(!is.na(selected$status) & selected$status == "ok")
        row[[if (family == "baseline") "baseline_complete" else "marginal_hw_complete"]] <- ok
        if (ok) {
          shifts <- selected$estimate[-1L] - selected$estimate[1L]
          prefix <- if (family == "baseline") "baseline" else "marginal_hw"
          row[[paste0(prefix, "_coefficient_sensitivity")]] <- mean(abs(shifts))
          row[[paste0(prefix, "_width_sensitivity")]] <- mean(abs(log(selected$width[-1L] / selected$width[1L])))
          row[[paste0(family, "_quarter_shift")]] <- shifts[1L]
          row[[paste0(family, "_four_shift")]] <- shifts[2L]
        }
      }
      methods <- paste0("marginal_", analysis, c("", "_regression_half", "_regression_double"))
      selected <- r[match(methods, r$method), ]
      row$marginal_regression_complete <- nrow(selected) == 3L && all(!is.na(selected$status) & selected$status == "ok")
      if (row$marginal_regression_complete) {
        shifts <- selected$estimate[-1L] - selected$estimate[1L]
        row$marginal_regression_coefficient_sensitivity <- mean(abs(shifts))
        row$marginal_regression_width_sensitivity <- mean(abs(log(selected$width[-1L] / selected$width[1L])))
        row$marginal_regression_half_shift <- shifts[1L]
        row$marginal_regression_double_shift <- shifts[2L]
      }
      row$marginal_complete <- row$marginal_hw_complete && row$marginal_regression_complete
      if (row$marginal_complete) {
        row$marginal_coefficient_sensitivity <- max(row$marginal_hw_coefficient_sensitivity, row$marginal_regression_coefficient_sensitivity)
        row$marginal_width_sensitivity <- max(row$marginal_hw_width_sensitivity, row$marginal_regression_width_sensitivity)
      }
      if (row$baseline_complete && row$marginal_complete) {
        row$coefficient_difference <- row$marginal_coefficient_sensitivity - row$baseline_coefficient_sensitivity
        row$width_difference <- row$marginal_width_sensitivity - row$baseline_width_sensitivity
      }
      triples[[length(triples) + 1L]] <- row
    }
  }
  triples <- do.call(rbind, triples)
  for (cell in prior_cells) for (analysis in c("direct", "cv3")) {
    z <- triples[triples$cell == cell & triples$analysis == analysis, ]
    good <- z[z$baseline_complete & z$marginal_complete, ]
    row <- data.frame(cell = cell, analysis = analysis, requested_triples = nrow(z),
      baseline_triples = sum(z$baseline_complete), marginal_triples = sum(z$marginal_complete),
      marginal_hw_triples = sum(z$marginal_hw_complete), marginal_regression_triples = sum(z$marginal_regression_complete),
      complete_pairs = nrow(good))
    for (metric in c("coefficient", "width")) {
      for (prefix in c("baseline", "marginal", "marginal_hw", "marginal_regression")) {
        field <- paste0(prefix, "_", metric, "_sensitivity")
        row[[paste0("mean_", field)]] <- if (nrow(good)) mean(good[[field]]) else NA_real_
      }
      field <- paste0(metric, "_difference")
      row[[field]] <- if (nrow(good)) mean(good[[field]]) else NA_real_
      row[[paste0(metric, "_mcse")]] <- se_mean(good[[field]])
    }
    summary[[length(summary) + 1L]] <- row
  }
  summary <- do.call(rbind, summary); rownames(summary) <- NULL
  direct <- summary[summary$analysis == "direct", ]
  overall <- mean(direct$coefficient_difference)
  mcse <- sqrt(sum(direct$coefficient_mcse^2)) / length(prior_cells)
  width <- mean(direct$width_difference)
  width_mcse <- sqrt(sum(direct$width_mcse^2)) / length(prior_cells)
  gate <- list(difference = overall, mcse = mcse,
    lower = overall - stats::qnorm(.975) * mcse, upper = overall + stats::qnorm(.975) * mcse,
    no_positive_cell = all(is.finite(direct$coefficient_difference) & direct$coefficient_difference <= 0),
    pass = is.finite(overall) && is.finite(mcse) && overall + stats::qnorm(.975) * mcse < 0 &&
      all(is.finite(direct$coefficient_difference) & direct$coefficient_difference <= 0),
    sensitivity_definition = "candidate_max_HW_and_regression_axes_vs_baseline_IW_on_common_3_plus_5_fits",
    secondary_width_difference = width, secondary_width_mcse = width_mcse,
    secondary_width_lower = width - stats::qnorm(.975) * width_mcse,
    secondary_width_upper = width + stats::qnorm(.975) * width_mcse)
  list(triples = triples, summary = summary, gate = gate)
}

calibration_screen <- function(records) {
  methods <- c("baseline_direct", "marginal_direct")
  panels <- lapply(cells, function(cell) {
    z <- records[records$cell == cell, ]; reps <- sort(unique(z$replicate))
    out <- list()
    for (metric in c("covered", "error")) {
      a <- matrix(NA_real_, length(reps), length(methods), dimnames = list(reps, methods))
      for (method in methods) {
        r <- z[z$method == method, ]; r <- r[match(reps, r$replicate), ]
        a[, method] <- ifelse(r$status == "ok", r[[metric]], NA_real_)
      }
      out[[metric]] <- a
    }
    out
  })
  names(panels) <- cells
  conditional <- vapply(panels, function(z) colMeans(z$covered, na.rm = TRUE), numeric(length(methods)))
  biases <- vapply(panels, function(z) colMeans(z$error, na.rm = TRUE), numeric(length(methods)))
  primary_distance <- rowMeans(abs(conditional[, primary_cells, drop = FALSE] - .95))
  contextual_bias <- rowMeans(abs(biases[, contextual_cells, drop = FALSE]))
  RNGkind("Mersenne-Twister", "Inversion", "Rejection"); set.seed(config$bootstrap_seed)
  bootstrap <- matrix(NA_real_, config$bootstrap_reps, 2L,
    dimnames = list(NULL, c("contextual_absolute_bias_difference", "primary_coverage_distance_difference")))
  for (b in seq_len(config$bootstrap_reps)) {
    coverage <- bias <- matrix(NA_real_, length(methods), length(cells), dimnames = list(methods, cells))
    for (cell in cells) {
      z <- panels[[cell]]; chosen <- sample.int(nrow(z$covered), nrow(z$covered), replace = TRUE)
      coverage[, cell] <- colMeans(z$covered[chosen, , drop = FALSE], na.rm = TRUE)
      bias[, cell] <- colMeans(z$error[chosen, , drop = FALSE], na.rm = TRUE)
    }
    bias_distance <- rowMeans(abs(bias[, contextual_cells, drop = FALSE]))
    coverage_distance <- rowMeans(abs(coverage[, primary_cells, drop = FALSE] - .95))
    bootstrap[b, ] <- c(bias_distance[2L] - bias_distance[1L], coverage_distance[2L] - coverage_distance[1L])
  }
  bounds <- apply(bootstrap, 2L, function(v) if (all(is.finite(v)))
    stats::quantile(v, probs = c(.025, .975), type = 7, names = FALSE) else c(NA_real_, NA_real_))
  contextual_pairs <- lapply(contextual_cells, function(cell) {
    z <- records[records$cell == cell, ]; a <- z[z$method == "marginal_direct", ]; b <- z[z$method == "baseline_direct", ]
    pair <- merge(a, b, by = "replicate", suffixes = c("_marginal", "_baseline"))
    delta <- pair$covered_marginal - pair$covered_baseline
    data.frame(cell = cell, pairs = nrow(pair), difference = mean(delta), mcse = se_mean(delta))
  })
  contextual_pairs <- do.call(rbind, contextual_pairs)
  delta <- mean(contextual_pairs$difference)
  delta_mcse <- sqrt(sum(contextual_pairs$mcse^2)) / length(contextual_cells)
  band <- function(group) all(is.finite(conditional[2L, group]) & conditional[2L, group] >= .925 & conditional[2L, group] <= .975)
  bias_pass <- function(group) all(is.finite(biases[2L, group]) & abs(biases[2L, group]) <= .02)
  contextual <- list(coverage_band_pass = band(contextual_cells), bias_bound_pass = bias_pass(contextual_cells),
    coverage_difference = delta, coverage_mcse = delta_mcse,
    coverage_lower = delta - stats::qnorm(.975) * delta_mcse,
    coverage_upper = delta + stats::qnorm(.975) * delta_mcse,
    no_negative_cell = all(is.finite(contextual_pairs$difference) & contextual_pairs$difference >= 0),
    mean_absolute_bias = as.list(contextual_bias),
    absolute_bias_difference = contextual_bias[2L] - contextual_bias[1L],
    absolute_bias_lower = bounds[1L, 1L], absolute_bias_upper = bounds[2L, 1L])
  contextual$pass <- contextual$coverage_band_pass && contextual$bias_bound_pass &&
    contextual$no_negative_cell && is.finite(contextual$coverage_lower) && contextual$coverage_lower > 0 &&
    is.finite(contextual$absolute_bias_difference) && contextual$absolute_bias_difference < 0
  ordinary <- list(coverage_band_pass = band(primary_cells), bias_bound_pass = bias_pass(primary_cells),
    mean_absolute_distance = as.list(primary_distance),
    distance_difference = primary_distance[2L] - primary_distance[1L],
    distance_lower = bounds[1L, 2L], distance_upper = bounds[2L, 2L],
    pass = band(primary_cells) && bias_pass(primary_cells))
  list(bootstrap = data.frame(bootstrap = seq_len(config$bootstrap_reps), bootstrap),
    contextual_pairs = contextual_pairs, contextual = contextual, ordinary = ordinary,
    conditional = setNames(lapply(seq_along(methods), function(k) setNames(as.list(conditional[k, ]), cells)), methods),
    bias = setNames(lapply(seq_along(methods), function(k) setNames(as.list(biases[k, ]), cells)), methods))
}

availability_screen <- function(fits, records) {
  rows <- lapply(cells, function(cell) {
    z <- fits[fits$cell == cell, ]; a <- z[z$fit_id == "marginal", ]; b <- z[z$fit_id == "baseline", ]
    pair <- merge(a, b, by = "replicate", suffixes = c("_marginal", "_baseline"))
    candidate_fail <- pair$fit_status_marginal != "ok"; baseline_fail <- pair$fit_status_baseline != "ok"
    delta <- as.numeric(candidate_fail) - as.numeric(baseline_fail)
    candidate_records <- records[records$cell == cell & records$method == "marginal_direct", ]
    data.frame(cell = cell, requests = nrow(pair), baseline_failures = sum(baseline_fail),
      marginal_failures = sum(candidate_fail), baseline_failure_rate = mean(baseline_fail),
      marginal_failure_rate = mean(candidate_fail), marginal_interval_failures = sum(candidate_records$status != "ok"),
      marginal_interval_failure_rate = mean(candidate_records$status != "ok"), difference = mean(delta), mcse = se_mean(delta))
  })
  rows <- do.call(rbind, rows)
  delta <- mean(rows$difference); mcse <- sqrt(sum(rows$mcse^2)) / nrow(rows)
  fit_rate <- all(rows$marginal_failure_rate <= .01)
  interval_rate <- all(rows$marginal_interval_failure_rate <= .01)
  no_increase <- all(rows$marginal_failures <= rows$baseline_failures)
  list(cells = rows, gate = list(difference = delta, mcse = mcse,
    lower = delta - stats::qnorm(.975) * mcse, upper = delta + stats::qnorm(.975) * mcse,
    every_fit_return_rate_at_least_99 = fit_rate, every_interval_return_rate_at_least_99 = interval_rate,
    no_observed_fit_failure_increase = no_increase, pass = fit_rate && interval_rate && no_increase))
}

all_records <- all_fits <- all_diagnostics <- list()
started <- Sys.time()
for (index in seq_along(cells)) {
  cell <- cells[index]
  reps <- config$primary_reps
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
cohorts <- length(cells) * config$primary_reps
expected_analyses <- 5L * cohorts + 12L * length(prior_cells) * config$prior_reps
expected_fits <- 3L * cohorts + 6L * length(prior_cells) * config$prior_reps
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
  utils::write.csv(calibration$contextual_pairs, file.path(out, "contextual-paired.csv"), row.names = FALSE, na = "")
  availability <- availability_screen(fits, records)
  utils::write.csv(availability$cells, file.path(out, "availability.csv"), row.names = FALSE, na = "")
  jsonlite::write_json(list(stage = mode, contextual = calibration$contextual, ordinary = calibration$ordinary,
    coefficient_prior_sensitivity = prior$gate, availability = availability$gate,
    conditional_coverage = calibration$conditional, bias = calibration$bias,
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
