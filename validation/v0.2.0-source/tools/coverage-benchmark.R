#!/usr/bin/env Rscript
# Prespecified held-out study; see validation/COVERAGE_PLAN.md.
args <- commandArgs(trailingOnly = TRUE)
mode <- if (length(args)) match.arg(args[1], c("pilot", "confirmatory", "clinical")) else "pilot"
file_arg <- grep("^--file=", commandArgs(), value = TRUE)
root <- normalizePath(file.path(dirname(sub("^--file=", "", file_arg[1])), ".."))
outdir <- if (length(args) >= 2L) args[2] else file.path(root, "validation", paste0("coverage-", mode))
if (dir.exists(outdir) && length(list.files(outdir, all.files = TRUE, no.. = TRUE))) {
  stop("Use a new or empty output directory; existing evidence is never overwritten.")
}
hash_sources <- function() {
  paths <- c(sort(list.files(file.path(root, "R"), "\\.R$", full.names = TRUE)),
    file.path(root, c("DESCRIPTION", "NAMESPACE", "tools/coverage-benchmark.R",
                      "tools/lab-fixture.R", "validation/COVERAGE_PLAN.md")))
  data.frame(file = substring(paths, nchar(root) + 2L),
    sha256 = vapply(paths, function(p) digest::digest(file = p, algo = "sha256"), character(1)),
    row.names = NULL)
}
hashes_before <- hash_sources()
pkgload::load_all(root, quiet = TRUE)
source(file.path(root, "tools", "lab-fixture.R"), local = TRUE)
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
workers <- suppressWarnings(as.integer(Sys.getenv("CATMECH_BENCHMARK_WORKERS", "4")))
if (length(workers) != 1L || is.na(workers) || workers < 1L || workers > 16L) stop("Invalid CATMECH_BENCHMARK_WORKERS (1 to 16).")
if (.Platform$OS.type == "windows") workers <- 1L
config <- list(mode = mode, workers = workers, m = 20L, chains = 4L, burnin = 200L, iterations = 400L,
  primary_reps = if (mode == "pilot") 5L else 1000L,
  gibbs_reps = if (mode == "pilot") 2L else 250L,
  clinical_reps = if (mode == "pilot") 3L else 250L,
  synthetic_n = 300L, clinical_n = 600L, centered_bootstrap_reps = 400L)
jsonlite::write_json(config, file.path(outdir, "settings.json"), pretty = TRUE, auto_unbox = TRUE)
synthetic <- data.frame(
  scenario = c("missing_only", "perfect_flags", "false_positive_flags", "category_error", "null_effect", "heavy_tails"),
  missing_rate = rep(.2, 6), flag_rate = c(0, rep(.1, 5)),
  corruption_probability = c(1, 1, .5, 1, .5, 1),
  category_error_rate = c(0, 0, 0, .1, 0, 0),
  error_distribution = c(rep("normal", 5), "heavy_tailed"),
  truth = c(rep(.8, 4), 0, .8), family = "synthetic", stringsAsFactors = FALSE)
clinical <- data.frame(scenario = c("clean", "flagged_errors", "missed_flags", "site_shift", "rounded_records"),
                       family = "clinical", stringsAsFactors = FALSE)
scenario_list <- list()
if (mode != "clinical") scenario_list <- lapply(seq_len(nrow(synthetic)), function(i) as.list(synthetic[i, ]))
if (mode != "confirmatory") scenario_list <- c(scenario_list, lapply(seq_len(nrow(clinical)), function(i) as.list(clinical[i, ])))
thresholds <- c(-Inf, 5.5, 7, Inf)
category_levels <- c("Low", "Medium", "High")
base_methods <- c("oracle", "naive_observed", "trusted_observed", "trusted_pmm", "trusted_normal", "model_ignore")

# All scientific analyses use the same precomputed covariate design. lm.fit is
# equivalent to lm here and avoids rebuilding a model.frame for every imputation.
fit_lab <- function(lab, prepared) {
  use <- is.finite(lab)
  x <- cbind(prepared$analysis_x[use, , drop = FALSE], lab_value = lab[use])
  y <- prepared$input$outcome[use]
  fit <- stats::lm.fit(x, y)
  p <- ncol(x); n <- nrow(x)
  if (fit$rank != p || n <= p) stop("Scientific analysis is rank deficient or has no residual degrees of freedom.")
  position <- match(p, fit$qr$pivot)
  variance_factor <- chol2inv(qr.R(fit$qr))[position, position]
  se <- sqrt(sum(fit$residuals^2) / (n - p) * variance_factor)
  list(estimate = unname(fit$coefficients[p]), se = se, df = n - p, n_used = n)
}
pool_fits <- function(fits) {
  q <- vapply(fits, function(x) x$estimate, numeric(1))
  u <- vapply(fits, function(x) x$se^2, numeric(1))
  if (length(unique(vapply(fits, function(x) x$n_used, numeric(1)))) != 1L) stop("Imputations have different analysis sample sizes.")
  z <- mice::pool.scalar(q, u, n = fits[[1]]$n_used, k = fits[[1]]$n_used - fits[[1]]$df)
  list(estimate = z$qbar, se = sqrt(z$t), df = z$df, n_used = fits[[1]]$n_used)
}
prepare_data <- function(scenario, scenario_index, replicate) {
  if (scenario$family == "synthetic") {
    seed <- as.integer((if (mode == "pilot") 7000000 else 20000000) + 100000 * scenario_index + replicate)
    d <- simulate_catmech_lab(n = config$synthetic_n, seed = seed, thresholds = thresholds,
      missing_rate = scenario$missing_rate, flag_rate = scenario$flag_rate,
      corruption_probability = scenario$corruption_probability,
      category_error_rate = scenario$category_error_rate, error_distribution = scenario$error_distribution)
    if (scenario$scenario == "null_effect") {
      set.seed(seed + 125000L)
      d$outcome <- stats::rnorm(nrow(d))
    }
    predictors <- c("age", "sex", "inflammation", "outcome")
    input <- d[, c(predictors, "lab_value", "lab_category", "outlier_flag")]
    actual_truth <- d$lab_true
    truth <- scenario$truth
    stream <- as.character(seed)
  } else {
    f <- make_lab_fixture(n = config$clinical_n, replicate = replicate,
      phase = if (mode == "pilot") "development" else "heldout", scenario = scenario$scenario)
    check_lab_fixture(f)
    predictors <- c("age", "sex", "inflammation", "SITEID", "outcome")
    input <- f$data[, predictors]
    input$lab_value <- f$data$LBSTRESN
    input$lab_category <- f$data$LBNRIND
    input$outlier_flag <- f$data$QCFLAG
    actual_truth <- f$truth$LAB_TRUE
    truth <- f$manifest$known_slope
    seed <- as.integer(30000000 + 100000 * scenario_index + replicate)
    stream <- paste0(if (mode == "pilot") "development" else "heldout", ":", replicate)
  }
  analysis_predictors <- setdiff(predictors, "outcome")
  analysis_x <- stats::model.matrix(stats::reformulate(analysis_predictors), data = input)
  model_x <- stats::model.matrix(stats::reformulate(predictors), data = input)[, -1L, drop = FALSE]
  pmm_x <- stats::model.matrix(stats::reformulate(c(predictors, "lab_category")), data = input)[, -1L, drop = FALSE]
  flagged <- as.character(input$outlier_flag) == "Yes"
  list(input = input, actual_truth = actual_truth, truth = truth, predictors = predictors,
       analysis_x = analysis_x, model_x = model_x, pmm_x = pmm_x, flagged = flagged,
       seed = seed, stream = stream)
}
extract_diagnostics <- function(result) {
  z <- result$model_diagnostics$summary
  if (is.null(z) || !is.data.frame(z) || !nrow(z)) return(c(rhat_max = NA_real_, ess_min = NA_real_, ess_tail_min = NA_real_))
  r <- if ("rhat" %in% names(z)) z$rhat else numeric()
  e <- if ("ess_bulk" %in% names(z)) z$ess_bulk else numeric()
  t <- if ("ess_tail" %in% names(z)) z$ess_tail else numeric()
  c(rhat_max = if (length(r) && any(!is.na(r))) max(r, na.rm = TRUE) else NA_real_,
    ess_min = if (length(e) && any(!is.na(e))) min(e, na.rm = TRUE) else NA_real_,
    ess_tail_min = if (length(t) && any(!is.na(t))) min(t, na.rm = TRUE) else NA_real_)
}
run_method <- function(method, p, scenario) {
  y <- p$input$lab_value
  if (method == "oracle") return(fit_lab(p$actual_truth, p))
  if (method == "naive_observed") return(fit_lab(y, p))
  if (method == "trusted_observed") {
    y[p$flagged] <- NA_real_
    return(fit_lab(y, p))
  }
  seed <- p$seed + 1000000L
  if (method %in% c("trusted_pmm", "trusted_normal")) {
    y[p$flagged] <- NA_real_; ry <- !is.na(y); wy <- !ry
    set.seed(seed)
    fits <- lapply(seq_len(config$m), function(j) {
      completed <- y
      if (any(wy)) completed[wy] <- if (method == "trusted_pmm") {
        mice::mice.impute.pmm(y, ry, p$pmm_x, wy = wy, donors = 5L, matchtype = 1L)
      } else mice::mice.impute.norm(y, ry, p$model_x, wy = wy)
      fit_lab(completed, p)
    })
    return(pool_fits(fits))
  }
  if (method == "legacy") {
    s <- catmech_spec("lab_value", "lab_category", category_levels, thresholds,
      flag = "outlier_flag", corruption_probability = scenario$corruption_probability,
      conflict_policy = "reconstruct", min_donors = 1L, n_donors = 5L, fallback = "error")
    imp <- catmech_impute(p$input, s, p$predictors, m = config$m, seed = seed)
  } else {
    category_mode <- switch(method, model_ignore = "ignore", model_exact = "exact",
                             model_noisy = "noisy", model_wrong_noise = "noisy")
    model_args <- list(data = p$input, value = "lab_value", predictors = p$predictors,
      flag = "outlier_flag", flag_yes = "Yes", flag_no = "No", category_mode = category_mode,
      m = config$m, seed = seed, burnin = config$burnin,
      iterations = config$iterations, chains = config$chains)
    if (category_mode != "ignore") {
      model_args$category <- "lab_category"; model_args$thresholds <- thresholds
      model_args$levels <- category_levels
    }
    if (category_mode == "noisy") {
      error <- if (method == "model_wrong_noise") .02 else .10
      H <- matrix(error / 2, 3, 3, dimnames = list(category_levels, category_levels)); diag(H) <- 1 - error
      model_args$confusion <- H
    }
    imp <- do.call(catmech_model, model_args)
  }
  result <- pool_fits(lapply(imp$completed, function(d) fit_lab(d$lab_value, p)))
  result$model_diagnostics <- imp$model_diagnostics
  result
}
methods_for <- function(scenario, replicate) {
  methods <- base_methods
  if (scenario$family == "synthetic") {
    methods <- c(methods, "legacy")
    if (replicate <= config$gibbs_reps || scenario$scenario == "category_error") methods <- c(methods, "model_exact")
    if (scenario$scenario == "category_error" && replicate <= config$gibbs_reps) {
      methods <- c(methods, "model_noisy", "model_wrong_noise")
    }
  }
  methods
}

run_replicate <- function(s, r) {
  scenario <- scenario_list[[s]]
  records <- list(); decompositions <- list(); idx <- 0L; didx <- 0L
    p <- prepare_data(scenario, s, r)
    keep <- !is.na(p$input$lab_value) & !p$flagged
    oracle <- fit_lab(p$actual_truth, p)
    kept_truth <- p$actual_truth; kept_truth[!keep] <- NA_real_
    kept_recorded <- p$input$lab_value; kept_recorded[!keep] <- NA_real_
    oracle_kept <- fit_lab(kept_truth, p); recorded_kept <- fit_lab(kept_recorded, p)
    didx <- didx + 1L
    decompositions[[didx]] <- data.frame(family = scenario$family, scenario = scenario$scenario, replicate = r,
      oracle_full = oracle$estimate, oracle_retained = oracle_kept$estimate, recorded_retained = recorded_kept$estimate,
      selection_difference = oracle_kept$estimate - oracle$estimate,
      residual_error_difference = recorded_kept$estimate - oracle_kept$estimate,
      retained_n = sum(keep), truth = p$truth)
    for (method in methods_for(scenario, r)) {
      warnings <- character()
      result <- tryCatch(withCallingHandlers(run_method(method, p, scenario),
        warning = function(w) { warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning") }),
        error = function(e) e)
      failed <- inherits(result, "error")
      status <- if (failed) "error" else "ok"
      message <- if (failed) conditionMessage(result) else ""
      if (inherits(result, "catmech_convergence_error")) status <- "convergence_failure"
      if (failed && method == "model_exact" && grepl("conflict|contradic|incompatible", message, ignore.case = TRUE)) status <- "refused"
      diagnostics <- extract_diagnostics(result)
      row <- data.frame(family = scenario$family, scenario = scenario$scenario, replicate = r,
        data_seed = p$seed, data_stream = p$stream, method = method, truth = p$truth,
        estimate = NA_real_, se = NA_real_, df = NA_real_, n_used = NA_real_,
        lower = NA_real_, upper = NA_real_, width = NA_real_, coverage = NA_real_,
        error = NA_real_, rhat_max = diagnostics["rhat_max"], ess_min = diagnostics["ess_min"], ess_tail_min = diagnostics["ess_tail_min"],
        status = status, message = message, warnings = paste(unique(warnings), collapse = " | "),
        stringsAsFactors = FALSE, row.names = NULL)
      if (!failed) {
        if (!is.finite(result$estimate) || !is.finite(result$se) || result$se < 0 || is.na(result$df) || result$df <= 0) {
          row$status <- "error"; row$message <- "Nonfinite analysis estimate or invalid uncertainty."
        } else {
          row$estimate <- result$estimate; row$se <- result$se; row$df <- result$df; row$n_used <- result$n_used
          half <- stats::qt(.975, result$df) * result$se
          row$lower <- result$estimate - half; row$upper <- result$estimate + half; row$width <- 2 * half
          row$coverage <- as.numeric(row$lower <= p$truth && row$upper >= p$truth)
          row$error <- result$estimate - p$truth
        }
      }
      idx <- idx + 1L; records[[idx]] <- row
    }
  list(records = records, decompositions = decompositions)
}
records <- list(); decompositions <- list()
start <- Sys.time()
for (s in seq_along(scenario_list)) {
  scenario <- scenario_list[[s]]
  reps <- if (scenario$family == "synthetic") config$primary_reps else config$clinical_reps
  chunks <- split(seq_len(reps), ceiling(seq_len(reps) / 50L))
  for (chunk in chunks) {
    # Every data set and method has its own explicit seed, so scheduling cannot
    # choose the scientific result. Windows uses the identical sequential path.
    batch <- parallel::mclapply(chunk, function(r) run_replicate(s, r),
      mc.cores = config$workers, mc.preschedule = TRUE, mc.set.seed = FALSE)
    if (any(vapply(batch, inherits, logical(1), what = "try-error"))) {
      stop("A replicate worker failed; earlier checkpoint retained. Do not promote this run.")
    }
    for (item in batch) {
      records <- c(records, item$records)
      decompositions <- c(decompositions, item$decompositions)
    }
    recent <- do.call(rbind, records)
    utils::write.csv(recent, file.path(outdir, "checkpoint-replicates.csv"), row.names = FALSE)
    cat(scenario$family, scenario$scenario, max(chunk), "/", reps, "elapsed",
      round(as.numeric(difftime(Sys.time(), start, units = "secs"))), "s; status counts: ",
      paste(names(table(recent$status)), table(recent$status), collapse = ", "), "\n")
    flush.console()
  }
}
if (!identical(hashes_before, hash_sources())) stop("Source/plan changed during execution. Checkpoint retained; this run is not verified.")
results <- do.call(rbind, records)
utils::write.csv(results, file.path(outdir, "replicates.csv"), row.names = FALSE)
utils::write.csv(do.call(rbind, decompositions), file.path(outdir, "bias-decomposition.csv"), row.names = FALSE)
mcse <- function(x) if (length(x) > 1L) stats::sd(x) / sqrt(length(x)) else NA_real_
wilson <- function(k, n) {
  if (!n) return(c(NA_real_, NA_real_))
  p <- k / n; z <- stats::qnorm(.975); denom <- 1 + z^2 / n
  centre <- (p + z^2 / (2 * n)) / denom
  half <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / denom
  c(max(0, centre - half), min(1, centre + half))
}
groups <- split(results, interaction(results$family, results$scenario, results$method, drop = TRUE))
summary_rows <- lapply(seq_along(groups), function(i) {
  g <- groups[[i]]; good <- g[g$status == "ok", ]; n <- nrow(good); requested <- nrow(g)
  coverage <- if (n) mean(good$coverage) else NA_real_; ci <- wilson(sum(good$coverage), n)
  rmse <- if (n) sqrt(mean(good$error^2)) else NA_real_
  bias <- if (n) mean(good$error) else NA_real_
  centered <- if (n) mean(good$lower - bias <= good$truth & good$upper - bias >= good$truth) else NA_real_
  centered_se <- NA_real_
  if (n > 1L) {
    set.seed(9500000L + i)
    b <- replicate(config$centered_bootstrap_reps, {
      z <- good[sample.int(n, n, replace = TRUE), ]; shift <- mean(z$error)
      mean(z$lower - shift <= z$truth & z$upper - shift >= z$truth)
    })
    centered_se <- stats::sd(b)
  }
  yielded <- sum(good$coverage) / requested; yield_ci <- wilson(sum(good$coverage), requested)
  data.frame(family = g$family[1], scenario = g$scenario[1], method = g$method[1], requested = requested,
    returned_intervals = n, refused = sum(g$status == "refused"), convergence_failures = sum(g$status == "convergence_failure"),
    other_errors = sum(g$status == "error"), interval_return_rate = n / requested,
    bias = bias, bias_mcse = mcse(good$error), rmse = rmse,
    rmse_mcse = if (n > 1 && rmse > 0) mcse(good$error^2) / (2 * rmse) else NA_real_,
    empirical_sd = if (n > 1) stats::sd(good$estimate) else NA_real_,
    mean_se = if (n) mean(good$se) else NA_real_, mean_se_mcse = mcse(good$se),
    mean_width = if (n) mean(good$width) else NA_real_, width_mcse = mcse(good$width),
    coverage_conditional = coverage, coverage_mcse = if (n) sqrt(coverage * (1 - coverage) / n) else NA_real_,
    coverage_wilson_lower = ci[1], coverage_wilson_upper = ci[2],
    covering_interval_per_request = yielded, covering_yield_mcse = sqrt(yielded * (1 - yielded) / requested),
    covering_yield_wilson_lower = yield_ci[1], covering_yield_wilson_upper = yield_ci[2],
    bias_centered_coverage_diagnostic = centered, centered_bootstrap_mcse = centered_se,
    approximate_nominal_criterion = if (n) ci[1] >= .925 && ci[2] <= .975 else NA)
})
summary_table <- do.call(rbind, summary_rows)
utils::write.csv(summary_table, file.path(outdir, "summary.csv"), row.names = FALSE)
status_rows <- do.call(rbind, lapply(groups, function(g) do.call(rbind, lapply(c("ok", "refused", "convergence_failure", "error"), function(status) {
  k <- sum(g$status == status); n <- nrow(g); ci <- wilson(k, n); rate <- k / n
  data.frame(family = g$family[1], scenario = g$scenario[1], method = g$method[1], status = status,
    count = k, requested = n, rate = rate, mcse = sqrt(rate * (1 - rate) / n), lower = ci[1], upper = ci[2])
}))))
utils::write.csv(status_rows, file.path(outdir, "request-status.csv"), row.names = FALSE)
paired <- list(); pi <- 0L
for (family in unique(results$family)) for (scenario in unique(results$scenario[results$family == family])) {
  g <- results[results$family == family & results$scenario == scenario & results$status == "ok", ]
  for (method in intersect(c("model_ignore", "model_exact", "model_noisy", "model_wrong_noise", "legacy"), unique(g$method))) {
    a <- g[g$method == method, ]
    for (baseline in c("trusted_pmm", "trusted_normal")) {
      b <- g[g$method == baseline, ]; both <- merge(a, b, by = "replicate", suffixes = c("_new", "_baseline"))
      de <- both$error_new - both$error_baseline; dm <- both$error_new^2 - both$error_baseline^2
      dc <- both$coverage_new - both$coverage_baseline; dw <- both$width_new - both$width_baseline
      pi <- pi + 1L
      paired[[pi]] <- data.frame(family = family, scenario = scenario, method = method, comparator = baseline,
        complete_pairs = nrow(both), error_difference = if (nrow(both)) mean(de) else NA_real_, error_mcse = mcse(de),
        mse_difference = if (nrow(both)) mean(dm) else NA_real_, mse_mcse = mcse(dm),
        coverage_difference = if (nrow(both)) mean(dc) else NA_real_, coverage_difference_mcse = mcse(dc),
        width_difference = if (nrow(both)) mean(dw) else NA_real_, width_difference_mcse = mcse(dw))
    }
  }
}
utils::write.csv(do.call(rbind, paired), file.path(outdir, "paired-comparisons.csv"), row.names = FALSE)
decomposition <- do.call(rbind, decompositions)
ds <- split(decomposition, interaction(decomposition$family, decomposition$scenario, drop = TRUE))
dsum <- do.call(rbind, lapply(ds, function(g) data.frame(family = g$family[1], scenario = g$scenario[1], replicates = nrow(g),
  selection_difference = mean(g$selection_difference), selection_mcse = mcse(g$selection_difference),
  residual_error_difference = mean(g$residual_error_difference), residual_error_mcse = mcse(g$residual_error_difference))))
utils::write.csv(dsum, file.path(outdir, "bias-decomposition-summary.csv"), row.names = FALSE)
utils::write.csv(hashes_before, file.path(outdir, "source-hashes.csv"), row.names = FALSE)
writeLines(c(paste("Completed UTC:", format(Sys.time(), tz = "UTC", usetz = TRUE)),
  paste("Elapsed seconds:", round(as.numeric(difftime(Sys.time(), start, units = "secs")), 1)),
  capture.output(sessionInfo())), file.path(outdir, "session-info.txt"))
fmt <- function(x, digits = 3L) ifelse(is.na(x), "NA", formatC(x, digits = digits, format = "f"))
rows <- vapply(seq_len(nrow(summary_table)), function(i) {
  z <- summary_table[i, ]
  paste0("| ", z$family, "/", z$scenario, " | ", z$method, " | ", z$returned_intervals, "/", z$requested,
    " | ", z$refused, "/", z$convergence_failures, "/", z$other_errors, " | ", fmt(z$bias), " (", fmt(z$bias_mcse), ")",
    " | ", fmt(z$rmse), " | ", fmt(z$coverage_conditional), " [", fmt(z$coverage_wilson_lower), ", ", fmt(z$coverage_wilson_upper), "]",
    " | ", fmt(z$mean_width), " | ", fmt(z$bias_centered_coverage_diagnostic), " |")
}, character(1))
writeLines(c(paste("# Coverage evaluation:", mode), "",
  "The table reports conditional coverage only among returned intervals. Refusal, nonconvergence and other errors remain in the requested denominator. Bias-centered coverage is a retrospective simulation diagnostic using hidden truth, not a corrected confidence interval.", "",
  "| Scenario | Method | Intervals/requested | Refused/MCMC/error | Bias (MCSE) | RMSE | Coverage [95% Wilson] | Mean width | Bias-centered diagnostic |",
  "|---|---|---:|---:|---:|---:|---:|---:|---:|", rows, "",
  "All MCSEs, interval yields, status rates and paired comparisons are in the adjacent CSV files. A model is not automatically validated by completing this script. Gibbs results use a smaller prespecified subset; primary and stress settings have different interpretation boundaries. Consult COVERAGE_PLAN.md before making a claim.", "",
  "The direct ignore-category model is an implementation of established normal multiple imputation after masking quality flags. Exact/noisy category models additionally assume the declared category likelihood. Missed corruption, inaccurate category likelihoods, model misspecification and incompatible analysis can defeat inference.", "",
  "source-hashes.csv verifies the files read before execution against the end of this run. settings.json records the fixed computational settings. The checkpoint was removed only after successful source-integrity verification and complete result generation."), file.path(outdir, "RESULTS.md"))
unlink(file.path(outdir, "checkpoint-replicates.csv"))
cat("Completed", mode, "evaluation:", nrow(results), "requests; status counts:",
    paste(names(table(results$status)), table(results$status), collapse = ", "), "\n")
cat("Output:", normalizePath(outdir), "\n")
