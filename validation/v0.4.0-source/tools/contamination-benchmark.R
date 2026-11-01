#!/usr/bin/env Rscript
# Prospective v0.4 study; settings locked in validation/CONTAMINATION_PLAN.md.
# Usage: Rscript tools/contamination-benchmark.R pilot|study [NEW_DIRECTORY] [PILOT_REPS]
args <- commandArgs(trailingOnly = TRUE)
mode <- if (length(args)) match.arg(args[1], c("pilot", "study")) else "pilot"
file_arg <- grep("^--file=", commandArgs(), value = TRUE)
root <- normalizePath(file.path(dirname(sub("^--file=", "", file_arg[1])), ".."))
outdir <- if (length(args) >= 2L) args[2] else file.path(root, "validation", paste0("contamination-", mode))
if (dir.exists(outdir) && length(list.files(outdir, all.files = TRUE, no.. = TRUE)))
  stop("Use a new or empty output directory; existing evidence is never overwritten.")
hash_sources <- function() {
  paths <- c(sort(list.files(file.path(root, "R"), "\\.R$", full.names = TRUE)),
    file.path(root, c("DESCRIPTION", "NAMESPACE", "tools/contamination-benchmark.R", "tools/contamination-fixtures.R",
                      "tools/lab-fixture.R", "validation/CONTAMINATION_PLAN.md")))
  data.frame(file = substring(paths, nchar(root) + 2L),
    sha256 = vapply(paths, function(p) digest::digest(file = p, algo = "sha256"), character(1)), row.names = NULL)
}
hashes_before <- hash_sources()
pkgload::load_all(root, quiet = TRUE)
source(file.path(root, "tools/lab-fixture.R"), local = TRUE)
source(file.path(root, "tools/contamination-fixtures.R"), local = TRUE)
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
workers <- suppressWarnings(as.integer(Sys.getenv("CATMECH_BENCHMARK_WORKERS", "14")))
if (length(workers) != 1L || is.na(workers) || workers < 1L || workers > 16L) stop("Invalid worker count (1 to 16).")
if (.Platform$OS.type == "windows") workers <- 1L
pilot_reps <- if (length(args) >= 3L) as.integer(args[3]) else 5L
config <- list(mode = mode, workers = workers, m = 20L, n = 600L,
  primary_reps = if (mode == "pilot") pilot_reps else 1000L,
  stress_reps = if (mode == "pilot") pilot_reps else 150L,
  holdout_reps = if (mode == "pilot") min(pilot_reps, 5L) else 40L,
  holdout_folds = 5L, holdout_m = 40L, centered_bootstrap_reps = 400L,
  fixture_master_seed = 410273L, phase = if (mode == "pilot") "development" else "heldout")
jsonlite::write_json(config, file.path(outdir, "settings.json"), pretty = TRUE, auto_unbox = TRUE)
scenario_list <- contamination_scenarios(config$primary_reps, config$stress_reps)
all_methods <- c("oracle", "gold_only", "v02_mask_normal", "reference_normal", "cole_mi", "raking", "contamination")
matching_scenarios <- c("missed_flags", "reference_05pct", "reference_10pct", "reference_50pct", "heavy_tail_target", "nonlinear_target")
analysis_formula <- outcome ~ lab_value + age + sex + inflammation + SITEID

fit_lab <- function(lab, p) {
  use <- is.finite(lab)
  x <- cbind(p$analysis_x[use, , drop = FALSE], lab_value = lab[use])
  fit <- stats::lm.fit(x, p$input$outcome[use]); k <- ncol(x); n <- nrow(x)
  if (fit$rank != k || n <= k) stop("Scientific analysis is rank deficient or has no residual degrees of freedom.")
  position <- match(k, fit$qr$pivot)
  variance_factor <- chol2inv(qr.R(fit$qr))[position, position]
  list(estimate = unname(fit$coefficients[k]), se = sqrt(sum(fit$residuals^2) / (n - k) * variance_factor), df = n - k, n_used = n)
}
pool_fits <- function(fits) {
  q <- vapply(fits, function(x) x$estimate, numeric(1)); u <- vapply(fits, function(x) x$se^2, numeric(1))
  z <- mice::pool.scalar(q, u, n = fits[[1]]$n_used, k = fits[[1]]$n_used - fits[[1]]$df)
  list(estimate = z$qbar, se = sqrt(z$t), df = z$df, n_used = fits[[1]]$n_used)
}
# The strongest simple established imputer: a validation-sample linear model given
# the recorded value and its flag, with a separate fit for records lacking a value.
# Flag terms need ten flagged calibration rows; fewer leave their slope unstable.
fit_cole <- function(p) {
  y <- p$input$gold_value; w <- p$input$lab_value
  ref <- !is.na(y); obs <- !is.na(w)
  yes <- as.numeric(p$input$flag == "Yes")
  w0 <- ifelse(obs, w, 0)
  x_obs <- cbind(p$model_x, recorded = w0)
  if (sum(yes[obs & ref]) >= 10L) x_obs <- cbind(x_obs, flag = yes, recorded_flag = w0 * yes)
  # Aliased columns are dropped first; a sparse flag term must not fail this comparator.
  full_rank <- function(x, ry) {
    q <- qr(cbind(1, x[ry, , drop = FALSE]))
    x[, setdiff(q$pivot[seq_len(q$rank)], 1L) - 1L, drop = FALSE]
  }
  x_obs <- full_rank(x_obs, obs & ref)
  x_mis <- full_rank(p$model_x, ref)
  set.seed(p$imputation_seed)
  pool_fits(lapply(seq_len(config$m), function(j) {
    completed <- y
    wy <- obs & !ref
    if (any(wy)) completed[wy] <- mice::mice.impute.norm(y, obs & ref, x_obs, wy = wy)
    wy <- !obs & !ref
    if (any(wy)) completed[wy] <- mice::mice.impute.norm(y, ref, x_mis, wy = wy)
    fit_lab(completed, p)
  }))
}
# Generalized raking (Lumley, Shaw and Dai 2011) on the design-based side. Calibrates
# Bernoulli sampling weights to full-cohort influence functions of the analysis fit
# to a regression-calibrated proxy. It needs no recording model, only the sampling
# probabilities; under nonignorable selection the analyst has only the average rate.
fit_raking <- function(p) {
  y <- p$input$gold_value; w <- p$input$lab_value; outcome <- p$input$outcome
  ref <- !is.na(y); obs <- !is.na(w); pi_i <- p$sampling_probability
  yes <- as.numeric(p$input$flag == "Yes")
  w0 <- ifelse(obs, w, 0)
  x_obs <- cbind(1, p$model_x, recorded = w0)
  if (sum(yes[obs & ref]) >= 10L) x_obs <- cbind(x_obs, flag = yes, recorded_flag = w0 * yes)
  x_mis <- cbind(1, p$model_x)
  proxy <- numeric(length(y))
  use <- which(obs & ref)
  beta <- stats::lm.wfit(x_obs[use, , drop = FALSE], y[use], 1 / pi_i[use])$coefficients
  beta[is.na(beta)] <- 0
  proxy[obs] <- as.vector(x_obs[obs, , drop = FALSE] %*% beta)
  use <- which(ref)
  beta <- stats::lm.wfit(x_mis[use, , drop = FALSE], y[use], 1 / pi_i[use])$coefficients
  proxy[!obs] <- as.vector(x_mis[!obs, , drop = FALSE] %*% beta)
  design <- cbind(p$analysis_x, lab_value = proxy)
  naive <- stats::lm.fit(design, outcome)
  h <- design * naive$residuals
  h <- t(backsolve(qr.R(naive$qr), t(h), transpose = TRUE))
  aux <- cbind(1, h)
  target <- colSums(aux)
  d <- 1 / pi_i[ref]; a <- aux[ref, , drop = FALSE]
  lambda <- numeric(ncol(aux))
  for (iter in seq_len(100L)) {
    weight <- d * exp(as.vector(a %*% lambda))
    gap <- target - colSums(a * weight)
    if (max(abs(gap)) < 1e-8 * max(1, abs(target))) break
    step <- solve(crossprod(a * sqrt(weight)), gap)
    size <- 1
    objective <- function(l) sum(d * exp(as.vector(a %*% l))) - sum(target * l)
    while (size > 1e-4 && objective(lambda + size * step) > objective(lambda)) size <- size / 2
    lambda <- lambda + size * as.vector(step)
    if (iter == 100L) stop("Raking calibration did not converge.")
  }
  weight <- d * exp(as.vector(a %*% lambda))
  x <- cbind(p$analysis_x[ref, , drop = FALSE], lab_value = y[ref])
  fit <- stats::lm.wfit(x, outcome[ref], weight)
  k <- ncol(x); position <- k
  inverse <- solve(crossprod(x * sqrt(weight)))
  psi <- (x * fit$residuals) %*% inverse
  b <- stats::lm.wfit(a, psi[, position], weight)$coefficients
  g <- psi[, position] - as.vector(a %*% b)
  phase_one <- sum(weight * psi[, position]^2)
  phase_two <- sum((1 - pi_i[ref]) * weight^2 * g^2)
  list(estimate = unname(fit$coefficients[position]), se = sqrt(phase_one + phase_two),
       df = sum(ref) - k, n_used = sum(ref))
}
fit_contamination <- function(p, s) {
  catmech_contamination(p$input, "lab_value", "gold_value", p$predictors,
    corruption_strata = s$corruption_strata, error_strata = s$error_strata,
    agreement_tol = s$agreement_tol, reference_sd = s$reference_sd,
    m = config$m, seed = p$imputation_seed)
}
corruption_scores <- function(imp, p, s, r) {
  table <- imp$contamination$corruption
  scored <- which(table$evidence == "posterior_model" & !is.na(table$posterior_corruption))
  if (!length(scored)) return(NULL)
  y <- p$corrupted[scored]; q <- table$posterior_corruption[scored]
  ranks <- rank(q)
  positives <- sum(y); negatives <- sum(!y)
  auc <- if (positives && negatives) (sum(ranks[y]) - positives * (positives + 1) / 2) / (positives * negatives) else NA_real_
  data.frame(family = s$family, scenario = s$scenario, replicate = r, rows = length(scored),
    actual_errors = positives, expected_errors = sum(q), brier = mean((q - y)^2),
    brier_prevalence = mean((mean(y) - y)^2), auc = auc,
    mean_probability_when_wrong = if (positives) mean(q[y]) else NA_real_,
    mean_probability_when_correct = if (negatives) mean(q[!y]) else NA_real_,
    completion_rmse = sqrt(mean((table$posterior_target - p$truth)^2)),
    recorded_rmse = sqrt(mean((p$input$lab_value[scored] - p$truth[scored])^2)))
}
run_method <- function(method, p, s, r) {
  y <- p$input$lab_value
  if (method == "oracle") return(list(fit = fit_lab(p$truth, p)))
  if (method == "gold_only") return(list(fit = fit_lab(p$input$gold_value, p)))
  if (method == "cole_mi") return(list(fit = fit_cole(p)))
  if (method == "raking") return(list(fit = fit_raking(p)))
  imp <- if (method %in% c("reference_normal", "reference_pmm", "reference_pmm_abb")) {
    catmech_reference(p$input, "lab_value", "gold_value", p$predictors, method = sub("reference_", "", method),
                      m = config$m, seed = p$imputation_seed)
  } else if (method == "v02_mask_normal") {
    catmech_model(p$input, "lab_value", p$predictors, flag = "flag", category_mode = "ignore",
                  m = config$m, seed = p$imputation_seed)
  } else fit_contamination(p, s)
  fit <- pool_fits(lapply(imp$completed, function(d) fit_lab(d$lab_value, p)))
  out <- list(fit = fit)
  if (method == "contamination") {
    out$scores <- corruption_scores(imp, p, s, r)
    if (isTRUE(s$holdout) && r <= config$holdout_reps) {
      h <- catmech_holdout_check(imp, folds = config$holdout_folds, m = config$holdout_m, seed = p$imputation_seed)
      row <- if (is.null(h$summary)) data.frame(predictive_coverage = NA_real_, ks_p_value = NA_real_,
        calibration_z = NA_real_, errors_observed = NA_real_, errors_expected = NA_real_) else
        h$summary[c("predictive_coverage", "ks_p_value", "calibration_z", "errors_observed", "errors_expected")]
      out$holdout <- cbind(data.frame(family = s$family, scenario = s$scenario, replicate = r, status = h$status,
        flags = paste(h$flags, collapse = "|"), failures = length(h$failures), stringsAsFactors = FALSE), row)
    }
  }
  out
}
run_replicate <- function(s, r) {
  p <- make_contamination_case(s, r, config$phase, config$fixture_master_seed, make_lab_fixture)
  records <- list(); scores <- NULL; holdout <- NULL
  reference_rows <- !is.na(p$input$gold_value)
  for (method in c(if (is.null(s$methods)) all_methods else s$methods,
                   if (s$scenario %in% matching_scenarios) c("reference_pmm", "reference_pmm_abb"))) {
    warnings <- character()
    result <- tryCatch(withCallingHandlers(run_method(method, p, s, r), warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning")
    }), error = identity)
    failed <- inherits(result, "error")
    status <- if (!failed) "ok" else if (inherits(result, "catmech_convergence_error")) "convergence_failure" else
      if (grepl("No calibration rows support|full column rank|more observations than coefficients|Fewer trusted/reference|lost design support", conditionMessage(result))) "refused" else "error"
    row <- data.frame(family = s$family, scenario = s$scenario, replicate = r, data_seed = p$selection_seed,
      data_stream = p$stream, imputation_seed = p$imputation_seed, method = method, truth = p$slope,
      available_n = sum(reference_rows), n = nrow(p$input), estimate = NA_real_, se = NA_real_, df = NA_real_, n_used = NA_real_,
      lower = NA_real_, upper = NA_real_, width = NA_real_, coverage = NA_real_, error = NA_real_,
      status = status, message = if (failed) conditionMessage(result) else "",
      warnings = paste(unique(warnings), collapse = " | "), stringsAsFactors = FALSE)
    if (!failed) {
      fit <- result$fit
      if (!is.finite(fit$estimate) || !is.finite(fit$se) || fit$se <= 0 || is.na(fit$df) || fit$df <= 0) {
        row$status <- "error"; row$message <- "Nonfinite estimate or degenerate uncertainty."
      } else {
        half <- stats::qt(.975, fit$df) * fit$se
        row$estimate <- fit$estimate; row$se <- fit$se; row$df <- fit$df; row$n_used <- fit$n_used
        row$lower <- fit$estimate - half; row$upper <- fit$estimate + half; row$width <- row$upper - row$lower
        row$coverage <- as.numeric(row$lower <= p$slope && row$upper >= p$slope); row$error <- fit$estimate - p$slope
      }
      scores <- result$scores
      holdout <- result$holdout
    }
    records[[length(records) + 1L]] <- row
  }
  list(records = do.call(rbind, records), scores = scores, holdout = holdout)
}
records <- list(); scores <- list(); holdouts <- list(); start <- Sys.time()
for (s in scenario_list) {
  for (chunk in split(seq_len(s$reps), ceiling(seq_len(s$reps) / 120L))) {
    batch <- parallel::mclapply(chunk, function(r) run_replicate(s, r), mc.cores = workers, mc.preschedule = TRUE, mc.set.seed = FALSE)
    if (any(vapply(batch, inherits, logical(1), what = "try-error"))) stop("A worker failed; checkpoint is not verified evidence.")
    records <- c(records, lapply(batch, `[[`, "records")); scores <- c(scores, lapply(batch, `[[`, "scores"))
    holdouts <- c(holdouts, lapply(batch, `[[`, "holdout"))
    saveRDS(list(records = records, scores = scores, holdouts = holdouts), file.path(outdir, "checkpoint.rds"))
    cat(s$family, s$scenario, max(chunk), "/", s$reps, "; elapsed", round(as.numeric(difftime(Sys.time(), start, units = "secs"))), "seconds\n"); flush.console()
  }
}
if (!identical(hashes_before, hash_sources())) stop("Source or plan changed during run. Checkpoint retained; evidence is not verified.")
results <- do.call(rbind, records)
score_table <- do.call(rbind, scores)
holdout_table <- do.call(rbind, holdouts)
utils::write.csv(results, file.path(outdir, "replicates.csv"), row.names = FALSE)
utils::write.csv(score_table, file.path(outdir, "corruption-scores.csv"), row.names = FALSE)
utils::write.csv(if (is.null(holdout_table)) data.frame() else holdout_table, file.path(outdir, "holdout-checks.csv"), row.names = FALSE)
mcse <- function(x) if (length(x) > 1L) stats::sd(x) / sqrt(length(x)) else NA_real_
wilson <- function(k, n) {
  if (!n) return(c(NA_real_, NA_real_))
  p <- k / n; z <- stats::qnorm(.975); denom <- 1 + z^2 / n
  centre <- (p + z^2 / (2 * n)) / denom; half <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / denom
  c(max(0, centre - half), min(1, centre + half))
}
groups <- split(results, interaction(results$family, results$scenario, results$method, drop = TRUE))
summary_table <- do.call(rbind, lapply(seq_along(groups), function(i) {
  g <- groups[[i]]; good <- g[g$status == "ok", ]; n <- nrow(good); requested <- nrow(g)
  coverage <- if (n) mean(good$coverage) else NA_real_; ci <- wilson(sum(good$coverage), n)
  rmse <- if (n) sqrt(mean(good$error^2)) else NA_real_; bias <- if (n) mean(good$error) else NA_real_
  centered <- if (n) mean(good$lower - bias <= good$truth & good$upper - bias >= good$truth) else NA_real_
  centered_se <- sd_mcse <- NA_real_
  if (n > 1L) {
    set.seed(9910000L + i)
    b <- replicate(config$centered_bootstrap_reps, {
      z <- good[sample.int(n, n, replace = TRUE), ]; shift <- mean(z$error)
      c(centered = mean(z$lower - shift <= z$truth & z$upper - shift >= z$truth), empirical_sd = stats::sd(z$estimate))
    })
    centered_se <- stats::sd(b[1, ]); sd_mcse <- stats::sd(b[2, ])
  }
  primary <- g$family[1] == "primary" && g$method[1] == "contamination"
  data.frame(family = g$family[1], scenario = g$scenario[1], method = g$method[1], requested = requested,
    returned_intervals = n, refused = sum(g$status == "refused"), convergence_failures = sum(g$status == "convergence_failure"),
    other_errors = sum(g$status == "error"), interval_return_rate = n / requested,
    warned = sum(nzchar(g$warnings)), mean_available_n = mean(g$available_n), bias = bias, bias_mcse = mcse(good$error), rmse = rmse,
    rmse_mcse = if (n > 1 && rmse > 0) mcse(good$error^2) / (2 * rmse) else NA_real_,
    empirical_sd = if (n > 1) stats::sd(good$estimate) else NA_real_, empirical_sd_bootstrap_mcse = sd_mcse,
    mean_se = if (n) mean(good$se) else NA_real_, mean_se_mcse = mcse(good$se), mean_width = if (n) mean(good$width) else NA_real_, width_mcse = mcse(good$width),
    coverage_conditional = coverage, coverage_mcse = if (n) sqrt(coverage * (1 - coverage) / n) else NA_real_,
    coverage_wilson_lower = ci[1], coverage_wilson_upper = ci[2], covering_interval_per_request = if (requested) sum(good$coverage) / requested else NA_real_,
    bias_centered_coverage_diagnostic = centered, centered_bootstrap_mcse = centered_se,
    primary_criterion_applicable = primary && requested >= 1000L,
    approximate_nominal_criterion = if (primary && n >= 1000L) ci[1] >= .925 && ci[2] <= .975 else NA)
}))
utils::write.csv(summary_table, file.path(outdir, "summary.csv"), row.names = FALSE)
status_rows <- do.call(rbind, lapply(groups, function(g) do.call(rbind, lapply(c("ok", "refused", "convergence_failure", "error"), function(status) {
  k <- sum(g$status == status); n <- nrow(g); ci <- wilson(k, n); rate <- k / n
  data.frame(family = g$family[1], scenario = g$scenario[1], method = g$method[1], status = status,
    count = k, requested = n, rate = rate, mcse = sqrt(rate * (1 - rate) / n), lower = ci[1], upper = ci[2])
}))))
utils::write.csv(status_rows, file.path(outdir, "request-status.csv"), row.names = FALSE)
paired <- list()
for (scenario in unique(results$scenario)) {
  g <- results[results$scenario == scenario & results$status == "ok", ]
  for (comparator in c("reference_normal", "cole_mi", "raking", "gold_only", "oracle", "v02_mask_normal")) {
    if (!all(c("contamination", comparator) %in% g$method)) next
    both <- merge(g[g$method == "contamination", ], g[g$method == comparator, ], by = "replicate", suffixes = c("_new", "_baseline"))
    de <- both$error_new - both$error_baseline; dm <- both$error_new^2 - both$error_baseline^2
    dc <- both$coverage_new - both$coverage_baseline; dw <- both$width_new - both$width_baseline
    paired[[length(paired) + 1L]] <- data.frame(family = g$family[1], scenario = scenario, method = "contamination", comparator = comparator,
      complete_pairs = nrow(both), error_difference = mean(de), error_mcse = mcse(de), mse_difference = mean(dm), mse_mcse = mcse(dm),
      coverage_difference = mean(dc), coverage_difference_mcse = mcse(dc), width_difference = mean(dw), width_difference_mcse = mcse(dw))
  }
}
utils::write.csv(if (length(paired)) do.call(rbind, paired) else data.frame(), file.path(outdir, "paired-comparisons.csv"), row.names = FALSE)
score_groups <- split(score_table, score_table$scenario)
score_summary <- do.call(rbind, lapply(score_groups, function(g) data.frame(family = g$family[1], scenario = g$scenario[1], replicates = nrow(g),
  mean_actual_errors = mean(g$actual_errors), mean_expected_errors = mean(g$expected_errors),
  mean_error_count_gap = mean(g$expected_errors - g$actual_errors), error_count_gap_mcse = mcse(g$expected_errors - g$actual_errors),
  brier = mean(g$brier), brier_prevalence = mean(g$brier_prevalence), auc = mean(g$auc, na.rm = TRUE),
  mean_probability_when_wrong = mean(g$mean_probability_when_wrong, na.rm = TRUE),
  mean_probability_when_correct = mean(g$mean_probability_when_correct, na.rm = TRUE),
  completion_rmse = mean(g$completion_rmse), recorded_rmse = mean(g$recorded_rmse))))
utils::write.csv(score_summary, file.path(outdir, "corruption-scores-summary.csv"), row.names = FALSE)
if (!is.null(holdout_table)) {
  holdout_groups <- split(holdout_table, holdout_table$scenario)
  utils::write.csv(do.call(rbind, lapply(holdout_groups, function(g) data.frame(family = g$family[1], scenario = g$scenario[1], replicates = nrow(g),
    detected = sum(g$status == "miscalibration_detected"), unresolved = sum(g$status == "unresolved_failed_folds"),
    detection_rate = mean(g$status == "miscalibration_detected"), mean_coverage = mean(g$predictive_coverage, na.rm = TRUE)))),
    file.path(outdir, "holdout-summary.csv"), row.names = FALSE)
}
utils::write.csv(hashes_before, file.path(outdir, "source-hashes.csv"), row.names = FALSE)
writeLines(c(paste("Completed UTC:", format(Sys.time(), tz = "UTC", usetz = TRUE)),
  paste("Elapsed seconds:", round(as.numeric(difftime(Sys.time(), start, units = "secs")), 1)), capture.output(sessionInfo())), file.path(outdir, "session-info.txt"))
fmt <- function(x) ifelse(is.na(x), "NA", formatC(x, digits = 3, format = "f"))
rows <- vapply(seq_len(nrow(summary_table)), function(i) {
  z <- summary_table[i, ]
  paste0("| ", z$family, "/", z$scenario, " | ", z$method, " | ", z$returned_intervals, "/", z$requested,
    " | ", fmt(z$bias), " (", fmt(z$bias_mcse), ") | ", fmt(z$coverage_conditional), " (", fmt(z$coverage_mcse), ") [",
    fmt(z$coverage_wilson_lower), ", ", fmt(z$coverage_wilson_upper), "] | ", fmt(z$mean_width), " (", fmt(z$width_mcse), ") |")
}, character(1))
writeLines(c(paste("# v0.4 contamination evaluation:", mode), "",
  "Prespecified settings, assumptions and interpretation boundaries are in ../CONTAMINATION_PLAN.md. Coverage is conditional on returned intervals; request-status.csv retains refusals and failures.", "",
  "| Setting | Method | Intervals/requested | Bias (MCSE) | Coverage (MCSE) [95% Wilson] | Width (MCSE) |", "|---|---|---:|---:|---:|---:|", rows, "",
  "The contamination model's correct-model cells use the additive normal error it assumes, so they alone cannot show robustness. The wrong_model cells are the deliberate counterweight. Summary, paired contrasts, corruption scores and hold-out checks are in adjacent CSVs. No result is made valid merely by this script completing."), file.path(outdir, "RESULTS.md"))
unlink(file.path(outdir, "checkpoint.rds"))
cat("Completed", mode, ":", nrow(results), "requests;", paste(names(table(results$status)), table(results$status), collapse = ", "), "\n")
cat("Output:", normalizePath(outdir), "\n")
