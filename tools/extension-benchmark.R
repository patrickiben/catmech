#!/usr/bin/env Rscript
# Prospective v0.3 study; settings locked in validation/EXTENSION_PLAN.md.
args <- commandArgs(trailingOnly = TRUE)
mode <- if (length(args)) match.arg(args[1], c("pilot", "distribution", "reference", "stress", "controls")) else "pilot"
file_arg <- grep("^--file=", commandArgs(), value = TRUE)
root <- normalizePath(file.path(dirname(sub("^--file=", "", file_arg[1])), ".."))
outdir <- if (length(args) >= 2L) args[2] else file.path(root, "validation", paste0("extension-", mode))
if (dir.exists(outdir) && length(list.files(outdir, all.files = TRUE, no.. = TRUE)))
  stop("Use a new or empty output directory; existing evidence is never overwritten.")
hash_sources <- function() {
  paths <- c(sort(list.files(file.path(root, "R"), "\\.R$", full.names = TRUE)),
    file.path(root, c("DESCRIPTION", "NAMESPACE", "tools/extension-benchmark.R", "tools/lab-fixture.R", "validation/EXTENSION_PLAN.md")))
  data.frame(file = substring(paths, nchar(root) + 2L),
    sha256 = vapply(paths, function(p) digest::digest(file = p, algo = "sha256"), character(1)), row.names = NULL)
}
hashes_before <- hash_sources()
pkgload::load_all(root, quiet = TRUE)
source(file.path(root, "tools/lab-fixture.R"), local = TRUE)
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
workers <- suppressWarnings(as.integer(Sys.getenv("CATMECH_BENCHMARK_WORKERS", "4")))
if (length(workers) != 1L || is.na(workers) || workers < 1L || workers > 16L) stop("Invalid worker count (1 to 16).")
if (.Platform$OS.type == "windows") workers <- 1L
config <- list(mode = mode, workers = workers, m = 20L, donors = 5L, synthetic_n = 300L,
  reference_n = 600L, primary_reps = if (mode == "pilot") 5L else 1000L,
  support_reps = if (mode == "pilot") 5L else 250L,
  noisy_gold_reps = if (mode == "pilot") 5L else 500L,
  control_reps = if (mode == "pilot") 5L else 20L, centered_bootstrap_reps = 400L,
  fixture_master_seed = 390173L)
jsonlite::write_json(config, file.path(outdir, "settings.json"), pretty = TRUE, auto_unbox = TRUE)
synthetic <- data.frame(scenario = c("missing_only", "perfect_flags", "false_positive_flags", "category_error", "null_effect", "heavy_tails", "heavy_category_error", "strong_mar"),
  missing_rate = c(rep(.2, 7), 0), flag_rate = c(0, rep(.1, 6), 0),
  corruption_probability = c(1, 1, .5, 1, .5, 1, 1, 1),
  category_error_rate = c(0, 0, 0, .1, 0, 0, .1, 0),
  error_distribution = c(rep("normal", 5), "heavy_tailed", "heavy_tailed", "normal"),
  truth = c(rep(.8, 4), 0, rep(.8, 3)), family = "distribution", stringsAsFactors = FALSE)
scenario_list <- list()
if (mode %in% c("pilot", "distribution")) scenario_list <- lapply(seq_len(nrow(synthetic)), function(i) {
  z <- as.list(synthetic[i, ]); z$index <- i; z$reps <- config$primary_reps; z
})
add_reference <- function(scenario, fixture_scenario = scenario, fraction = .25, selection = "random", family = "reference", reps = config$primary_reps) {
  scenario_list[[length(scenario_list) + 1L]] <<- list(scenario = scenario, fixture_scenario = fixture_scenario,
    fraction = fraction, selection = selection, family = family, reps = reps)
}
if (mode %in% c("pilot", "reference")) for (s in c("clean", "flagged_errors", "missed_flags", "site_shift", "rounded_records")) add_reference(s)
if (mode %in% c("pilot", "stress")) {
  for (s in c("missed_flags", "site_shift")) for (f in c(.10, .50))
    add_reference(paste0(s, "_gold", f * 100), s, f, family = "support", reps = config$support_reps)
  add_reference("nonignorable_selection", "clean", selection = "mnar", family = "stress")
  add_reference("noisy_reference", "missed_flags", selection = "noisy", family = "stress", reps = config$noisy_gold_reps)
}
if (mode %in% c("pilot", "controls")) {
  add_reference("all_reference", "missed_flags", fraction = 1, family = "control", reps = config$control_reps)
  add_reference("no_reference", "missed_flags", fraction = 0, family = "control", reps = config$control_reps)
}
thresholds <- c(-Inf, 5.5, 7, Inf)
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
  if (length(unique(vapply(fits, function(x) x$n_used, numeric(1)))) != 1L) stop("Analysis sample sizes differ.")
  z <- mice::pool.scalar(q, u, n = fits[[1]]$n_used, k = fits[[1]]$n_used - fits[[1]]$df)
  list(estimate = z$qbar, se = sqrt(z$t), df = z$df, n_used = fits[[1]]$n_used)
}
prepare_data <- function(s, r) {
  if (s$family == "distribution") {
    seed <- as.integer((if (mode == "pilot") 50000000 else 60000000) + 100000 * s$index + r)
    d <- simulate_catmech_lab(n = config$synthetic_n, seed = seed, thresholds = thresholds,
      missing_rate = s$missing_rate, flag_rate = s$flag_rate, corruption_probability = s$corruption_probability,
      category_error_rate = s$category_error_rate, error_distribution = s$error_distribution)
    if (s$scenario == "null_effect") { set.seed(seed + 125000L); d$outcome <- stats::rnorm(nrow(d)) }
    if (s$scenario == "strong_mar") {
      set.seed(seed + 42000L); d$lab_value[stats::runif(nrow(d)) < stats::plogis(3 * (d$inflammation - 2))] <- NA_real_
    }
    predictors <- c("age", "sex", "inflammation", "outcome")
    input <- d[, c(predictors, "lab_value", "lab_category", "outlier_flag")]
    actual_truth <- d$lab_true; truth <- s$truth; stream <- as.character(seed)
    imputation_seed <- seed + 1000000L
  } else {
    f <- make_lab_fixture(n = config$reference_n, replicate = r,
      phase = if (mode == "pilot") "development" else "heldout", seed = config$fixture_master_seed,
      scenario = s$fixture_scenario)
    check_lab_fixture(f)
    predictors <- c("age", "sex", "inflammation", "SITEID", "outcome")
    input <- f$data[, predictors]; input$lab_value <- f$data$LBSTRESN
    input$lab_category <- f$data$LBNRIND; input$outlier_flag <- f$data$QCFLAG
    actual_truth <- f$truth$LAB_TRUE; truth <- f$manifest$known_slope
    seed <- as.integer((if (mode == "pilot") 70100000 else 80100000) + r)
    set.seed(seed); selection_u <- stats::runif(nrow(input)); probability <- rep(s$fraction, nrow(input))
    if (s$selection == "mnar") {
      mx <- stats::model.matrix(stats::reformulate(predictors), data = input)
      hidden_residual <- as.numeric(scale(stats::lm.fit(mx, actual_truth)$residuals))
      intercept <- stats::uniroot(function(a) mean(stats::plogis(a + 2 * hidden_residual)) - .25, c(-30, 30))$root
      probability <- stats::plogis(intercept + 2 * hidden_residual)
    }
    gold <- actual_truth
    if (s$selection == "noisy") { set.seed(86100000L + r); gold <- gold + stats::rnorm(nrow(input), 0, .08) }
    gold[selection_u >= probability] <- NA_real_; input$gold_value <- gold
    stream <- paste0("master", config$fixture_master_seed, ":", if (mode == "pilot") "development" else "heldout", ":", r)
    imputation_seed <- as.integer((if (mode == "pilot") 75000000 else 90000000) + r)
  }
  analysis_predictors <- setdiff(predictors, "outcome")
  list(input = input, actual_truth = actual_truth, truth = truth, predictors = predictors,
    analysis_x = stats::model.matrix(stats::reformulate(analysis_predictors), data = input),
    analysis_formula = stats::reformulate(c("lab_value", analysis_predictors), response = "outcome"),
    model_x = stats::model.matrix(stats::reformulate(predictors), data = input)[, -1L, drop = FALSE],
    category_x = stats::model.matrix(stats::reformulate(c(predictors, "lab_category")), data = input)[, -1L, drop = FALSE],
    flagged = as.character(input$outlier_flag) == "Yes", seed = seed, stream = stream, imputation_seed = imputation_seed)
}
methods_for <- function(s) {
  if (s$family == "control") return(c("reference_normal", "reference_pmm"))
  if (s$family == "distribution") return(c("oracle", "trusted_observed", "normal_no_cat", "normal_cat", "pmm_no_cat", "pmm_cat",
    "mice_normal_no_cat", "mice_normal_cat", "mice_pmm_no_cat", "mice_pmm_cat"))
  c("oracle", "v02_mask_normal", "gold_only", "reference_normal", "reference_pmm", "mice_reference_normal", "mice_reference_pmm")
}
run_method <- function(method, p, s) {
  y <- p$input$lab_value
  if (method == "oracle") return(fit_lab(p$actual_truth, p))
  if (method == "gold_only") return(fit_lab(p$input$gold_value, p))
  if (method == "trusted_observed") { y[p$flagged] <- NA_real_; return(fit_lab(y, p)) }
  has_category <- grepl("_cat$", method) && !grepl("_no_cat$", method)
  predictors <- c(p$predictors, if (has_category) "lab_category")
  if (startsWith(method, "mice_")) {
    if (grepl("reference", method)) y <- p$input$gold_value else y[p$flagged] <- NA_real_
    ry <- !is.na(y); wy <- !ry; x <- if (has_category) p$category_x else p$model_x
    set.seed(p$imputation_seed)
    fits <- lapply(seq_len(config$m), function(j) {
      completed <- y
      completed[wy] <- if (grepl("pmm", method)) mice::mice.impute.pmm(y, ry, x, wy = wy, donors = config$donors, matchtype = 1L, ridge = 1e-5)
        else mice::mice.impute.norm(y, ry, x, wy = wy)
      fit_lab(completed, p)
    })
    return(pool_fits(fits))
  }
  if (startsWith(method, "reference_")) {
    imp <- catmech_reference(p$input, value = "lab_value", reference = "gold_value", predictors = predictors,
      method = if (method == "reference_pmm") "pmm" else "normal", donors = config$donors, m = config$m, seed = p$imputation_seed)
  } else if (startsWith(method, "pmm_")) {
    imp <- catmech_pmm(p$input, value = "lab_value", predictors = predictors, flag = "outlier_flag",
      donors = config$donors, m = config$m, seed = p$imputation_seed)
  } else {
    imp <- catmech_model(p$input, value = "lab_value", predictors = predictors, flag = "outlier_flag",
      category_mode = "ignore", m = config$m, seed = p$imputation_seed)
  }
  if (s$family == "control" && s$scenario == "all_reference") {
    for (d in imp$completed) if (!identical(d$lab_value, p$actual_truth)) stop("All-reference target differs from truth.")
    pooled <- catmech_pool_lm(imp, p$analysis_formula); pooled <- pooled[pooled$term == "lab_value", ]
    expected <- fit_lab(p$actual_truth, p); half <- stats::qt(.975, expected$df) * expected$se
    delta <- max(abs(c(pooled$estimate - expected$estimate, pooled$std.error - expected$se,
      pooled$df - expected$df, pooled$lower - (expected$estimate - half), pooled$upper - (expected$estimate + half))))
    if (!is.finite(delta) || delta > 1e-10) stop("All-reference complete-data estimate/SE/df/CI mismatch.")
    expected$control_max_abs_difference <- delta
    return(expected)
  }
  pool_fits(lapply(imp$completed, function(d) fit_lab(d$lab_value, p)))
}
run_replicate <- function(s, r) {
  p <- prepare_data(s, r); records <- list(); decomposition <- NULL
  keep <- if (s$family == "distribution") !is.na(p$input$lab_value) & !p$flagged else !is.na(p$input$gold_value)
  recorded <- if (s$family == "distribution") p$input$lab_value else p$input$gold_value
  if (sum(keep) > ncol(p$analysis_x) + 1L) {
    full <- fit_lab(p$actual_truth, p); retained_truth <- p$actual_truth; retained_truth[!keep] <- NA_real_
    recorded[!keep] <- NA_real_; true_kept <- fit_lab(retained_truth, p); recorded_kept <- fit_lab(recorded, p)
    decomposition <- data.frame(family = s$family, scenario = s$scenario, replicate = r,
      oracle_full = full$estimate, oracle_retained = true_kept$estimate, recorded_retained = recorded_kept$estimate,
      selection_difference = true_kept$estimate - full$estimate,
      residual_error_difference = recorded_kept$estimate - true_kept$estimate, retained_n = sum(keep), truth = p$truth)
  }
  for (method in methods_for(s)) {
    warnings <- character()
    result <- tryCatch(withCallingHandlers(run_method(method, p, s), warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning")
    }), error = identity)
    failed <- inherits(result, "error"); status <- if (failed) "error" else "ok"
    message <- if (failed) conditionMessage(result) else ""
    if (failed && s$scenario == "no_reference" && grepl("independently observed reference", message, fixed = TRUE)) status <- "refused"
    row <- data.frame(family = s$family, scenario = s$scenario, replicate = r, data_seed = p$seed,
      data_stream = p$stream, imputation_seed = p$imputation_seed, method = method, truth = p$truth,
      available_n = sum(keep), n = nrow(p$input), estimate = NA_real_, se = NA_real_, df = NA_real_, n_used = NA_real_,
      lower = NA_real_, upper = NA_real_, width = NA_real_, coverage = NA_real_, error = NA_real_,
      control_max_abs_difference = NA_real_, status = status, message = message,
      warnings = paste(unique(warnings), collapse = " | "), stringsAsFactors = FALSE)
    if (!failed) {
      if (!is.finite(result$estimate) || !is.finite(result$se) || result$se <= 0 || is.na(result$df) || result$df <= 0) {
        row$status <- "error"; row$message <- "Nonfinite estimate or degenerate uncertainty."
      } else {
        row$estimate <- result$estimate; row$se <- result$se; row$df <- result$df; row$n_used <- result$n_used
        half <- stats::qt(.975, result$df) * result$se
        row$lower <- result$estimate - half; row$upper <- result$estimate + half; row$width <- 2 * half
        row$coverage <- as.numeric(row$lower <= p$truth && row$upper >= p$truth); row$error <- result$estimate - p$truth
        if (!is.null(result$control_max_abs_difference)) row$control_max_abs_difference <- result$control_max_abs_difference
      }
    }
    records[[length(records) + 1L]] <- row
  }
  list(records = do.call(rbind, records), decomposition = decomposition)
}
records <- list(); decompositions <- list(); start <- Sys.time()
for (s in scenario_list) {
  chunks <- split(seq_len(s$reps), ceiling(seq_len(s$reps) / 200L))
  for (chunk in chunks) {
    batch <- parallel::mclapply(chunk, function(r) run_replicate(s, r), mc.cores = workers, mc.preschedule = TRUE, mc.set.seed = FALSE)
    if (any(vapply(batch, inherits, logical(1), what = "try-error"))) stop("A worker failed; checkpoint is not verified evidence.")
    records <- c(records, lapply(batch, `[[`, "records")); decompositions <- c(decompositions, lapply(batch, `[[`, "decomposition"))
    # RDS checkpoints avoid repeated large CSV formatting costs. Final CSVs are portable.
    saveRDS(list(records = records, decompositions = decompositions), file.path(outdir, "checkpoint.rds"))
    cat(s$family, s$scenario, max(chunk), "/", s$reps, "; elapsed", round(as.numeric(difftime(Sys.time(), start, units = "secs"))), "seconds\n"); flush.console()
  }
}
if (!identical(hashes_before, hash_sources())) stop("Source or plan changed during run. Checkpoint retained; evidence is not verified.")
results <- do.call(rbind, records); decomposition <- do.call(rbind, decompositions)
utils::write.csv(results, file.path(outdir, "replicates.csv"), row.names = FALSE)
utils::write.csv(decomposition, file.path(outdir, "bias-decomposition.csv"), row.names = FALSE)
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
    set.seed(9900000L + i)
    b <- replicate(config$centered_bootstrap_reps, {
      z <- good[sample.int(n, n, replace = TRUE), ]; shift <- mean(z$error)
      c(centered = mean(z$lower - shift <= z$truth & z$upper - shift >= z$truth), empirical_sd = stats::sd(z$estimate))
    })
    centered_se <- stats::sd(b[1, ]); sd_mcse <- stats::sd(b[2, ])
  }
  yielded <- sum(good$coverage) / requested; yield_ci <- wilson(sum(good$coverage), requested)
  primary <- (g$family[1] == "distribution" && g$scenario[1] %in% c("missing_only", "perfect_flags", "false_positive_flags", "category_error", "null_effect", "strong_mar") && g$method[1] %in% c("normal_no_cat", "mice_normal_no_cat")) ||
    (g$family[1] == "reference" && g$method[1] %in% c("reference_normal", "mice_reference_normal"))
  data.frame(family = g$family[1], scenario = g$scenario[1], method = g$method[1], requested = requested,
    returned_intervals = n, refused = sum(g$status == "refused"), other_errors = sum(g$status == "error"), interval_return_rate = n / requested,
    mean_available_n = mean(g$available_n), bias = bias, bias_mcse = mcse(good$error), rmse = rmse,
    rmse_mcse = if (n > 1 && rmse > 0) mcse(good$error^2) / (2 * rmse) else NA_real_,
    empirical_sd = if (n > 1) stats::sd(good$estimate) else NA_real_, empirical_sd_bootstrap_mcse = sd_mcse,
    mean_se = if (n) mean(good$se) else NA_real_, mean_se_mcse = mcse(good$se), mean_width = if (n) mean(good$width) else NA_real_, width_mcse = mcse(good$width),
    coverage_conditional = coverage, coverage_mcse = if (n) sqrt(coverage * (1 - coverage) / n) else NA_real_,
    coverage_wilson_lower = ci[1], coverage_wilson_upper = ci[2], covering_interval_per_request = yielded,
    covering_yield_mcse = sqrt(yielded * (1 - yielded) / requested), covering_yield_wilson_lower = yield_ci[1], covering_yield_wilson_upper = yield_ci[2],
    bias_centered_coverage_diagnostic = centered, centered_bootstrap_mcse = centered_se,
    primary_criterion_applicable = primary && requested >= 1000L,
    approximate_nominal_criterion = if (primary && n >= 1000L) ci[1] >= .925 && ci[2] <= .975 else NA)
}))
utils::write.csv(summary_table, file.path(outdir, "summary.csv"), row.names = FALSE)
status_rows <- do.call(rbind, lapply(groups, function(g) do.call(rbind, lapply(c("ok", "refused", "error"), function(status) {
  k <- sum(g$status == status); n <- nrow(g); ci <- wilson(k, n); rate <- k / n
  data.frame(family = g$family[1], scenario = g$scenario[1], method = g$method[1], status = status,
    count = k, requested = n, rate = rate, mcse = sqrt(rate * (1 - rate) / n), lower = ci[1], upper = ci[2])
}))))
utils::write.csv(status_rows, file.path(outdir, "request-status.csv"), row.names = FALSE)
paired <- list()
for (family in unique(results$family)) for (scenario in unique(results$scenario[results$family == family])) {
  g <- results[results$family == family & results$scenario == scenario & results$status == "ok", ]
  comparisons <- if (family == "distribution") list(c("pmm_no_cat", "normal_no_cat"), c("pmm_cat", "normal_cat"),
    c("normal_cat", "normal_no_cat"), c("pmm_cat", "pmm_no_cat"), c("normal_no_cat", "mice_normal_no_cat"),
    c("normal_cat", "mice_normal_cat"), c("pmm_no_cat", "mice_pmm_no_cat"), c("pmm_cat", "mice_pmm_cat")) else
    list(c("reference_normal", "mice_reference_normal"), c("reference_pmm", "mice_reference_pmm"),
      c("reference_pmm", "reference_normal"), c("reference_normal", "gold_only"), c("reference_pmm", "gold_only"),
      c("reference_normal", "v02_mask_normal"), c("reference_pmm", "v02_mask_normal"))
  for (pair in comparisons) {
    if (!all(pair %in% g$method)) next
    both <- merge(g[g$method == pair[1], ], g[g$method == pair[2], ], by = "replicate", suffixes = c("_new", "_baseline"))
    de <- both$error_new - both$error_baseline; dm <- both$error_new^2 - both$error_baseline^2
    dc <- both$coverage_new - both$coverage_baseline; dw <- both$width_new - both$width_baseline
    paired[[length(paired) + 1L]] <- data.frame(family = family, scenario = scenario, method = pair[1], comparator = pair[2],
      complete_pairs = nrow(both), error_difference = mean(de), error_mcse = mcse(de), mse_difference = mean(dm), mse_mcse = mcse(dm),
      coverage_difference = mean(dc), coverage_difference_mcse = mcse(dc), width_difference = mean(dw), width_difference_mcse = mcse(dw))
  }
}
utils::write.csv(if (length(paired)) do.call(rbind, paired) else data.frame(), file.path(outdir, "paired-comparisons.csv"), row.names = FALSE)
ds <- split(decomposition, interaction(decomposition$family, decomposition$scenario, drop = TRUE))
dsum <- do.call(rbind, lapply(ds, function(g) data.frame(family = g$family[1], scenario = g$scenario[1], replicates = nrow(g),
  selection_difference = mean(g$selection_difference), selection_mcse = mcse(g$selection_difference),
  residual_error_difference = mean(g$residual_error_difference), residual_error_mcse = mcse(g$residual_error_difference))))
utils::write.csv(dsum, file.path(outdir, "bias-decomposition-summary.csv"), row.names = FALSE)
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
writeLines(c(paste("# v0.3 extension evaluation:", mode), "",
  "Prespecified settings, assumptions and interpretation boundaries are in ../EXTENSION_PLAN.md. Coverage is conditional on returned intervals; request-status.csv retains failures. Same-reference results across recording scenarios reuse cohorts and reference subsets, so they are not independent replications.", "",
  "| Setting | Method | Intervals/requested | Bias (MCSE) | Coverage (MCSE) [95% Wilson] | Width (MCSE) |", "|---|---|---:|---:|---:|---:|", rows, "",
  "Normal/PMM comparisons must hold category information constant. Reference workflows obtain additional measured truth; they do not identify missed errors without that information. Reference-only and independent mice comparators use the same gold measurements. Nonignorable selection, noisy gold and sparse donor support are deliberate failure probes.", "",
  "Summary, paired contrasts, bias decomposition and all requests are in adjacent CSVs. Bias-centered coverage is retrospective and uses hidden truth; it is not a real-data correction. Source hashes match before and after execution. No result is made valid merely by this script completing."), file.path(outdir, "RESULTS.md"))
unlink(file.path(outdir, "checkpoint.rds"))
cat("Completed", mode, ":", nrow(results), "requests;", paste(names(table(results$status)), table(results$status), collapse = ", "), "\n")
cat("Output:", normalizePath(outdir), "\n")
