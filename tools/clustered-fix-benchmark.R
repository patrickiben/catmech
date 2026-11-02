#!/usr/bin/env Rscript
# Prospective test of replacements for the resampling in clustered imputation; settings locked in validation/CLUSTERED_FIX_PLAN.md.
# Every interval comes from the package's own catmech_reference() and catmech_pool_lm(). Candidates replace only catmech_reference_cluster_fill().
# Usage: Rscript tools/clustered-fix-benchmark.R pilot|study [NEW_DIRECTORY] [PILOT_REPS]
args <- commandArgs(trailingOnly = TRUE)
mode <- if (length(args)) match.arg(args[1], c("pilot", "study")) else "pilot"
file_arg <- grep("^--file=", commandArgs(), value = TRUE)
root <- normalizePath(file.path(dirname(sub("^--file=", "", file_arg[1])), ".."))
outdir <- if (length(args) >= 2L) args[2] else file.path(root, "validation", paste0("clustered-fix-", mode))
if (dir.exists(outdir) && length(list.files(outdir, all.files = TRUE, no.. = TRUE)))
  stop("Use a new or empty output directory; existing evidence is never overwritten.")
hash_sources <- function() {
  paths <- c(sort(list.files(file.path(root, "R"), "\\.R$", full.names = TRUE)),
    file.path(root, c("DESCRIPTION", "NAMESPACE", "tools/clustered-fix-benchmark.R", "tools/clustered-fix-candidates.R", "tools/clustered-fix-score.py", "validation/CLUSTERED_FIX_PLAN.md", "validation/real-audit/rwanda-pool.csv")))
  data.frame(file = substring(paths, nchar(root) + 2L),
    sha256 = vapply(paths, function(p) digest::digest(file = p, algo = "sha256"), character(1)), row.names = NULL)
}
hashes_before <- hash_sources()
pkgload::load_all(root, quiet = TRUE)
source(file.path(root, "tools/clustered-fix-candidates.R"), local = TRUE)
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
workers <- suppressWarnings(as.integer(Sys.getenv("CATMECH_BENCHMARK_WORKERS", "14")))
if (length(workers) != 1L || is.na(workers) || workers < 1L || workers > 16L) stop("Invalid worker count (1 to 16).")
pilot_reps <- if (length(args) >= 3L) as.integer(args[3]) else 3L
pilot <- mode == "pilot"
config <- list(mode = mode, workers = workers, m = 20L, fractions = c(.25, .50), reps = if (pilot) pilot_reps else 500L, ratio_resamples = 2000L, ratio_seed = 1500000000L,
  held_out = if (pilot) list(cohort = 1420100000, selection = 1430100000, imputation = 1440000000) else list(cohort = 860100000, selection = 870100000, imputation = 880000000),
  development = if (pilot) list(cohort = 1450100000, selection = 1460100000, imputation = 1470000000) else list(cohort = 400100000, selection = 410100000, imputation = 420000000),
  simulation = if (pilot) list(cohort = 1400100000, imputation = 1410000000) else list(cohort = 1300100000, imputation = 1310000000),
  simulation_parameters = list(b0 = 2, b1 = .6, size_log_mean = 3, size_log_sd = .5, x_noise_sd = 3, cluster_variance = 20, row_variance = 80, stratum_sd = 2))
jsonlite::write_json(config, file.path(outdir, "settings.json"), pretty = TRUE, auto_unbox = TRUE)

rwanda <- utils::read.csv(file.path(root, "validation/real-audit/rwanda-pool.csv"), stringsAsFactors = FALSE)
rwanda <- rwanda[rwanda$facility_type == "health_center", ]
add <- function(...) { parts <- list(...); out <- parts[[1]]; for (p in parts[-1]) out <- out + p; out }
HELD_OUT <- list(
  lbw = list(w = with(rwanda, add(lbw_notpremature_h, lbw_premature_h)), r = rwanda$lbw_r, x = rwanda$hf_deliveries_h),
  resus = list(w = rwanda$newborns_nocryres_succ_h, r = rwanda$nocry_res_success_r, x = rwanda$hf_deliveries_h),
  stillbirths = list(w = with(rwanda, add(sbmacerated_h, sbfresh_h)), r = with(rwanda, add(sbmacerated_r, sbfresh_r)), x = rwanda$hf_deliveries_h),
  deliveries = list(w = rwanda$hf_deliveries_h, r = rwanda$hf_deliveries_r, x = rwanda$ancnewreg_h),
  anc1_early = list(w = rwanda$anc1vis1trim_h, r = rwanda$anc1vis_16weeks_r, x = rwanda$ancnewreg_h))
DEVELOPMENT <- list(
  pnc = list(w = rwanda$pnc_newborn_24h_h, r = rwanda$pnc_newborn_24h_r, x = rwanda$hf_deliveries_h),
  syphilis = list(w = rwanda$anc_syphilis_test_h, r = rwanda$anc_syphilis_test_r, x = rwanda$ancnewreg_h),
  iron = list(w = rwanda$womenreciron_h, r = rwanda$womenreciron_r, x = rwanda$ancnewreg_h),
  anc4 = list(w = rwanda$anc_4visits_h, r = rwanda$anc_4visits_r, x = rwanda$ancnewreg_h))
UNUSED <- list(
  ancnewreg = list(w = rwanda$ancnewreg_h, r = rwanda$ancnewreg_r, x = rwanda$hf_deliveries_h),
  livebirths = list(w = rwanda$livebirths_h, r = rwanda$livebirths_r, x = rwanda$ancnewreg_h))
real_analysis <- function(id, values, index, level, seeds) {
  keep <- stats::complete.cases(values$w, values$r, values$x)
  pool <- data.frame(unit = paste0("f", rwanda$facility_id[keep]), stratum = factor(paste0("c", rwanda$hosp_id[keep])),
                     lab_value = values$w[keep], lab_true = values$r[keep], covariate = values$x[keep], stringsAsFactors = FALSE)
  x <- stats::model.matrix(~ covariate + stratum, data = pool)
  list(id = id, index = index, level = level, seeds = seeds, pool = pool, rows_by_unit = split(seq_len(nrow(pool)), pool$unit),
       units_by_stratum = lapply(split(pool$unit, pool$stratum), unique), truth = unname(stats::lm.fit(x, pool$lab_true)$coefficients["covariate"]),
       recorded_slope = unname(stats::lm.fit(x, pool$lab_value)$coefficients["covariate"]))
}
real <- c(Map(real_analysis, names(HELD_OUT), HELD_OUT, seq_along(HELD_OUT), "held_out", list(config$held_out)),
          Map(real_analysis, names(DEVELOPMENT), DEVELOPMENT, seq_along(DEVELOPMENT), "development", list(config$development)),
          Map(real_analysis, names(UNUSED), UNUSED, length(HELD_OUT) + seq_along(UNUSED), "unused", list(config$held_out)))

DESIGNS <- list(rwanda = c(14L, 9L, 8L, 13L, 5L, 9L, 11L), pairs = rep(8L, 10L), none = 70L, none_210 = 210L)
simulate_cell <- function(cell, seed, fraction) {
  set.seed(seed)
  a <- config$simulation_parameters
  sizes <- DESIGNS[[cell$design]]; rows <- cell$rows; g <- sum(sizes); n <- g * rows
  stratified <- length(sizes) > 1L
  id <- rep(seq_len(g), each = rows)
  stratum_of <- rep(seq_along(sizes), sizes); row_stratum <- stratum_of[id]
  size <- rep(stats::rlnorm(g, a$size_log_mean, a$size_log_sd), each = rows)
  x <- size + stats::rnorm(n, 0, a$x_noise_sd)
  if (cell$errors == "rare") {
    level <- if (stratified) stats::runif(length(sizes), 0, .2)[row_stratum] else 0
    frailty <- rep(stats::rgamma(g, .5, .5), each = rows)
    truth <- as.numeric(stats::rpois(n, pmax(.01, .1 + .02 * x + level) * frailty))
    slip <- stats::runif(n) < .2
    recorded <- pmax(0, truth + slip * sample(c(-1, 1), n, replace = TRUE))
    target <- .02
  } else {
    level <- if (stratified) stats::rnorm(length(sizes), 0, a$stratum_sd)[row_stratum] else 0
    residual <- if (cell$errors == "row_only") stats::rnorm(n, 0, sqrt(a$cluster_variance + a$row_variance)) else
      rep(stats::rnorm(g, 0, sqrt(a$cluster_variance)), each = rows) + stats::rnorm(n, 0, sqrt(a$row_variance))
    truth <- a$b0 + a$b1 * x + level + residual
    recorded <- switch(cell$errors,
      shift_low = truth + rep(stats::rnorm(g, 0, 3), each = rows) + stats::rnorm(n, 0, 6),
      shift = truth + rep(stats::rnorm(g, 0, 5), each = rows) + stats::rnorm(n, 0, 3),
      under = (1 - rep(stats::rbeta(g, 2, 6), each = rows)) * truth + stats::rnorm(n, 0, 2),
      outlier = truth + rep((stats::runif(g) < .1) * stats::rnorm(g, 0, 20), each = rows) + stats::rnorm(n, 0, 3),
      row_only = truth + stats::rnorm(n, 0, 6),
      stop("Unknown errors"))
    target <- a$b1
  }
  chosen <- logical(g); start <- 0L
  for (nh in sizes) {
    chosen[start + sample.int(nh, max(1L, as.integer(round(fraction * nh))))] <- TRUE
    start <- start + nh
  }
  audited <- chosen[id]
  gold <- truth; gold[!audited] <- NA_real_
  input <- data.frame(covariate = x, cluster_id = id, lab_value = recorded, gold_value = gold)
  if (stratified) input$stratum <- factor(paste0("s", row_stratum))
  list(input = input, truth = truth, target = target, audited = audited, strata = if (stratified) "stratum" else NULL,
       formula = if (stratified) lab_value ~ covariate + stratum else lab_value ~ covariate,
       predictors = if (stratified) c("covariate", "stratum") else "covariate")
}
real_case <- function(z, fraction_index, r) {
  set.seed(as.integer(z$seeds$cohort + z$index * 1e6 + r))
  drawn <- lapply(z$units_by_stratum, function(u) u[sample.int(length(u), length(u), replace = TRUE)])
  set.seed(as.integer(z$seeds$selection + z$index * 1e6 + fraction_index * 1e5 + r))
  fraction <- config$fractions[fraction_index]
  picked <- lapply(drawn, function(u) seq_along(u) %in% sample.int(length(u), max(1L, as.integer(round(fraction * length(u))))))
  parts <- list(); position <- 0L
  for (s in names(drawn)) for (j in seq_along(drawn[[s]])) {
    position <- position + 1L
    rows <- z$rows_by_unit[[drawn[[s]][j]]]
    parts[[position]] <- data.frame(row = rows, cluster_id = position, audited = picked[[s]][j])
  }
  index <- do.call(rbind, parts)
  d <- z$pool[index$row, ]; rownames(d) <- NULL
  gold <- d$lab_true; gold[!index$audited] <- NA_real_
  list(input = data.frame(covariate = d$covariate, stratum = d$stratum, cluster_id = index$cluster_id, lab_value = d$lab_value, gold_value = gold),
       truth = d$lab_true, target = z$truth, formula = lab_value ~ covariate + stratum, predictors = c("covariate", "stratum"), strata = "stratum", audited = index$audited)
}

namespace <- asNamespace("catmech")
package_fill <- get("catmech_reference_cluster_fill", envir = namespace)
set_fill <- function(f) {
  unlockBinding("catmech_reference_cluster_fill", namespace)
  assign("catmech_reference_cluster_fill", f, envir = namespace)
  lockBinding("catmech_reference_cluster_fill", namespace)
}
fills <- list(linear_mi_clustered = package_fill, no_strata = package_fill, rao_wu = rao_wu_fill, jackknife = jackknife_fill, random_intercept = random_intercept_fill)
methods <- c("oracle", "linear_mi", names(fills))
run_method <- function(method, p, seed) {
  input <- p$input
  if (method == "oracle") input$gold_value <- p$truth
  fit <- if (method %in% c("oracle", "linear_mi")) {
    catmech_reference(input, "lab_value", "gold_value", p$predictors, method = "normal", m = config$m, seed = seed, use_recorded = TRUE)
  } else {
    set_fill(fills[[method]])
    on.exit(set_fill(package_fill), add = TRUE)
    catmech_reference(input, "lab_value", "gold_value", p$predictors, method = "normal", m = config$m, seed = seed, use_recorded = TRUE,
                      cluster = "cluster_id", cluster_strata = if (method == "no_strata") NULL else p$strata)
  }
  pooled <- catmech_pool_lm(fit, p$formula, cluster = "cluster_id")
  pooled[pooled$term == "covariate", ]
}
run_case <- function(s, p, r, seed) {
  oracle_estimate <- NA_real_; records <- list()
  for (method in methods) {
    warnings <- character()
    result <- tryCatch(withCallingHandlers(run_method(method, p, seed), warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning") }), error = identity)
    failed <- inherits(result, "error")
    status <- if (!failed) "ok" else if (grepl("Fewer than five|lost design support|full column rank|more observations than coefficients|constant predictor|no residual variance|can be left out", conditionMessage(result))) "refused" else "error"
    row <- data.frame(level = s$level, analysis = s$analysis, scenario = s$scenario, fraction = s$fraction, replicate = r, imputation_seed = seed, method = method, truth = p$target,
      audited_rows = sum(p$audited), audited_clusters = length(unique(p$input$cluster_id[p$audited])), n = nrow(p$input), estimate = NA_real_, se = NA_real_, df = NA_real_, lower = NA_real_, upper = NA_real_,
      width = NA_real_, coverage = NA_real_, error = NA_real_, recovery_error = NA_real_, within_variance = NA_real_, between_variance = NA_real_, status = status, message = if (failed) conditionMessage(result) else "",
      warnings = paste(unique(warnings), collapse = " | "), stringsAsFactors = FALSE)
    if (!failed) {
      if (method == "oracle") oracle_estimate <- result$estimate
      row$estimate <- result$estimate; row$se <- result$std.error; row$df <- result$df; row$lower <- result$lower; row$upper <- result$upper
      row$width <- result$upper - result$lower; row$coverage <- as.numeric(result$lower <= p$target && result$upper >= p$target)
      row$error <- result$estimate - p$target; row$recovery_error <- result$estimate - oracle_estimate
      row$within_variance <- result$within_variance; row$between_variance <- result$between_variance
    }
    records[[length(records) + 1L]] <- row
  }
  do.call(rbind, records)
}

scenarios <- list()
for (z in real) for (j in seq_along(config$fractions)) scenarios[[length(scenarios) + 1L]] <- list(kind = "real", z = z, fraction_index = j, level = z$level, analysis = z$id, fraction = config$fractions[j],
  scenario = paste0(z$level, ":", z$id, ":audit_", 100 * config$fractions[j], "pct"))
grid <- expand.grid(design = c("rwanda", "pairs", "none"), errors = c("shift_low", "shift", "under", "outlier", "rare", "row_only"), fraction = config$fractions, stringsAsFactors = FALSE)
grid$rows <- 3L
grid <- rbind(grid, data.frame(design = "none_210", errors = "row_only", fraction = config$fractions, rows = 1L))
grid$label <- ifelse(grid$rows == 1L, "control_one_row", paste0(grid$design, "_", grid$errors))
for (i in seq_len(nrow(grid))) scenarios[[length(scenarios) + 1L]] <- list(kind = "simulation", cell = as.list(grid[i, ]), index = i, level = "simulation", analysis = grid$label[i], fraction = grid$fraction[i],
  scenario = paste0("simulation:", grid$label[i], ":audit_", 100 * grid$fraction[i], "pct"))
run_replicate <- function(s, r) {
  if (s$kind == "real") {
    p <- real_case(s$z, s$fraction_index, r)
    run_case(s, p, r, as.integer(s$z$seeds$imputation + s$z$index * 1e6 + s$fraction_index * 1e5 + r))
  } else {
    p <- simulate_cell(s$cell, as.integer(config$simulation$cohort + s$index * 1e5 + r), s$fraction)
    run_case(s, p, r, as.integer(config$simulation$imputation + s$index * 1e5 + r))
  }
}
pool_rows <- do.call(rbind, lapply(real, function(z) {
  e <- z$pool$lab_value - z$pool$lab_true
  data.frame(level = z$level, analysis = z$id, rows = nrow(z$pool), units = length(z$rows_by_unit), target_coefficient = z$truth, recorded_coefficient = z$recorded_slope, wrong_rows = sum(abs(e) > 1e-8))
}))
utils::write.csv(pool_rows, file.path(outdir, "pool-targets.csv"), row.names = FALSE)
utils::write.csv(grid, file.path(outdir, "simulation-grid.csv"), row.names = FALSE)

records <- list(); start <- Sys.time()
for (s in scenarios) {
  for (chunk in split(seq_len(config$reps), ceiling(seq_len(config$reps) / 140L))) {
    batch <- parallel::mclapply(chunk, function(r) run_replicate(s, r), mc.cores = workers, mc.preschedule = TRUE, mc.set.seed = FALSE)
    if (any(vapply(batch, inherits, logical(1), what = "try-error"))) stop("A worker failed; checkpoint is not verified evidence.")
    records <- c(records, batch)
  }
  saveRDS(records, file.path(outdir, "checkpoint.rds"))
  cat(s$scenario, "; elapsed", round(as.numeric(difftime(Sys.time(), start, units = "secs"))), "seconds\n"); flush.console()
}
if (!identical(hashes_before, hash_sources())) stop("Source, plan or pool changed during run. Checkpoint retained; evidence is not verified.")
results <- do.call(rbind, records)
out <- gzfile(file.path(outdir, "replicates.csv.gz"), "wt")
utils::write.csv(results, out, row.names = FALSE)
close(out)
wilson <- function(k, n) {
  if (!n) return(c(NA_real_, NA_real_))
  p <- k / n; z <- stats::qnorm(.975); denom <- 1 + z^2 / n
  centre <- (p + z^2 / (2 * n)) / denom; half <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / denom
  c(max(0, centre - half), min(1, centre + half))
}
summary_table <- do.call(rbind, lapply(split(results, interaction(results$scenario, results$method, drop = TRUE)), function(g) {
  good <- g[g$status == "ok", ]; n <- nrow(good); ci <- wilson(sum(good$coverage), n)
  data.frame(level = g$level[1], analysis = g$analysis[1], scenario = g$scenario[1], fraction = g$fraction[1], method = g$method[1], requested = nrow(g), returned_intervals = n,
    refused = sum(g$status == "refused"), other_errors = sum(g$status == "error"), warned = sum(nzchar(g$warnings)), mean_audited_clusters = mean(g$audited_clusters), target = g$truth[1],
    recovery_rmse = if (n) sqrt(mean(good$recovery_error^2)) else NA_real_, rmse = if (n) sqrt(mean(good$error^2)) else NA_real_, bias = if (n) mean(good$error) else NA_real_,
    empirical_sd = if (n > 1) stats::sd(good$estimate) else NA_real_, mean_se = if (n) mean(good$se) else NA_real_, mean_width = if (n) mean(good$width) else NA_real_, median_width = if (n) stats::median(good$width) else NA_real_,
    mean_between = if (n) mean(good$between_variance) else NA_real_, median_between = if (n) stats::median(good$between_variance) else NA_real_, median_within = if (n) stats::median(good$within_variance) else NA_real_,
    wide_requests = if (n) sum(good$width > 5 * stats::median(good$width)) else NA_integer_, smallest_df = if (n) min(good$df) else NA_real_,
    coverage = if (n) mean(good$coverage) else NA_real_, coverage_wilson_lower = ci[1], coverage_wilson_upper = ci[2])
}))
rownames(summary_table) <- NULL
utils::write.csv(summary_table, file.path(outdir, "summary.csv"), row.names = FALSE)
set.seed(config$ratio_seed)
ratio_rows <- list()
for (s in scenarios) for (method in names(fills)) {
  g <- results[results$scenario == s$scenario & results$status == "ok", ]
  both <- merge(g[g$method == method, c("replicate", "recovery_error", "error", "coverage")], g[g$method == "linear_mi", c("replicate", "recovery_error", "error", "coverage")], by = "replicate", suffixes = c("_a", "_b"))
  if (!nrow(both)) next
  draw <- function(a, b) {
    resampled <- replicate(config$ratio_resamples, { i <- sample.int(length(a), replace = TRUE); sqrt(sum(a[i]^2) / sum(b[i]^2)) })
    c(sqrt(sum(a^2) / sum(b^2)), stats::quantile(resampled, c(.025, .975), names = FALSE, na.rm = TRUE))
  }
  recovery <- draw(both$recovery_error_a, both$recovery_error_b); pool <- draw(both$error_a, both$error_b)
  ratio_rows[[length(ratio_rows) + 1L]] <- data.frame(level = s$level, analysis = s$analysis, scenario = s$scenario, fraction = s$fraction, method = method, complete_pairs = nrow(both),
    recovery_ratio = recovery[1], recovery_lower = recovery[2], recovery_upper = recovery[3], pool_ratio = pool[1], pool_lower = pool[2], pool_upper = pool[3],
    coverage_difference = mean(both$coverage_a - both$coverage_b), coverage_difference_mcse = stats::sd(both$coverage_a - both$coverage_b) / sqrt(nrow(both)))
}
utils::write.csv(do.call(rbind, ratio_rows), file.path(outdir, "ratios.csv"), row.names = FALSE)
utils::write.csv(hashes_before, file.path(outdir, "source-hashes.csv"), row.names = FALSE)
writeLines(c(paste("Completed UTC:", format(Sys.time(), tz = "UTC", usetz = TRUE)),
  paste("Elapsed seconds:", round(as.numeric(difftime(Sys.time(), start, units = "secs")), 1)), capture.output(sessionInfo())), file.path(outdir, "session-info.txt"))
unlink(file.path(outdir, "checkpoint.rds"))
cat("Completed", mode, ":", nrow(results), "requests;", paste(names(table(results$status)), table(results$status), collapse = ", "), "\n")
cat("Output:", normalizePath(outdir), "\n")
