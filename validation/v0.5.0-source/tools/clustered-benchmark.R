#!/usr/bin/env Rscript
# Prospective test of imputation for clustered audits; settings locked in validation/CLUSTERED_PLAN.md.
# Every interval comes from the package's own catmech_reference() and catmech_pool_lm().
# Usage: Rscript tools/clustered-benchmark.R pilot|study [NEW_DIRECTORY] [PILOT_REPS]
args <- commandArgs(trailingOnly = TRUE)
mode <- if (length(args)) match.arg(args[1], c("pilot", "study")) else "pilot"
file_arg <- grep("^--file=", commandArgs(), value = TRUE)
root <- normalizePath(file.path(dirname(sub("^--file=", "", file_arg[1])), ".."))
outdir <- if (length(args) >= 2L) args[2] else file.path(root, "validation", paste0("clustered-", mode))
if (dir.exists(outdir) && length(list.files(outdir, all.files = TRUE, no.. = TRUE)))
  stop("Use a new or empty output directory; existing evidence is never overwritten.")
hash_sources <- function() {
  paths <- c(sort(list.files(file.path(root, "R"), "\\.R$", full.names = TRUE)),
    file.path(root, c("DESCRIPTION", "NAMESPACE", "tools/clustered-benchmark.R", "validation/CLUSTERED_PLAN.md", "validation/real-audit/rwanda-pool.csv")))
  data.frame(file = substring(paths, nchar(root) + 2L),
    sha256 = vapply(paths, function(p) digest::digest(file = p, algo = "sha256"), character(1)), row.names = NULL)
}
hashes_before <- hash_sources()
pkgload::load_all(root, quiet = TRUE)
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
workers <- suppressWarnings(as.integer(Sys.getenv("CATMECH_BENCHMARK_WORKERS", "14")))
if (length(workers) != 1L || is.na(workers) || workers < 1L || workers > 16L) stop("Invalid worker count (1 to 16).")
pilot_reps <- if (length(args) >= 3L) as.integer(args[3]) else 3L
config <- list(mode = mode, workers = workers, m = 20L, fractions = c(.25, .50), reps = if (mode == "pilot") pilot_reps else 500L, ratio_resamples = 2000L, ratio_seed = 930000000L,
  held_out = list(cohort = if (mode == "pilot") 960100000 else 860100000, selection = if (mode == "pilot") 970100000 else 870100000, imputation = if (mode == "pilot") 980000000 else 880000000),
  development = list(cohort = 400100000, selection = 410100000, imputation = 420000000),
  simulation = list(cohort = if (mode == "pilot") 990100000 else 900100000, imputation = if (mode == "pilot") 1000000000 else 910000000),
  simulation_parameters = list(b0 = 2, b1 = .6, size_log_mean = 3, size_log_sd = .5, x_noise_sd = 3, residual_variance = 100, wrong_rate = .4, rate_beta = c(1, 1.5), error_sd = 6, tracking_slope = .25))
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
          Map(real_analysis, names(DEVELOPMENT), DEVELOPMENT, seq_along(DEVELOPMENT), "development", list(config$development)))

simulate_cell <- function(cell, seed, fraction) {
  set.seed(seed)
  a <- config$simulation_parameters
  g <- cell$clusters; rows <- cell$rows
  id <- rep(seq_len(g), each = rows)
  size <- rep(stats::rlnorm(g, a$size_log_mean, a$size_log_sd), each = rows)
  x <- size + stats::rnorm(g * rows, 0, a$x_noise_sd)
  u <- rep(stats::rnorm(g, 0, sqrt(cell$icc * a$residual_variance)), each = rows)
  truth <- a$b0 + a$b1 * x + u + stats::rnorm(g * rows, 0, sqrt((1 - cell$icc) * a$residual_variance))
  rate <- if (cell$errors == "iid") rep(a$wrong_rate, g * rows) else rep(stats::rbeta(g, a$rate_beta[1], a$rate_beta[2]), each = rows)
  wrong <- stats::runif(g * rows) < rate
  shift <- if (cell$errors == "tracking") a$tracking_slope * (x - mean(x)) else 0
  recorded <- truth + wrong * (shift + stats::rnorm(g * rows, 0, a$error_sd))
  chosen <- seq_len(g) %in% sample.int(g, max(1L, as.integer(round(fraction * g))))
  audited <- chosen[id]
  gold <- truth; gold[!audited] <- NA_real_
  list(input = data.frame(covariate = x, cluster_id = id, lab_value = recorded, gold_value = gold), truth = truth, target = a$b1, formula = lab_value ~ covariate,
       predictors = "covariate", strata = NULL, audited = audited)
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
methods <- c("oracle", "recorded", "linear_mi", "linear_mi_clustered")
run_method <- function(method, p, seed) {
  input <- p$input
  if (method == "oracle") input$gold_value <- p$truth
  if (method == "recorded") input$gold_value <- input$lab_value
  fit <- if (method == "linear_mi_clustered") {
    catmech_reference(input, "lab_value", "gold_value", p$predictors, method = "normal", m = config$m, seed = seed, use_recorded = TRUE, cluster = "cluster_id", cluster_strata = p$strata)
  } else catmech_reference(input, "lab_value", "gold_value", p$predictors, method = "normal", m = config$m, seed = seed, use_recorded = TRUE)
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
    status <- if (!failed) "ok" else if (grepl("Fewer than five|lost design support|full column rank|more observations than coefficients|constant predictor", conditionMessage(result))) "refused" else "error"
    row <- data.frame(level = s$level, analysis = s$analysis, scenario = s$scenario, fraction = s$fraction, replicate = r, imputation_seed = seed, method = method, truth = p$target,
      audited_rows = sum(p$audited), audited_clusters = length(unique(p$input$cluster_id[p$audited])), n = nrow(p$input), estimate = NA_real_, se = NA_real_, df = NA_real_, lower = NA_real_, upper = NA_real_,
      width = NA_real_, coverage = NA_real_, error = NA_real_, recovery_error = NA_real_, status = status, message = if (failed) conditionMessage(result) else "",
      warnings = paste(unique(warnings), collapse = " | "), stringsAsFactors = FALSE)
    if (!failed) {
      if (method == "oracle") oracle_estimate <- result$estimate
      row$estimate <- result$estimate; row$se <- result$std.error; row$df <- result$df; row$lower <- result$lower; row$upper <- result$upper
      row$width <- result$upper - result$lower; row$coverage <- as.numeric(result$lower <= p$target && result$upper >= p$target)
      row$error <- result$estimate - p$target; row$recovery_error <- result$estimate - oracle_estimate
    }
    records[[length(records) + 1L]] <- row
  }
  do.call(rbind, records)
}

scenarios <- list()
for (z in real) for (j in seq_along(config$fractions)) scenarios[[length(scenarios) + 1L]] <- list(kind = "real", z = z, fraction_index = j, level = z$level, analysis = z$id, fraction = config$fractions[j],
  scenario = paste0(z$level, ":", z$id, ":audit_", 100 * config$fractions[j], "pct"))
grid <- expand.grid(clusters = c(30L, 70L, 150L), icc = c(0, .2, .5), errors = c("iid", "persistent", "tracking"), fraction = config$fractions, stringsAsFactors = FALSE)
grid$rows <- 3L
grid <- rbind(grid, data.frame(clusters = 210L, icc = 0, errors = "iid", fraction = config$fractions, rows = 1L))
for (i in seq_len(nrow(grid))) {
  label <- with(grid[i, ], if (rows == 1L) "control_one_row" else sprintf("g%d_icc%02d_%s", clusters, round(100 * icc), errors))
  scenarios[[length(scenarios) + 1L]] <- list(kind = "simulation", cell = as.list(grid[i, ]), index = i, level = "simulation", analysis = label, fraction = grid$fraction[i],
    scenario = paste0("simulation:", label, ":audit_", 100 * grid$fraction[i], "pct"))
}
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
utils::write.csv(results, file.path(outdir, "replicates.csv"), row.names = FALSE)
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
    empirical_sd = if (n > 1) stats::sd(good$estimate) else NA_real_, mean_se = if (n) mean(good$se) else NA_real_, mean_width = if (n) mean(good$width) else NA_real_,
    coverage = if (n) mean(good$coverage) else NA_real_, coverage_wilson_lower = ci[1], coverage_wilson_upper = ci[2])
}))
rownames(summary_table) <- NULL
utils::write.csv(summary_table, file.path(outdir, "summary.csv"), row.names = FALSE)
set.seed(config$ratio_seed)
ratio_rows <- list()
for (s in scenarios) {
  g <- results[results$scenario == s$scenario & results$status == "ok", ]
  both <- merge(g[g$method == "linear_mi_clustered", c("replicate", "recovery_error", "error", "coverage")], g[g$method == "linear_mi", c("replicate", "recovery_error", "error", "coverage")], by = "replicate", suffixes = c("_a", "_b"))
  if (!nrow(both)) next
  draw <- function(a, b) {
    resampled <- replicate(config$ratio_resamples, { i <- sample.int(length(a), replace = TRUE); sqrt(sum(a[i]^2) / sum(b[i]^2)) })
    c(sqrt(sum(a^2) / sum(b^2)), stats::quantile(resampled, c(.025, .975), names = FALSE))
  }
  recovery <- draw(both$recovery_error_a, both$recovery_error_b); pool <- draw(both$error_a, both$error_b)
  ratio_rows[[length(ratio_rows) + 1L]] <- data.frame(level = s$level, analysis = s$analysis, scenario = s$scenario, fraction = s$fraction, complete_pairs = nrow(both),
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
