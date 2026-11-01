#!/usr/bin/env Rscript
# The frozen study runner kept scores and hold-out checks in a variable that later
# methods overwrote, so the six cells that also ran reference PMM lost those two
# auxiliary records. This script recovers them with the frozen seeds and refuses to
# write anything unless every contamination estimate equals the main study's exactly.
# Usage: Rscript tools/contamination-supplement.R [NEW_DIRECTORY]
args <- commandArgs(trailingOnly = TRUE)
file_arg <- grep("^--file=", commandArgs(), value = TRUE)
root <- normalizePath(file.path(dirname(sub("^--file=", "", file_arg[1])), ".."))
outdir <- if (length(args)) args[1] else file.path(root, "validation", "contamination-supplement")
if (dir.exists(outdir) && length(list.files(outdir, all.files = TRUE, no.. = TRUE)))
  stop("Use a new or empty output directory; existing evidence is never overwritten.")
pkgload::load_all(root, quiet = TRUE)
source(file.path(root, "tools/lab-fixture.R"), local = TRUE)
source(file.path(root, "tools/contamination-fixtures.R"), local = TRUE)
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
workers <- suppressWarnings(as.integer(Sys.getenv("CATMECH_BENCHMARK_WORKERS", "14")))
study <- file.path(root, "validation", "contamination-study")
settings <- jsonlite::read_json(file.path(study, "settings.json"))
main <- utils::read.csv(file.path(study, "replicates.csv"), stringsAsFactors = FALSE)
main <- main[main$method == "contamination", ]
hashes <- utils::read.csv(file.path(study, "source-hashes.csv"), stringsAsFactors = FALSE)
for (i in seq_len(nrow(hashes)))
  if (digest::digest(file = file.path(root, hashes$file[i]), algo = "sha256") != hashes$sha256[i])
    stop("Current source differs from the study's frozen source: ", hashes$file[i])
missing_cells <- c("missed_flags", "reference_05pct", "reference_10pct", "reference_50pct", "heavy_tail_target", "nonlinear_target")
scenarios <- contamination_scenarios(settings$primary_reps, settings$stress_reps)
scenarios <- scenarios[vapply(scenarios, function(s) s$scenario %in% missing_cells, logical(1))]
fit_slope <- function(lab, p) {
  x <- cbind(p$analysis_x, lab_value = lab)
  stats::lm.fit(x, p$input$outcome)$coefficients[[ncol(x)]]
}
run_replicate <- function(s, r) {
  p <- make_contamination_case(s, r, settings$phase, settings$fixture_master_seed, make_lab_fixture)
  fit <- tryCatch(suppressWarnings(catmech_contamination(p$input, "lab_value", "gold_value", p$predictors,
    corruption_strata = s$corruption_strata, error_strata = s$error_strata, agreement_tol = s$agreement_tol,
    reference_sd = s$reference_sd, m = settings$m, seed = p$imputation_seed)), error = identity)
  saved <- main[main$scenario == s$scenario & main$replicate == r, ]
  if (nrow(saved) != 1L) stop("Missing main-study record.")
  if (inherits(fit, "error")) {
    if (saved$status == "ok") stop("A returned main-study interval did not reproduce: ", s$scenario, " ", r)
    return(NULL)
  }
  if (saved$status != "ok") stop("A refused main-study request returned a fit: ", s$scenario, " ", r)
  q <- vapply(fit$completed, function(d) fit_slope(d$lab_value, p), numeric(1))
  if (abs(mean(q) - saved$estimate) > 1e-9) stop("Estimate did not reproduce: ", s$scenario, " ", r)
  table <- fit$contamination$corruption
  scored <- which(table$evidence == "posterior_model" & !is.na(table$posterior_corruption))
  y <- p$corrupted[scored]; probability <- table$posterior_corruption[scored]
  ranks <- rank(probability)
  positives <- sum(y); negatives <- sum(!y)
  score <- data.frame(family = s$family, scenario = s$scenario, replicate = r, rows = length(scored),
    actual_errors = positives, expected_errors = sum(probability), brier = mean((probability - y)^2),
    brier_prevalence = mean((mean(y) - y)^2),
    auc = if (positives && negatives) (sum(ranks[y]) - positives * (positives + 1) / 2) / (positives * negatives) else NA_real_,
    mean_probability_when_wrong = if (positives) mean(probability[y]) else NA_real_,
    mean_probability_when_correct = if (negatives) mean(probability[!y]) else NA_real_,
    completion_rmse = sqrt(mean((table$posterior_target - p$truth)^2)),
    recorded_rmse = sqrt(mean((p$input$lab_value[scored] - p$truth[scored])^2)))
  holdout <- NULL
  if (isTRUE(s$holdout) && r <= settings$holdout_reps) {
    h <- catmech_holdout_check(fit, folds = settings$holdout_folds, m = settings$holdout_m, seed = p$imputation_seed)
    row <- if (is.null(h$summary)) data.frame(predictive_coverage = NA_real_, ks_p_value = NA_real_, calibration_z = NA_real_,
      errors_observed = NA_real_, errors_expected = NA_real_) else
      h$summary[c("predictive_coverage", "ks_p_value", "calibration_z", "errors_observed", "errors_expected")]
    holdout <- cbind(data.frame(family = s$family, scenario = s$scenario, replicate = r, status = h$status,
      flags = paste(h$flags, collapse = "|"), failures = length(h$failures), stringsAsFactors = FALSE), row)
  }
  list(score = score, holdout = holdout)
}
start <- Sys.time()
results <- list()
for (s in scenarios) {
  batch <- parallel::mclapply(seq_len(s$reps), function(r) run_replicate(s, r), mc.cores = workers, mc.preschedule = TRUE, mc.set.seed = FALSE)
  failed <- vapply(batch, inherits, logical(1), what = "try-error")
  if (any(failed)) stop(as.character(batch[[which(failed)[1]]]))
  results <- c(results, batch)
  cat(s$scenario, "done; elapsed", round(as.numeric(difftime(Sys.time(), start, units = "secs"))), "seconds\n"); flush.console()
}
score_table <- do.call(rbind, lapply(results, `[[`, "score"))
holdout_table <- do.call(rbind, lapply(results, `[[`, "holdout"))
utils::write.csv(score_table, file.path(outdir, "corruption-scores.csv"), row.names = FALSE)
utils::write.csv(holdout_table, file.path(outdir, "holdout-checks.csv"), row.names = FALSE)
mcse <- function(x) if (length(x) > 1L) stats::sd(x) / sqrt(length(x)) else NA_real_
utils::write.csv(do.call(rbind, lapply(split(score_table, score_table$scenario), function(g) data.frame(family = g$family[1], scenario = g$scenario[1],
  replicates = nrow(g), mean_actual_errors = mean(g$actual_errors), mean_expected_errors = mean(g$expected_errors),
  mean_error_count_gap = mean(g$expected_errors - g$actual_errors), error_count_gap_mcse = mcse(g$expected_errors - g$actual_errors),
  brier = mean(g$brier), brier_prevalence = mean(g$brier_prevalence), auc = mean(g$auc, na.rm = TRUE),
  mean_probability_when_wrong = mean(g$mean_probability_when_wrong, na.rm = TRUE),
  mean_probability_when_correct = mean(g$mean_probability_when_correct, na.rm = TRUE),
  completion_rmse = mean(g$completion_rmse), recorded_rmse = mean(g$recorded_rmse)))),
  file.path(outdir, "corruption-scores-summary.csv"), row.names = FALSE)
utils::write.csv(do.call(rbind, lapply(split(holdout_table, holdout_table$scenario), function(g) data.frame(family = g$family[1], scenario = g$scenario[1], replicates = nrow(g),
  detected = sum(g$status == "miscalibration_detected"), unresolved = sum(g$status == "unresolved_failed_folds"),
  detection_rate = mean(g$status == "miscalibration_detected"), mean_coverage = mean(g$predictive_coverage, na.rm = TRUE)))),
  file.path(outdir, "holdout-summary.csv"), row.names = FALSE)
writeLines(c(paste("Recovered", nrow(score_table), "corruption-score records and", nrow(holdout_table), "hold-out records."),
  "Every contamination estimate matched the frozen main study to within 1e-9, and every refused request was refused again.",
  paste("Elapsed seconds:", round(as.numeric(difftime(Sys.time(), start, units = "secs")), 1))), file.path(outdir, "README.md"))
cat(readLines(file.path(outdir, "README.md")), sep = "\n")
