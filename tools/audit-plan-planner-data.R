#!/usr/bin/env Rscript
# Bundles the audit-plan study's numbers and exact test vectors for a browser page, and compares
# pool-level k with the frozen study's widths and RMSE. Reads saved CSV files and runs no new study
# cohort; the vectors use their own small synthetic data.
# Usage: Rscript tools/audit-plan-planner-data.R [STUDY_DIRECTORY]
args <- commandArgs(trailingOnly = TRUE)
file_arg <- grep("^--file=", commandArgs(), value = TRUE)
root <- normalizePath(file.path(dirname(sub("^--file=", "", file_arg[1])), ".."))
folder <- if (length(args)) args[1] else file.path(root, "validation", "audit-plan-study")
output <- file.path(folder, "planner-data.json")
comparison_output <- file.path(folder, "frozen-comparison.csv")
if (file.exists(output) || file.exists(comparison_output)) stop("Output already exists; existing evidence is never overwritten.")
pkgload::load_all(root, quiet = TRUE)
read <- function(name) utils::read.csv(file.path(folder, name), stringsAsFactors = FALSE)
summary_table <- read("summary.csv")
k_table <- read("k-summary.csv")
pool <- read("pool-k.csv")
frozen <- utils::read.csv(file.path(root, "validation/real-data-study/summary.csv"), stringsAsFactors = FALSE)

frozen_block <- function(name) {
  scenarios <- c("0.05" = "audit_05pct", "0.1" = "audit_10pct", "0.25" = "audit_25pct", "0.5" = "audit_50pct")
  entries <- list()
  for (fraction in names(scenarios)) {
    width <- function(method) frozen$mean_width[frozen$family == name & frozen$scenario == scenarios[[fraction]] & frozen$method == method]
    if (length(width("cole_mi")))
      entries[[length(entries) + 1L]] <- list(fraction = as.numeric(fraction), linear_mi_mean_width = width("cole_mi"),
        reference_normal_mean_width = width("reference_normal"), oracle_mean_width = width("oracle"),
        cohorts = frozen$requested[frozen$family == name & frozen$scenario == scenarios[[fraction]] & frozen$method == "oracle"])
  }
  list(source = "validation/real-data-study/summary.csv",
       note = "Mean interval widths from the frozen study, whose audits were random Bernoulli samples. Linear MI there is the method labeled cole_mi.",
       by_fraction = entries)
}
spread_block <- function(name) {
  rows <- k_table[k_table$cell == name, ]
  lapply(split(rows, rows$pilot_records), function(r) as.list(r[c("pilot_records", "plans_returned", "k_mean", "k_median", "k_q05", "k_q25", "k_q75", "k_q95",
    "interval_covers_pool", "median_interval_lower", "median_interval_upper", "discard_k_median", "discard_k_q05", "discard_k_q95")]))
}
accuracy_block <- function(name) {
  rows <- summary_table[summary_table$cell == name, c("pilot_records", "target_fraction", "compared", "median_error", "p90_error", "interval_coverage",
                                                     "predictive_coverage", "floor_median_error", "retest_median_error")]
  unname(lapply(split(rows, seq_len(nrow(rows))), as.list))
}
labels <- c(ctgov = "ClinicalTrials.gov: enrollment against participant-flow Started, log scale; slope of log duration on log size",
            nhanes = "NHANES 2017-18: self-reported against measured weight, kg; slope of systolic pressure on weight")
datasets <- lapply(stats::setNames(pool$dataset, pool$dataset), function(name) {
  row <- pool[pool$dataset == name, ]
  list(label = labels[[name]], pool_rows = row$pool_rows, pool_k = row$k, pool_discard_k = row$discard_k, pool_classical = row$classical,
       pool_discard_classical = row$discard_classical, target_coefficient = row$target_coefficient,
       frozen_study = frozen_block(name), pilot_estimates = spread_block(name), planner_accuracy = accuracy_block(name))
})

# What the frozen study's widths and RMSE imply for k, against the pool values from this study.
comparison <- do.call(rbind, lapply(pool$dataset, function(name) {
  do.call(rbind, lapply(list(c(.05, "audit_05pct"), c(.10, "audit_10pct"), c(.25, "audit_25pct"), c(.50, "audit_50pct")), function(cell) {
    fraction <- as.numeric(cell[1])
    field <- function(method, column) frozen[[column]][frozen$family == name & frozen$scenario == cell[2] & frozen$method == method]
    if (!length(field("cole_mi", "mean_width"))) return(NULL)
    odds <- (1 - fraction) / fraction
    row <- pool[pool$dataset == name, ]
    width_ratio <- field("cole_mi", "mean_width") / field("oracle", "mean_width")
    discard_ratio <- field("reference_normal", "mean_width") / field("oracle", "mean_width")
    rmse_ratio <- field("cole_mi", "rmse") / field("oracle", "rmse")
    data.frame(dataset = name, fraction = fraction, cohorts = field("oracle", "requested"), linear_mi_mean_width = field("cole_mi", "mean_width"),
      reference_normal_mean_width = field("reference_normal", "mean_width"), oracle_mean_width = field("oracle", "mean_width"),
      linear_mi_rmse = field("cole_mi", "rmse"), oracle_rmse = field("oracle", "rmse"),
      k_implied_by_width = (width_ratio^2 - 1) / odds, k_implied_by_rmse = (rmse_ratio^2 - 1) / odds, discard_k_implied_by_width = (discard_ratio^2 - 1) / odds,
      pool_k = row$k, pool_classical = row$classical, pool_discard_k = row$discard_k,
      frozen_width_ratio = width_ratio, pool_k_width_ratio = sqrt(1 + row$k * odds), pool_classical_width_ratio = sqrt(1 + row$classical * odds),
      frozen_discard_ratio = discard_ratio, pool_discard_ratio = sqrt(1 + row$discard_k * odds),
      discard_understatement = 1 - sqrt(1 + row$discard_k * odds) / discard_ratio)
  }))
}))
utils::write.csv(comparison, comparison_output, row.names = FALSE)

# Test vectors: small numeric pilots, written with four decimals so the inputs are exact in JSON.
vector <- function(name, description, seed, n, pilot, make_recorded) {
  restore <- preserve_rng()
  on.exit(restore(), add = TRUE)
  set.seed(seed)
  x1 <- round(stats::rnorm(n), 4)
  x2 <- as.numeric(stats::runif(n) < .4)
  truth <- round(1 + .5 * x1 + .3 * x2 + stats::rnorm(n, 0, .8), 4)
  outcome <- round(2 + .7 * truth + .2 * x1 - .4 * x2 + stats::rnorm(n, 0, .6), 4)
  recorded <- round(make_recorded(truth), 4)
  reference <- rep(NA_real_, n)
  rows <- sort(sample.int(n, pilot))
  reference[rows] <- truth[rows]
  frame <- data.frame(outcome = outcome, recorded = recorded, reference = reference, x1 = x1, x2 = x2)
  fractions <- c(.05, .10, .25, .50, .75)
  plan <- suppressWarnings(catmech_audit_plan(frame, "recorded", "reference", c("x1", "x2", "outcome"), outcome ~ recorded + x1 + x2, "recorded",
                                              fractions = fractions, bootstrap = 500L, m = 20L, seed = seed))
  list(name = name, description = description, formula = "outcome ~ recorded + x1 + x2", term = "recorded", target_ratio = 1.25,
       outcome = outcome, recorded = recorded, reference = reference, covariates = cbind(x1, x2), r_seed = seed,
       expected = list(k = plan$k$estimate[1], discard_k = plan$k$estimate[2], classical = plan$k$classical[1],
         table = unname(lapply(split(plan$table[c("fraction", "records", "additional_records", "relative_width", "discard_relative_width")], seq_along(fractions)), as.list)),
         answer = as.list(plan$answer[c("target_ratio", "fraction", "records", "smallest_grid_fraction")]),
         flags = I(plan$flags),
         r_bootstrap_only = list(note = "These depend on R's random stream, so a port should not expect them to match exactly.",
           k_lower = plan$k$lower[1], k_upper = plan$k$upper[1], discard_k_lower = plan$k$lower[2], discard_k_upper = plan$k$upper[2])))
}
vectors <- list(
  vector("noisy_recording", "Every recorded value has small normal noise.", 810001L, 60L, 30L, function(t) t + stats::rnorm(length(t), 0, .3)),
  vector("exact_or_different", "Most recorded values are exact; about 20% are off by a moderate amount.", 810002L, 80L, 40L,
         function(t) t + (stats::runif(length(t)) < .2) * stats::rnorm(length(t), .2, .5)),
  vector("missing_recorded", "Noisy recording, with about 15% of recorded values missing.", 810003L, 70L, 45L, function(t) {
    w <- t + stats::rnorm(length(t), 0, .4)
    w[stats::runif(length(t)) < .15] <- NA_real_
    w
  }))

k_calculation <- c(
  "Inputs: n records with outcome O, numeric covariates Z, recorded value W (may be missing) and reference Y on the pilot only. The analysis is O ~ y + Z with term y, so its design row is D(y) = [1, y, Z]. p is its column count.",
  "Step 1, imputation with W: regress Y on G = [1, Z, O, W] by least squares over pilot rows with W. tau2 = RSS / (rows - columns). If every such row has |W - Y| <= 1e-8 * max(1, |Y|), use mean W and tau2 = 0 instead.",
  "Step 2, imputation without W: regress Y on G0 = [1, Z, O] over all pilot rows the same way. Records with W use the step 1 model; records without W use this one. When discarding recorded values, every record uses this one.",
  "Step 3, for each record: mu is its fitted value from its model and t is its model's tau2.",
  "Step 4, analysis fit: y_i is Y on pilot rows and mu elsewhere; e_i is 0 on pilot rows and t elsewhere. M is the sum of D(y_i)' D(y_i), plus the sum of e_i added to the [term, term] cell. Solve M b = sum of D(y_i)' O_i. C is the inverse of M.",
  "Step 5: s2 = (sum of (O_i - D(y_i) b)^2 + b_term^2 * sum of e_i) / (n - p). w = n * C[term, ]. a = w_term. c = -b_term.",
  "Step 6, over every record including pilot rows, using mu rather than Y: u_i = D(mu_i) w and v_i = O_i - D(mu_i) b. var_i = t_i * (u_i * c + v_i * a)^2 + 2 * (a * c)^2 * t_i^2.",
  "Step 7: k = mean(var_i) / (s2 * a).",
  "Plan table: relative_width(f) = sqrt(1 + k * (1 - f) / f). records(f) = min(n, ceil(f * n - 1e-8)). Fraction reaching target ratio r: k / (r^2 - 1 + k).",
  "Bootstrap interval: resample pilot rows with replacement, redo steps 1 to 7 with the resampled rows in the regressions only, using RSS / max(rows - 2 * columns, 1) for tau2, then take limits from the draws. R uses bias-corrected and accelerated limits; see catmech_audit_bca in R/audit-plan.R.")

bundle <- list(schema = "catmech_audit_plan_planner_data_1", package_version = as.character(utils::packageVersion("catmech")),
  verdict = read("verdict.csv")$verdict, verdict_rule = "See validation/AUDIT_PLAN_PLAN.md.", datasets = datasets,
  js_test_vectors = vectors, k_calculation = k_calculation,
  limits = c("Widths only. The planner says nothing about bias or coverage.",
             "Assumes a random pilot and linear MI as the analysis. It warns about gross errors but does not repair them.",
             "These are two public datasets, neither a trial's source-data verification, in cohorts of 600."))
jsonlite::write_json(bundle, output, pretty = TRUE, auto_unbox = TRUE, digits = NA, na = "null", matrix = "rowmajor")
cat("Wrote", output, "\n")
