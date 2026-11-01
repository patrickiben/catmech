#!/usr/bin/env Rscript
# Post-hoc explanation for validation/REAL_AUDIT_SYNTHESIS.md, added after viewing the study
# results and not part of the plan. Rebuilds the study's cohorts from their seeds, checks them
# against the saved estimates, then computes swap-in: reference values on audited rows and
# recorded values everywhere else. Plain OLS only; no imputation model is fitted.
# Usage: Rscript tools/real-audit-posthoc.R [NEW_DIRECTORY]
args <- commandArgs(trailingOnly = TRUE)
file_arg <- grep("^--file=", commandArgs(), value = TRUE)
root <- normalizePath(file.path(dirname(sub("^--file=", "", file_arg[1])), ".."))
outdir <- if (length(args)) args[1] else file.path(root, "validation", "real-audit-posthoc")
if (dir.exists(outdir) && length(list.files(outdir, all.files = TRUE, no.. = TRUE)))
  stop("Use a new or empty output directory; existing evidence is never overwritten.")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
study <- file.path(root, "validation", "real-audit-study")
config <- jsonlite::read_json(file.path(study, "settings.json"), simplifyVector = TRUE)
saved <- utils::read.csv(file.path(study, "replicates.csv"), stringsAsFactors = FALSE)
saved <- saved[saved$method %in% c("oracle", "recorded", "linear_mi", "contamination"), c("scenario", "replicate", "method", "estimate", "status", "truth")]

rwanda <- utils::read.csv(file.path(root, "validation/real-audit/rwanda-pool.csv"), stringsAsFactors = FALSE)
rwanda <- rwanda[rwanda$facility_type == "health_center", ]
effect <- utils::read.csv(file.path(root, "validation/real-audit/effect-size-pool.csv"), stringsAsFactors = FALSE)
effect <- effect[effect$info %in% 0:1, ]
RWANDA <- list(pnc = c(w = "pnc_newborn_24h_h", r = "pnc_newborn_24h_r", x = "hf_deliveries_h"),
  syphilis = c(w = "anc_syphilis_test_h", r = "anc_syphilis_test_r", x = "ancnewreg_h"),
  iron = c(w = "womenreciron_h", r = "womenreciron_r", x = "ancnewreg_h"),
  anc4 = c(w = "anc_4visits_h", r = "anc_4visits_r", x = "ancnewreg_h"))
month_pool <- function(cols) {
  d <- rwanda[stats::complete.cases(rwanda[cols]), ]
  data.frame(unit = paste0("f", d$facility_id), stratum = paste0("c", d$hosp_id), lab_value = d[[cols[["w"]]]],
             lab_true = d[[cols[["r"]]]], covariate = d[[cols[["x"]]]], stringsAsFactors = FALSE)
}
quarter_pool <- function(cols) {
  d <- rwanda[c("facility_id", "hosp_id", cols)]
  whole <- vapply(split(stats::complete.cases(d), d$facility_id), all, logical(1))
  d <- d[d$facility_id %in% names(whole)[whole], ]
  totals <- stats::aggregate(d[cols], by = list(facility_id = d$facility_id, hosp_id = d$hosp_id), FUN = sum)
  data.frame(unit = paste0("f", totals$facility_id), stratum = paste0("c", totals$hosp_id), lab_value = totals[[cols[["w"]]]],
             lab_true = totals[[cols[["r"]]]], covariate = totals[[cols[["x"]]]], stringsAsFactors = FALSE)
}
effect_pool <- function(family) {
  d <- effect[effect$family == family, ]
  counts <- table(d$meta)
  d <- d[d$meta %in% names(counts)[counts >= 2L], ]
  data.frame(unit = paste0("e", d$row), stratum = d$meta, lab_value = d$recorded, lab_true = d$reference,
             covariate = log(d$n), stringsAsFactors = FALSE)
}
analyses <- c(
  lapply(names(RWANDA), function(a) list(id = a, level = "rwanda_month", pool = month_pool(RWANDA[[a]]), reps = config$month_reps)),
  lapply(names(RWANDA), function(a) list(id = a, level = "rwanda_quarter", pool = quarter_pool(RWANDA[[a]]), reps = config$quarter_reps)),
  list(list(id = "effect_smd", level = "effect_size", pool = effect_pool("smd"), reps = config$effect_reps),
       list(id = "effect_r", level = "effect_size", pool = effect_pool("r"), reps = config$effect_reps)))
for (i in seq_along(analyses)) {
  z <- analyses[[i]]
  z$index <- i
  z$pool$stratum <- factor(z$pool$stratum, levels = sort(unique(z$pool$stratum)))
  z$rows_by_unit <- split(seq_len(nrow(z$pool)), z$pool$unit)
  z$units_by_stratum <- lapply(split(z$pool$unit, z$pool$stratum), unique)
  analyses[[i]] <- z
}
slope <- function(y, x) unname(stats::lm.fit(x, y)$coefficients["covariate"])
make_case <- function(z, fraction_index, r) {
  set.seed(as.integer(config$cohort_base + z$index * 1e6 + r))
  drawn <- lapply(z$units_by_stratum, function(u) u[sample.int(length(u), length(u), replace = TRUE)])
  set.seed(as.integer(config$selection_base + z$index * 1e6 + fraction_index * 1e5 + r))
  fraction <- config$fractions[fraction_index]
  picked <- lapply(drawn, function(u) seq_along(u) %in% sample.int(length(u), max(1L, as.integer(round(fraction * length(u))))))
  parts <- list(); position <- 0L
  for (s in names(drawn)) for (j in seq_along(drawn[[s]])) {
    position <- position + 1L
    rows <- z$rows_by_unit[[drawn[[s]][j]]]
    parts[[position]] <- data.frame(row = rows, audited = picked[[s]][j])
  }
  index <- do.call(rbind, parts)
  d <- z$pool[index$row, ]
  list(d = d, audited = index$audited, x = stats::model.matrix(~ covariate + stratum, data = d))
}
rmse <- function(x) sqrt(mean(x^2))
rows <- list(); checked <- 0L
for (z in analyses) for (j in seq_along(config$fractions)) {
  scenario <- paste0(z$level, ":", z$id, ":audit_", 100 * config$fractions[j], "pct")
  g <- saved[saved$scenario == scenario, ]
  out <- parallel::mclapply(seq_len(z$reps), function(r) {
    p <- make_case(z, j, r)
    oracle <- slope(p$d$lab_true, p$x); recorded <- slope(p$d$lab_value, p$x)
    swap <- p$d$lab_value; swap[p$audited] <- p$d$lab_true[p$audited]
    data.frame(replicate = r, oracle = oracle, recorded = recorded, swap_in = slope(swap, p$x))
  }, mc.cores = 12L, mc.set.seed = FALSE)
  out <- do.call(rbind, out)
  for (method in c("oracle", "recorded")) {
    s <- g[g$method == method, ]
    m <- out[[method]][match(s$replicate, out$replicate)]
    if (max(abs(m - s$estimate)) > 1e-7) stop("Rebuilt cohorts do not match the study: ", scenario, " ", method)
    checked <- checked + nrow(s)
  }
  pick <- function(method) { s <- g[g$method == method & g$status == "ok", ]; s$estimate[match(out$replicate, s$replicate)] }
  contamination <- pick("contamination"); linear <- pick("linear_mi")
  both <- !is.na(contamination)
  truth <- g$truth[1]
  rows[[length(rows) + 1L]] <- data.frame(level = z$level, analysis = z$id, scenario = scenario, fraction = config$fractions[j], cohorts = z$reps,
    recorded_recovery_rmse = rmse(out$recorded - out$oracle), swap_in_recovery_rmse = rmse(out$swap_in - out$oracle),
    swap_in_bias = mean(out$swap_in) - truth, recorded_bias = mean(out$recorded) - truth,
    linear_mi_recovery_rmse = rmse((linear - out$oracle)[!is.na(linear)]),
    pairs_with_contamination = sum(both), contamination_recovery_rmse = rmse((contamination - out$oracle)[both]),
    swap_in_recovery_rmse_same_cohorts = rmse((out$swap_in - out$oracle)[both]),
    contamination_bias = mean(contamination[both]) - truth, swap_in_bias_same_cohorts = mean(out$swap_in[both]) - truth,
    correlation_contamination_swap_in = stats::cor(contamination[both], out$swap_in[both]),
    rms_gap_contamination_swap_in = rmse((contamination - out$swap_in)[both]))
}
summary_table <- do.call(rbind, rows)
utils::write.csv(summary_table, file.path(outdir, "swap-in-summary.csv"), row.names = FALSE)
writeLines(c("Post-hoc analysis, added after viewing the study results. Not part of validation/REAL_AUDIT_PLAN.md.",
  paste0("Rebuilt cohorts matched ", checked, " saved oracle and recorded estimates to within 1e-7."),
  "Swap-in: reference values on audited rows, recorded values elsewhere, then plain OLS. No interval is computed."),
  file.path(outdir, "README.md"))
cat(readLines(file.path(outdir, "README.md")), sep = "\n")
