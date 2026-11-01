#!/usr/bin/env Rscript
# Post-hoc checks added after independent review of validation/REAL_AUDIT_SYNTHESIS.md. Not part of
# validation/REAL_AUDIT_PLAN.md, and none of these results enter its decision rules. Each table tests a claim the
# reviewers questioned: the estimand, what drives the slope gaps, whether a verdict survives fresh seeds, the
# stated mechanism, and what limits the intervals. Scratch estimators here are not package features.
# Usage: Rscript tools/real-audit-review-checks.R [NEW_DIRECTORY] [SCALE]
args <- commandArgs(trailingOnly = TRUE)
file_arg <- grep("^--file=", commandArgs(), value = TRUE)
root <- normalizePath(file.path(dirname(sub("^--file=", "", file_arg[1])), ".."))
outdir <- if (length(args)) args[1] else file.path(root, "validation", "real-audit-review")
scale <- if (length(args) >= 2L) as.numeric(args[2]) else 1
if (dir.exists(outdir) && length(list.files(outdir, all.files = TRUE, no.. = TRUE)))
  stop("Use a new or empty output directory; existing evidence is never overwritten.")
pkgload::load_all(root, quiet = TRUE)
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
cores <- 14L
reps <- function(n) max(4L, as.integer(round(n * scale)))
start <- Sys.time()
note <- function(...) { cat(..., "; elapsed", round(as.numeric(difftime(Sys.time(), start, units = "secs"))), "seconds\n"); flush.console() }

study <- file.path(root, "validation", "real-audit-study")
config <- jsonlite::read_json(file.path(study, "settings.json"), simplifyVector = TRUE)
saved <- utils::read.csv(file.path(study, "replicates.csv"), stringsAsFactors = FALSE)
rwanda <- utils::read.csv(file.path(root, "validation/real-audit/rwanda-pool.csv"), stringsAsFactors = FALSE)
rwanda <- rwanda[rwanda$facility_type == "health_center", ]
effect <- utils::read.csv(file.path(root, "validation/real-audit/effect-size-pool.csv"), stringsAsFactors = FALSE)
effect <- effect[effect$info %in% 0:1, ]
RWANDA <- list(pnc = c(w = "pnc_newborn_24h_h", r = "pnc_newborn_24h_r", x = "hf_deliveries_h"),
  syphilis = c(w = "anc_syphilis_test_h", r = "anc_syphilis_test_r", x = "ancnewreg_h"),
  iron = c(w = "womenreciron_h", r = "womenreciron_r", x = "ancnewreg_h"),
  anc4 = c(w = "anc_4visits_h", r = "anc_4visits_r", x = "ancnewreg_h"))
RECOUNT <- c(pnc = "hf_deliveries_r", syphilis = "ancnewreg_r", iron = "ancnewreg_r", anc4 = "ancnewreg_r")
month_pool <- function(id) {
  cols <- RWANDA[[id]]
  d <- rwanda[stats::complete.cases(rwanda[cols]), ]
  data.frame(unit = paste0("f", d$facility_id), stratum = paste0("c", d$hosp_id), lab_value = d[[cols[["w"]]]], lab_true = d[[cols[["r"]]]],
             covariate = d[[cols[["x"]]]], covariate_recount = d[[RECOUNT[[id]]]], stringsAsFactors = FALSE)
}
quarter_pool <- function(id) {
  cols <- RWANDA[[id]]
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
  data.frame(unit = paste0("e", d$row), stratum = d$meta, lab_value = d$recorded, lab_true = d$reference, covariate = log(d$n), stringsAsFactors = FALSE)
}
prepare <- function(z, index) {
  z$index <- index
  z$pool$stratum <- factor(z$pool$stratum, levels = sort(unique(z$pool$stratum)))
  z$rows_by_unit <- split(seq_len(nrow(z$pool)), z$pool$unit)
  z$units_by_stratum <- lapply(split(z$pool$unit, z$pool$stratum), unique)
  x <- stats::model.matrix(~ covariate + stratum, data = z$pool)
  z$truth <- unname(stats::lm.fit(x, z$pool$lab_true)$coefficients["covariate"])
  z$recorded_slope <- unname(stats::lm.fit(x, z$pool$lab_value)$coefficients["covariate"])
  z
}
analyses <- c(
  lapply(names(RWANDA), function(a) list(id = a, level = "rwanda_month", pool = month_pool(a), tol = 0)),
  lapply(names(RWANDA), function(a) list(id = a, level = "rwanda_quarter", pool = quarter_pool(a), tol = 0)),
  list(list(id = "effect_smd", level = "effect_size", pool = effect_pool("smd"), tol = .049),
       list(id = "effect_r", level = "effect_size", pool = effect_pool("r"), tol = .025)))
analyses <- Map(prepare, analyses, seq_along(analyses))
names(analyses) <- vapply(analyses, function(z) paste0(z$level, ":", z$id), character(1))
FOUR <- paste0("rwanda_month:", names(RWANDA))

make_case <- function(z, fraction, cohort_seed, selection_seed, resample = TRUE, pool = z$pool) {
  if (resample) {
    set.seed(cohort_seed)
    drawn <- lapply(z$units_by_stratum, function(u) u[sample.int(length(u), length(u), replace = TRUE)])
  } else drawn <- z$units_by_stratum
  set.seed(selection_seed)
  picked <- lapply(drawn, function(u) {
    size <- max(1L, as.integer(round(fraction * length(u))))
    list(chosen = seq_along(u) %in% sample.int(length(u), size), probability = size / length(u))
  })
  parts <- list(); position <- 0L
  for (s in names(drawn)) for (j in seq_along(drawn[[s]])) {
    position <- position + 1L
    rows <- z$rows_by_unit[[drawn[[s]][j]]]
    parts[[position]] <- data.frame(row = rows, cluster = position, audited = picked[[s]]$chosen[j], probability = picked[[s]]$probability, unit = drawn[[s]][j])
  }
  index <- do.call(rbind, parts)
  d <- pool[index$row, ]; rownames(d) <- NULL
  gold <- d$lab_true; gold[!index$audited] <- NA_real_
  input <- data.frame(covariate = d$covariate, stratum = d$stratum, lab_value = d$lab_value, gold_value = gold)
  x <- stats::model.matrix(~ covariate + stratum, data = input)
  list(input = input, truth = d$lab_true, target = z$truth, predictors = c("covariate", "stratum"), tol = z$tol, audited = index$audited,
       sampling_probability = index$probability, cluster = index$cluster, unit = index$unit, analysis_x = x, term_index = match("covariate", colnames(x)))
}
study_case <- function(z, fraction_index, r) {
  make_case(z, config$fractions[fraction_index], as.integer(config$cohort_base + z$index * 1e6 + r),
            as.integer(config$selection_base + z$index * 1e6 + fraction_index * 1e5 + r))
}
ols_term <- function(y, p, x = p$analysis_x, cluster = p$cluster) {
  use <- is.finite(y); x <- x[use, , drop = FALSE]; k <- ncol(x); n <- nrow(x)
  fit <- stats::lm.fit(x, y[use])
  if (fit$rank != k || n <= k) stop("Rank deficient analysis.")
  position <- match(p$term_index, fit$qr$pivot)
  inverse <- chol2inv(qr.R(fit$qr))
  scores <- as.vector(x[, fit$qr$pivot, drop = FALSE] %*% inverse[position, ]) * fit$residuals
  sums <- rowsum(scores, cluster[use], reorder = FALSE)
  g <- nrow(sums)
  list(estimate = unname(fit$coefficients[p$term_index]), se = sqrt(sum(fit$residuals^2) / (n - k) * inverse[position, position]), df = n - k, n_used = n,
       se_hc1 = sqrt(sum(scores^2) * n / (n - k)), se_cluster = sqrt(sum(sums^2) * (g / (g - 1)) * ((n - 1) / (n - k))), df_cluster = g - 1, clusters = g)
}
pool_fits <- function(fits) {
  q <- vapply(fits, function(x) x$estimate, numeric(1)); u_c <- vapply(fits, function(x) x$se_cluster^2, numeric(1))
  z_c <- mice::pool.scalar(q, u_c, n = fits[[1]]$clusters, k = 1)
  list(estimate = z_c$qbar, se_cluster = sqrt(z_c$t), df_cluster = z_c$df)
}
covers <- function(fit, target) as.numeric(abs(fit$estimate - target) <= stats::qt(.975, fit$df_cluster) * fit$se_cluster)
run_linear <- function(p, seed) pool_fits(lapply(catmech_reference(p$input, "lab_value", "gold_value", p$predictors, method = "normal", m = 20L, seed = seed, use_recorded = TRUE)$completed, function(d) ols_term(d$lab_value, p)))
run_model <- function(p, seed, input = p$input, ...) pool_fits(lapply(suppressWarnings(catmech_contamination(input, "lab_value", "gold_value", p$predictors, agreement_tol = p$tol, m = 20L, seed = seed, ...))$completed, function(d) ols_term(d$lab_value, p)))
rmse <- function(x) sqrt(mean(x^2))
ratio_limits <- function(a, b, seed) {
  set.seed(seed)
  draws <- replicate(2000L, { i <- sample.int(length(a), replace = TRUE); sqrt(sum(a[i]^2) / sum(b[i]^2)) })
  c(ratio = sqrt(sum(a^2) / sum(b^2)), lower = unname(stats::quantile(draws, .025)), upper = unname(stats::quantile(draws, .975)))
}
apply_jobs <- function(n, f) {
  out <- parallel::mclapply(seq_len(n), f, mc.cores = cores, mc.preschedule = TRUE, mc.set.seed = FALSE)
  failed <- vapply(out, inherits, logical(1), what = "try-error")
  if (any(failed)) stop(as.character(out[[which(failed)[1]]]))
  out
}
write_table <- function(x, name) utils::write.csv(x, file.path(outdir, name), row.names = FALSE)

# 1. What each slope is, and what a recount of both sides would give.
set.seed(20261102)
estimand <- do.call(rbind, lapply(FOUR, function(key) {
  z <- analyses[[key]]; d <- z$pool
  slopes <- function(rows) {
    a <- d[rows, ]
    x <- stats::model.matrix(~ covariate + stratum, data = a)
    k <- stats::complete.cases(a$covariate_recount)
    xr <- stats::model.matrix(~ covariate_recount + stratum, data = a[k, ])
    c(reported = unname(stats::lm.fit(x, a$lab_value)$coefficients["covariate"]), one_item = unname(stats::lm.fit(x, a$lab_true)$coefficients["covariate"]),
      both_sides = unname(stats::lm.fit(xr, a$lab_true[k])$coefficients["covariate_recount"]))
  }
  full <- slopes(seq_len(nrow(d)))
  draws <- t(replicate(reps(2000), {
    picked <- unlist(lapply(z$units_by_stratum, function(u) u[sample.int(length(u), length(u), replace = TRUE)]), use.names = FALSE)
    slopes(unlist(z$rows_by_unit[picked], use.names = FALSE))
  }))
  gap_one <- draws[, "reported"] - draws[, "one_item"]; gap_both <- draws[, "reported"] - draws[, "both_sides"]
  data.frame(analysis = z$id, rows = nrow(d), reported = full[["reported"]], one_item_target = full[["one_item"]], both_sides_recount = full[["both_sides"]],
    gap_one_item = full[["reported"]] - full[["one_item"]], gap_one_item_percent = 100 * (full[["reported"]] / full[["one_item"]] - 1),
    gap_one_item_lower = unname(stats::quantile(gap_one, .025)), gap_one_item_upper = unname(stats::quantile(gap_one, .975)),
    gap_both_sides = full[["reported"]] - full[["both_sides"]], gap_both_sides_percent = 100 * (full[["reported"]] / full[["both_sides"]] - 1),
    gap_both_sides_lower = unname(stats::quantile(gap_both, .025)), gap_both_sides_upper = unname(stats::quantile(gap_both, .975)),
    covariate_rows_wrong = sum(abs(d$covariate - d$covariate_recount) > 1e-8, na.rm = TRUE), correlation_target_error_covariate_error = stats::cor(d$lab_value - d$lab_true, d$covariate - d$covariate_recount, use = "complete.obs"))
}))
write_table(estimand, "estimand-slopes.csv"); note("estimand")

# 2. What drives each gap: copied values, a few facilities, and recounts of zero.
gap_of <- function(d) {
  x <- stats::model.matrix(~ covariate + stratum, data = droplevels(d))
  unname(stats::lm.fit(x, d$lab_value)$coefficients["covariate"] - stats::lm.fit(x, d$lab_true)$coefficients["covariate"])
}
tops <- list()
structure_rows <- do.call(rbind, lapply(FOUR, function(key) {
  z <- analyses[[key]]; d <- z$pool; e <- d$lab_value - d$lab_true; wrong <- abs(e) > 1e-8
  full <- gap_of(d)
  units <- names(z$rows_by_unit)
  without <- vapply(units, function(u) gap_of(d[d$unit != u, ]), numeric(1))
  ranked <- units[order(-sign(full) * (full - without))]
  top <- ranked[1:3]
  tops[[z$id]] <<- top
  keep <- !d$unit %in% top
  zero <- d$lab_true == 0 & d$lab_value >= 10
  reduced <- d[!zero, ]
  xr <- stats::model.matrix(~ covariate + stratum, data = droplevels(reduced))
  data.frame(analysis = z$id, rows = nrow(d), wrong_rows = sum(wrong), reported_equals_covariate = sum(abs(d$lab_value - d$covariate) < 1e-8),
    wrong_rows_copying_covariate = sum(wrong & abs(d$lab_value - d$covariate) < 1e-8), top3 = paste(top, collapse = " "), top3_rows = sum(!keep), gap = full, gap_without_top3 = gap_of(d[keep, ]),
    correlation = stats::cor(e, d$covariate), correlation_without_top3 = stats::cor(e[keep], d$covariate[keep]),
    within_catchment_correlation = stats::cor(stats::lm.fit(stats::model.matrix(~ stratum, d), e)$residuals, stats::lm.fit(stats::model.matrix(~ stratum, d), d$covariate)$residuals),
    zero_recount_rows = sum(zero), zero_recount_facilities = length(unique(d$unit[zero])), one_item_without_zero = unname(stats::lm.fit(xr, reduced$lab_true)$coefficients["covariate"]),
    reported_without_zero = unname(stats::lm.fit(xr, reduced$lab_value)$coefficients["covariate"]), gap_without_zero = gap_of(reduced))
}))
write_table(structure_rows, "gap-structure.csv"); note("gap structure")

# 3. Do single audits correct the bias, or does the average hide cancellation? Study cohorts at 25%.
conditional <- do.call(rbind, lapply(FOUR, function(key) {
  z <- analyses[[key]]; top <- tops[[z$id]]
  flags <- unlist(apply_jobs(config$month_reps, function(r) { p <- study_case(z, 1L, r); any(p$audited & p$unit %in% top) }))
  scenario <- paste0(key, ":audit_25pct")
  do.call(rbind, lapply(c("recorded", "linear_mi", "contamination", "audited_only"), function(m) {
    g <- saved[saved$scenario == scenario & saved$method == m & saved$status == "ok", ]
    hit <- flags[g$replicate]
    data.frame(analysis = z$id, method = m, cohorts_none_audited = sum(!hit), mean_error_none_audited = mean(g$error[!hit]), cohorts_any_audited = sum(hit), mean_error_any_audited = mean(g$error[hit]),
      coverage_none_audited = mean(g$coverage_cluster[!hit]), coverage_any_audited = mean(g$coverage_cluster[hit]), mean_width_cluster = mean(g$width_cluster))
  }))
}))
write_table(conditional, "conditional-on-key-facilities.csv"); note("conditional")

# 4. Does Rule 1 survive fresh seeds? No MCMC, so many cohorts are cheap. Also the sensitivity to zero recounts.
extend <- function(z, fraction_index, n, label = z$id, offset = 0) {
  rows <- apply_jobs(n, function(r) {
    base <- z$index * 1e6 + offset + r
    p <- make_case(z, config$fractions[fraction_index], as.integer(700100000 + base), as.integer(710100000 + base + fraction_index * 1e5))
    oracle <- ols_term(p$truth, p)$estimate; recorded <- ols_term(p$input$lab_value, p)$estimate
    linear <- tryCatch(run_linear(p, as.integer(720000000 + base + fraction_index * 1e5))$estimate, error = function(e) NA_real_)
    c(oracle = oracle, recorded = recorded, linear = linear)
  })
  m <- do.call(rbind, rows); m <- m[is.finite(m[, "linear"]), , drop = FALSE]
  recovery <- ratio_limits(m[, "linear"] - m[, "oracle"], m[, "recorded"] - m[, "oracle"], 730000000 + z$index * 1000 + fraction_index)
  pool <- ratio_limits(m[, "linear"] - z$truth, m[, "recorded"] - z$truth, 740000000 + z$index * 1000 + fraction_index)
  data.frame(analysis = label, fraction = config$fractions[fraction_index], cohorts = nrow(m), recovery_ratio = recovery[["ratio"]], recovery_lower = recovery[["lower"]], recovery_upper = recovery[["upper"]],
    pool_ratio = pool[["ratio"]], pool_lower = pool[["lower"]], pool_upper = pool[["upper"]])
}
extension <- do.call(rbind, lapply(FOUR, function(key) rbind(extend(analyses[[key]], 1L, reps(2000)), extend(analyses[[key]], 2L, reps(2000)))))
z <- analyses[["rwanda_month:syphilis"]]
z$pool <- z$pool[!(z$pool$lab_true == 0 & z$pool$lab_value >= 10), ]
z$rows_by_unit <- split(seq_len(nrow(z$pool)), z$pool$unit)
z$units_by_stratum <- lapply(split(z$pool$unit, z$pool$stratum), unique)
x <- stats::model.matrix(~ covariate + stratum, data = z$pool)
z$truth <- unname(stats::lm.fit(x, z$pool$lab_true)$coefficients["covariate"])
extension <- rbind(extension, extend(z, 1L, reps(1000), "syphilis_without_zero_recounts", 500000))
write_table(extension, "rule1-extension.csv"); note("rule 1 extension")

# 5. How often would checking errors against the covariate on audited rows flag anything? Study cohorts at 25%.
detection <- do.call(rbind, lapply(FOUR, function(key) {
  z <- analyses[[key]]
  res <- do.call(rbind, apply_jobs(config$month_reps, function(r) {
    p <- study_case(z, 1L, r)
    e <- p$input$lab_value - p$truth; e[!p$audited] <- NA
    fit <- ols_term(e, p)
    c(slope = fit$estimate, p = 2 * stats::pt(-abs(fit$estimate / fit$se_cluster), fit$df_cluster))
  }))
  data.frame(analysis = z$id, cohorts = nrow(res), pool_gap = z$recorded_slope - z$truth, mean_slope = mean(res[, "slope"]), share_significant = mean(res[, "p"] < .05), share_wrong_sign = mean(sign(res[, "slope"]) != sign(z$recorded_slope - z$truth)))
}))
write_table(detection, "audited-rows-detection.csv"); note("detection")

# 6. What the oracle's default shortfall is made of. Study cohorts.
oracle_intervals <- do.call(rbind, lapply(FOUR, function(key) {
  z <- analyses[[key]]
  res <- do.call(rbind, apply_jobs(config$month_reps, function(r) {
    p <- study_case(z, 1L, r); f <- ols_term(p$truth, p)
    c(default = abs(f$estimate - z$truth) <= stats::qt(.975, f$df) * f$se, hc1 = abs(f$estimate - z$truth) <= stats::qt(.975, f$df) * f$se_hc1, cr1 = abs(f$estimate - z$truth) <= stats::qt(.975, f$df_cluster) * f$se_cluster)
  }))
  data.frame(analysis = z$id, default = mean(res[, "default"]), hc1_unequal_variance = mean(res[, "hc1"]), cr1_clustered = mean(res[, "cr1"]))
}))
write_table(oracle_intervals, "oracle-interval-types.csv"); note("oracle intervals")

# 7. Raking with the same small-sample factor as every other method, from saved results.
raking <- saved[saved$method == "raking" & saved$status == "ok", ]
k <- raking$audited_rows - raking$df; g <- raking$df_cluster + 1
factor <- sqrt((g / (g - 1)) * ((raking$audited_rows - 1) / (raking$audited_rows - k)))
half <- stats::qt(.975, raking$df_cluster) * raking$se_cluster * factor
raking$coverage_factor <- as.numeric(abs(raking$estimate - raking$truth) <= half); raking$width_factor <- 2 * half; raking$multiplier <- factor
like <- do.call(rbind, lapply(split(raking, raking$scenario), function(a) data.frame(scenario = a$scenario[1], multiplier_mean = mean(a$multiplier), coverage_saved = mean(a$coverage_cluster), coverage_with_factor = mean(a$coverage_factor), width_saved = mean(a$width_cluster), width_with_factor = mean(a$width_factor))))
write_table(like, "raking-like-with-like.csv"); note("raking factor")

# 8. A design-weighted difference estimator with the raking proxy: does the raking column speak for design-based correction?
difference <- do.call(rbind, lapply(analyses, function(z) do.call(rbind, lapply(1:2, function(j) {
  n <- switch(z$level, rwanda_month = config$month_reps, rwanda_quarter = config$quarter_reps, effect_size = config$effect_reps)
  res <- do.call(rbind, apply_jobs(n, function(r) {
    p <- study_case(z, j, r)
    x <- p$analysis_x; y <- p$input$gold_value; w <- p$input$lab_value; ref <- p$audited; d <- 1 / p$sampling_probability[ref]
    xo <- cbind(x, recorded = w)
    b <- stats::lm.wfit(xo[ref, , drop = FALSE], y[ref], d)$coefficients; b[is.na(b)] <- 0
    proxy <- as.vector(xo %*% b)
    full <- stats::lm.fit(x, proxy)$coefficients[p$term_index]
    xa <- x[ref, , drop = FALSE]
    correction <- tryCatch(solve(crossprod(xa, xa * d), crossprod(xa, (y[ref] - proxy[ref]) * d))[p$term_index], error = function(e) NA_real_)
    c(estimate = unname(full) + correction, oracle = ols_term(p$truth, p)$estimate)
  }))
  res <- res[is.finite(res[, "estimate"]), , drop = FALSE]
  scenario <- paste0(z$level, ":", z$id, ":audit_", 100 * config$fractions[j], "pct")
  look <- function(m) { a <- saved[saved$scenario == scenario & saved$method == m & saved$status == "ok", ]; sqrt(mean(a$recovery_error^2)) }
  data.frame(scenario = scenario, cohorts = nrow(res), difference_recovery_rmse = rmse(res[, "estimate"] - res[, "oracle"]), linear_mi_recovery_rmse = look("linear_mi"), raking_recovery_rmse = look("raking"))
}))))
write_table(difference, "difference-estimator.csv"); note("difference estimator")

# 9. Mechanism: shuffle errors within catchment before each cohort, which removes any link to the covariate and any copying.
mechanism <- do.call(rbind, lapply(FOUR, function(key) {
  z <- analyses[[key]]
  res <- do.call(rbind, apply_jobs(reps(200), function(r) {
    set.seed(as.integer(750000000 + z$index * 1e6 + r))
    pool <- z$pool
    e <- pool$lab_value - pool$lab_true
    for (s in levels(pool$stratum)) { i <- which(pool$stratum == s); e[i] <- e[i][sample.int(length(i))] }
    pool$lab_value <- pool$lab_true + e
    p <- make_case(z, .25, as.integer(760000000 + z$index * 1e6 + r), as.integer(770000000 + z$index * 1e6 + r), pool = pool)
    oracle <- ols_term(p$truth, p); recorded <- ols_term(p$input$lab_value, p)
    linear <- tryCatch(run_linear(p, as.integer(780000000 + z$index * 1e6 + r)), error = function(e) NULL)
    model <- tryCatch(run_model(p, as.integer(790000000 + z$index * 1e6 + r)), error = function(e) NULL)
    c(oracle_estimate = oracle$estimate, oracle_cover = covers(oracle, z$truth), recorded_estimate = recorded$estimate,
      linear_estimate = if (is.null(linear)) NA else linear$estimate, linear_cover = if (is.null(linear)) NA else covers(linear, z$truth),
      model_estimate = if (is.null(model)) NA else model$estimate, model_cover = if (is.null(model)) NA else covers(model, z$truth), model_se = if (is.null(model)) NA else model$se_cluster)
  }))
  ok <- is.finite(res[, "model_estimate"]); li <- is.finite(res[, "linear_estimate"])
  data.frame(analysis = z$id, cohorts = nrow(res), recorded_bias = mean(res[, "recorded_estimate"]) - z$truth, oracle_coverage = mean(res[, "oracle_cover"]),
    linear_returned = sum(li), linear_bias = mean(res[li, "linear_estimate"]) - z$truth, linear_coverage = mean(res[li, "linear_cover"]), linear_recovery_rmse = rmse(res[li, "linear_estimate"] - res[li, "oracle_estimate"]),
    model_returned = sum(ok), model_bias = mean(res[ok, "model_estimate"]) - z$truth, model_coverage = mean(res[ok, "model_cover"]), model_recovery_rmse = rmse(res[ok, "model_estimate"] - res[ok, "oracle_estimate"]),
    model_se_over_sd = mean(res[ok, "model_se"]) / stats::sd(res[ok, "model_estimate"]))
}))
write_table(mechanism, "mechanism-shuffled-errors.csv"); note("mechanism")

# 10. A configuration a careful user could try: separate error means by facility-size group. Study cohorts at 25%, paired with the default fit.
size_strata <- do.call(rbind, lapply(FOUR, function(key) {
  z <- analyses[[key]]
  size <- tapply(z$pool$covariate, z$pool$unit, mean)
  groups <- cut(size, stats::quantile(size, c(0, 1/3, 2/3, 1)), include.lowest = TRUE, labels = c("small", "middle", "large"))
  names(groups) <- names(size)
  scenario <- paste0(key, ":audit_25pct")
  n <- reps(150)
  res <- do.call(rbind, apply_jobs(n, function(r) {
    p <- study_case(z, 1L, r)
    input <- p$input; input$size_group <- factor(groups[p$unit], levels = c("small", "middle", "large"))
    fit <- tryCatch(run_model(p, as.integer(config$imputation_base + z$index * 1e6 + 1e5 + r), input = input, error_strata = "size_group"), error = function(e) NULL)
    c(replicate = r, estimate = if (is.null(fit)) NA else fit$estimate, cover = if (is.null(fit)) NA else covers(fit, z$truth), oracle = ols_term(p$truth, p)$estimate)
  }))
  g <- saved[saved$scenario == scenario & saved$method == "contamination" & saved$status == "ok", c("replicate", "estimate", "coverage_cluster", "recovery_error")]
  l <- saved[saved$scenario == scenario & saved$method == "linear_mi" & saved$status == "ok", c("replicate", "recovery_error", "error", "coverage_cluster")]
  both <- merge(merge(as.data.frame(res), g, by = "replicate", suffixes = c("", "_default")), l, by = "replicate", suffixes = c("", "_linear"))
  both <- both[is.finite(both$estimate), ]
  data.frame(analysis = z$id, requested = n, returned = sum(is.finite(res[, "estimate"])), paired = nrow(both), strata_recovery_rmse = rmse(both$estimate - both$oracle), default_recovery_rmse = rmse(both$recovery_error),
    linear_recovery_rmse = rmse(both$recovery_error_linear), strata_bias = mean(both$estimate) - z$truth, default_bias = mean(both$estimate_default) - z$truth, linear_bias = mean(both$error),
    strata_coverage = mean(both$cover), default_coverage = mean(both$coverage_cluster), linear_coverage = mean(both$coverage_cluster_linear))
}))
write_table(size_strata, "error-strata-by-size.csv"); note("error strata")

# 11. What limits linear MI's intervals: imputations that treat rows as independent. Facility bootstrap in each imputation; study cohorts at 25%.
cluster_mi <- do.call(rbind, lapply(FOUR, function(key) {
  z <- analyses[[key]]
  scenario <- paste0(key, ":audit_25pct")
  n <- reps(200)
  res <- do.call(rbind, apply_jobs(n, function(r) {
    p <- study_case(z, 1L, r)
    set.seed(as.integer(800000000 + z$index * 1e6 + r))
    x <- cbind(p$analysis_x, recorded = p$input$lab_value)
    y <- p$input$gold_value; au <- p$audited
    audited_clusters <- split(unique(p$cluster[au]), as.character(p$input$stratum[match(unique(p$cluster[au]), p$cluster)]))
    fits <- lapply(seq_len(20L), function(j) {
      pick <- unlist(lapply(audited_clusters, function(u) if (length(u) == 1L) u else u[sample.int(length(u), replace = TRUE)]), use.names = FALSE)
      rows <- unlist(lapply(pick, function(cl) which(p$cluster == cl)), use.names = FALSE)
      fit <- stats::lm.fit(x[rows, , drop = FALSE], y[rows])
      b <- fit$coefficients; b[is.na(b)] <- 0
      s2 <- sum(fit$residuals^2) / max(1, length(rows) - fit$rank)
      completed <- y
      completed[!au] <- as.vector(x[!au, , drop = FALSE] %*% b) + stats::rnorm(sum(!au), 0, sqrt(s2))
      ols_term(completed, p)
    })
    f <- pool_fits(fits)
    c(replicate = r, estimate = f$estimate, cover = covers(f, z$truth), se = f$se_cluster, oracle = ols_term(p$truth, p)$estimate)
  }))
  l <- saved[saved$scenario == scenario & saved$method == "linear_mi" & saved$status == "ok", c("replicate", "recovery_error", "coverage_cluster", "se_cluster", "estimate")]
  o <- saved[saved$scenario == scenario & saved$method == "oracle", c("replicate", "coverage_cluster")]
  both <- merge(merge(as.data.frame(res), l, by = "replicate", suffixes = c("", "_saved")), o, by = "replicate", suffixes = c("", "_oracle"))
  data.frame(analysis = z$id, cohorts = nrow(both), facility_bootstrap_coverage = mean(both$cover), package_coverage = mean(both$coverage_cluster), oracle_coverage = mean(both$coverage_cluster_oracle),
    facility_bootstrap_recovery_rmse = rmse(both$estimate - both$oracle), package_recovery_rmse = rmse(both$recovery_error),
    facility_bootstrap_se_over_sd = mean(both$se) / stats::sd(both$estimate), package_se_over_sd = mean(both$se_cluster) / stats::sd(both$estimate_saved))
}))
write_table(cluster_mi, "linear-mi-facility-bootstrap.csv"); note("facility bootstrap MI")

# 12. Hold the pool fixed and redraw only the audit, where recorded's error is simply its gap.
fixed_swap <- do.call(rbind, lapply(analyses, function(z) do.call(rbind, lapply(1:2, function(j) {
  x <- stats::model.matrix(~ covariate + stratum, data = z$pool)
  oracle <- z$truth
  res <- unlist(apply_jobs(reps(2000), function(r) {
    p <- make_case(z, config$fractions[j], 1L, as.integer(810000000 + z$index * 1e6 + j * 1e5 + r), resample = FALSE)
    swap <- p$input$lab_value; swap[p$audited] <- p$truth[p$audited]
    ols_term(swap, p)$estimate - oracle
  }))
  data.frame(level = z$level, analysis = z$id, fraction = config$fractions[j], recorded_absolute_gap = abs(z$recorded_slope - z$truth), swap_in_recovery_rmse = rmse(res), swap_over_recorded = rmse(res) / abs(z$recorded_slope - z$truth))
}))))
write_table(fixed_swap, "fixed-pool-swap-in.csv"); note("fixed pool swap-in")
fixed_methods <- do.call(rbind, lapply(FOUR, function(key) {
  z <- analyses[[key]]
  res <- do.call(rbind, apply_jobs(reps(200), function(r) {
    p <- make_case(z, .25, 1L, as.integer(820000000 + z$index * 1e6 + r), resample = FALSE)
    linear <- tryCatch(run_linear(p, as.integer(830000000 + z$index * 1e6 + r))$estimate, error = function(e) NA_real_)
    model <- tryCatch(run_model(p, as.integer(840000000 + z$index * 1e6 + r))$estimate, error = function(e) NA_real_)
    c(linear = linear, model = model)
  }))
  li <- is.finite(res[, "linear"]); both <- li & is.finite(res[, "model"])
  gap <- abs(z$recorded_slope - z$truth)
  data.frame(analysis = z$id, audits = nrow(res), recorded_absolute_gap = gap, linear_recovery_rmse = rmse(res[li, "linear"] - z$truth), linear_over_recorded = rmse(res[li, "linear"] - z$truth) / gap,
    pairs = sum(both), model_recovery_rmse = rmse(res[both, "model"] - z$truth), model_over_linear = rmse(res[both, "model"] - z$truth) / rmse(res[both, "linear"] - z$truth))
}))
write_table(fixed_methods, "fixed-pool-methods.csv"); note("fixed pool methods")

writeLines(c("Post-hoc checks added after independent review of validation/REAL_AUDIT_SYNTHESIS.md. They are not part of validation/REAL_AUDIT_PLAN.md and enter none of its decision rules.",
  "Scratch estimators here, including the facility-bootstrap imputation and the difference estimator, are not package features. Elapsed seconds: ", round(as.numeric(difftime(Sys.time(), start, units = "secs")), 1)), file.path(outdir, "README.md"))
cat("Output:", normalizePath(outdir), "\n")
