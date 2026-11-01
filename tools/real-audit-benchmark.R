#!/usr/bin/env Rscript
# Prospective test on real audit data; settings locked in validation/REAL_AUDIT_PLAN.md.
# Usage: Rscript tools/real-audit-benchmark.R pilot|study [NEW_DIRECTORY] [PILOT_REPS]
args <- commandArgs(trailingOnly = TRUE)
mode <- if (length(args)) match.arg(args[1], c("pilot", "study")) else "pilot"
file_arg <- grep("^--file=", commandArgs(), value = TRUE)
root <- normalizePath(file.path(dirname(sub("^--file=", "", file_arg[1])), ".."))
outdir <- if (length(args) >= 2L) args[2] else file.path(root, "validation", paste0("real-audit-", mode))
if (dir.exists(outdir) && length(list.files(outdir, all.files = TRUE, no.. = TRUE)))
  stop("Use a new or empty output directory; existing evidence is never overwritten.")
hash_sources <- function() {
  paths <- c(sort(list.files(file.path(root, "R"), "\\.R$", full.names = TRUE)),
    file.path(root, c("DESCRIPTION", "NAMESPACE", "tools/real-audit-benchmark.R", "validation/REAL_AUDIT_PLAN.md",
                      "validation/real-audit/rwanda-pool.csv", "validation/real-audit/effect-size-pool.csv")))
  data.frame(file = substring(paths, nchar(root) + 2L),
    sha256 = vapply(paths, function(p) digest::digest(file = p, algo = "sha256"), character(1)), row.names = NULL)
}
hashes_before <- hash_sources()
pkgload::load_all(root, quiet = TRUE)
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
workers <- suppressWarnings(as.integer(Sys.getenv("CATMECH_BENCHMARK_WORKERS", "14")))
if (length(workers) != 1L || is.na(workers) || workers < 1L || workers > 16L) stop("Invalid worker count (1 to 16).")
if (.Platform$OS.type == "windows") workers <- 1L
pilot_reps <- if (length(args) >= 3L) as.integer(args[3]) else 3L
config <- list(mode = mode, workers = workers, m = 20L, fractions = c(.25, .50),
  month_reps = if (mode == "pilot") pilot_reps else 500L, quarter_reps = if (mode == "pilot") pilot_reps else 200L,
  effect_reps = if (mode == "pilot") pilot_reps else 500L, ratio_resamples = 2000L, ratio_seed = 430000000L,
  cohort_base = if (mode == "pilot") 500100000 else 400100000,
  selection_base = if (mode == "pilot") 510100000 else 410100000,
  imputation_base = if (mode == "pilot") 520000000 else 420000000)
jsonlite::write_json(config, file.path(outdir, "settings.json"), pretty = TRUE, auto_unbox = TRUE)

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
  lapply(names(RWANDA), function(a) list(id = a, level = "rwanda_month", pool = month_pool(RWANDA[[a]]), tol = 0, reps = config$month_reps)),
  lapply(names(RWANDA), function(a) list(id = a, level = "rwanda_quarter", pool = quarter_pool(RWANDA[[a]]), tol = 0, reps = config$quarter_reps)),
  list(list(id = "effect_smd", level = "effect_size", pool = effect_pool("smd"), tol = .049, reps = config$effect_reps),
       list(id = "effect_r", level = "effect_size", pool = effect_pool("r"), tol = .025, reps = config$effect_reps)))
for (i in seq_along(analyses)) {
  z <- analyses[[i]]
  z$index <- i
  z$strata <- sort(unique(z$pool$stratum))
  z$pool$stratum <- factor(z$pool$stratum, levels = z$strata)
  z$rows_by_unit <- split(seq_len(nrow(z$pool)), z$pool$unit)
  z$units_by_stratum <- lapply(split(z$pool$unit, z$pool$stratum), unique)
  x <- stats::model.matrix(~ covariate + stratum, data = z$pool)
  z$truth <- unname(stats::lm.fit(x, z$pool$lab_true)$coefficients["covariate"])
  z$recorded_coefficient <- unname(stats::lm.fit(x, z$pool$lab_value)$coefficients["covariate"])
  analyses[[i]] <- z
}

ols_term <- function(y, p) {
  use <- is.finite(y)
  x <- p$analysis_x[use, , drop = FALSE]; k <- ncol(x); n <- nrow(x)
  fit <- stats::lm.fit(x, y[use])
  if (fit$rank != k || n <= k) stop("Scientific analysis is rank deficient or has no residual degrees of freedom.")
  position <- match(p$term_index, fit$qr$pivot)
  inverse <- chol2inv(qr.R(fit$qr))
  scores <- as.vector(x[, fit$qr$pivot, drop = FALSE] %*% inverse[position, ]) * fit$residuals
  sums <- rowsum(scores, p$cluster[use], reorder = FALSE)
  g <- nrow(sums)
  list(estimate = unname(fit$coefficients[p$term_index]), se = sqrt(sum(fit$residuals^2) / (n - k) * inverse[position, position]),
       df = n - k, n_used = n, se_cluster = sqrt(sum(sums^2) * (g / (g - 1)) * ((n - 1) / (n - k))), df_cluster = g - 1, clusters = g)
}
pool_fits <- function(fits) {
  q <- vapply(fits, function(x) x$estimate, numeric(1))
  u <- vapply(fits, function(x) x$se^2, numeric(1)); u_c <- vapply(fits, function(x) x$se_cluster^2, numeric(1))
  first <- fits[[1]]
  z <- mice::pool.scalar(q, u, n = first$n_used, k = first$n_used - first$df)
  z_c <- mice::pool.scalar(q, u_c, n = first$clusters, k = 1)
  list(estimate = z$qbar, se = sqrt(z$t), df = z$df, n_used = first$n_used, se_cluster = sqrt(z_c$t), df_cluster = z_c$df, clusters = first$clusters)
}
fit_raking <- function(p) {
  y <- p$input$gold_value; w <- p$input$lab_value
  ref <- !is.na(y); pi_i <- p$sampling_probability
  x_obs <- cbind(p$analysis_x, recorded = w)
  use <- which(ref)
  beta <- stats::lm.wfit(x_obs[use, , drop = FALSE], y[use], 1 / pi_i[use])$coefficients
  beta[is.na(beta)] <- 0
  proxy <- as.vector(x_obs %*% beta)
  naive <- stats::lm.fit(p$analysis_x, proxy)
  position <- match(p$term_index, naive$qr$pivot)
  influence <- as.vector(p$analysis_x[, naive$qr$pivot, drop = FALSE] %*% chol2inv(qr.R(naive$qr))[position, ]) * naive$residuals
  aux <- cbind(1, influence * nrow(p$analysis_x))
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
  x <- p$analysis_x[ref, , drop = FALSE]
  fit <- stats::lm.wfit(x, y[ref], weight)
  k <- p$term_index
  psi <- (x * fit$residuals) %*% solve(crossprod(x * sqrt(weight)))
  b <- stats::lm.wfit(a, psi[, k], weight)$coefficients
  b[is.na(b)] <- 0
  g <- psi[, k] - as.vector(a %*% b)
  cluster <- p$cluster[ref]
  score <- rowsum(weight * psi[, k], cluster, reorder = FALSE)
  mean_weight <- rowsum(weight, cluster, reorder = FALSE) / rowsum(rep(1, length(weight)), cluster, reorder = FALSE)
  residual <- rowsum(weight * g, cluster, reorder = FALSE)
  cluster_pi <- rowsum(pi_i[ref], cluster, reorder = FALSE) / rowsum(rep(1, length(weight)), cluster, reorder = FALSE)
  list(estimate = unname(fit$coefficients[k]), se = sqrt(sum(weight * psi[, k]^2) + sum((1 - pi_i[ref]) * weight^2 * g^2)),
       df = sum(ref) - ncol(x), n_used = sum(ref), se_cluster = sqrt(sum(score^2 / mean_weight) + sum((1 - cluster_pi) * residual^2)),
       df_cluster = nrow(score) - 1, clusters = nrow(score))
}
make_case <- function(z, fraction_index, r) {
  set.seed(as.integer(config$cohort_base + z$index * 1e6 + r))
  drawn <- lapply(z$units_by_stratum, function(u) u[sample.int(length(u), length(u), replace = TRUE)])
  selection_seed <- as.integer(config$selection_base + z$index * 1e6 + fraction_index * 1e5 + r)
  set.seed(selection_seed)
  fraction <- config$fractions[fraction_index]
  picked <- lapply(drawn, function(u) {
    size <- max(1L, as.integer(round(fraction * length(u))))
    list(chosen = seq_along(u) %in% sample.int(length(u), size), probability = size / length(u))
  })
  parts <- list(); position <- 0L
  for (s in names(drawn)) for (j in seq_along(drawn[[s]])) {
    position <- position + 1L
    rows <- z$rows_by_unit[[drawn[[s]][j]]]
    parts[[position]] <- data.frame(row = rows, cluster = position, audited = picked[[s]]$chosen[j], probability = picked[[s]]$probability)
  }
  index <- do.call(rbind, parts)
  d <- z$pool[index$row, ]
  rownames(d) <- NULL
  gold <- d$lab_true
  gold[!index$audited] <- NA_real_
  input <- data.frame(covariate = d$covariate, stratum = d$stratum, lab_value = d$lab_value, gold_value = gold)
  x <- stats::model.matrix(~ covariate + stratum, data = input)
  list(input = input, truth = d$lab_true, target = z$truth, predictors = c("covariate", "stratum"), tol = z$tol,
       corrupted = abs(d$lab_value - d$lab_true) > z$tol + 1e-8, sampling_probability = index$probability, cluster = index$cluster,
       analysis_x = x, term_index = match("covariate", colnames(x)), selection_seed = selection_seed,
       imputation_seed = as.integer(config$imputation_base + z$index * 1e6 + fraction_index * 1e5 + r))
}
corruption_scores <- function(imp, p, s, r) {
  table <- imp$contamination$corruption
  scored <- which(table$evidence == "posterior_model" & !is.na(table$posterior_corruption))
  if (!length(scored)) return(NULL)
  y <- p$corrupted[scored]; q <- table$posterior_corruption[scored]
  ranks <- rank(q)
  positives <- sum(y); negatives <- sum(!y)
  data.frame(level = s$level, analysis = s$analysis, scenario = s$scenario, replicate = r, rows = length(scored),
    actual_errors = positives, expected_errors = sum(q), brier = mean((q - y)^2), brier_prevalence = mean((mean(y) - y)^2),
    auc = if (positives && negatives) (sum(ranks[y]) - positives * (positives + 1) / 2) / (positives * negatives) else NA_real_,
    completion_rmse = sqrt(mean((table$posterior_target[scored] - p$truth[scored])^2)),
    recorded_rmse = sqrt(mean((p$input$lab_value[scored] - p$truth[scored])^2)))
}
methods <- c("oracle", "recorded", "audited_only", "linear_mi", "raking", "contamination")
run_method <- function(method, p, s, r) {
  if (method == "oracle") return(list(fit = ols_term(p$truth, p)))
  if (method == "recorded") return(list(fit = ols_term(p$input$lab_value, p)))
  if (method == "audited_only") return(list(fit = ols_term(p$input$gold_value, p)))
  if (method == "raking") return(list(fit = fit_raking(p)))
  imp <- if (method == "linear_mi") {
    catmech_reference(p$input, "lab_value", "gold_value", p$predictors, method = "normal", m = config$m, seed = p$imputation_seed, use_recorded = TRUE)
  } else {
    catmech_contamination(p$input, "lab_value", "gold_value", p$predictors, agreement_tol = p$tol, m = config$m, seed = p$imputation_seed)
  }
  out <- list(fit = pool_fits(lapply(imp$completed, function(d) ols_term(d$lab_value, p))))
  if (method == "contamination") out$scores <- corruption_scores(imp, p, s, r)
  out
}
run_replicate <- function(s, r) {
  z <- analyses[[s$index]]
  p <- make_case(z, s$fraction_index, r)
  records <- list(); scores <- NULL; oracle_estimate <- NA_real_
  for (method in methods) {
    warnings <- character()
    result <- tryCatch(withCallingHandlers(run_method(method, p, s, r), warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning")
    }), error = identity)
    failed <- inherits(result, "error")
    status <- if (!failed) "ok" else if (inherits(result, "catmech_convergence_error")) "convergence_failure" else
      if (grepl("No calibration rows support|full column rank|more observations than coefficients|Fewer trusted/reference|lost design support|rank deficient|constant predictor", conditionMessage(result))) "refused" else "error"
    row <- data.frame(level = s$level, analysis = s$analysis, scenario = s$scenario, fraction = s$fraction, replicate = r,
      selection_seed = p$selection_seed, imputation_seed = p$imputation_seed, method = method, truth = p$target,
      audited_rows = sum(!is.na(p$input$gold_value)), n = nrow(p$input), estimate = NA_real_, se = NA_real_, df = NA_real_,
      lower = NA_real_, upper = NA_real_, width = NA_real_, coverage = NA_real_, se_cluster = NA_real_, df_cluster = NA_real_,
      lower_cluster = NA_real_, upper_cluster = NA_real_, width_cluster = NA_real_, coverage_cluster = NA_real_,
      error = NA_real_, recovery_error = NA_real_, status = status, message = if (failed) conditionMessage(result) else "",
      warnings = paste(unique(warnings), collapse = " | "), stringsAsFactors = FALSE)
    if (!failed) {
      fit <- result$fit
      if (method == "oracle") oracle_estimate <- fit$estimate
      if (!is.finite(fit$estimate) || !is.finite(fit$se) || fit$se <= 0 || is.na(fit$df) || fit$df <= 0 ||
          !is.finite(fit$se_cluster) || fit$se_cluster <= 0 || is.na(fit$df_cluster) || fit$df_cluster <= 0) {
        row$status <- "error"; row$message <- "Nonfinite estimate or degenerate uncertainty."
      } else {
        half <- stats::qt(.975, fit$df) * fit$se
        half_cluster <- stats::qt(.975, fit$df_cluster) * fit$se_cluster
        row$estimate <- fit$estimate; row$se <- fit$se; row$df <- fit$df
        row$lower <- fit$estimate - half; row$upper <- fit$estimate + half; row$width <- 2 * half
        row$coverage <- as.numeric(row$lower <= p$target && row$upper >= p$target)
        row$se_cluster <- fit$se_cluster; row$df_cluster <- fit$df_cluster
        row$lower_cluster <- fit$estimate - half_cluster; row$upper_cluster <- fit$estimate + half_cluster; row$width_cluster <- 2 * half_cluster
        row$coverage_cluster <- as.numeric(row$lower_cluster <= p$target && row$upper_cluster >= p$target)
        row$error <- fit$estimate - p$target; row$recovery_error <- fit$estimate - oracle_estimate
      }
      if (!is.null(result$scores)) scores <- result$scores
    }
    records[[length(records) + 1L]] <- row
  }
  list(records = do.call(rbind, records), scores = scores)
}

pool_rows <- do.call(rbind, lapply(analyses, function(z) {
  d <- z$pool
  x <- stats::model.matrix(~ covariate + stratum, data = d)
  fit <- stats::lm.fit(x, d$lab_true)
  error <- d$lab_value - d$lab_true
  wrong <- abs(error) > z$tol + 1e-8
  data.frame(level = z$level, analysis = z$id, rows = nrow(d), units = length(z$rows_by_unit), strata = length(z$strata),
    target_coefficient = z$truth, recorded_coefficient = z$recorded_coefficient, recorded_minus_target = z$recorded_coefficient - z$truth,
    agreeing_rows = sum(!wrong), wrong_rows = sum(wrong), recorded_above_reference = sum(error > z$tol + 1e-8), recorded_below_reference = sum(error < -z$tol - 1e-8),
    error_sd = stats::sd(error), residual_sd = sqrt(sum(fit$residuals^2) / (nrow(d) - ncol(x))), error_to_residual = stats::sd(error) / sqrt(sum(fit$residuals^2) / (nrow(d) - ncol(x))))
}))
utils::write.csv(pool_rows, file.path(outdir, "pool-targets.csv"), row.names = FALSE)

scenario_list <- list()
for (i in seq_along(analyses)) for (j in seq_along(config$fractions)) {
  z <- analyses[[i]]
  scenario_list[[length(scenario_list) + 1L]] <- list(index = i, fraction_index = j, level = z$level, analysis = z$id, reps = z$reps,
    fraction = config$fractions[j], scenario = paste0(z$level, ":", z$id, ":audit_", 100 * config$fractions[j], "pct"))
}
records <- list(); scores <- list(); start <- Sys.time()
for (s in scenario_list) {
  for (chunk in split(seq_len(s$reps), ceiling(seq_len(s$reps) / 140L))) {
    batch <- parallel::mclapply(chunk, function(r) run_replicate(s, r), mc.cores = workers, mc.preschedule = TRUE, mc.set.seed = FALSE)
    if (any(vapply(batch, inherits, logical(1), what = "try-error"))) stop("A worker failed; checkpoint is not verified evidence.")
    records <- c(records, lapply(batch, `[[`, "records")); scores <- c(scores, lapply(batch, `[[`, "scores"))
    saveRDS(list(records = records, scores = scores), file.path(outdir, "checkpoint.rds"))
    cat(s$scenario, max(chunk), "/", s$reps, "; elapsed", round(as.numeric(difftime(Sys.time(), start, units = "secs"))), "seconds\n"); flush.console()
  }
}
if (!identical(hashes_before, hash_sources())) stop("Source, plan or pool changed during run. Checkpoint retained; evidence is not verified.")
results <- do.call(rbind, records)
score_table <- do.call(rbind, scores)
utils::write.csv(results, file.path(outdir, "replicates.csv"), row.names = FALSE)
utils::write.csv(score_table, file.path(outdir, "corruption-scores.csv"), row.names = FALSE)
mcse <- function(x) if (length(x) > 1L) stats::sd(x) / sqrt(length(x)) else NA_real_
wilson <- function(k, n) {
  if (!n) return(c(NA_real_, NA_real_))
  p <- k / n; z <- stats::qnorm(.975); denom <- 1 + z^2 / n
  centre <- (p + z^2 / (2 * n)) / denom; half <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / denom
  c(max(0, centre - half), min(1, centre + half))
}
groups <- split(results, interaction(results$scenario, results$method, drop = TRUE))
summary_table <- do.call(rbind, lapply(groups, function(g) {
  good <- g[g$status == "ok", ]; n <- nrow(good); requested <- nrow(g)
  ci <- wilson(sum(good$coverage), n); ci_cluster <- wilson(sum(good$coverage_cluster), n)
  data.frame(level = g$level[1], analysis = g$analysis[1], scenario = g$scenario[1], fraction = g$fraction[1], method = g$method[1],
    requested = requested, returned_intervals = n, refused = sum(g$status == "refused"), convergence_failures = sum(g$status == "convergence_failure"),
    other_errors = sum(g$status == "error"), interval_return_rate = n / requested, warned = sum(nzchar(g$warnings)),
    mean_audited_rows = mean(g$audited_rows), mean_rows = mean(g$n), target = g$truth[1],
    recovery_rmse = if (n) sqrt(mean(good$recovery_error^2)) else NA_real_, recovery_bias = if (n) mean(good$recovery_error) else NA_real_,
    bias = if (n) mean(good$error) else NA_real_, bias_mcse = mcse(good$error), rmse = if (n) sqrt(mean(good$error^2)) else NA_real_,
    empirical_sd = if (n > 1) stats::sd(good$estimate) else NA_real_,
    mean_se = if (n) mean(good$se) else NA_real_, mean_width = if (n) mean(good$width) else NA_real_, coverage = if (n) mean(good$coverage) else NA_real_,
    coverage_wilson_lower = ci[1], coverage_wilson_upper = ci[2],
    mean_se_cluster = if (n) mean(good$se_cluster) else NA_real_, mean_width_cluster = if (n) mean(good$width_cluster) else NA_real_,
    coverage_cluster = if (n) mean(good$coverage_cluster) else NA_real_, coverage_cluster_wilson_lower = ci_cluster[1], coverage_cluster_wilson_upper = ci_cluster[2])
}))
rownames(summary_table) <- NULL
utils::write.csv(summary_table, file.path(outdir, "summary.csv"), row.names = FALSE)
set.seed(config$ratio_seed)
pairs <- list(c("linear_mi", "recorded"), c("contamination", "linear_mi"), c("contamination", "recorded"), c("raking", "recorded"),
              c("raking", "linear_mi"), c("audited_only", "recorded"), c("contamination", "raking"))
ratio_rows <- list()
for (s in scenario_list) {
  g <- results[results$scenario == s$scenario & results$status == "ok", ]
  for (pair in pairs) {
    both <- merge(g[g$method == pair[1], c("replicate", "recovery_error")], g[g$method == pair[2], c("replicate", "recovery_error")],
                  by = "replicate", suffixes = c("_a", "_b"))
    if (!nrow(both)) next
    a <- both$recovery_error_a^2; b <- both$recovery_error_b^2
    resampled <- replicate(config$ratio_resamples, { i <- sample.int(length(a), replace = TRUE); sqrt(sum(a[i]) / sum(b[i])) })
    limits <- stats::quantile(resampled, c(.025, .975), names = FALSE, na.rm = TRUE)
    ratio_rows[[length(ratio_rows) + 1L]] <- data.frame(level = s$level, analysis = s$analysis, scenario = s$scenario, fraction = s$fraction,
      method = pair[1], comparator = pair[2], complete_pairs = nrow(both), recovery_rmse_ratio = sqrt(sum(a) / sum(b)), lower = limits[1], upper = limits[2])
  }
}
utils::write.csv(do.call(rbind, ratio_rows), file.path(outdir, "recovery-ratios.csv"), row.names = FALSE)
utils::write.csv(hashes_before, file.path(outdir, "source-hashes.csv"), row.names = FALSE)
writeLines(c(paste("Completed UTC:", format(Sys.time(), tz = "UTC", usetz = TRUE)),
  paste("Elapsed seconds:", round(as.numeric(difftime(Sys.time(), start, units = "secs")), 1)), capture.output(sessionInfo())), file.path(outdir, "session-info.txt"))
unlink(file.path(outdir, "checkpoint.rds"))
cat("Completed", mode, ":", nrow(results), "requests;", paste(names(table(results$status)), table(results$status), collapse = ", "), "\n")
cat("Output:", normalizePath(outdir), "\n")
