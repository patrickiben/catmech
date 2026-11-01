#!/usr/bin/env Rscript
# Prospective real-data study; settings locked in validation/REAL_DATA_PLAN.md.
# Usage: Rscript tools/real-data-benchmark.R pilot|study [NEW_DIRECTORY] [PILOT_REPS]
args <- commandArgs(trailingOnly = TRUE)
mode <- if (length(args)) match.arg(args[1], c("pilot", "study")) else "pilot"
file_arg <- grep("^--file=", commandArgs(), value = TRUE)
root <- normalizePath(file.path(dirname(sub("^--file=", "", file_arg[1])), ".."))
outdir <- if (length(args) >= 2L) args[2] else file.path(root, "validation", paste0("real-data-", mode))
if (dir.exists(outdir) && length(list.files(outdir, all.files = TRUE, no.. = TRUE)))
  stop("Use a new or empty output directory; existing evidence is never overwritten.")
hash_sources <- function() {
  paths <- c(sort(list.files(file.path(root, "R"), "\\.R$", full.names = TRUE)),
    file.path(root, c("DESCRIPTION", "NAMESPACE", "tools/real-data-benchmark.R", "validation/REAL_DATA_PLAN.md",
                      "validation/real-data/ctgov-pool.csv", "validation/real-data/nhanes-pool.csv")))
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
config <- list(mode = mode, workers = workers, m = 20L, n = 600L,
  primary_reps = if (mode == "pilot") pilot_reps else 1000L, other_reps = if (mode == "pilot") pilot_reps else 200L,
  holdout_reps = if (mode == "pilot") min(pilot_reps, 2L) else 40L, holdout_folds = 5L, holdout_m = 40L,
  split = if (mode == "pilot") "development" else "evaluation",
  cohort_base = if (mode == "pilot") 300100000 else 200100000,
  selection_base = if (mode == "pilot") 310100000 else 210100000,
  imputation_base = if (mode == "pilot") 320000000 else 220000000, nhanes_offset = 5000000)
jsonlite::write_json(config, file.path(outdir, "settings.json"), pretty = TRUE, auto_unbox = TRUE)

load_pool <- function(dataset) {
  if (dataset == "ctgov") {
    d <- utils::read.csv(file.path(root, "validation/real-data/ctgov-pool.csv"), stringsAsFactors = FALSE)
    d <- d[d$split == config$split, ]
    data.frame(id = d$nct_id, outcome = log(d$duration_months), lab_true = log(d$reference_n), lab_value = log(d$recorded_n),
      phase = factor(d$phase, levels = c("not_applicable", "early", "phase2", "phase3_4")),
      sponsor = factor(d$sponsor, levels = c("other", "industry")), log_sites = log(d$sites), start_year = d$start_year)
  } else {
    d <- utils::read.csv(file.path(root, "validation/real-data/nhanes-pool.csv"), stringsAsFactors = FALSE)
    d <- d[d$split == config$split, ]
    data.frame(id = d$seqn, outcome = d$sbp, lab_true = d$reference_kg, lab_value = d$recorded_kg, age = d$age,
      sex = factor(d$sex, levels = c("female", "male")),
      race = factor(d$race, levels = c("white", "black", "mexican_american", "other_hispanic", "other")), height = d$height)
  }
}
datasets <- list(
  ctgov = list(pool = load_pool("ctgov"), covariates = c("phase", "sponsor", "log_sites", "start_year"), offset = 0,
               tol = 0, strata = "sponsor", audit_by = c("phase", "sponsor")),
  nhanes = list(pool = load_pool("nhanes"), covariates = c("age", "sex", "race", "height"), offset = config$nhanes_offset,
                tol = 1, strata = "sex", audit_by = c("sex", "race")))
for (name in names(datasets)) {
  z <- datasets[[name]]
  fit <- stats::lm(stats::reformulate(c(z$covariates, "lab_true"), response = "outcome"), data = z$pool)
  datasets[[name]]$truth <- unname(stats::coef(fit)["lab_true"])
  datasets[[name]]$truth_se <- unname(sqrt(stats::vcov(fit)["lab_true", "lab_true"]))
  datasets[[name]]$pool_rows <- nrow(z$pool)
}
utils::write.csv(do.call(rbind, lapply(names(datasets), function(name) with(datasets[[name]], data.frame(dataset = name,
  split = config$split, pool_rows = pool_rows, target_coefficient = truth, target_se = truth_se)))),
  file.path(outdir, "pool-targets.csv"), row.names = FALSE)
pool_audits <- lapply(names(datasets), function(name) {
  z <- datasets[[name]]
  catmech_recording_audit(z$pool, "lab_value", "lab_true", by = z$audit_by, agreement_tol = z$tol)
})
names(pool_audits) <- names(datasets)
utils::write.csv(do.call(rbind, lapply(names(pool_audits), function(name) cbind(dataset = name, status = pool_audits[[name]]$status, pool_audits[[name]]$overall))),
  file.path(outdir, "pool-audit-overall.csv"), row.names = FALSE)
utils::write.csv(do.call(rbind, lapply(names(pool_audits), function(name) cbind(dataset = name, pool_audits[[name]]$by))),
  file.path(outdir, "pool-audit-by-group.csv"), row.names = FALSE)
utils::write.csv(do.call(rbind, lapply(names(pool_audits), function(name) cbind(dataset = name, pool_audits[[name]]$heterogeneity))),
  file.path(outdir, "pool-audit-heterogeneity.csv"), row.names = FALSE)
cells <- function(dataset) {
  fractions <- c(audit_25pct = .25, audit_10pct = .10, audit_50pct = .50, if (dataset == "ctgov") c(audit_05pct = .05))
  out <- lapply(names(fractions), function(name) list(dataset = dataset, cell = name, fraction = fractions[[name]], targeted = FALSE,
    reps = if (name == "audit_25pct") config$primary_reps else config$other_reps))
  c(out, list(list(dataset = dataset, cell = "targeted_audit", fraction = NA_real_, targeted = TRUE, reps = config$other_reps)))
}
scenario_list <- c(cells("ctgov"), cells("nhanes"))
matching <- c("audit_25pct", "audit_10pct", "audit_50pct")
all_methods <- c("oracle", "naive_recorded", "gold_only", "reference_normal", "cole_mi", "raking", "contamination", "contamination_strata")

make_case <- function(s, r) {
  z <- datasets[[s$dataset]]
  set.seed(as.integer(config$cohort_base + z$offset + r))
  d <- z$pool[sample.int(nrow(z$pool), config$n, replace = TRUE), ]
  rownames(d) <- NULL
  selection_seed <- as.integer(config$selection_base + z$offset + r)
  set.seed(selection_seed)
  u <- stats::runif(nrow(d))
  probability <- rep(s$fraction, nrow(d))
  if (s$targeted) {
    centre <- (d$lab_value - mean(d$lab_value, na.rm = TRUE)) / stats::sd(d$lab_value, na.rm = TRUE)
    probability <- ifelse(is.na(d$lab_value), .15, ifelse(abs(centre) > 1.25, .60, .15))
  }
  gold <- d$lab_true
  gold[u >= probability] <- NA_real_
  predictors <- c(z$covariates, "outcome")
  input <- d[predictors]
  input$lab_value <- d$lab_value
  input$flag <- factor(rep("No", nrow(d)), levels = c("No", "Yes"))
  input$gold_value <- gold
  list(input = input, truth = d$lab_true, target = z$truth, predictors = predictors, tol = z$tol, strata = z$strata,
       corrupted = !is.na(d$lab_value) & abs(d$lab_value - d$lab_true) > z$tol + 1e-8,
       sampling_probability = probability,
       analysis_x = stats::model.matrix(stats::reformulate(z$covariates), data = input),
       model_x = stats::model.matrix(stats::reformulate(predictors), data = input)[, -1L, drop = FALSE],
       selection_seed = selection_seed, imputation_seed = as.integer(config$imputation_base + z$offset + r),
       stream = paste0(s$dataset, ":", config$split, ":", r))
}
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
fit_cole <- function(p) {
  y <- p$input$gold_value; w <- p$input$lab_value
  ref <- !is.na(y); obs <- !is.na(w)
  w0 <- ifelse(obs, w, 0)
  full_rank <- function(x, ry) {
    q <- qr(cbind(1, x[ry, , drop = FALSE]))
    x[, setdiff(q$pivot[seq_len(q$rank)], 1L) - 1L, drop = FALSE]
  }
  x_obs <- full_rank(cbind(p$model_x, recorded = w0), obs & ref)
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
fit_raking <- function(p) {
  y <- p$input$gold_value; w <- p$input$lab_value; outcome <- p$input$outcome
  ref <- !is.na(y); obs <- !is.na(w); pi_i <- p$sampling_probability
  x_obs <- cbind(1, p$model_x, recorded = ifelse(obs, w, 0))
  x_mis <- cbind(1, p$model_x)
  proxy <- numeric(length(y))
  use <- which(obs & ref)
  beta <- stats::lm.wfit(x_obs[use, , drop = FALSE], y[use], 1 / pi_i[use])$coefficients
  beta[is.na(beta)] <- 0
  proxy[obs] <- as.vector(x_obs[obs, , drop = FALSE] %*% beta)
  if (any(!obs)) {
    use <- which(ref)
    beta <- stats::lm.wfit(x_mis[use, , drop = FALSE], y[use], 1 / pi_i[use])$coefficients
    beta[is.na(beta)] <- 0
    proxy[!obs] <- as.vector(x_mis[!obs, , drop = FALSE] %*% beta)
  }
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
  k <- ncol(x)
  psi <- (x * fit$residuals) %*% solve(crossprod(x * sqrt(weight)))
  b <- stats::lm.wfit(a, psi[, k], weight)$coefficients
  g <- psi[, k] - as.vector(a %*% b)
  list(estimate = unname(fit$coefficients[k]), se = sqrt(sum(weight * psi[, k]^2) + sum((1 - pi_i[ref]) * weight^2 * g^2)),
       df = sum(ref) - k, n_used = sum(ref))
}
corruption_scores <- function(imp, p, s, r, method) {
  table <- imp$contamination$corruption
  scored <- which(table$evidence == "posterior_model" & !is.na(table$posterior_corruption))
  if (!length(scored)) return(NULL)
  y <- p$corrupted[scored]; q <- table$posterior_corruption[scored]
  ranks <- rank(q)
  positives <- sum(y); negatives <- sum(!y)
  data.frame(family = s$dataset, scenario = s$cell, method = method, replicate = r, rows = length(scored),
    actual_errors = positives, expected_errors = sum(q), brier = mean((q - y)^2), brier_prevalence = mean((mean(y) - y)^2),
    auc = if (positives && negatives) (sum(ranks[y]) - positives * (positives + 1) / 2) / (positives * negatives) else NA_real_,
    completion_rmse = sqrt(mean((table$posterior_target - p$truth)^2)),
    recorded_rmse = sqrt(mean((p$input$lab_value[scored] - p$truth[scored])^2)))
}
run_method <- function(method, p, s, r) {
  w <- p$input$lab_value
  if (method == "oracle") return(list(fit = fit_lab(p$truth, p)))
  if (method == "gold_only") return(list(fit = fit_lab(p$input$gold_value, p)))
  if (method == "cole_mi") return(list(fit = fit_cole(p)))
  if (method == "raking") return(list(fit = fit_raking(p)))
  if (method == "naive_recorded" && !anyNA(w)) return(list(fit = fit_lab(w, p)))
  contamination_call <- function(strata, tol) catmech_contamination(p$input, "lab_value", "gold_value", p$predictors,
    corruption_strata = strata, agreement_tol = tol, m = config$m, seed = p$imputation_seed)
  imp <- switch(method,
    naive_recorded = catmech_model(p$input, "lab_value", p$predictors, category_mode = "ignore", m = config$m, seed = p$imputation_seed),
    reference_normal = catmech_reference(p$input, "lab_value", "gold_value", p$predictors, method = "normal", m = config$m, seed = p$imputation_seed),
    reference_pmm = catmech_reference(p$input, "lab_value", "gold_value", p$predictors, method = "pmm", m = config$m, seed = p$imputation_seed),
    reference_pmm_abb = catmech_reference(p$input, "lab_value", "gold_value", p$predictors, method = "pmm_abb", m = config$m, seed = p$imputation_seed),
    contamination = contamination_call(NULL, p$tol),
    contamination_strata = contamination_call(p$strata, p$tol),
    contamination_no_tol = contamination_call(NULL, 0))
  out <- list(fit = pool_fits(lapply(imp$completed, function(d) fit_lab(d$lab_value, p))))
  if (startsWith(method, "contamination")) {
    out$scores <- corruption_scores(imp, p, s, r, method)
    if (method == "contamination" && s$cell == "audit_25pct" && r <= config$holdout_reps) {
      h <- catmech_holdout_check(imp, folds = config$holdout_folds, m = config$holdout_m, seed = p$imputation_seed)
      row <- if (is.null(h$summary)) data.frame(predictive_coverage = NA_real_, ks_p_value = NA_real_, calibration_z = NA_real_,
        errors_observed = NA_real_, errors_expected = NA_real_, brier_model = NA_real_, brier_stratum_rate = NA_real_, rmse_model = NA_real_, rmse_recorded = NA_real_) else
        h$summary[c("predictive_coverage", "ks_p_value", "calibration_z", "errors_observed", "errors_expected", "brier_model", "brier_stratum_rate", "rmse_model", "rmse_recorded")]
      out$holdout <- cbind(data.frame(family = s$dataset, scenario = s$cell, replicate = r, status = h$status,
        flags = paste(h$flags, collapse = "|"), failures = length(h$failures), stringsAsFactors = FALSE), row)
    }
  }
  out
}
run_replicate <- function(s, r) {
  p <- make_case(s, r)
  records <- list(); scores <- list(); holdout <- NULL
  methods <- c(all_methods, if (s$cell %in% matching) c("reference_pmm", "reference_pmm_abb"),
               if (s$dataset == "nhanes" && s$cell == "audit_25pct") "contamination_no_tol")
  for (method in methods) {
    warnings <- character()
    result <- tryCatch(withCallingHandlers(run_method(method, p, s, r), warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning")
    }), error = identity)
    failed <- inherits(result, "error")
    status <- if (!failed) "ok" else if (inherits(result, "catmech_convergence_error")) "convergence_failure" else
      if (grepl("No calibration rows support|full column rank|more observations than coefficients|Fewer trusted/reference|lost design support|rank deficient|constant predictor", conditionMessage(result))) "refused" else "error"
    row <- data.frame(family = s$dataset, scenario = s$cell, replicate = r, data_seed = p$selection_seed,
      data_stream = p$stream, imputation_seed = p$imputation_seed, method = method, truth = p$target,
      available_n = sum(!is.na(p$input$gold_value)), n = nrow(p$input), estimate = NA_real_, se = NA_real_, df = NA_real_, n_used = NA_real_,
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
        row$coverage <- as.numeric(row$lower <= p$target && row$upper >= p$target); row$error <- fit$estimate - p$target
      }
      if (!is.null(result$scores)) scores[[length(scores) + 1L]] <- result$scores
      if (!is.null(result$holdout)) holdout <- result$holdout
    }
    records[[length(records) + 1L]] <- row
  }
  list(records = do.call(rbind, records), scores = do.call(rbind, scores), holdout = holdout)
}
records <- list(); scores <- list(); holdouts <- list(); start <- Sys.time()
for (s in scenario_list) {
  for (chunk in split(seq_len(s$reps), ceiling(seq_len(s$reps) / 140L))) {
    batch <- parallel::mclapply(chunk, function(r) run_replicate(s, r), mc.cores = workers, mc.preschedule = TRUE, mc.set.seed = FALSE)
    if (any(vapply(batch, inherits, logical(1), what = "try-error"))) stop("A worker failed; checkpoint is not verified evidence.")
    records <- c(records, lapply(batch, `[[`, "records")); scores <- c(scores, lapply(batch, `[[`, "scores"))
    holdouts <- c(holdouts, lapply(batch, `[[`, "holdout"))
    saveRDS(list(records = records, scores = scores, holdouts = holdouts), file.path(outdir, "checkpoint.rds"))
    cat(s$dataset, s$cell, max(chunk), "/", s$reps, "; elapsed", round(as.numeric(difftime(Sys.time(), start, units = "secs"))), "seconds\n"); flush.console()
  }
}
if (!identical(hashes_before, hash_sources())) stop("Source, plan or pool changed during run. Checkpoint retained; evidence is not verified.")
results <- do.call(rbind, records)
score_table <- do.call(rbind, scores)
holdout_table <- do.call(rbind, holdouts)
utils::write.csv(results, file.path(outdir, "replicates.csv"), row.names = FALSE)
utils::write.csv(score_table, file.path(outdir, "corruption-scores.csv"), row.names = FALSE)
utils::write.csv(holdout_table, file.path(outdir, "holdout-checks.csv"), row.names = FALSE)
mcse <- function(x) if (length(x) > 1L) stats::sd(x) / sqrt(length(x)) else NA_real_
wilson <- function(k, n) {
  if (!n) return(c(NA_real_, NA_real_))
  p <- k / n; z <- stats::qnorm(.975); denom <- 1 + z^2 / n
  centre <- (p + z^2 / (2 * n)) / denom; half <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / denom
  c(max(0, centre - half), min(1, centre + half))
}
groups <- split(results, interaction(results$family, results$scenario, results$method, drop = TRUE))
summary_table <- do.call(rbind, lapply(groups, function(g) {
  good <- g[g$status == "ok", ]; n <- nrow(good); requested <- nrow(g)
  coverage <- if (n) mean(good$coverage) else NA_real_; ci <- wilson(sum(good$coverage), n)
  rmse <- if (n) sqrt(mean(good$error^2)) else NA_real_
  data.frame(family = g$family[1], scenario = g$scenario[1], method = g$method[1], requested = requested,
    returned_intervals = n, refused = sum(g$status == "refused"), convergence_failures = sum(g$status == "convergence_failure"),
    other_errors = sum(g$status == "error"), interval_return_rate = n / requested, warned = sum(nzchar(g$warnings)),
    mean_available_n = mean(g$available_n), target = g$truth[1], bias = if (n) mean(good$error) else NA_real_, bias_mcse = mcse(good$error), rmse = rmse,
    empirical_sd = if (n > 1) stats::sd(good$estimate) else NA_real_, mean_se = if (n) mean(good$se) else NA_real_, mean_se_mcse = mcse(good$se),
    mean_width = if (n) mean(good$width) else NA_real_, width_mcse = mcse(good$width),
    coverage_conditional = coverage, coverage_mcse = if (n) sqrt(coverage * (1 - coverage) / n) else NA_real_,
    coverage_wilson_lower = ci[1], coverage_wilson_upper = ci[2], covering_interval_per_request = sum(good$coverage) / requested)
}))
rownames(summary_table) <- NULL
utils::write.csv(summary_table, file.path(outdir, "summary.csv"), row.names = FALSE)
status_rows <- do.call(rbind, lapply(groups, function(g) do.call(rbind, lapply(c("ok", "refused", "convergence_failure", "error"), function(status) {
  k <- sum(g$status == status); n <- nrow(g); ci <- wilson(k, n); rate <- k / n
  data.frame(family = g$family[1], scenario = g$scenario[1], method = g$method[1], status = status,
    count = k, requested = n, rate = rate, mcse = sqrt(rate * (1 - rate) / n), lower = ci[1], upper = ci[2])
}))))
utils::write.csv(status_rows, file.path(outdir, "request-status.csv"), row.names = FALSE)
paired <- list()
for (s in scenario_list) {
  g <- results[results$family == s$dataset & results$scenario == s$cell & results$status == "ok", ]
  for (comparator in c("oracle", "naive_recorded", "reference_normal", "cole_mi", "raking", "gold_only")) {
    both <- merge(g[g$method == "contamination", ], g[g$method == comparator, ], by = "replicate", suffixes = c("_new", "_baseline"))
    if (!nrow(both)) next
    de <- both$error_new - both$error_baseline; dm <- both$error_new^2 - both$error_baseline^2
    dc <- both$coverage_new - both$coverage_baseline; dw <- both$width_new - both$width_baseline
    paired[[length(paired) + 1L]] <- data.frame(family = s$dataset, scenario = s$cell, method = "contamination", comparator = comparator,
      complete_pairs = nrow(both), error_difference = mean(de), error_mcse = mcse(de), mse_difference = mean(dm), mse_mcse = mcse(dm),
      coverage_difference = mean(dc), coverage_difference_mcse = mcse(dc), width_difference = mean(dw), width_difference_mcse = mcse(dw))
  }
}
utils::write.csv(do.call(rbind, paired), file.path(outdir, "paired-comparisons.csv"), row.names = FALSE)
utils::write.csv(hashes_before, file.path(outdir, "source-hashes.csv"), row.names = FALSE)
writeLines(c(paste("Completed UTC:", format(Sys.time(), tz = "UTC", usetz = TRUE)),
  paste("Elapsed seconds:", round(as.numeric(difftime(Sys.time(), start, units = "secs")), 1)), capture.output(sessionInfo())), file.path(outdir, "session-info.txt"))
unlink(file.path(outdir, "checkpoint.rds"))
cat("Completed", mode, ":", nrow(results), "requests;", paste(names(table(results$status)), table(results$status), collapse = ", "), "\n")
cat("Output:", normalizePath(outdir), "\n")
