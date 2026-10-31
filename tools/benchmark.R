#!/usr/bin/env Rscript
# Exploratory, paired simulation. See validation/BENCHMARK_PLAN.md first.
args <- commandArgs(trailingOnly = TRUE)
parse_count <- function(i, default) {
  val <- if (length(args) >= i) suppressWarnings(as.numeric(args[i])) else default
  if (length(val) != 1L || !is.finite(val) || val < 2 || val != floor(val)) {
    stop("Replicate and imputation counts must be integers of at least 2.")
  }
  as.integer(val)
}
reps <- parse_count(1L, 100L)
m <- parse_count(2L, 10L)
file_arg <- grep("^--file=", commandArgs(), value = TRUE)
root <- normalizePath(file.path(dirname(sub("^--file=", "", file_arg[1])), ".."))
outdir <- if (length(args) >= 3L) args[3] else file.path(root, "validation", "benchmark")
if (!requireNamespace("pkgload", quietly = TRUE)) stop("Install pkgload to run this source benchmark.")
if (dir.exists(outdir) && length(list.files(outdir, all.files = TRUE, no.. = TRUE))) {
  stop("Output directory must be new or empty; supply a new third argument for a fresh run.")
}
source_paths <- function() c(sort(list.files(file.path(root, "R"), pattern = "\\.R$", full.names = TRUE)),
  file.path(root, c("DESCRIPTION", "NAMESPACE")), file.path(root, "tools", "benchmark.R"), file.path(root, "validation", "BENCHMARK_PLAN.md"))
hash_sources <- function() {
  paths <- source_paths()
  data.frame(file = substring(paths, nchar(root) + 2L),
    sha256 = vapply(paths, function(f) digest::digest(file = f, algo = "sha256"), character(1)),
    row.names = NULL)
}
source_hashes_before <- hash_sources()
pkgload::load_all(root, quiet = TRUE)
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

truth <- 0.8
predictors <- c("age", "sex", "inflammation", "outcome")
input_names <- c("patient_id", predictors, "lab_value", "lab_category", "outlier_flag")
analysis_formula <- outcome ~ lab_value + age + sex + inflammation
thresholds <- c(-Inf, 5.5, 7, Inf)
scenarios <- data.frame(
  scenario = c("missing_only", "perfect_flags", "false_positive_flags", "category_error", "heavy_tails"),
  missing_rate = rep(0.2, 5), flag_rate = c(0, rep(0.1, 4)),
  corruption_probability = c(1, 1, 0.5, 1, 1),
  category_error_rate = c(0, 0, 0, 0.1, 0),
  error_distribution = c(rep("normal", 4), "heavy_tailed"),
  stringsAsFactors = FALSE
)
methods <- c("catmech", "ordinary_pmm", "trusted_pmm", "naive_observed", "trusted_observed", "oracle")

fit_one <- function(dat) {
  fit <- stats::lm(analysis_formula, data = dat)
  coef <- summary(fit)$coefficients
  if (!"lab_value" %in% rownames(coef)) stop("Analysis coefficient not estimable.")
  list(estimate = unname(coef["lab_value", "Estimate"]),
       se = unname(coef["lab_value", "Std. Error"]),
       df = stats::df.residual(fit), n_used = stats::nobs(fit))
}
fit_completed <- function(completed) {
  fits <- lapply(completed, fit_one)
  n_used <- vapply(fits, function(x) x$n_used, numeric(1))
  if (length(unique(n_used)) != 1L) stop("Different analysis sizes across imputations.")
  pooled <- mice::pool.scalar(
    Q = vapply(fits, function(x) x$estimate, numeric(1)),
    U = vapply(fits, function(x) x$se^2, numeric(1)),
    n = n_used[1], k = n_used[1] - fits[[1]]$df
  )
  list(estimate = pooled$qbar, se = sqrt(pooled$t), df = pooled$df,
       n_used = n_used[1])
}
run_pmm <- function(dat, trusted, seed) {
  use <- dat[, c("lab_value", predictors, "lab_category", "outlier_flag"), drop = FALSE]
  if (trusted) use$lab_value[as.character(use$outlier_flag) == "Yes"] <- NA_real_
  meth <- stats::setNames(rep("", ncol(use)), names(use)); meth["lab_value"] <- "pmm"
  pred <- matrix(0L, ncol(use), ncol(use), dimnames = list(names(use), names(use)))
  pred["lab_value", c(predictors, "lab_category")] <- 1L
  if (!trusted) pred["lab_value", "outlier_flag"] <- 1L
  imp <- mice::mice(use, method = meth, predictorMatrix = pred,
                    m = m, maxit = 1L, donors = 5L, matchtype = 1L,
                    seed = seed, printFlag = FALSE)
  if (!is.null(imp$loggedEvents)) {
    warning(paste("mice logged events:", paste(unique(imp$loggedEvents$out), collapse = "; ")),
            call. = FALSE)
  }
  fit_completed(lapply(seq_len(m), function(j) mice::complete(imp, j)))
}
run_method <- function(method, dat, scenario, seed) {
  input <- dat[, input_names, drop = FALSE]
  if (method == "oracle") {
    input$lab_value <- dat$lab_true
    return(fit_one(input))
  }
  if (method == "naive_observed") return(fit_one(input))
  if (method == "trusted_observed") {
    input$lab_value[as.character(input$outlier_flag) == "Yes"] <- NA_real_
    return(fit_one(input))
  }
  if (method == "ordinary_pmm") return(run_pmm(input, FALSE, seed))
  if (method == "trusted_pmm") return(run_pmm(input, TRUE, seed))
  spec <- catmech_spec(
    value = "lab_value", category = "lab_category",
    levels = c("Low", "Medium", "High"), thresholds = thresholds,
    flag = "outlier_flag", flag_yes = "Yes", flag_no = "No",
    corruption_probability = scenario$corruption_probability,
    conflict_policy = "reconstruct", min_donors = 1L, n_donors = 5L,
    fallback = "error"
  )
  imp <- catmech_impute(input, spec, predictors = predictors, m = m, seed = seed)
  if (!identical(input, dat[, input_names, drop = FALSE])) stop("Input changed during imputation.")
  for (completed in imp$completed) {
    if (any(!is.finite(completed$lab_value))) stop("Nonfinite completed target.")
    derived <- cut(completed$lab_value, thresholds, labels = c("Low", "Medium", "High"),
                   right = TRUE, include.lowest = TRUE)
    # The released specification declares right-closed intervals.
    if (any(as.character(derived) != as.character(completed$lab_category))) {
      stop("Completed target violates its exact category.")
    }
  }
  fit_completed(imp$completed)
}

records <- vector("list", nrow(scenarios) * reps * length(methods))
idx <- 0L
start <- Sys.time()
for (s in seq_len(nrow(scenarios))) {
  scenario <- scenarios[s, , drop = FALSE]
  for (r in seq_len(reps)) {
    data_seed <- as.integer(100000L + 10000L * s + r)
    dat <- simulate_catmech_lab(
      n = 300, seed = data_seed, thresholds = thresholds,
      corruption_probability = scenario$corruption_probability,
      missing_rate = scenario$missing_rate, flag_rate = scenario$flag_rate,
      category_error_rate = scenario$category_error_rate,
      error_distribution = scenario$error_distribution
    )
    for (method in methods) {
      idx <- idx + 1L
      warnings <- character()
      result <- tryCatch(withCallingHandlers(
        run_method(method, dat, scenario, data_seed + 1000000L),
        warning = function(w) {
          warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning")
        }), error = function(e) e)
      failed <- inherits(result, "error")
      row <- data.frame(scenario = scenario$scenario, replicate = r, data_seed = data_seed,
                        method = method, estimate = NA_real_, se = NA_real_, df = NA_real_,
                        n_used = NA_real_, lower = NA_real_, upper = NA_real_,
                        coverage = NA_real_, error = NA_real_, status = "error",
                        message = if (failed) conditionMessage(result) else "",
                        warnings = paste(unique(warnings), collapse = " | "),
                        stringsAsFactors = FALSE)
      if (!failed) {
        if (!is.finite(result$estimate) || !is.finite(result$se) || result$se < 0 ||
            is.na(result$df) || result$df <= 0) {
          row$message <- "Nonfinite estimate/standard error or invalid degrees of freedom."
        } else {
          row$estimate <- result$estimate; row$se <- result$se
          row$df <- result$df; row$n_used <- result$n_used
          critical <- stats::qt(0.975, df = result$df)
          row$lower <- result$estimate - critical * result$se
          row$upper <- result$estimate + critical * result$se
          row$coverage <- as.numeric(row$lower <= truth && row$upper >= truth)
          row$error <- result$estimate - truth; row$status <- "ok"
        }
      }
      records[[idx]] <- row
    }
    if (r == 1L || r %% 20L == 0L || r == reps) {
      cat(scenario$scenario, r, "/", reps, "elapsed", round(as.numeric(difftime(Sys.time(), start, units = "secs"))), "s\n")
      flush.console()
    }
  }
}
if (!identical(source_hashes_before, hash_sources())) {
  stop("Source or benchmark plan changed during execution; rerun against frozen files. Results were not written.")
}
results <- do.call(rbind, records)
utils::write.csv(results, file.path(outdir, "replicates.csv"), row.names = FALSE)
mcse_mean <- function(x) if (length(x) > 1L) stats::sd(x) / sqrt(length(x)) else NA_real_
wilson_interval <- function(successes, n) {
  if (!n) return(c(NA_real_, NA_real_))
  z <- stats::qnorm(0.975); p <- successes / n; denom <- 1 + z^2 / n
  centre <- (p + z^2 / (2 * n)) / denom
  half <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / denom
  c(max(0, centre - half), min(1, centre + half))
}
groups <- split(results, interaction(results$scenario, results$method, drop = TRUE))
summary_rows <- lapply(groups, function(g) {
  good <- g[g$status == "ok", , drop = FALSE]
  n <- nrow(good); fail <- 1 - n / nrow(g)
  rmse <- if (n) sqrt(mean(good$error^2)) else NA_real_
  coverage <- if (n) mean(good$coverage) else NA_real_
  coverage_interval <- wilson_interval(sum(good$coverage), n)
  failure_interval <- wilson_interval(nrow(g) - n, nrow(g))
  width <- good$upper - good$lower
  data.frame(scenario = g$scenario[1], method = g$method[1], requested = nrow(g),
             successful = n, failures = nrow(g) - n, failure_rate = fail,
             failure_mcse = sqrt(fail * (1 - fail) / nrow(g)),
             failure_ci_lower = failure_interval[1], failure_ci_upper = failure_interval[2],
             bias = if (n) mean(good$error) else NA_real_, bias_mcse = mcse_mean(good$error),
             rmse = rmse, rmse_mcse = if (n > 1L && rmse > 0) mcse_mean(good$error^2) / (2 * rmse) else NA_real_,
             mean_se = if (n) mean(good$se) else NA_real_, mean_se_mcse = mcse_mean(good$se),
             coverage = coverage, coverage_mcse = if (n) sqrt(coverage * (1 - coverage) / n) else NA_real_,
             coverage_ci_lower = coverage_interval[1], coverage_ci_upper = coverage_interval[2],
             mean_interval_width = if (n) mean(width) else NA_real_, interval_width_mcse = mcse_mean(width),
             estimand_truth = truth, metrics_conditioned_on_success = TRUE)
})
summary_table <- do.call(rbind, summary_rows)
rownames(summary_table) <- NULL
utils::write.csv(summary_table, file.path(outdir, "summary.csv"), row.names = FALSE)
paired <- list(); pidx <- 0L
for (scenario in scenarios$scenario) {
  ref <- results[results$scenario == scenario & results$method == "catmech" & results$status == "ok", ]
  for (method in setdiff(methods, "catmech")) {
    other <- results[results$scenario == scenario & results$method == method & results$status == "ok", ]
    both <- merge(ref, other, by = "replicate", suffixes = c("_cat", "_other"))
    d_error <- both$error_cat - both$error_other
    d_mse <- both$error_cat^2 - both$error_other^2
    d_cover <- both$coverage_cat - both$coverage_other
    pidx <- pidx + 1L
    paired[[pidx]] <- data.frame(scenario = scenario, comparator = method, complete_pairs = nrow(both),
      error_difference = if (nrow(both)) mean(d_error) else NA_real_, error_difference_mcse = mcse_mean(d_error),
      mse_difference = if (nrow(both)) mean(d_mse) else NA_real_, mse_difference_mcse = mcse_mean(d_mse),
      coverage_difference = if (nrow(both)) mean(d_cover) else NA_real_, coverage_difference_mcse = mcse_mean(d_cover))
  }
}
utils::write.csv(do.call(rbind, paired), file.path(outdir, "paired_comparisons.csv"), row.names = FALSE)
utils::write.csv(scenarios, file.path(outdir, "scenarios.csv"), row.names = FALSE)
writeLines(c(paste("Run completed:", format(Sys.time(), tz = "UTC", usetz = TRUE)),
             paste("Replicates per scenario:", reps), paste("Imputations:", m),
             paste("Elapsed seconds:", round(as.numeric(difftime(Sys.time(), start, units = "secs")), 1)),
             "Metrics are exploratory and conditional on successful replicates.",
             capture.output(sessionInfo())), file.path(outdir, "session-info.txt"))
utils::write.csv(source_hashes_before, file.path(outdir, "source-hashes.csv"), row.names = FALSE)
fmt <- function(x, digits = 3L) ifelse(is.na(x), "NA", formatC(x, digits = digits, format = "f"))
rows <- vapply(seq_len(nrow(summary_table)), function(i) {
  z <- summary_table[i, ]
  paste0("| ", z$scenario, " | ", z$method, " | ", z$successful, "/", z$requested,
    " | ", fmt(z$bias), " (", fmt(z$bias_mcse), ") | ", fmt(z$rmse), " (", fmt(z$rmse_mcse), ")",
    " | ", fmt(z$coverage, 2), " (", fmt(z$coverage_mcse), "); [", fmt(z$coverage_ci_lower), ", ",
    fmt(z$coverage_ci_upper), "] | ", fmt(z$mean_interval_width), " (", fmt(z$interval_width_mcse), ") |")
}, character(1))
pair_table <- do.call(rbind, paired)
trusted_pairs <- pair_table[pair_table$comparator == "trusted_pmm", ]
paired_rows <- vapply(seq_len(nrow(trusted_pairs)), function(i) {
  z <- trusted_pairs[i, ]
  paste0("| ", z$scenario, " | ", z$complete_pairs, " | ", fmt(z$mse_difference, 5),
    " | ", fmt(z$mse_difference_mcse, 5), " | ", fmt(z$coverage_difference, 2), " | ",
    fmt(z$coverage_difference_mcse), " |")
}, character(1))
lines <- c("# Exploratory benchmark results", "",
  paste0(reps, " data sets per scenario, 300 records per data set, ", m, " imputations. ",
    "The true analysis coefficient is 0.8. All results below are conditional on a successful fit; failures are retained in replicates.csv."), "",
  paste0(sum(nzchar(results$warnings)), " fits emitted recorded warnings; ", sum(results$status != "ok"), " fits failed. In missing_only, ordinary PMM logs removal of the constant all-No flag; this does not change its predictor information."), "",
  "Each parenthesized value is a Monte Carlo standard error. Coverage also includes a 95% Wilson interval for its simulation probability. A zero plug-in MCSE at 0% or 100% is not certainty.", "",
  "| Scenario | Method | Successful/requested | Bias (MCSE) | RMSE (MCSE) | Coverage (MCSE); Wilson interval | Mean interval width (MCSE) |",
  "|---|---|---:|---:|---:|---:|---:|", rows, "",
  "## Paired comparison with the strongest matched PMM baseline", "",
  "The baseline treats all flagged measurements as missing and supplies the same covariates, outcome and category to standard mice PMM. Negative squared-error differences favor catmech. Differences smaller than two MCSEs do not support a superiority claim.", "",
  "| Scenario | Pairs | Mean squared-error difference | MCSE | Coverage difference | MCSE |",
  "|---|---:|---:|---:|---:|---:|", paired_rows, "",
  "## Interpretation and boundary", "",
  "The perfect-flags scenario uses q = 1. The uncertain-flags scenario uses q = 0.5, equal to the generative corruption fraction but not a calibrated conditional posterior. The highest category has no finite upper bound: a severely corrupted value may remain category-compatible and be retained by the sensitivity mixture. Exact category consistency alone therefore cannot validate inference.", "",
  "The category-error scenario deliberately violates the exact-category assumption. Forcing incorrect categories can damage both bias and coverage. The heavy-tail scenario changes the latent residual distribution to a standardized Student t with three degrees of freedom.", "",
  "These comparisons test a specific coefficient and five synthetic settings. They neither prove proper multiple imputation nor validate confidence intervals for new data or analysis models. Inspect the trusted PMM comparison before attributing a gain over naive observed-data analysis to a new method.", "",
  "## Artifacts", "",
  "replicates.csv records every fit, warning and failure; summary.csv includes all MCSEs, failure rates and Wilson intervals; paired_comparisons.csv includes every comparator; scenarios.csv records the data-generating settings; source-hashes.csv and session-info.txt identify the code and environment. The script rejects source edits during a run.")
writeLines(lines, file.path(outdir, "RESULTS.md"))
print(summary_table[, c("scenario", "method", "successful", "bias", "rmse", "coverage", "coverage_mcse")], row.names = FALSE)
cat("Saved exploratory results to", normalizePath(outdir), "\n")
