#' Pool a prespecified linear analysis of completed datasets
#'
#' Uses Rubin's scalar rules with the complete-data residual degrees of freedom.
#' A result is conditional on the imputation/analysis models and their assumptions.
#' Failed convergence is rejected. The legacy retain/replace sensitivity engine
#' requires explicit opt-in and is labeled as sensitivity output.
#' @param result A catmech_result from [catmech_model()], [catmech_pmm()],
#'   [catmech_reference()], [catmech_contamination()] or [catmech_impute()].
#' @param formula A prespecified formula for [stats::lm()]. Missing analysis
#'   values are rejected rather than silently dropping records.
#' @param conf_level Confidence level between zero and one.
#' @param allow_sensitivity Explicitly allow pooled calculations for the legacy
#'   sensitivity engine or experimental partial-reference PMM. Results retain a
#'   sensitivity-only label; this does not establish their inferential validity.
#' @param cluster Optional name of a complete column in the completed data identifying
#'   clusters. Each completed analysis then uses a cluster-robust variance, pooled by
#'   Rubin's rules with complete-data degrees of freedom of clusters minus one. A fit
#'   from `catmech_reference(cluster = )` uses its own column by default. With few
#'   clusters, or a few influential ones, these intervals can be too narrow. A coefficient
#'   determined by a single cluster, such as a site type found at one site, draws a warning,
#'   because its interval can be far too narrow. This variance repairs the analysis step only:
#'   linear imputation without `cluster` still undercovered on real items after it.
#' @param cluster_variance With `cluster`, `"cv3"` (the default) for the delete-one-cluster
#'   jackknife, or `"cr1s"` for the form with factor (G / (G - 1)) ((n - 1) / (n - p)), the
#'   default through version 0.6.0. A planned test chose `"cv3"`, narrowly, as documented
#'   in `validation/VARIANCE_MI_SYNTHESIS.md`. It gives wider intervals than `"cr1s"`.
#' @details When every audited row agreed and recorded values were kept, the interpretation
#'   ends in `_recorded_kept`, and a warning repeats that those records are treated as correct.
#'   When the workflow declares that no target cells were imputed, all
#'   completed analysis frames must agree and ordinary complete-data t intervals
#'   are returned. A zero estimated between-imputation variance alone does not
#'   activate this rule. Predictive matching and reference-subset results retain
#'   distinct interpretation labels; none guarantees coverage under misspecification.
#'   Partial-reference PMM requires explicit sensitivity opt-in because held-out
#'   studies found material undercoverage even with random accurate references.
#'   This package policy does not assert that every PMM application is invalid.
#'   Only univariate linear models with at least one named coefficient are supported.
#'
#'   When the audited value is an analysis predictor, the analysis outcome must be among the
#'   imputation `predictors`; otherwise its coefficient is biased toward zero and a warning
#'   says so. When the audited value is the outcome, the analysis predictors belong in the
#'   imputation model instead.
#' @return A data.frame with term, estimate, std.error, df, lower, upper,
#'   within_variance, between_variance, m, n and interpretation, plus clusters
#'   when a cluster-robust variance was used.
#' @examples
#' set.seed(123)
#' n <- 40
#' df <- data.frame(
#'   lab = stats::rnorm(n, 10, 2),
#'   age = stats::rnorm(n, 50, 10),
#'   outcome = stats::rnorm(n, 5, 1)
#' )
#' df$lab_ref <- df$lab
#' df$lab_ref[11:n] <- NA
#' # lab is a predictor below, so the outcome joins the imputation model.
#' res <- catmech_reference(df, value = "lab", reference = "lab_ref",
#'                          predictors = c("age", "outcome"), m = 2, seed = 1)
#' pooled <- catmech_pool_lm(res, outcome ~ lab + age)
#' print(pooled)
#' @export
catmech_pool_lm <- function(result, formula, conf_level = 0.95,
                            allow_sensitivity = FALSE, cluster = NULL,
                            cluster_variance = c("cv3", "cr1s")) {
  cluster_variance <- match.arg(cluster_variance)
  if (!inherits(result, "catmech_result")) stop("Expected a catmech_result.", call. = FALSE)
  if (!inherits(formula, "formula")) stop("formula must be an R formula.", call. = FALSE)
  if (!is.numeric(conf_level) || length(conf_level) != 1L ||
      !is.finite(conf_level) || conf_level <= 0 || conf_level >= 1)
    stop("conf_level must be strictly between zero and one.", call. = FALSE)
  if (!is.logical(allow_sensitivity) || length(allow_sensitivity) != 1L || is.na(allow_sensitivity))
    stop("allow_sensitivity must be TRUE or FALSE.", call. = FALSE)
  engine <- result$manifest$engine
  labels <- c(gaussian_model = "conditional_model_based",
    trusted_pmm = "conditional_predictive_matching",
    reference_normal = "conditional_reference_model",
    reference_pmm = "conditional_reference_matching",
    reference_pmm_abb = "conditional_reference_matching",
    reference_normal_recorded = "conditional_reference_model",
    reference_joint_cluster = "conditional_joint_cluster_model",
    reference_marginal_cluster = "conditional_marginal_cluster_model",
    reference_contamination = "conditional_contamination_model")
  recognized <- is.character(engine) && length(engine) == 1L && !is.na(engine) && engine %in% names(labels)
  if (!is.null(engine) && !recognized)
    stop("Unrecognized workflow engine; pooled inference is refused.", call. = FALSE)
  if (!recognized && !allow_sensitivity)
    stop("This is a retain/replace sensitivity result. Set allow_sensitivity=TRUE only for explicitly labeled sensitivity calculations.", call. = FALSE)
  if (recognized && !isTRUE(result$model_diagnostics$convergence_pass))
    stop("Model convergence checks did not pass; pooled inference is refused.", call. = FALSE)
  if (length(result$completed) < 2L) stop("At least two completed datasets are required.", call. = FALSE)
  fits <- lapply(result$completed, function(d) stats::lm(formula, data = d,
                                                       na.action = stats::na.fail))
  if (any(vapply(fits, function(f) {
    coefficients <- stats::coef(f)
    !is.null(dim(coefficients)) || !length(coefficients) || is.null(names(coefficients))
  }, logical(1))))
    stop("Only univariate linear models with at least one named coefficient are supported.", call. = FALSE)
  terms <- names(stats::coef(fits[[1]]))
  ns <- vapply(fits, stats::nobs, numeric(1))
  dfs <- vapply(fits, stats::df.residual, numeric(1))
  if (length(unique(ns)) != 1L || length(unique(dfs)) != 1L || dfs[1] <= 0)
    stop("Analysis sizes and positive residual degrees of freedom must agree across completions.", call. = FALSE)
  for (fit in fits) {
    if (!identical(names(stats::coef(fit)), terms) || any(!is.finite(stats::coef(fit))) ||
        any(!is.finite(stats::vcov(fit))))
      stop("Every coefficient must be estimable in every completed analysis.", call. = FALSE)
  }
  # Missingness is workflow metadata, not inferred from coincidentally equal
  # finite-imputation estimates (which can occur with small donor pools).
  no_imputation <- FALSE
  if (recognized) {
    if (!is.null(result$manifest$imputed_rows)) {
      no_imputation <- is.numeric(result$manifest$imputed_rows) &&
        length(result$manifest$imputed_rows) == 0L
    } else if (identical(engine, "gaussian_model")) {
      requested <- result$diagnostics$requested
      no_imputation <- is.numeric(requested) && length(requested) == length(fits) &&
        all(is.finite(requested) & requested == 0)
    }
  }
  if (no_imputation && !all(vapply(fits, function(f) {
    identical(stats::model.frame(f), stats::model.frame(fits[[1]]))
  }, logical(1))))
    stop("No-imputation metadata conflicts with the completed analysis data.", call. = FALSE)
  reference_pmm_sensitivity <- identical(engine, "reference_pmm") && !no_imputation
  if (reference_pmm_sensitivity && !allow_sensitivity)
    stop("Partial-reference PMM is experimental and showed undercoverage in validation. Set allow_sensitivity=TRUE only for sensitivity calculations.", call. = FALSE)
  catmech_analysis_role_checks(result, formula, engine, "catmech_pool_lm", no_imputation = no_imputation)
  cluster_name <- if (is.null(cluster)) result[["spec"]][["cluster"]] else cluster
  robust <- NULL
  if (!is.null(cluster_name)) {
    if (!is.character(cluster_name) || length(cluster_name) != 1L || is.na(cluster_name) ||
        !cluster_name %in% names(result$completed[[1]]))
      stop("cluster must name a column in the completed data.", call. = FALSE)
    ids <- result$completed[[1]][[cluster_name]]
    if (anyNA(ids)) stop("cluster must not contain missing values.", call. = FALSE)
    if (any(!nzchar(trimws(as.character(ids))))) stop("cluster must not contain blank values.", call. = FALSE)
    clusters <- length(unique(ids))
    if (clusters < 2L) stop("A cluster-robust variance needs at least two clusters.", call. = FALSE)
    robust <- lapply(fits, if (cluster_variance == "cv3") catmech_cluster_jackknife else catmech_cluster_vcov, ids = ids)
    alone <- catmech_cluster_alone(fits[[1]], ids)
    if (length(alone))
      warning("A single cluster determines these coefficients, so their cluster-robust intervals can be far too narrow: ",
              paste(alone, collapse = ", "), ".", call. = FALSE)
  }
  recorded_kept <- identical(result[["spec"]][["recorded_fit"]], "recorded_kept_all_calibration_rows_agree")
  if (recorded_kept)
    warning("Every audited row agreed, so unaudited recorded values were kept, and these intervals treat them as correct.", call. = FALSE)
  output <- lapply(terms, function(term) {
    Q <- vapply(fits, function(f) unname(stats::coef(f)[term]), numeric(1))
    U <- if (is.null(robust)) vapply(fits, function(f) unname(stats::vcov(f)[term, term]), numeric(1)) else
      vapply(robust, function(v) unname(v[term, term]), numeric(1))
    pooled <- if (!is.null(robust)) {
      if (no_imputation) list(qbar = Q[1], t = U[1], ubar = U[1], b = 0, df = clusters - 1) else
        mice::pool.scalar(Q, U, n = clusters, k = 1)
    } else if (no_imputation) list(qbar = Q[1], t = U[1], ubar = U[1],
                                      b = 0, df = dfs[1]) else
      mice::pool.scalar(Q, U, n = ns[1], k = ns[1] - dfs[1])
    if (!is.finite(pooled$qbar) || !is.finite(pooled$t) || pooled$t <= 0 ||
        is.na(pooled$df) || pooled$df <= 0)
      stop("Pooled uncertainty is degenerate or nonfinite; a confidence interval cannot be reported.", call. = FALSE)
    if (!is.null(robust) && pooled$t <= (sqrt(.Machine$double.eps) * abs(pooled$qbar))^2)
      stop("Pooled uncertainty is degenerate or nonfinite; a confidence interval cannot be reported.", call. = FALSE)
    se <- sqrt(pooled$t)
    critical <- stats::qt((1 + conf_level)/2, pooled$df)
    row <- data.frame(term = term, estimate = pooled$qbar, std.error = se, df = pooled$df,
      lower = pooled$qbar - critical*se, upper = pooled$qbar + critical*se,
      within_variance = pooled$ubar, between_variance = pooled$b,
      m = length(fits), n = ns[1],
      interpretation = if (recognized && !reference_pmm_sensitivity) unname(labels[engine]) else "sensitivity_only",
      stringsAsFactors = FALSE)
    if (!is.null(robust)) {
      row$clusters <- clusters
      row$cluster_variance <- cluster_variance
      row$interpretation <- paste0(row$interpretation, "_cluster_robust")
    }
    if (recorded_kept) row$interpretation <- paste0(row$interpretation, "_recorded_kept")
    row
  })
  res <- do.call(rbind, output)
  structure(res, class = c("catmech_pooled_lm", "data.frame"))
}

#' Print a pooled linear model result
#'
#' @param x A `catmech_pooled_lm` object.
#' @param digits Number of decimal digits to print (default 4).
#' @param ... Unused.
#' @return `x`, invisibly.
#' @export
print.catmech_pooled_lm <- function(x, digits = 4L, ...) {
  m_val <- if (!is.null(x$m[1L])) x$m[1L] else length(x)
  n_val <- if (!is.null(x$n[1L])) x$n[1L] else NA_integer_
  cat("-- Catmech Pooled Linear Model --\n")
  cat(sprintf("Completions (m): %d | Observations (n): %d\n\n", m_val, n_val))

  tbl <- as.data.frame(x)
  num_cols <- intersect(names(tbl), c("estimate", "std.error", "df", "lower", "upper",
                                      "within_variance", "between_variance"))
  for (col in num_cols) {
    tbl[[col]] <- ifelse(is.finite(tbl[[col]]), round(tbl[[col]], digits), tbl[[col]])
  }
  disp_cols <- intersect(c("term", "estimate", "std.error", "df", "lower", "upper", "interpretation"), names(tbl))
  print(tbl[, disp_cols, drop = FALSE], row.names = FALSE)
  invisible(x)
}

#' Pool a prespecified generalized linear analysis of completed datasets
#'
#' Fits a generalized linear model ([stats::glm()]) across each completed dataset
#' and pools parameter estimates and standard errors using Rubin's rules with
#' Barnard & Rubin (1999) small-sample degrees of freedom.
#'
#' @param result A `catmech_result` from [catmech_model()], [catmech_pmm()],
#'   [catmech_reference()], [catmech_contamination()] or [catmech_impute()].
#' @param formula A prespecified formula for [stats::glm()].
#' @param family A description of the error distribution and link function to be used
#'   in the model. Defaults to `stats::binomial()` (for logistic regression).
#'   Can be a character string naming a family function or the result of a call to a family function.
#' @param conf_level Confidence level between zero and one (default 0.95).
#' @param exponentiate Logical: if `TRUE`, adds exponentiated estimates and confidence
#'   intervals (`estimate_exp`, `lower_exp`, `upper_exp`), commonly interpreted as odds
#'   ratios (for logistic regression) or incidence rate ratios (for Poisson regression).
#' @param allow_sensitivity Explicitly allow pooled calculations for sensitivity output.
#' @param cluster Optional name of a column identifying clusters for cluster-robust inference.
#'   Each completed analysis then uses the CR1S sandwich, factor (G / (G - 1)) ((n - 1) / (n - p)),
#'   built from the unscaled inverse information so that it does not depend on an estimated
#'   dispersion, and pooling uses complete-data degrees of freedom of clusters minus one.
#' @param ... Additional arguments passed to [stats::glm()] (e.g. `control = list(maxit = 50)`).
#' @details When the audited value is an analysis predictor, include the analysis outcome
#'   among the imputation `predictors`. Imputing a predictor without the outcome biases its
#'   coefficient toward zero; a warning says so. In a check recorded in
#'   `validation/REVIEW_FIXES_2026-11-07.md`, a logistic log odds ratio of 0.80 came out at
#'   0.61, with 65\% coverage of nominal 95\% intervals, when the outcome was left out.
#' @return A data frame of class `c("catmech_pooled_glm", "data.frame")` containing:
#'   \item{term}{Model coefficient names.}
#'   \item{estimate}{Pooled log-scale parameter estimate.}
#'   \item{std.error}{Total pooled standard error.}
#'   \item{df}{Barnard & Rubin adjusted degrees of freedom.}
#'   \item{lower, upper}{Confidence interval bounds on the link scale.}
#'   \item{estimate_exp, lower_exp, upper_exp}{Exponentiated parameters and bounds (if `exponentiate = TRUE`).}
#'   \item{within_variance}{Average within-imputation variance.}
#'   \item{between_variance}{Between-imputation variance.}
#'   \item{m, n}{Number of completions and total observations.}
#'   \item{family, link}{The GLM family and link function used.}
#'   \item{interpretation}{Inference label.}
#' @examples
#' set.seed(123)
#' n <- 50
#' df <- data.frame(
#'   lab = stats::rnorm(n, 10, 2),
#'   age = stats::rnorm(n, 50, 10),
#'   y = stats::rbinom(n, 1, 0.4)
#' )
#' df$lab_ref <- df$lab
#' df$lab_ref[11:n] <- NA
#' # lab is a predictor below, so the outcome y joins the imputation model.
#' res <- catmech_reference(df, value = "lab", reference = "lab_ref",
#'                          predictors = c("age", "y"), m = 2, seed = 1)
#' pooled_glm <- catmech_pool_glm(res, y ~ lab + age, family = stats::binomial(),
#'                                exponentiate = TRUE)
#' print(pooled_glm)
#' @export
catmech_pool_glm <- function(result, formula, family = stats::binomial(),
                             conf_level = 0.95, exponentiate = FALSE,
                             allow_sensitivity = FALSE, cluster = NULL, ...) {
  if (!inherits(result, "catmech_result")) stop("Expected a catmech_result.", call. = FALSE)
  if (!inherits(formula, "formula")) stop("formula must be an R formula.", call. = FALSE)
  if (is.character(family)) family <- get(family, mode = "function", envir = parent.frame())()
  if (is.function(family)) family <- family()
  if (!inherits(family, "family")) stop("family must be a valid glm family object or name.", call. = FALSE)
  if (!is.numeric(conf_level) || length(conf_level) != 1L ||
      !is.finite(conf_level) || conf_level <= 0 || conf_level >= 1)
    stop("conf_level must be strictly between zero and one.", call. = FALSE)
  if (!is.logical(exponentiate) || length(exponentiate) != 1L || is.na(exponentiate))
    stop("exponentiate must be TRUE or FALSE.", call. = FALSE)
  if (!is.logical(allow_sensitivity) || length(allow_sensitivity) != 1L || is.na(allow_sensitivity))
    stop("allow_sensitivity must be TRUE or FALSE.", call. = FALSE)

  engine <- result$manifest$engine
  labels <- c(gaussian_model = "conditional_model_based",
    trusted_pmm = "conditional_predictive_matching",
    reference_normal = "conditional_reference_model",
    reference_pmm = "conditional_reference_matching",
    reference_pmm_abb = "conditional_reference_matching",
    reference_normal_recorded = "conditional_reference_model",
    reference_joint_cluster = "conditional_joint_cluster_model",
    reference_marginal_cluster = "conditional_marginal_cluster_model",
    reference_contamination = "conditional_contamination_model")
  recognized <- is.character(engine) && length(engine) == 1L && !is.na(engine) && engine %in% names(labels)
  if (!is.null(engine) && !recognized)
    stop("Unrecognized workflow engine; pooled inference is refused.", call. = FALSE)
  if (!recognized && !allow_sensitivity)
    stop("This is a retain/replace sensitivity result. Set allow_sensitivity=TRUE only for explicitly labeled sensitivity calculations.", call. = FALSE)
  if (recognized && !isTRUE(result$model_diagnostics$convergence_pass))
    stop("Model convergence checks did not pass; pooled inference is refused.", call. = FALSE)
  if (length(result$completed) < 2L) stop("At least two completed datasets are required.", call. = FALSE)

  fits <- lapply(result$completed, function(d) stats::glm(formula, data = d, family = family,
                                                          na.action = stats::na.fail, ...))
  if (any(vapply(fits, function(f) !isTRUE(f$converged), logical(1L))))
    stop("GLM estimation failed convergence in at least one completed dataset.", call. = FALSE)

  terms <- names(stats::coef(fits[[1L]]))
  ns <- vapply(fits, stats::nobs, numeric(1L))
  dfs <- vapply(fits, stats::df.residual, numeric(1L))
  if (length(unique(ns)) != 1L || length(unique(dfs)) != 1L || dfs[1L] <= 0)
    stop("Analysis sizes and positive residual degrees of freedom must agree across completions.", call. = FALSE)

  for (fit in fits) {
    if (!identical(names(stats::coef(fit)), terms) || any(!is.finite(stats::coef(fit))) ||
        any(!is.finite(stats::vcov(fit))))
      stop("Every coefficient must be estimable in every completed analysis.", call. = FALSE)
  }

  no_imputation <- FALSE
  if (recognized) {
    if (!is.null(result$manifest$imputed_rows)) {
      no_imputation <- is.numeric(result$manifest$imputed_rows) &&
        length(result$manifest$imputed_rows) == 0L
    } else if (identical(engine, "gaussian_model")) {
      requested <- result$diagnostics$requested
      no_imputation <- is.numeric(requested) && length(requested) == length(fits) &&
        all(is.finite(requested) & requested == 0)
    }
  }

  reference_pmm_sensitivity <- identical(engine, "reference_pmm") && !no_imputation
  if (reference_pmm_sensitivity && !allow_sensitivity)
    stop("Partial-reference PMM is experimental and showed undercoverage in validation. Set allow_sensitivity=TRUE only for sensitivity calculations.", call. = FALSE)
  catmech_analysis_role_checks(result, formula, engine, "catmech_pool_glm", no_imputation = no_imputation)

  cluster_name <- if (is.null(cluster)) result[["spec"]][["cluster"]] else cluster
  robust <- NULL
  if (!is.null(cluster_name)) {
    if (!is.character(cluster_name) || length(cluster_name) != 1L || is.na(cluster_name) ||
        !cluster_name %in% names(result$completed[[1L]]))
      stop("cluster must name a column in the completed data.", call. = FALSE)
    ids <- result$completed[[1L]][[cluster_name]]
    if (anyNA(ids)) stop("cluster must not contain missing values.", call. = FALSE)
    if (any(!nzchar(trimws(as.character(ids))))) stop("cluster must not contain blank values.", call. = FALSE)
    clusters <- length(unique(ids))
    if (clusters < 2L) stop("A cluster-robust variance needs at least two clusters.", call. = FALSE)
    robust <- lapply(fits, catmech_glm_cluster_vcov, ids = ids)
  }

  recorded_kept <- identical(result[["spec"]][["recorded_fit"]], "recorded_kept_all_calibration_rows_agree")
  if (recorded_kept)
    warning("Every audited row agreed, so unaudited recorded values were kept, and these intervals treat them as correct.", call. = FALSE)

  output <- lapply(terms, function(term) {
    Q <- vapply(fits, function(f) unname(stats::coef(f)[term]), numeric(1L))
    U <- if (is.null(robust)) vapply(fits, function(f) unname(stats::vcov(f)[term, term]), numeric(1L)) else
      vapply(robust, function(v) unname(v[term, term]), numeric(1L))
    pooled <- if (!is.null(robust)) {
      if (no_imputation) list(qbar = Q[1L], t = U[1L], ubar = U[1L], b = 0, df = clusters - 1) else
        mice::pool.scalar(Q, U, n = clusters, k = 1)
    } else if (no_imputation) {
      list(qbar = Q[1L], t = U[1L], ubar = U[1L], b = 0, df = dfs[1L])
    } else {
      mice::pool.scalar(Q, U, n = ns[1L], k = ns[1L] - dfs[1L])
    }
    if (!is.finite(pooled$qbar) || !is.finite(pooled$t) || pooled$t <= 0 ||
        is.na(pooled$df) || pooled$df <= 0)
      stop("Pooled uncertainty is degenerate or nonfinite; a confidence interval cannot be reported.", call. = FALSE)
    se <- sqrt(pooled$t)
    critical <- stats::qt((1 + conf_level) / 2, pooled$df)
    lo <- pooled$qbar - critical * se
    hi <- pooled$qbar + critical * se
    row <- data.frame(
      term = term, estimate = pooled$qbar, std.error = se, df = pooled$df,
      lower = lo, upper = hi,
      within_variance = pooled$ubar, between_variance = pooled$b,
      m = length(fits), n = ns[1L],
      family = family$family, link = family$link,
      interpretation = if (recognized && !reference_pmm_sensitivity) unname(labels[engine]) else "sensitivity_only",
      stringsAsFactors = FALSE
    )
    if (exponentiate) {
      row$estimate_exp <- exp(pooled$qbar)
      row$lower_exp <- exp(lo)
      row$upper_exp <- exp(hi)
    }
    if (!is.null(robust)) {
      row$clusters <- clusters
      row$interpretation <- paste0(row$interpretation, "_cluster_robust")
    }
    if (recorded_kept) row$interpretation <- paste0(row$interpretation, "_recorded_kept")
    row
  })
  res <- do.call(rbind, output)
  structure(res, class = c("catmech_pooled_glm", "data.frame"))
}

#' Print a pooled GLM result
#'
#' @param x A `catmech_pooled_glm` object.
#' @param digits Number of decimal digits to print (default 4).
#' @param ... Unused.
#' @return `x`, invisibly.
#' @export
print.catmech_pooled_glm <- function(x, digits = 4L, ...) {
  fam <- if (!is.null(x$family[1L])) x$family[1L] else "glm"
  lnk <- if (!is.null(x$link[1L])) x$link[1L] else ""
  m_val <- if (!is.null(x$m[1L])) x$m[1L] else length(x)
  n_val <- if (!is.null(x$n[1L])) x$n[1L] else NA_integer_
  cat(sprintf("-- Catmech Pooled GLM: %s (%s link) --\n", fam, lnk))
  cat(sprintf("Completions (m): %d | Observations (n): %d\n\n", m_val, n_val))

  tbl <- as.data.frame(x)
  num_cols <- intersect(names(tbl), c("estimate", "std.error", "df", "lower", "upper",
                                      "estimate_exp", "lower_exp", "upper_exp",
                                      "within_variance", "between_variance"))
  for (col in num_cols) {
    tbl[[col]] <- ifelse(is.finite(tbl[[col]]), round(tbl[[col]], digits), tbl[[col]])
  }
  disp_cols <- c("term", "estimate", "std.error", "df", "lower", "upper")
  if ("estimate_exp" %in% names(tbl)) {
    disp_cols <- c(disp_cols, "estimate_exp", "lower_exp", "upper_exp")
  }
  print(tbl[, disp_cols, drop = FALSE], row.names = FALSE)
  invisible(x)
}

#' Pool a prespecified Cox proportional hazards analysis of completed datasets
#'
#' Fits a Cox proportional hazards model ([survival::coxph()]) across each completed dataset
#' and pools parameter estimates and standard errors using Rubin's rules with
#' Barnard & Rubin (1999) small-sample degrees of freedom based on event counts.
#'
#' @param result A `catmech_result` from [catmech_model()], [catmech_pmm()],
#'   [catmech_reference()], [catmech_contamination()] or [catmech_impute()].
#' @param formula A prespecified formula for [survival::coxph()] (e.g. `Surv(time, status) ~ lab_value + age`).
#' @param conf_level Confidence level between zero and one (default 0.95).
#' @param exponentiate Logical: if `TRUE` (default), adds exponentiated estimates and confidence
#'   intervals (`estimate_exp`, `lower_exp`, `upper_exp`), representing Hazard Ratios (HR).
#' @param allow_sensitivity Explicitly allow pooled calculations for sensitivity output.
#' @param cluster Optional name of a column identifying clusters for cluster-robust inference.
#' @param ... Additional arguments passed to [survival::coxph()] (e.g. `ties = "efron"`).
#' @details When the audited value is an analysis predictor, the imputation model needs the
#'   outcome: include the event indicator and the Nelson-Aalen cumulative hazard, which
#'   [mice::nelsonaalen()] computes, among the imputation `predictors` (White and Royston,
#'   2009, <doi:10.1002/sim.3618>). A warning is given when the event indicator is missing.
#'   The complete-data degrees of freedom are events minus coefficients, the same default
#'   that `mice::pool()` uses for Cox fits.
#' @return A data frame of class `c("catmech_pooled_coxph", "data.frame")` containing:
#'   \item{term}{Model coefficient names.}
#'   \item{estimate}{Pooled log hazard ratio parameter estimate.}
#'   \item{std.error}{Total pooled standard error.}
#'   \item{df}{Barnard & Rubin adjusted degrees of freedom.}
#'   \item{lower, upper}{Confidence interval bounds on the log hazard scale.}
#'   \item{estimate_exp, lower_exp, upper_exp}{Exponentiated parameters and bounds (Hazard Ratios) if `exponentiate = TRUE`.}
#'   \item{within_variance}{Average within-imputation variance.}
#'   \item{between_variance}{Between-imputation variance.}
#'   \item{m, n, nevent}{Number of completions, total observations, and events.}
#'   \item{interpretation}{Inference label.}
#' @examples
#' set.seed(123)
#' n <- 50
#' df <- data.frame(
#'   time = stats::rexp(n, 0.1),
#'   status = stats::rbinom(n, 1, 0.5),
#'   lab = stats::rnorm(n, 10, 2),
#'   age = stats::rnorm(n, 50, 10)
#' )
#' df$lab_ref <- df$lab
#' df$lab_ref[11:n] <- NA
#' # lab is a predictor below, so the event indicator and the Nelson-Aalen
#' # cumulative hazard join the imputation model.
#' df$hazard <- mice::nelsonaalen(df, time, status)
#' res <- catmech_reference(df, value = "lab", reference = "lab_ref",
#'                          predictors = c("age", "status", "hazard"), m = 2, seed = 1)
#' pooled_cox <- catmech_pool_coxph(res, survival::Surv(time, status) ~ lab + age,
#'                                  exponentiate = TRUE)
#' print(pooled_cox)
#' @export
catmech_pool_coxph <- function(result, formula, conf_level = 0.95,
                               exponentiate = TRUE, allow_sensitivity = FALSE,
                               cluster = NULL, ...) {
  if (!requireNamespace("survival", quietly = TRUE))
    stop("Package 'survival' is required for catmech_pool_coxph().", call. = FALSE)
  if (!inherits(result, "catmech_result")) stop("Expected a catmech_result.", call. = FALSE)
  if (!inherits(formula, "formula")) stop("formula must be an R formula.", call. = FALSE)
  if (!is.numeric(conf_level) || length(conf_level) != 1L ||
      !is.finite(conf_level) || conf_level <= 0 || conf_level >= 1)
    stop("conf_level must be strictly between zero and one.", call. = FALSE)
  if (!is.logical(exponentiate) || length(exponentiate) != 1L || is.na(exponentiate))
    stop("exponentiate must be TRUE or FALSE.", call. = FALSE)
  if (!is.logical(allow_sensitivity) || length(allow_sensitivity) != 1L || is.na(allow_sensitivity))
    stop("allow_sensitivity must be TRUE or FALSE.", call. = FALSE)

  engine <- result$manifest$engine
  labels <- c(gaussian_model = "conditional_model_based",
    trusted_pmm = "conditional_predictive_matching",
    reference_normal = "conditional_reference_model",
    reference_pmm = "conditional_reference_matching",
    reference_pmm_abb = "conditional_reference_matching",
    reference_normal_recorded = "conditional_reference_model",
    reference_joint_cluster = "conditional_joint_cluster_model",
    reference_marginal_cluster = "conditional_marginal_cluster_model",
    reference_contamination = "conditional_contamination_model")
  recognized <- is.character(engine) && length(engine) == 1L && !is.na(engine) && engine %in% names(labels)
  if (!is.null(engine) && !recognized)
    stop("Unrecognized workflow engine; pooled inference is refused.", call. = FALSE)
  if (!recognized && !allow_sensitivity)
    stop("This is a retain/replace sensitivity result. Set allow_sensitivity=TRUE only for explicitly labeled sensitivity calculations.", call. = FALSE)
  if (recognized && !isTRUE(result$model_diagnostics$convergence_pass))
    stop("Model convergence checks did not pass; pooled inference is refused.", call. = FALSE)
  if (length(result$completed) < 2L) stop("At least two completed datasets are required.", call. = FALSE)

  cluster_name <- if (is.null(cluster)) result[["spec"]][["cluster"]] else cluster
  if (!is.null(cluster_name)) {
    if (!is.character(cluster_name) || length(cluster_name) != 1L || is.na(cluster_name) ||
        !cluster_name %in% names(result$completed[[1L]]))
      stop("cluster must name a column in the completed data.", call. = FALSE)
    ids <- result$completed[[1L]][[cluster_name]]
    if (anyNA(ids)) stop("cluster must not contain missing values.", call. = FALSE)
    if (any(!nzchar(trimws(as.character(ids))))) stop("cluster must not contain blank values.", call. = FALSE)
    clusters <- length(unique(ids))
    if (clusters < 2L) stop("A cluster-robust variance needs at least two clusters.", call. = FALSE)
  }

  fits <- lapply(result$completed, function(d) {
    if (!is.null(cluster_name)) {
      cl_call <- as.call(c(list(quote(survival::coxph), formula = formula, data = quote(d),
                                cluster = as.name(cluster_name), na.action = quote(stats::na.fail)),
                           list(...)))
      eval(cl_call)
    } else {
      cl_call <- as.call(c(list(quote(survival::coxph), formula = formula, data = quote(d),
                                na.action = quote(stats::na.fail)),
                           list(...)))
      eval(cl_call)
    }
  })

  terms <- names(stats::coef(fits[[1L]]))
  ns <- vapply(fits, function(f) f$n, numeric(1L))
  nevents <- vapply(fits, function(f) f$nevent, numeric(1L))
  if (length(unique(ns)) != 1L || length(unique(nevents)) != 1L)
    stop("Analysis sample size and event counts must agree across completions.", call. = FALSE)

  p <- length(terms)
  df_com <- nevents[1L] - p
  if (df_com <= 0)
    stop("Cox model requires more events than parameters.", call. = FALSE)

  for (fit in fits) {
    if (!identical(names(stats::coef(fit)), terms) || any(!is.finite(stats::coef(fit))) ||
        any(!is.finite(stats::vcov(fit))))
      stop("Every coefficient must be estimable in every completed analysis.", call. = FALSE)
  }

  no_imputation <- FALSE
  if (recognized) {
    if (!is.null(result$manifest$imputed_rows)) {
      no_imputation <- is.numeric(result$manifest$imputed_rows) &&
        length(result$manifest$imputed_rows) == 0L
    } else if (identical(engine, "gaussian_model")) {
      requested <- result$diagnostics$requested
      no_imputation <- is.numeric(requested) && length(requested) == length(fits) &&
        all(is.finite(requested) & requested == 0)
    }
  }

  reference_pmm_sensitivity <- identical(engine, "reference_pmm") && !no_imputation
  if (reference_pmm_sensitivity && !allow_sensitivity)
    stop("Partial-reference PMM is experimental and showed undercoverage in validation. Set allow_sensitivity=TRUE only for sensitivity calculations.", call. = FALSE)
  catmech_analysis_role_checks(result, formula, engine, "catmech_pool_coxph", survival = TRUE,
                               no_imputation = no_imputation)

  recorded_kept <- identical(result[["spec"]][["recorded_fit"]], "recorded_kept_all_calibration_rows_agree")
  if (recorded_kept)
    warning("Every audited row agreed, so unaudited recorded values were kept, and these intervals treat them as correct.", call. = FALSE)

  output <- lapply(terms, function(term) {
    Q <- vapply(fits, function(f) unname(stats::coef(f)[term]), numeric(1L))
    U <- vapply(fits, function(f) unname(stats::vcov(f)[term, term]), numeric(1L))
    pooled <- if (!is.null(cluster_name)) {
      if (no_imputation) list(qbar = Q[1L], t = U[1L], ubar = U[1L], b = 0, df = clusters - 1) else
        mice::pool.scalar(Q, U, n = clusters, k = 1)
    } else if (no_imputation) {
      list(qbar = Q[1L], t = U[1L], ubar = U[1L], b = 0, df = df_com)
    } else {
      mice::pool.scalar(Q, U, n = nevents[1L], k = p)
    }
    if (!is.finite(pooled$qbar) || !is.finite(pooled$t) || pooled$t <= 0 ||
        is.na(pooled$df) || pooled$df <= 0)
      stop("Pooled uncertainty is degenerate or nonfinite; a confidence interval cannot be reported.", call. = FALSE)
    se <- sqrt(pooled$t)
    critical <- stats::qt((1 + conf_level) / 2, pooled$df)
    lo <- pooled$qbar - critical * se
    hi <- pooled$qbar + critical * se
    row <- data.frame(
      term = term, estimate = pooled$qbar, std.error = se, df = pooled$df,
      lower = lo, upper = hi,
      within_variance = pooled$ubar, between_variance = pooled$b,
      m = length(fits), n = ns[1L], nevent = nevents[1L],
      interpretation = if (recognized && !reference_pmm_sensitivity) unname(labels[engine]) else "sensitivity_only",
      stringsAsFactors = FALSE
    )
    if (exponentiate) {
      row$estimate_exp <- exp(pooled$qbar)
      row$lower_exp <- exp(lo)
      row$upper_exp <- exp(hi)
    }
    if (!is.null(cluster_name)) {
      row$clusters <- clusters
      row$interpretation <- paste0(row$interpretation, "_cluster_robust")
    }
    if (recorded_kept) row$interpretation <- paste0(row$interpretation, "_recorded_kept")
    row
  })
  res <- do.call(rbind, output)
  structure(res, class = c("catmech_pooled_coxph", "data.frame"))
}

#' Print a pooled Cox PH result
#'
#' @param x A `catmech_pooled_coxph` object.
#' @param digits Number of decimal digits to print (default 4).
#' @param ... Unused.
#' @return `x`, invisibly.
#' @export
print.catmech_pooled_coxph <- function(x, digits = 4L, ...) {
  m_val <- if (!is.null(x$m[1L])) x$m[1L] else length(x)
  n_val <- if (!is.null(x$n[1L])) x$n[1L] else NA_integer_
  ev_val <- if (!is.null(x$nevent[1L])) x$nevent[1L] else NA_integer_
  cat(sprintf("-- Catmech Pooled Cox PH Model --\n"))
  cat(sprintf("Completions (m): %d | Observations (n): %d | Events: %d\n\n", m_val, n_val, ev_val))

  tbl <- as.data.frame(x)
  num_cols <- intersect(names(tbl), c("estimate", "std.error", "df", "lower", "upper",
                                      "estimate_exp", "lower_exp", "upper_exp",
                                      "within_variance", "between_variance"))
  for (col in num_cols) {
    tbl[[col]] <- ifelse(is.finite(tbl[[col]]), round(tbl[[col]], digits), tbl[[col]])
  }
  disp_cols <- c("term", "estimate", "std.error", "df", "lower", "upper")
  if ("estimate_exp" %in% names(tbl)) {
    disp_cols <- c(disp_cols, "estimate_exp", "lower_exp", "upper_exp")
  }
  print(tbl[, disp_cols, drop = FALSE], row.names = FALSE)
  invisible(x)
}

plot_pooled_forest <- function(x, terms = NULL, intercept = FALSE, exponentiate = NULL,
                               col_point = "#1f77b4", col_ci = "#1f77b4", col_ref = "darkred",
                               lwd = 2, pch = 19, cex = 1.2, main = NULL, xlab = NULL, ...) {
  tbl <- as.data.frame(x)
  if (!intercept && nrow(tbl) > 1L) {
    tbl <- tbl[!tbl$term %in% c("(Intercept)", "Intercept"), , drop = FALSE]
  }
  if (!is.null(terms)) {
    tbl <- tbl[tbl$term %in% terms, , drop = FALSE]
  }
  if (nrow(tbl) == 0L) {
    message("No terms to plot.")
    return(invisible(x))
  }

  has_exp <- "estimate_exp" %in% names(tbl) && all(is.finite(tbl$estimate_exp))
  use_exp <- if (is.null(exponentiate)) has_exp else isTRUE(exponentiate)

  if (use_exp && !has_exp) {
    tbl$estimate_exp <- exp(tbl$estimate)
    tbl$lower_exp <- exp(tbl$lower)
    tbl$upper_exp <- exp(tbl$upper)
  }

  vals <- if (use_exp) tbl$estimate_exp else tbl$estimate
  lows <- if (use_exp) tbl$lower_exp else tbl$lower
  highs <- if (use_exp) tbl$upper_exp else tbl$upper
  ref_val <- if (use_exp) 1 else 0

  k <- nrow(tbl)
  y_pos <- rev(seq_len(k))

  xlim <- range(c(lows, highs, ref_val), na.rm = TRUE)
  x_padding <- max(diff(xlim) * 0.15, 0.05)
  xlim <- c(xlim[1L] - x_padding, xlim[2L] + x_padding)

  if (is.null(xlab)) {
    xlab <- if (use_exp) {
      if (inherits(x, "catmech_pooled_coxph")) "Hazard Ratio (95% CI)" else "Odds Ratio (95% CI)"
    } else {
      "Estimate (95% CI)"
    }
  }

  if (is.null(main)) {
    fam <- if (!is.null(tbl$family[1L])) paste0(" (", tbl$family[1L], ")") else ""
    main <- paste0("Catmech Pooled Forest Plot", fam)
  }

  old_mar <- graphics::par("mar")
  on.exit(graphics::par(mar = old_mar))
  max_len <- max(nchar(as.character(tbl$term)))
  graphics::par(mar = c(5, max(4.5, min(10, max_len * 0.5)), 4, 6))

  graphics::plot(
    vals, y_pos,
    type = "n",
    xlim = xlim, ylim = c(0.5, k + 0.5),
    xlab = xlab, ylab = "",
    yaxt = "n", main = main,
    panel.first = graphics::grid(col = "gray90", lty = "dotted"),
    ...
  )

  graphics::abline(v = ref_val, lty = 2, col = col_ref, lwd = 1.5)

  for (i in seq_len(k)) {
    graphics::segments(lows[i], y_pos[i], highs[i], y_pos[i], col = col_ci, lwd = lwd)
    tick_h <- 0.12
    graphics::segments(lows[i], y_pos[i] - tick_h, lows[i], y_pos[i] + tick_h, col = col_ci, lwd = lwd)
    graphics::segments(highs[i], y_pos[i] - tick_h, highs[i], y_pos[i] + tick_h, col = col_ci, lwd = lwd)
  }

  graphics::points(vals, y_pos, pch = pch, col = col_point, cex = cex)
  graphics::axis(2, at = y_pos, labels = tbl$term, las = 1, cex.axis = 0.9)

  num_labels <- sprintf("%.2f [%.2f, %.2f]", vals, lows, highs)
  graphics::mtext(num_labels, side = 4, at = y_pos, las = 1, line = 0.5, cex = 0.8, col = "gray30")

  invisible(x)
}

#' Plot pooled GLM coefficients (Forest Plot)
#'
#' @param x A `catmech_pooled_glm` object.
#' @param terms Character vector of terms to plot (default `NULL` for all).
#' @param intercept Logical: include intercept (default `FALSE`).
#' @param exponentiate Logical: exponentiate estimates (default `TRUE` if available).
#' @param col_point Point color (default `"#1f77b4"`).
#' @param col_ci CI color (default `"#1f77b4"`).
#' @param col_ref Null reference line color (default `"darkred"`).
#' @param lwd Line width (default 2).
#' @param pch Point symbol (default 19).
#' @param cex Point size (default 1.2).
#' @param main Plot title.
#' @param xlab X-axis label.
#' @param ... Additional arguments passed to [graphics::plot()].
#' @return `x`, invisibly.
#' @examples
#' set.seed(123)
#' n <- 40
#' df <- data.frame(
#'   lab = stats::rnorm(n, 10, 2),
#'   age = stats::rnorm(n, 50, 10),
#'   y = stats::rbinom(n, 1, 0.4)
#' )
#' df$lab_ref <- df$lab
#' df$lab_ref[11:n] <- NA
#' # lab is a predictor below, so the outcome y joins the imputation model.
#' res <- catmech_reference(df, value = "lab", reference = "lab_ref",
#'                          predictors = c("age", "y"), m = 2, seed = 1)
#' pooled <- catmech_pool_glm(res, y ~ lab + age, family = stats::binomial())
#' plot(pooled)
#' @export
plot.catmech_pooled_glm <- function(x, terms = NULL, intercept = FALSE, exponentiate = NULL,
                                   col_point = "#1f77b4", col_ci = "#1f77b4", col_ref = "darkred",
                                   lwd = 2, pch = 19, cex = 1.2, main = NULL, xlab = NULL, ...) {
  plot_pooled_forest(x, terms = terms, intercept = intercept, exponentiate = exponentiate,
                     col_point = col_point, col_ci = col_ci, col_ref = col_ref,
                     lwd = lwd, pch = pch, cex = cex, main = main, xlab = xlab, ...)
}

#' Plot pooled Cox PH hazard ratios (Forest Plot)
#'
#' @param x A `catmech_pooled_coxph` object.
#' @param terms Character vector of terms to plot (default `NULL` for all).
#' @param intercept Logical: include intercept (default `FALSE`).
#' @param exponentiate Logical: exponentiate estimates (default `TRUE` if available).
#' @param col_point Point color (default `"#1f77b4"`).
#' @param col_ci CI color (default `"#1f77b4"`).
#' @param col_ref Null reference line color (default `"darkred"`).
#' @param lwd Line width (default 2).
#' @param pch Point symbol (default 19).
#' @param cex Point size (default 1.2).
#' @param main Plot title.
#' @param xlab X-axis label.
#' @param ... Additional arguments passed to [graphics::plot()].
#' @return `x`, invisibly.
#' @examples
#' set.seed(123)
#' n <- 40
#' df <- data.frame(
#'   time = stats::rexp(n, 0.1), status = stats::rbinom(n, 1, 0.5),
#'   lab = stats::rnorm(n, 10, 2), age = stats::rnorm(n, 50, 10)
#' )
#' df$lab_ref <- df$lab
#' df$lab_ref[11:n] <- NA
#' df$hazard <- mice::nelsonaalen(df, time, status)
#' res <- catmech_reference(df, value = "lab", reference = "lab_ref",
#'                          predictors = c("age", "status", "hazard"), m = 2, seed = 1)
#' pooled <- catmech_pool_coxph(res, survival::Surv(time, status) ~ lab + age)
#' plot(pooled)
#' @export
plot.catmech_pooled_coxph <- function(x, terms = NULL, intercept = FALSE, exponentiate = NULL,
                                     col_point = "#1f77b4", col_ci = "#1f77b4", col_ref = "darkred",
                                     lwd = 2, pch = 19, cex = 1.2, main = NULL, xlab = NULL, ...) {
  plot_pooled_forest(x, terms = terms, intercept = intercept, exponentiate = exponentiate,
                     col_point = col_point, col_ci = col_ci, col_ref = col_ref,
                     lwd = lwd, pch = pch, cex = cex, main = main, xlab = xlab, ...)
}

#' Plot pooled linear model coefficients (Forest Plot)
#'
#' @param x A `catmech_pooled_lm` object.
#' @param terms Character vector of terms to plot (default `NULL` for all).
#' @param intercept Logical: include intercept (default `FALSE`).
#' @param exponentiate Logical: exponentiate estimates (default `FALSE`).
#' @param col_point Point color (default `"#1f77b4"`).
#' @param col_ci CI color (default `"#1f77b4"`).
#' @param col_ref Null reference line color (default `"darkred"`).
#' @param lwd Line width (default 2).
#' @param pch Point symbol (default 19).
#' @param cex Point size (default 1.2).
#' @param main Plot title.
#' @param xlab X-axis label.
#' @param ... Additional arguments passed to [graphics::plot()].
#' @return `x`, invisibly.
#' @examples
#' set.seed(123)
#' n <- 40
#' df <- data.frame(
#'   lab = stats::rnorm(n, 10, 2),
#'   age = stats::rnorm(n, 50, 10),
#'   y = stats::rnorm(n, 5, 1)
#' )
#' df$lab_ref <- df$lab
#' df$lab_ref[11:n] <- NA
#' res <- catmech_reference(df, value = "lab", reference = "lab_ref",
#'                          predictors = c("age", "y"), m = 2, seed = 1)
#' pooled <- catmech_pool_lm(res, y ~ lab + age)
#' plot(pooled)
#' @export
plot.catmech_pooled_lm <- function(x, terms = NULL, intercept = FALSE, exponentiate = NULL,
                                  col_point = "#1f77b4", col_ci = "#1f77b4", col_ref = "darkred",
                                  lwd = 2, pch = 19, cex = 1.2, main = NULL, xlab = NULL, ...) {
  plot_pooled_forest(x, terms = terms, intercept = intercept, exponentiate = exponentiate,
                     col_point = col_point, col_ci = col_ci, col_ref = col_ref,
                     lwd = lwd, pch = pch, cex = cex, main = main, xlab = xlab, ...)
}

#' Turn a pooled catmech model into a tidy data frame
#'
#' @param x An object to be tidied.
#' @param ... Additional arguments.
#' @return A data.frame with standardized model columns.
#' @examples
#' set.seed(123)
#' n <- 40
#' df <- data.frame(
#'   lab = stats::rnorm(n, 10, 2),
#'   age = stats::rnorm(n, 50, 10),
#'   y = stats::rnorm(n, 5, 1)
#' )
#' df$lab_ref <- df$lab
#' df$lab_ref[11:n] <- NA
#' res <- catmech_reference(df, value = "lab", reference = "lab_ref",
#'                          predictors = c("age", "y"), m = 2, seed = 1)
#' pooled <- catmech_pool_lm(res, y ~ lab + age)
#' tidy(pooled)
#' @export
tidy <- function(x, ...) UseMethod("tidy")

#' Summarize a pooled catmech model with glance statistics
#'
#' @param x An object to be summarized.
#' @param ... Additional arguments.
#' @return A one-row data.frame.
#' @examples
#' set.seed(123)
#' n <- 40
#' df <- data.frame(
#'   lab = stats::rnorm(n, 10, 2),
#'   age = stats::rnorm(n, 50, 10),
#'   y = stats::rnorm(n, 5, 1)
#' )
#' df$lab_ref <- df$lab
#' df$lab_ref[11:n] <- NA
#' res <- catmech_reference(df, value = "lab", reference = "lab_ref",
#'                          predictors = c("age", "y"), m = 2, seed = 1)
#' pooled <- catmech_pool_lm(res, y ~ lab + age)
#' glance(pooled)
#' @export
glance <- function(x, ...) UseMethod("glance")

tidy_pooled_model <- function(x, conf.int = TRUE, conf.level = 0.95, exponentiate = FALSE, ...) {
  tbl <- as.data.frame(x)
  statistic <- tbl$estimate / tbl$std.error
  p_val <- if ("df" %in% names(tbl) && any(is.finite(tbl$df))) {
    2 * stats::pt(-abs(statistic), df = tbl$df)
  } else {
    2 * (1 - stats::pnorm(abs(statistic)))
  }
  df_val <- if ("df" %in% names(tbl)) tbl$df else NA_real_

  res <- data.frame(
    term = tbl$term,
    estimate = tbl$estimate,
    std.error = tbl$std.error,
    statistic = statistic,
    df = df_val,
    p.value = p_val,
    stringsAsFactors = FALSE
  )

  if (isTRUE(conf.int)) {
    if (isTRUE(conf.level == 0.95) && "lower" %in% names(tbl) && "upper" %in% names(tbl)) {
      res$conf.low <- tbl$lower
      res$conf.high <- tbl$upper
    } else {
      crit <- if (any(is.finite(df_val))) stats::qt((1 + conf.level) / 2, df_val) else stats::qnorm((1 + conf.level) / 2)
      res$conf.low <- tbl$estimate - crit * tbl$std.error
      res$conf.high <- tbl$estimate + crit * tbl$std.error
    }
  }

  if (isTRUE(exponentiate)) {
    res$estimate <- exp(res$estimate)
    if (isTRUE(conf.int)) {
      res$conf.low <- exp(res$conf.low)
      res$conf.high <- exp(res$conf.high)
    }
    res$std.error <- res$estimate * tbl$std.error
  }

  res
}

#' @export
tidy.catmech_pooled_glm <- function(x, conf.int = TRUE, conf.level = 0.95, exponentiate = FALSE, ...) {
  tidy_pooled_model(x, conf.int = conf.int, conf.level = conf.level, exponentiate = exponentiate, ...)
}

#' @export
tidy.catmech_pooled_coxph <- function(x, conf.int = TRUE, conf.level = 0.95, exponentiate = FALSE, ...) {
  tidy_pooled_model(x, conf.int = conf.int, conf.level = conf.level, exponentiate = exponentiate, ...)
}

#' @export
tidy.catmech_pooled_lm <- function(x, conf.int = TRUE, conf.level = 0.95, exponentiate = FALSE, ...) {
  tidy_pooled_model(x, conf.int = conf.int, conf.level = conf.level, exponentiate = exponentiate, ...)
}

#' @export
glance.catmech_pooled_glm <- function(x, ...) {
  data.frame(
    nobs = if (!is.null(x$n[1L])) x$n[1L] else NA_integer_,
    m = if (!is.null(x$m[1L])) x$m[1L] else NA_integer_,
    family = if (!is.null(x$family[1L])) x$family[1L] else NA_character_,
    link = if (!is.null(x$link[1L])) x$link[1L] else NA_character_,
    stringsAsFactors = FALSE
  )
}

#' @export
glance.catmech_pooled_coxph <- function(x, ...) {
  data.frame(
    nobs = if (!is.null(x$n[1L])) x$n[1L] else NA_integer_,
    n_events = if (!is.null(x$nevent[1L])) x$nevent[1L] else NA_integer_,
    m = if (!is.null(x$m[1L])) x$m[1L] else NA_integer_,
    stringsAsFactors = FALSE
  )
}

#' @export
glance.catmech_pooled_lm <- function(x, ...) {
  data.frame(
    nobs = if (!is.null(x$n[1L])) x$n[1L] else NA_integer_,
    m = if (!is.null(x$m[1L])) x$m[1L] else NA_integer_,
    stringsAsFactors = FALSE
  )
}

# Checks on the role of the audited value in a pooled analysis.
# 1. When the audited value is an analysis predictor, the analysis outcome belongs in the
#    imputation model. Imputing a predictor without the outcome biases its coefficient
#    toward zero: in validation/REVIEW_FIXES_2026-11-07.md, a logistic slope of 0.80 came
#    out at 0.61 with 65% coverage of nominal 95% intervals. For a Cox model the event
#    indicator and the Nelson-Aalen cumulative hazard play the outcome's part (White and
#    Royston 2009).
# 2. The clustered imputation options were not validated with the audited value as a predictor.
catmech_analysis_role_checks <- function(result, formula, engine, fun, survival = FALSE,
                                         no_imputation = FALSE) {
  spec <- result[["spec"]]
  value <- if (is.list(spec)) spec[["value"]] else NULL
  if (!is.character(value) || length(value) != 1L || is.na(value) || length(formula) != 3L)
    return(invisible(NULL))
  rhs <- tryCatch(all.vars(stats::delete.response(stats::terms(formula, data = result$completed[[1L]]))),
                  error = function(e) character())
  if (!value %in% rhs) return(invisible(NULL))
  # Engines keep their imputation predictors in spec, except catmech_model() and the
  # legacy engine, which keep them in the manifest.
  used <- if ("predictors" %in% names(spec)) spec[["predictors"]] else result[["manifest"]][["predictors"]]
  known <- "predictors" %in% names(spec) || "predictors" %in% names(result[["manifest"]])
  if (known && !no_imputation) {
    lhs <- formula[[2L]]
    outcome <- if (survival) catmech_surv_event_vars(lhs) else all.vars(lhs)
    missing <- setdiff(outcome, c(used, value))
    if (length(missing)) {
      if (survival)
        warning(sprintf("The audited value '%s' is a predictor in this analysis, but the event indicator '%s' was not an imputation predictor. Imputing a predictor without the outcome biases its hazard ratio toward 1. Refit with the event indicator and the Nelson-Aalen cumulative hazard (mice::nelsonaalen()) in `predictors`.",
          value, paste(missing, collapse = "', '")), call. = FALSE)
      else
        warning(sprintf("The audited value '%s' is a predictor in this analysis, but the outcome %s %s not an imputation predictor. Imputing a predictor without the outcome biases its coefficient toward zero. Refit with the outcome in `predictors`.",
          value, paste0("'", missing, "'", collapse = ", "), if (length(missing) == 1L) "was" else "were"),
          call. = FALSE)
    }
  }
  if (!is.null(spec[["cluster"]])) {
    if (identical(engine, "reference_marginal_cluster"))
      warning(if (identical(fun, "catmech_pool_lm"))
        "Joint clustered imputation has not been validated here with the audited value as an analysis predictor in catmech_pool_lm(); for direct population marginal regression when the outcome was included in predictors, use catmech_marginal_summary(result, formula = ...). See docs/JOINT_MARGINAL_METHOD.md." else
        sprintf("Joint clustered imputation has not been validated here with the audited value as an analysis predictor in %s(). See docs/JOINT_MARGINAL_METHOD.md.", fun),
        call. = FALSE)
    else if (identical(engine, "reference_joint_cluster"))
      warning("Joint clustered imputation has not been validated here with the audited value as an analysis predictor. The earlier moment model showed slope attenuation in that role; joint sampling alone does not establish compatible inference. See docs/JOINT_CLUSTER_METHOD.md.", call. = FALSE)
    else
      warning(sprintf("catmech_reference(cluster = ) gives an unaudited cluster a new intercept, which biased slopes toward zero in simulations when the audited value was a predictor. Prefer a fit without cluster, pooled with %s(cluster = ). See validation/VARIANCE_MI_SYNTHESIS.md.", fun), call. = FALSE)
  }
  invisible(NULL)
}

# Variables of the event indicator in a Surv() response: the event argument of
# Surv(time, event), Surv(start, stop, event) or Surv(time, event == k), matched by name or
# position. Other responses, such as a precomputed Surv column, cannot be checked.
catmech_surv_event_vars <- function(lhs) {
  if (!is.call(lhs)) return(character())
  fun <- lhs[[1L]]
  is_surv <- identical(fun, as.name("Surv")) ||
    (is.call(fun) && length(fun) == 3L && identical(fun[[3L]], as.name("Surv")))
  if (!is_surv) return(character())
  matched <- tryCatch(match.call(survival::Surv, lhs), error = function(e) NULL)
  if (is.null(matched)) return(character())
  event <- if (!is.null(matched$event)) matched$event else matched$time2
  if (is.null(event)) character() else all.vars(event)
}

catmech_glm_cluster_vcov <- function(fit, ids) {
  x <- stats::model.matrix(fit)
  # Aliased columns are dropped, matching summary(fit)$cov.unscaled below.
  x <- x[, !is.na(stats::coef(fit)), drop = FALSE]
  n <- nrow(x)
  k <- ncol(x)
  w <- fit$weights
  res_working <- stats::residuals(fit, type = "working")
  # Working weights times working residuals are the score contributions without the
  # dispersion, so the bread is the unscaled inverse information. vcov(fit) would carry
  # the estimated dispersion twice into the sandwich for gaussian, Gamma and quasi
  # families; binomial and poisson, with dispersion one, are unchanged.
  scores <- x * (w * res_working)
  sums <- rowsum(scores, as.character(ids), reorder = FALSE)
  g <- nrow(sums)
  bread <- summary(fit)$cov.unscaled
  meat <- crossprod(sums)
  out <- bread %*% meat %*% bread * (g / (g - 1)) * ((n - 1) / (n - k))
  dimnames(out) <- list(colnames(x), colnames(x))
  out
}

# Design and (X'X)^-1 of a linear fit over its estimable coefficients. Aliased columns
# (NA coefficients) are dropped: the pivoted QR of a rank-deficient fit would otherwise
# give a bread in the wrong order with a near-zero pivot. Full-rank fits are unchanged.
catmech_lm_bread <- function(fit) {
  x <- stats::model.matrix(fit)
  keep <- !is.na(stats::coef(fit))
  if (all(keep)) return(list(x = x, bread = chol2inv(qr.R(fit$qr))))
  x <- x[, keep, drop = FALSE]
  list(x = x, bread = chol2inv(qr.R(qr(x))))
}

# Delete-one-cluster jackknife covariance of one fitted linear model (CV3), by the closed
# form for each deletion. A singular block, from a cluster that alone determines a
# coefficient, uses a pseudoinverse.
catmech_cluster_jackknife <- function(fit, ids) {
  design <- catmech_lm_bread(fit)
  x <- design$x
  e <- stats::residuals(fit)
  bread <- design$bread
  units <- split(seq_len(nrow(x)), as.character(ids))
  out <- matrix(0, ncol(x), ncol(x), dimnames = list(colnames(x), colnames(x)))
  for (u in units) {
    xg <- x[u, , drop = FALSE]
    block <- eigen(diag(length(u)) - xg %*% bread %*% t(xg), symmetric = TRUE)
    keep <- block$values > 1e-10
    inverse <- block$vectors[, keep, drop = FALSE] %*% (t(block$vectors[, keep, drop = FALSE]) / block$values[keep])
    out <- out + tcrossprod(bread %*% crossprod(xg, inverse %*% e[u]))
  }
  g <- length(units)
  out * (g - 1) / g
}

# Coefficients that a single cluster determines: deleting that cluster leaves them inestimable.
# A zero eigenvalue v of the block gives the null direction (X'X)^-1 X_g' v of the design without it.
catmech_cluster_alone <- function(fit, ids) {
  design <- catmech_lm_bread(fit)
  x <- design$x
  bread <- design$bread
  flagged <- logical(ncol(x))
  for (u in split(seq_len(nrow(x)), as.character(ids))) {
    xg <- x[u, , drop = FALSE]
    block <- eigen(diag(length(u)) - xg %*% bread %*% t(xg), symmetric = TRUE)
    zero <- block$values <= 1e-10
    if (!any(zero)) next
    direction <- bread %*% crossprod(xg, block$vectors[, zero, drop = FALSE])
    flagged <- flagged | apply(abs(direction), 1L, max) > 1e-8 * max(abs(direction))
  }
  colnames(x)[flagged]
}

# Cluster-robust covariance of one fitted linear model with factor (G / (G - 1)) ((n - 1) / (n - p)).
catmech_cluster_vcov <- function(fit, ids) {
  design <- catmech_lm_bread(fit)
  x <- design$x
  n <- nrow(x); k <- ncol(x)
  sums <- rowsum(x * stats::residuals(fit), as.character(ids), reorder = FALSE)
  g <- nrow(sums)
  bread <- design$bread
  out <- bread %*% crossprod(sums) %*% bread * (g / (g - 1)) * ((n - 1) / (n - k))
  dimnames(out) <- list(colnames(x), colnames(x))
  out
}

#' Inspect a finite set of prespecified imputation assumptions
#'
#' Inspired by assumption-path audits: keep failed and unresolved scenarios
#' visible, and do not select the most favorable result. The interval envelope
#' is descriptive; it is not a confidence interval over arbitrary misspecification.
#' It covers only the supplied finite set and inherits each member's limitations.
#' @param results A uniquely named list of catmech_result objects or error
#'   conditions captured with tryCatch(..., error = identity). Supply every
#'   prespecified scenario, including failures.
#' @param formula A prespecified linear analysis formula, identical across scenarios.
#'   When `method = "marginal_summary"`, `formula` is optional and defaults to the
#'   population regression of the true measurement on the continuous predictors;
#'   it can also specify a custom response and covariate subset (e.g. `outcome ~ recorded + x`).
#'   When `method = "joint_summary"`, `formula` is not used.
#' @param term Name of the coefficient to compare.
#' @param reference Reference value, usually zero, for interval-direction reporting.
#' @param conf_level Conditional confidence level passed to [catmech_pool_lm()] or
#'   [catmech_marginal_summary()] / [catmech_joint_summary()].
#' @param expected_scenarios Optional character vector naming the complete
#'   prespecified set. Omitted members are recorded as failures; unplanned extras
#'   are rejected. Without this argument only supplied scenarios can be checked.
#' @param cluster,cluster_variance Passed to [catmech_pool_lm()], so that fits made with
#'   and without `catmech_reference(cluster = )` can be compared under the same variance.
#'   Without `cluster`, a fit with its own cluster column uses a cluster-robust variance and
#'   a fit without one a model-based variance, and their widths are not comparable.
#' @param method Estimation method: `"pool_lm"` (default) for multiple-imputation
#'   pooling of linear fits, `"pool_glm"` for generalized linear model pooling,
#'   `"pool_coxph"` for Cox proportional hazards model pooling via [catmech_pool_coxph()],
#'   `"marginal_summary"` for population marginal inference from
#'   [catmech_marginal_summary()], or `"joint_summary"` for fixed-effect summaries from
#'   [catmech_joint_summary()].
#' @param family Error distribution and link function for `method = "pool_glm"`.
#'   Defaults to `stats::binomial()`.
#' @return A list with table, status, complete, successful_envelope and
#'   interpretation. An incomplete scenario set always has an unresolved status.
#' @examples
#' set.seed(123)
#' df <- data.frame(
#'   lab = stats::rnorm(30, 10, 2),
#'   age = stats::rnorm(30, 50, 10),
#'   y = stats::rnorm(30, 5, 1)
#' )
#' df$ref <- df$lab
#' df$ref[11:30] <- NA
#' res1 <- catmech_reference(df, value = "lab", reference = "ref",
#'                           predictors = c("age", "y"), m = 2, seed = 1)
#' res2 <- catmech_reference(df, value = "lab", reference = "ref",
#'                           predictors = c("age", "y"), use_recorded = TRUE, m = 2, seed = 2)
#' sens <- catmech_sensitivity(
#'   list(standard = res1, recorded = res2),
#'   formula = y ~ lab + age, term = "lab"
#' )
#' print(sens)
#' @export
catmech_sensitivity <- function(results, formula = NULL, term, reference = 0,
                                conf_level = 0.95, expected_scenarios = NULL,
                                cluster = NULL, cluster_variance = c("cv3", "cr1s"),
                                method = c("pool_lm", "pool_glm", "pool_coxph", "marginal_summary", "joint_summary"),
                                family = stats::binomial()) {
  cluster_variance <- match.arg(cluster_variance)
  method <- match.arg(method)
  if (!is.list(results) || !length(results) || is.null(names(results)) ||
      anyNA(names(results)) || any(!nzchar(names(results))) || anyDuplicated(names(results)))
    stop("results must be a nonempty list with unique scenario names.", call. = FALSE)
  if (!is.null(expected_scenarios)) {
    if (!is.character(expected_scenarios) || !length(expected_scenarios) ||
        anyNA(expected_scenarios) || any(!nzchar(expected_scenarios)) ||
        anyDuplicated(expected_scenarios)) stop("expected_scenarios must be unique nonempty names.", call. = FALSE)
    if (length(setdiff(names(results), expected_scenarios)))
      stop("results includes unplanned scenarios outside expected_scenarios.", call. = FALSE)
    for (name in setdiff(expected_scenarios, names(results)))
      results[[name]] <- simpleError("Prespecified scenario was not supplied.")
    results <- results[expected_scenarios]
  }
  valid <- vapply(results, inherits, logical(1), "catmech_result")
  if (any(valid)) {
    original <- results[[which(valid)[1]]]$original
    if (!all(vapply(results[valid], function(fit) identical(fit$original, original), logical(1))))
      stop("All scenarios must use the same original data; different cohorts are not an assumption comparison.", call. = FALSE)
  }
  if (!is.character(term) || length(term) != 1L || is.na(term) || !nzchar(term))
    stop("term must name one coefficient.", call. = FALSE)
  if (!is.numeric(reference) || length(reference) != 1L || !is.finite(reference))
    stop("reference must be one finite number.", call. = FALSE)
  # An invalid shared analysis setting is not a scenario-level failure.
  if (method %in% c("pool_lm", "pool_glm", "pool_coxph")) {
    if (!inherits(formula, "formula")) stop("formula must be an R formula.", call. = FALSE)
  } else if (method == "marginal_summary") {
    if (!is.null(formula) && !inherits(formula, "formula"))
      stop("formula must be an R formula or NULL.", call. = FALSE)
  }
  if (!is.numeric(conf_level) || length(conf_level) != 1L || !is.finite(conf_level) ||
      conf_level <= 0 || conf_level >= 1) stop("conf_level must be between zero and one.", call. = FALSE)
  records <- lapply(names(results), function(name) {
    fit <- results[[name]]
    summary <- if (inherits(fit, "error")) fit else tryCatch({
      if (method == "pool_lm") {
        catmech_pool_lm(fit, formula, conf_level, allow_sensitivity = TRUE, cluster = cluster,
                        cluster_variance = cluster_variance)
      } else if (method == "pool_glm") {
        catmech_pool_glm(fit, formula, family = family, conf_level = conf_level,
                         allow_sensitivity = TRUE, cluster = cluster)
      } else if (method == "pool_coxph") {
        catmech_pool_coxph(fit, formula, conf_level = conf_level,
                           allow_sensitivity = TRUE, cluster = cluster)
      } else if (method == "marginal_summary") {
        catmech_marginal_summary(fit, level = conf_level, formula = formula)
      } else {
        catmech_joint_summary(fit, level = conf_level)
      }
    }, error = identity)
    if (!inherits(summary, "error") && !term %in% summary$term)
      summary <- simpleError(paste("Coefficient not found:", term))
    if (inherits(summary, "error")) return(data.frame(scenario = name, status = "failed",
      estimate = NA_real_, std.error = NA_real_, lower = NA_real_, upper = NA_real_,
      direction = "unresolved", interpretation = "not_available",
      message = conditionMessage(summary), stringsAsFactors = FALSE))
    row <- summary[summary$term == term, , drop = FALSE]
    se <- if ("std.error" %in% names(row)) row$std.error else row$posterior_sd
    direction <- if (row$lower > reference) "above" else if (row$upper < reference) "below" else "includes_reference"
    data.frame(scenario = name, status = "ok", estimate = row$estimate,
      std.error = se, lower = row$lower, upper = row$upper,
      direction = direction, interpretation = row$interpretation, message = "",
      stringsAsFactors = FALSE)
  })
  table <- do.call(rbind, records)
  good <- table$status == "ok"
  complete <- all(good)
  status <- if (!complete) "unresolved_incomplete_set" else if (length(unique(table$direction)) > 1L)
    "assumption_sensitive" else if (table$direction[1] == "includes_reference")
      "unresolved_in_all_scenarios" else "direction_stable_over_supplied_scenarios"
  envelope <- if (any(good)) c(lower = min(table$lower[good]), upper = max(table$upper[good])) else
    c(lower = NA_real_, upper = NA_real_)
  structure(list(table = table, status = status, complete = complete, expected_scenarios = names(results), successful_envelope = envelope,
    interpretation = "Descriptive finite-scenario envelope; not a calibrated confidence interval for unknown assumptions."),
    class = c("catmech_sensitivity", "list"))
}

#' Print a catmech sensitivity analysis result
#'
#' @param x A `catmech_sensitivity` object.
#' @param digits Number of decimal digits to print (default 4).
#' @param ... Unused.
#' @return `x`, invisibly.
#' @export
print.catmech_sensitivity <- function(x, digits = 4L, ...) {
  cat("-- Catmech Sensitivity Analysis --\n")
  cat("Status:              ", x$status, " (complete: ", x$complete, ")\n", sep = "")
  env_str <- if (all(is.finite(x$successful_envelope))) {
    sprintf("[%s, %s]", round(x$successful_envelope[1L], digits), round(x$successful_envelope[2L], digits))
  } else "[NA, NA]"
  cat("Successful envelope: ", env_str, "\n", sep = "")
  cat("Scenarios:           ", nrow(x$table), " evaluated\n\n", sep = "")

  tbl <- x$table
  num_cols <- intersect(names(tbl), c("estimate", "std.error", "lower", "upper"))
  for (col in num_cols) {
    tbl[[col]] <- ifelse(is.finite(tbl[[col]]), round(tbl[[col]], digits), tbl[[col]])
  }
  display_cols <- intersect(c("scenario", "status", "estimate", "std.error", "lower", "upper", "direction"), names(tbl))
  print(tbl[, display_cols, drop = FALSE], row.names = FALSE)
  invisible(x)
}

#' Print a catmech result
#'
#' Displays a compact summary of a fitted `catmech_result`, including the engine,
#' target variable, reference subset, predictors, clustering, and convergence status.
#'
#' @param x A `catmech_result` object.
#' @param ... Additional arguments passed to print methods (currently ignored).
#' @return `x`, invisibly.
#' @export
print.catmech_result <- function(x, ...) {
  manifest <- x$manifest
  spec <- x$spec
  engine <- if (is.list(manifest) && !is.null(manifest$engine)) manifest$engine else "catmech_result"
  cat("<catmech_result: ", engine, ">\n", sep = "")
  if (is.list(spec)) {
    if (!is.null(spec$value)) cat("  Target value:      ", spec$value, "\n", sep = "")
    if (!is.null(spec$reference) && is.data.frame(x$original) && spec$reference %in% names(x$original)) {
      ref_obs <- sum(!is.na(x$original[[spec$reference]]))
      ref_total <- nrow(x$original)
      cat("  Reference subset:  ", spec$reference, " (", ref_obs, "/", ref_total, " rows, ",
          sprintf("%.1f%%", 100 * ref_obs / max(1L, ref_total)), ")\n", sep = "")
    }
    preds <- if (!is.null(spec$predictors)) spec$predictors else
      if (is.list(manifest) && !is.null(manifest$predictors)) manifest$predictors else NULL
    if (!is.null(preds) && length(preds)) {
      pred_str <- if (length(preds) <= 5L) paste(preds, collapse = ", ") else
        paste0(paste(preds[1:5], collapse = ", "), ", ... (", length(preds), " total)")
      cat("  Predictors:        ", pred_str, "\n", sep = "")
    }
    if (!is.null(spec$cluster) && is.data.frame(x$original) && spec$cluster %in% names(x$original)) {
      n_clust <- length(unique(x$original[[spec$cluster]]))
      clust_model_str <- if (!is.null(spec$cluster_model)) paste0(" [model: ", spec$cluster_model, "]") else ""
      cat("  Cluster:           ", spec$cluster, " (", n_clust, " clusters)", clust_model_str, "\n", sep = "")
    }
  }
  if (!is.null(x$completed)) {
    cat("  Completions (m):   ", length(x$completed), "\n", sep = "")
  }
  if (is.list(manifest) && !is.null(manifest$imputed_rows) && is.data.frame(x$original)) {
    cat("  Imputed rows:      ", length(manifest$imputed_rows), "/", nrow(x$original), "\n", sep = "")
  }
  if (!is.null(x$model_diagnostics) && is.list(x$model_diagnostics)) {
    diag <- x$model_diagnostics
    if (!is.null(diag$convergence_pass)) {
      status_str <- if (isTRUE(diag$convergence_pass)) "passed" else "failed"
      iters_str <- if (!is.null(diag$iterations_used)) paste0(" (", diag$iterations_used, " iterations)") else ""
      cat("  Model convergence: ", status_str, iters_str, "\n", sep = "")
    }
  }
  invisible(x)
}

#' Summary of a catmech result
#'
#' @param object A `catmech_result` object.
#' @param ... Additional arguments passed to summary methods (currently ignored).
#' @return An object of class `summary.catmech_result`.
#' @export
summary.catmech_result <- function(object, ...) {
  out <- list(
    manifest = object$manifest,
    spec = object$spec,
    n_rows = if (is.data.frame(object$original)) nrow(object$original) else NA_integer_,
    m = length(object$completed),
    imputed_count = if (is.list(object$manifest) && !is.null(object$manifest$imputed_rows))
      length(object$manifest$imputed_rows) else NA_integer_
  )
  if (!is.null(object$model_diagnostics) && is.list(object$model_diagnostics)) {
    diag <- object$model_diagnostics
    out$convergence_pass <- diag$convergence_pass
    out$iterations_used <- diag$iterations_used
    if (is.data.frame(diag$summary)) {
      out$parameter_summary <- diag$summary
    }
  }
  class(out) <- "summary.catmech_result"
  out
}

#' Print a summary of a catmech result
#'
#' @param x An object of class `summary.catmech_result`.
#' @param ... Additional arguments (currently ignored).
#' @return `x`, invisibly.
#' @export
print.summary.catmech_result <- function(x, ...) {
  engine <- if (is.list(x$manifest) && !is.null(x$manifest$engine)) x$manifest$engine else "catmech_result"
  cat("<summary.catmech_result: ", engine, ">\n", sep = "")
  if (is.list(x$spec)) {
    if (!is.null(x$spec$value)) cat("  Target value:    ", x$spec$value, "\n", sep = "")
    if (!is.null(x$spec$reference)) cat("  Reference:       ", x$spec$reference, "\n", sep = "")
    preds <- if (!is.null(x$spec$predictors)) x$spec$predictors else
      if (is.list(x$manifest) && !is.null(x$manifest$predictors)) x$manifest$predictors else NULL
    if (!is.null(preds) && length(preds))
      cat("  Predictors:      ", paste(preds, collapse = ", "), "\n", sep = "")
    if (!is.null(x$spec$cluster))
      cat("  Cluster:         ", x$spec$cluster, if (!is.null(x$spec$cluster_model)) paste0(" (", x$spec$cluster_model, ")") else "", "\n", sep = "")
  }
  cat("  Observations:    ", x$n_rows, "\n", sep = "")
  cat("  Completions (m): ", x$m, "\n", sep = "")
  cat("  Imputed rows:    ", x$imputed_count, "\n", sep = "")
  if (!is.null(x$convergence_pass)) {
    cat("  Convergence:     ", if (isTRUE(x$convergence_pass)) "passed" else "failed",
        if (!is.null(x$iterations_used)) paste0(" (", x$iterations_used, " iterations)") else "", "\n", sep = "")
  }
  if (!is.null(x$parameter_summary) && is.data.frame(x$parameter_summary)) {
    cat("\nParameter diagnostics (first 8 monitors):\n")
    preview <- utils::head(x$parameter_summary, 8L)
    print(preview, row.names = FALSE)
    if (nrow(x$parameter_summary) > 8L)
      cat("... (", nrow(x$parameter_summary) - 8L, " additional monitors omitted)\n", sep = "")
  }
  invisible(x)
}

#' Print model diagnostics
#'
#' Displays a compact summary of MCMC sampler diagnostics without dumping raw traces.
#'
#' @param x An object of class `catmech_model_diagnostics`.
#' @param ... Additional arguments (currently ignored).
#' @return `x`, invisibly.
#' @export
print.catmech_model_diagnostics <- function(x, ...) {
  engine <- if (!is.null(x$engine)) x$engine else "sampler"
  cat("<catmech_model_diagnostics: ", engine, ">\n", sep = "")
  status <- if (!is.null(x$status)) x$status else "unknown"
  cat("  Status:          ", status,
      if (isTRUE(x$convergence_pass)) " (convergence checks passed)" else " (checks failed)", "\n", sep = "")
  if (is.array(x$traces)) {
    d <- dim(x$traces)
    if (length(d) == 3L) {
      cat("  MCMC traces:     ", d[1L], " iterations, ", d[2L], " chains, ", d[3L], " parameters\n", sep = "")
    } else if (length(d) == 2L) {
      cat("  Draws:           ", d[1L], " draws, ", d[2L], " parameters\n", sep = "")
    }
  }
  if (is.data.frame(x$summary)) {
    cat("  Summary (first 8 monitors):\n")
    preview <- utils::head(x$summary, 8L)
    print(preview, row.names = FALSE)
    if (nrow(x$summary) > 8L)
      cat("  ... (", nrow(x$summary) - 8L, " additional monitors in $summary)\n", sep = "")
  }
  invisible(x)
}
