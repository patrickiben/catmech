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
#' @return A data.frame with term, estimate, std.error, df, lower, upper,
#'   within_variance, between_variance, m, n and interpretation, plus clusters
#'   when a cluster-robust variance was used.
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
  if (!is.null(result[["spec"]][["cluster"]]) && !is.null(result[["spec"]][["value"]]) &&
      result[["spec"]][["value"]] %in% all.vars(stats::delete.response(stats::terms(formula)))) {
    if (identical(engine, "reference_joint_cluster"))
      warning("Joint clustered imputation has not been validated here with the audited value as an analysis predictor. The earlier moment model showed slope attenuation in that role; joint sampling alone does not establish compatible inference. See docs/JOINT_CLUSTER_METHOD.md.", call. = FALSE)
    else
      warning("catmech_reference(cluster = ) gives an unaudited cluster a new intercept, which biased slopes toward zero in simulations when the audited value was a predictor. Prefer a fit without cluster, pooled with catmech_pool_lm(cluster = ). See validation/VARIANCE_MI_SYNTHESIS.md.", call. = FALSE)
  }
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
  do.call(rbind, output)
}

# Delete-one-cluster jackknife covariance of one fitted linear model (CV3), by the closed
# form for each deletion. A singular block, from a cluster that alone determines a
# coefficient, uses a pseudoinverse.
catmech_cluster_jackknife <- function(fit, ids) {
  x <- stats::model.matrix(fit)
  e <- stats::residuals(fit)
  bread <- chol2inv(qr.R(fit$qr))
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
  x <- stats::model.matrix(fit)
  bread <- chol2inv(qr.R(fit$qr))
  flagged <- logical(ncol(x))
  for (u in split(seq_len(nrow(x)), as.character(ids))) {
    xg <- x[u, , drop = FALSE]
    block <- eigen(diag(length(u)) - xg %*% bread %*% t(xg), symmetric = TRUE)
    zero <- block$values <= 1e-10
    if (!any(zero)) next
    direction <- bread %*% crossprod(xg, block$vectors[, zero, drop = FALSE])
    flagged <- flagged | apply(abs(direction), 1L, max) > 1e-8 * max(abs(direction))
  }
  names(stats::coef(fit))[flagged]
}

# Cluster-robust covariance of one fitted linear model with factor (G / (G - 1)) ((n - 1) / (n - p)).
catmech_cluster_vcov <- function(fit, ids) {
  x <- stats::model.matrix(fit)
  n <- nrow(x); k <- ncol(x)
  sums <- rowsum(x * stats::residuals(fit), as.character(ids), reorder = FALSE)
  g <- nrow(sums)
  bread <- chol2inv(qr.R(fit$qr))
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
#' @param term Name of the coefficient to compare.
#' @param reference Reference value, usually zero, for interval-direction reporting.
#' @param conf_level Conditional confidence level passed to [catmech_pool_lm()].
#' @param expected_scenarios Optional character vector naming the complete
#'   prespecified set. Omitted members are recorded as failures; unplanned extras
#'   are rejected. Without this argument only supplied scenarios can be checked.
#' @param cluster,cluster_variance Passed to [catmech_pool_lm()], so that fits made with
#'   and without `catmech_reference(cluster = )` can be compared under the same variance.
#'   Without `cluster`, a fit with its own cluster column uses a cluster-robust variance and
#'   a fit without one a model-based variance, and their widths are not comparable.
#' @return A list with table, status, complete, successful_envelope and
#'   interpretation. An incomplete scenario set always has an unresolved status.
#' @export
catmech_sensitivity <- function(results, formula, term, reference = 0,
                                conf_level = 0.95, expected_scenarios = NULL,
                                cluster = NULL, cluster_variance = c("cv3", "cr1s")) {
  cluster_variance <- match.arg(cluster_variance)
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
  if (!inherits(formula,"formula")) stop("formula must be an R formula.", call.=FALSE)
  if (!is.numeric(conf_level) || length(conf_level)!=1L || !is.finite(conf_level) ||
      conf_level<=0 || conf_level>=1) stop("conf_level must be between zero and one.",call.=FALSE)
  records <- lapply(names(results), function(name) {
    fit <- results[[name]]
    summary <- if (inherits(fit, "error")) fit else tryCatch(
      catmech_pool_lm(fit, formula, conf_level, allow_sensitivity = TRUE, cluster = cluster,
                      cluster_variance = cluster_variance), error = identity)
    if (!inherits(summary, "error") && !term %in% summary$term)
      summary <- simpleError(paste("Coefficient not found:", term))
    if (inherits(summary, "error")) return(data.frame(scenario = name, status = "failed",
      estimate = NA_real_, std.error = NA_real_, lower = NA_real_, upper = NA_real_,
      direction = "unresolved", interpretation = "not_available",
      message = conditionMessage(summary), stringsAsFactors=FALSE))
    row <- summary[summary$term == term, , drop=FALSE]
    direction <- if (row$lower > reference) "above" else if (row$upper < reference) "below" else "includes_reference"
    data.frame(scenario = name, status = "ok", estimate = row$estimate,
      std.error = row$std.error, lower = row$lower, upper = row$upper,
      direction = direction, interpretation = row$interpretation, message = "",
      stringsAsFactors=FALSE)
  })
  table <- do.call(rbind, records)
  good <- table$status == "ok"
  complete <- all(good)
  status <- if (!complete) "unresolved_incomplete_set" else if (length(unique(table$direction)) > 1L)
    "assumption_sensitive" else if (table$direction[1] == "includes_reference")
      "unresolved_in_all_scenarios" else "direction_stable_over_supplied_scenarios"
  envelope <- if (any(good)) c(lower=min(table$lower[good]),upper=max(table$upper[good])) else
    c(lower=NA_real_,upper=NA_real_)
  list(table=table, status=status, complete=complete, expected_scenarios=names(results), successful_envelope=envelope,
    interpretation="Descriptive finite-scenario envelope; not a calibrated confidence interval for unknown assumptions.")
}
