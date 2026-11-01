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
#' @details When the workflow declares that no target cells were imputed, all
#'   completed analysis frames must agree and ordinary complete-data t intervals
#'   are returned. A zero estimated between-imputation variance alone does not
#'   activate this rule. Predictive matching and reference-subset results retain
#'   distinct interpretation labels; none guarantees coverage under misspecification.
#'   Partial-reference PMM requires explicit sensitivity opt-in because held-out
#'   studies found material undercoverage even with random accurate references.
#'   This package policy does not assert that every PMM application is invalid.
#'   Only univariate linear models with at least one named coefficient are supported.
#' @return A data.frame with term, estimate, std.error, df, lower, upper,
#'   within_variance, between_variance, m, n and interpretation.
#' @export
catmech_pool_lm <- function(result, formula, conf_level = 0.95,
                            allow_sensitivity = FALSE) {
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
  output <- lapply(terms, function(term) {
    Q <- vapply(fits, function(f) unname(stats::coef(f)[term]), numeric(1))
    U <- vapply(fits, function(f) unname(stats::vcov(f)[term, term]), numeric(1))
    pooled <- if (no_imputation) list(qbar = Q[1], t = U[1], ubar = U[1],
                                      b = 0, df = dfs[1]) else
      mice::pool.scalar(Q, U, n = ns[1], k = ns[1] - dfs[1])
    if (!is.finite(pooled$qbar) || !is.finite(pooled$t) || pooled$t <= 0 ||
        is.na(pooled$df) || pooled$df <= 0)
      stop("Pooled uncertainty is degenerate or nonfinite; a confidence interval cannot be reported.", call. = FALSE)
    se <- sqrt(pooled$t)
    critical <- stats::qt((1 + conf_level)/2, pooled$df)
    data.frame(term = term, estimate = pooled$qbar, std.error = se, df = pooled$df,
      lower = pooled$qbar - critical*se, upper = pooled$qbar + critical*se,
      within_variance = pooled$ubar, between_variance = pooled$b,
      m = length(fits), n = ns[1],
      interpretation = if (recognized && !reference_pmm_sensitivity) unname(labels[engine]) else "sensitivity_only",
      stringsAsFactors = FALSE)
  })
  do.call(rbind, output)
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
#' @return A list with table, status, complete, successful_envelope and
#'   interpretation. An incomplete scenario set always has an unresolved status.
#' @export
catmech_sensitivity <- function(results, formula, term, reference = 0,
                                conf_level = 0.95, expected_scenarios = NULL) {
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
      catmech_pool_lm(fit, formula, conf_level, allow_sensitivity = TRUE), error = identity)
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
