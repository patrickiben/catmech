#' Summarize the population coefficients of a fitted joint cluster model
#'
#' Summarizes the true-measurement fixed effects from every retained posterior
#' draw, on the original response and predictor scales. These are model-based
#' Bayesian credible intervals, conditional on the joint model and its priors.
#' They are not cluster-robust confidence intervals or summaries of a different
#' downstream analysis.
#'
#' @param result A result from [catmech_reference()] with
#'   `cluster_model = "joint"` and a fitted, converged posterior. The all-reference
#'   shortcut has no parameter posterior and is refused without refitting.
#' @param level Central posterior credible interval probability, strictly between
#'   zero and one.
#'
#' @details
#' The estimand is the population fixed-effect vector for the true measurement
#' in the original imputation design. Under the model's zero-mean random effects
#' independent of the covariates, these coefficients describe the marginal mean
#' of the true measurement as well as the mean conditional on a zero site effect.
#' They are not the least-squares coefficients of the realized completed cohort.
#' The function has no formula argument: it cannot change predictors, add an
#' offset or interaction, or use the repaired measurement as a predictor.
#'
#' All retained iterations from all chains are used, independently of the number
#' of completed datasets. The mean, sample standard deviation and equal-tailed
#' quantiles (R quantile type 7) summarize those draws. There is no additional
#' Rubin pooling, cluster sandwich variance, or Student t multiplier.
#'
#' Fixed effects are transformed jointly, including the contribution of predictor
#' centering to the original-scale intercept. Rank-normalized R-hat and bulk and
#' tail effective sample sizes are recomputed for these coefficients. R-hat must
#' be at most 1.05 and both effective sample sizes at least `100 * chains`.
#' Failure raises a `catmech_convergence_error` carrying `model_diagnostics` for
#' the transformed coefficient draws. The `posterior` package supplies Monte
#' Carlo errors for the mean and interval endpoints. These quantify simulation
#' error, not uncertainty about model adequacy or the population coefficient.
#' Passing diagnostics does not establish frequentist interval coverage.
#'
#' @return A data frame with one row per original design coefficient: `term`,
#'   posterior mean `estimate`, `posterior_sd`, `lower`, `upper`, `level`,
#'   `rhat`, `ess_bulk`, `ess_tail`, `mcse_mean`, `mcse_lower`, `mcse_upper`,
#'   `draws`, `iterations`, `chains`, and `interpretation`.
#' @export
catmech_joint_summary <- function(result, level = .95) {
  if (!is.list(result) || !inherits(result, "catmech_result") ||
      !is.list(result$manifest) ||
      !identical(result$manifest$engine, "reference_joint_cluster"))
    stop("Expected a result from cluster_model='joint'.", call. = FALSE)
  if (!is.numeric(level) || is.complex(level) || length(level) != 1L || !is.finite(level) ||
      level <= 0 || level >= 1)
    stop("level must be strictly between zero and one.", call. = FALSE)
  spec <- result$spec
  diagnostics <- result$model_diagnostics
  if (!is.list(spec) || !identical(spec$cluster_model, "joint") ||
      !is.list(diagnostics))
    stop("The joint model specification and diagnostics are missing or malformed.", call. = FALSE)
  if (identical(spec$recorded_fit, "not_fitted_all_reference_observed") ||
      identical(diagnostics$engine, "reference_complete"))
    stop("The all-reference shortcut did not fit a joint parameter posterior; no coefficient summary is available.", call. = FALSE)
  if (!isTRUE(diagnostics$convergence_pass))
    catmech_joint_summary_failure("Joint model convergence checks did not pass.", diagnostics)
  if (!identical(diagnostics$engine, "pan_joint_cluster"))
    stop("The result does not contain joint posterior sampler diagnostics.", call. = FALSE)
  traces <- diagnostics$traces
  dimensions <- dim(traces)
  if (!is.numeric(traces) || is.complex(traces) || length(dimensions) != 3L ||
      dimensions[1L] < 100L || dimensions[2L] < 4L || dimensions[3L] < 1L ||
      any(!is.finite(traces)))
    stop("Joint posterior traces must be a finite iterations-by-chains-by-parameters array with at least 100 iterations and four chains.", call. = FALSE)
  parameters <- dimnames(traces)[[3L]]
  if (is.null(parameters) || length(parameters) != dimensions[3L] ||
      anyNA(parameters) || any(!nzchar(parameters)) || anyDuplicated(parameters))
    stop("Joint posterior parameter names must be unique and nonempty.", call. = FALSE)
  control <- catmech_joint_cluster_control(spec$joint_control)
  if (!is.numeric(diagnostics$iterations_used) || is.complex(diagnostics$iterations_used) ||
      length(diagnostics$iterations_used) != 1L ||
      !is.finite(diagnostics$iterations_used) ||
      diagnostics$iterations_used != dimensions[1L] ||
      dimensions[1L] < control$iterations || dimensions[1L] > control$max_iterations ||
      dimensions[2L] != control$chains)
    stop("Joint posterior trace dimensions conflict with the retained sampler controls.", call. = FALSE)
  if (!is.data.frame(result$original) || !is.character(spec$predictors) ||
      anyNA(spec$predictors) || anyDuplicated(spec$predictors) ||
      any(!spec$predictors %in% names(result$original)))
    stop("The original imputation design is missing or malformed.", call. = FALSE)
  x <- catmech_model_matrix(result$original, spec$predictors)
  terms <- colnames(x)
  beta_t <- paste0("beta_T:", terms)
  beta_w <- paste0("beta_W:", terms)
  covariance <- c("sigma_TT", "sigma_TW", "sigma_WW", "psi_TT", "psi_TW", "psi_WW")
  if (!all(c(beta_t, beta_w, covariance) %in% parameters) ||
      !setequal(parameters[startsWith(parameters, "beta_T:")], beta_t) ||
      !setequal(parameters[startsWith(parameters, "beta_W:")], beta_w))
    stop("Joint coefficient traces do not match the original imputation design.", call. = FALSE)
  for (prefix in c("sigma", "psi")) {
    tt <- traces[, , paste0(prefix, "_TT")]
    tw <- traces[, , paste0(prefix, "_TW")]
    ww <- traces[, , paste0(prefix, "_WW")]
    if (any(tt <= 0 | ww <= 0))
      stop("Joint covariance traces must be positive definite.", call. = FALSE)
    correlation <- tw / sqrt(tt) / sqrt(ww)
    if (any(!is.finite(correlation)) || any(abs(correlation) >= 1))
      stop("Joint covariance traces must be positive definite.", call. = FALSE)
  }
  transform <- spec$joint_transform
  if (!is.list(transform) || !identical(transform, diagnostics$transform))
    stop("Joint transformations are missing or disagree with the sampler metadata.", call. = FALSE)
  extract <- function(name, expected, positive = FALSE) {
    z <- transform[[name]]
    if (!is.list(z) || is.null(names(z)) || anyNA(names(z)) ||
        anyDuplicated(names(z)) || !setequal(names(z), expected) ||
        !all(vapply(z, function(v) is.numeric(v) && !is.complex(v) && length(v) == 1L &&
          is.finite(v) && (!positive || v > 0), logical(1))))
      stop("Malformed joint transformation: ", name, ".", call. = FALSE)
    stats::setNames(vapply(z[expected], as.numeric, numeric(1)), expected)
  }
  response_center <- extract("response_center", c("T", "W"))
  response_scale <- extract("response_scale", c("T", "W"), positive = TRUE)
  predictor_center <- extract("predictor_center", terms)
  predictor_scale <- extract("predictor_scale", terms, positive = TRUE)
  if (!identical(terms[1L], "(Intercept)") || predictor_center[1L] != 0 ||
      predictor_scale[1L] != 1)
    stop("Joint transformations require the original unscaled intercept column.", call. = FALSE)
  coefficients <- traces[, , beta_t, drop = FALSE]
  dimnames(coefficients)[[3L]] <- terms
  for (j in seq_along(terms))
    coefficients[, , j] <- coefficients[, , j] * response_scale["T"] / predictor_scale[j]
  coefficients[, , 1L] <- coefficients[, , 1L] + response_center["T"]
  if (length(terms) > 1L) {
    for (j in 2L:length(terms))
      coefficients[, , 1L] <- coefficients[, , 1L] - predictor_center[j] * coefficients[, , j]
  }
  if (any(!is.finite(coefficients)))
    stop("Original-scale joint coefficient draws are nonfinite.", call. = FALSE)
  if ("beta_T_original:(Intercept)" %in% parameters &&
      any(abs(traces[, , "beta_T_original:(Intercept)"] - coefficients[, , 1L]) >
        sqrt(.Machine$double.eps) * (1 + abs(coefficients[, , 1L]))))
    stop("The original-scale intercept monitor conflicts with the coefficient transformations.", call. = FALSE)
  checked <- catmech_model_chain_diagnostics(coefficients)
  checked$engine <- "joint_population_fixed_effects"
  checked$level <- level
  if (!checked$convergence_pass)
    catmech_joint_summary_failure(paste0("Original-scale joint coefficients failed convergence checks ",
      "(R-hat <= 1.05 and bulk/tail ESS >= ", 100L * dimensions[2L], " required)."), checked)
  probabilities <- c((1 - level) / 2, (1 + level) / 2)
  output <- lapply(seq_along(terms), function(j) {
    samples <- coefficients[, , j]
    interval <- stats::quantile(samples, probabilities, names = FALSE, type = 7L)
    mcse <- tryCatch(c(posterior::mcse_mean(samples),
      posterior::mcse_quantile(samples, probabilities, names = FALSE)),
      error = function(e) rep(NA_real_, 3L))
    data.frame(term = terms[j], estimate = mean(samples), posterior_sd = stats::sd(as.vector(samples)),
      lower = interval[1L], upper = interval[2L], level = level,
      rhat = checked$summary$rhat[j], ess_bulk = checked$summary$ess_bulk[j],
      ess_tail = checked$summary$ess_tail[j], mcse_mean = mcse[1L],
      mcse_lower = mcse[2L], mcse_upper = mcse[3L],
      draws = prod(dimensions[1:2]), iterations = dimensions[1L], chains = dimensions[2L],
      interpretation = "conditional_joint_model_credible_interval", stringsAsFactors = FALSE)
  })
  output <- do.call(rbind, output)
  if (any(!is.finite(as.matrix(output[c("estimate", "posterior_sd", "lower", "upper",
      "mcse_mean", "mcse_lower", "mcse_upper")]))) || any(output$posterior_sd <= 0)) {
    checked$status <- "failed"
    checked$convergence_pass <- FALSE
    checked$coefficient_summary <- output
    catmech_joint_summary_failure("Finite posterior summaries and Monte Carlo errors could not be computed at this level.", checked)
  }
  output
}

catmech_joint_summary_failure <- function(message, diagnostics) {
  stop(structure(list(message = message, call = NULL, model_diagnostics = diagnostics),
    class = c("catmech_convergence_error", "error", "condition")))
}
