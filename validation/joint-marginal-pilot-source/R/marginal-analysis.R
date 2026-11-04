#' Summarize population marginal coefficients from a joint covariate model
#'
#' Summarizes the population regression of the true measurement on the original
#' continuous predictors from a fitted multivariate Gaussian cluster model.
#' Predictor means and covariances are estimated jointly with the measurement
#' model, allowing site effects in the measurement and predictors to be correlated.
#'
#' @param result A converged result from [catmech_reference()] with
#'   `cluster_model = "joint_marginal"`. An all-reference shortcut contains no
#'   fitted parameter posterior and cannot be summarized by this function.
#' @param level Central posterior credible interval probability, strictly between
#'   zero and one.
#' @details
#' The model is `(T,W,X) = mu + b_site + e`, with independent Gaussian site and
#' row vectors having unrestricted covariance matrices `Psi` and `Sigma`.
#' For every posterior draw, let `V = Psi + Sigma`. The marginal slopes are
#' `solve(V[X,X], V[X,T])` and the intercept is `mu[T] - mu[X]' beta`, transformed
#' back to the original units. This is the population linear regression for an
#' individual sampled under the model, not a within-site coefficient or the OLS
#' slope of the completed cohort. Site size must be noninformative. The Gaussian
#' covariate model is an additional assumption; categorical, incomplete and
#' deterministically related predictors are not supported by this engine.
#'
#' The covariance priors are hierarchical Huang--Wand priors with degrees of
#' freedom parameter two and marginal half-t(2) standard deviations. The scale
#' is the square root of `joint_control$prior_scale` in empirically standardized
#' units. The mean prior is flat. Credible intervals remain conditional on these
#' distributional and reference-selection assumptions. They are not
#' cluster-robust intervals and do not establish robustness to arbitrary
#' contextual effects, nonlinear associations or informative cluster sizes.
#'
#' Every retained draw from every chain is used independently of `m`. Equal-tailed
#' quantiles use R quantile type 7. Original-unit coefficient diagnostics are
#' recomputed: R-hat at most 1.05 and both bulk/tail ESS at least `100 * chains`.
#' The reported Monte Carlo errors concern simulation error only.
#' @return A data frame with `term`, `estimate`, `posterior_sd`, `lower`, `upper`,
#'   `level`, `rhat`, `ess_bulk`, `ess_tail`, `mcse_mean`, `mcse_lower`,
#'   `mcse_upper`, `draws`, `iterations`, `chains`, and `interpretation`.
#' @references Huang, A. and Wand, M. P. (2013). Simple marginally noninformative
#'   prior distributions for covariance matrices. Bayesian Analysis, 8, 439--452.
#'   \doi{10.1214/13-BA815}.
#' @export
catmech_marginal_summary <- function(result, level = .95) {
  if (!inherits(result, "catmech_result") || !is.list(result) ||
      !is.list(result$manifest) ||
      !identical(result$manifest$engine, "reference_marginal_cluster"))
    stop("Expected a result from cluster_model='joint_marginal'.", call. = FALSE)
  if (!is.numeric(level) || is.complex(level) || length(level) != 1L ||
      !is.finite(level) || level <= 0 || level >= 1)
    stop("level must be strictly between zero and one.", call. = FALSE)
  spec <- result$spec
  d <- result$model_diagnostics
  if (!is.list(spec) || !identical(spec$cluster_model, "joint_marginal") || !is.list(d))
    stop("The marginal model specification or diagnostics are missing.", call. = FALSE)
  if (identical(spec$recorded_fit, "not_fitted_all_reference_observed") ||
      identical(d$engine, "reference_complete"))
    stop("The all-reference shortcut did not fit a marginal parameter posterior.", call. = FALSE)
  if (!isTRUE(d$convergence_pass))
    catmech_joint_summary_failure("Marginal model convergence checks did not pass.", d)
  if (!identical(d$engine, "hw_marginal_cluster"))
    stop("The result does not contain marginal sampler diagnostics.", call. = FALSE)
  control <- catmech_joint_cluster_control(spec$joint_control)
  tr <- d$traces
  dims <- dim(tr)
  if (!is.numeric(tr) || is.complex(tr) || length(dims) != 3L ||
      dims[1L] < control$iterations || dims[1L] > control$max_iterations ||
      dims[2L] != control$chains || !identical(as.integer(d$iterations_used), dims[1L]) ||
      any(!is.finite(tr)))
    stop("Malformed marginal posterior traces or inconsistent sampler controls.", call. = FALSE)
  names <- dimnames(tr)[[3L]]
  if (is.null(names) || anyNA(names) || any(!nzchar(names)) || anyDuplicated(names))
    stop("Marginal posterior trace names must be unique and nonempty.", call. = FALSE)
  predictors <- spec$predictors
  if (!is.character(predictors) || anyNA(predictors) || anyDuplicated(predictors) ||
      "(Intercept)" %in% predictors || !is.data.frame(result$original) ||
      any(!predictors %in% names(result$original)) ||
      !all(vapply(result$original[predictors], function(x)
        is.numeric(x) && !is.complex(x) && is.null(dim(x)) && all(is.finite(x)), logical(1))))
    stop("The original continuous predictor design is missing or malformed.", call. = FALSE)
  terms <- c("(Intercept)", predictors)
  parameters <- paste0("marginal:", terms)
  internal <- c("T", "W", if (length(predictors)) paste0("X", seq_along(predictors)))
  k <- length(internal)
  pairs <- which(lower.tri(matrix(0, k, k), diag = TRUE), arr.ind = TRUE)
  suffix <- paste(internal[pairs[, 1L]], internal[pairs[, 2L]], sep = ":")
  required <- c(paste0("mu:", internal), paste0("sigma:", suffix),
    paste0("psi:", suffix), parameters)
  if (!all(required %in% names) || !setequal(names[startsWith(names, "marginal:")], parameters))
    stop("Marginal traces do not match the original predictor specification.", call. = FALSE)
  transform <- spec$joint_transform
  if (!is.list(transform) || !identical(transform, d$transform) ||
      !identical(as.character(unlist(transform$predictor_names)), unname(predictors)))
    stop("Marginal transformations are missing or inconsistent.", call. = FALSE)
  for (field in c("response_center", "response_scale")) {
    z <- transform[[field]]
    if (!is.list(z) || !identical(names(z), internal) ||
        !all(vapply(z, function(x) is.numeric(x) && !is.complex(x) &&
          length(x) == 1L && is.finite(x), logical(1))) ||
        (field == "response_scale" && any(unlist(z) <= 0)))
      stop("Malformed marginal transformation: ", field, ".", call. = FALSE)
  }
  # Reconstruct the estimand from the actual model draws, checking each covariance.
  coefficients <- tr[, , parameters, drop = FALSE]
  dimnames(coefficients)[[3L]] <- terms
  flat <- matrix(tr, nrow = prod(dims[1:2]), dimnames = list(NULL, names))
  rebuilt <- matrix(NA_real_, nrow(flat), length(terms))
  if (k <= 3L) {
    # Sylvester's criterion on correlation matrices checks the same positive
    # definiteness as Cholesky, with scale-free arithmetic for the common case.
    for (prefix in c("sigma:", "psi:")) {
      variances <- flat[, paste0(prefix, internal, ":", internal), drop = FALSE]
      if (any(variances <= 0))
        stop("Marginal covariance traces must be positive definite.", call. = FALSE)
      r12 <- flat[, paste0(prefix, "W:T")] / sqrt(variances[, 1L]) / sqrt(variances[, 2L])
      valid <- is.finite(r12) & abs(r12) < 1
      if (k == 3L) {
        r13 <- flat[, paste0(prefix, "X1:T")] / sqrt(variances[, 1L]) / sqrt(variances[, 3L])
        r23 <- flat[, paste0(prefix, "X1:W")] / sqrt(variances[, 2L]) / sqrt(variances[, 3L])
        determinant <- 1 + 2 * r12 * r13 * r23 - r12^2 - r13^2 - r23^2
        valid <- valid & is.finite(determinant) & determinant > 0 &
          abs(r13) < 1 & abs(r23) < 1
      }
      if (!all(valid))
        stop("Marginal covariance traces must be positive definite.", call. = FALSE)
    }
    center <- unlist(transform$response_center)
    scale <- unlist(transform$response_scale)
    rebuilt[, 1L] <- center["T"] + scale["T"] * flat[, "mu:T"]
    if (k == 3L) {
      rebuilt[, 2L] <- (flat[, "sigma:X1:T"] + flat[, "psi:X1:T"]) /
        (flat[, "sigma:X1:X1"] + flat[, "psi:X1:X1"]) * scale["T"] / scale["X1"]
      rebuilt[, 1L] <- rebuilt[, 1L] - rebuilt[, 2L] *
        (center["X1"] + scale["X1"] * flat[, "mu:X1"])
    }
  } else {
    for (i in seq_len(nrow(flat))) {
      matrices <- lapply(c("sigma:", "psi:"), function(prefix) {
        z <- matrix(0, k, k)
        z[pairs] <- flat[i, paste0(prefix, suffix)]
        z <- z + t(z) - diag(diag(z), k)
        if (inherits(try(chol(z), silent = TRUE), "try-error"))
          stop("Marginal covariance traces must be positive definite.", call. = FALSE)
        z
      })
      rebuilt[i, ] <- catmech_marginal_coefficients(flat[i, paste0("mu:", internal)],
        matrices[[1L]], matrices[[2L]], transform, predictors)
    }
  }
  if (any(!is.finite(rebuilt)) || any(abs(rebuilt - matrix(coefficients, nrow(flat))) >
      sqrt(.Machine$double.eps) * (1 + abs(rebuilt))))
    stop("Marginal coefficient traces conflict with the model parameters.", call. = FALSE)
  checked <- catmech_model_chain_diagnostics(coefficients)
  checked$engine <- "joint_population_marginal_coefficients"
  checked$level <- level
  if (!checked$convergence_pass)
    catmech_joint_summary_failure("Marginal coefficient draws failed convergence checks.", checked)
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
      draws = prod(dims[1:2]), iterations = dims[1L], chains = dims[2L],
      interpretation = "conditional_joint_marginal_credible_interval", stringsAsFactors = FALSE)
  })
  output <- do.call(rbind, output)
  if (any(!is.finite(as.matrix(output[c("estimate", "posterior_sd", "lower", "upper",
      "mcse_mean", "mcse_lower", "mcse_upper")]))) || any(output$posterior_sd <= 0))
    catmech_joint_summary_failure("Finite marginal summaries or Monte Carlo errors could not be computed.", checked)
  output
}
