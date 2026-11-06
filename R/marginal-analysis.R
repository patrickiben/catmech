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
#' @param formula Optional two-sided additive linear formula with an intercept
#'   (for example `outcome ~ recorded + x`) naming a response and covariates from
#'   the audited value (`spec$value`, mapped to the true measurement `T`) and the
#'   complete numeric predictors in `spec$predictors`. Use this when the audited
#'   measurement is an analysis predictor and the scientific outcome was included
#'   in `predictors` during [catmech_reference()]. Cannot be combined with
#'   `response` or `covariates`.
#' @param response Optional name of the population regression response, either
#'   `spec$value` or one of `spec$predictors`. Defaults to `spec$value`.
#' @param covariates Optional character vector of population regression covariates
#'   chosen from `c(spec$value, spec$predictors)` excluding `response`. When
#'   `response` is supplied and `covariates` is `NULL`, all remaining modeled
#'   variables in `c(spec$value, spec$predictors)` are used.
#' @details
#' The model is `(T,W,X) = mu + b_site + e`, with independent Gaussian site and
#' row vectors having unrestricted covariance matrices `Psi` and `Sigma`.
#' For every posterior draw, let `V = Psi + Sigma`. By default, the marginal slopes
#' are `solve(V[X,X], V[X,T])` and the intercept is `mu[T] - mu[X]' beta`,
#' transformed back to the original units. When `formula`, `response` or
#' `covariates` selects another modeled response `Y` and covariate set `P` from
#' `(T, X)` (excluding the noisy proxy `W`), the marginal slopes are
#' `solve(V[P,P], V[P,Y])` and the intercept is `mu[Y] - mu[P]' beta` in original
#' units. This allows the audited measurement `T` to serve as an analysis predictor
#' (`Y ~ T + X`) without the unaudited-cluster slope attenuation of the moment
#' model. This is the population linear regression for an
#' individual sampled under the model, not a within-site coefficient or the OLS
#' slope of the completed cohort. Site size must be noninformative. The Gaussian
#' covariate model is an additional assumption; categorical, incomplete and
#' deterministically related predictors are not supported by this engine.
#'
#' Each site/row covariance is factored into predictor covariance, regression
#' of the two measurements on predictors, and conditional measurement covariance.
#' The two covariance blocks have Huang--Wand priors with half-t(2) marginal
#' standard deviations and scale `sqrt(joint_control$prior_scale)` in empirically
#' standardized units. Those properties do not apply to the reconstructed full
#' covariance. Regression matrices have proper zero-mean matrix-normal priors,
#' with row covariance `regression_scale^2 * I` (default scale ten) and column
#' covariance equal to the conditional measurement covariance. The population
#' mean prior is flat. Credible intervals remain conditional on these
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
catmech_marginal_summary <- function(result, level = .95, formula = NULL,
                                     response = NULL, covariates = NULL) {
  if (!inherits(result, "catmech_result") || !is.list(result) ||
      !is.list(result$manifest) ||
      !identical(result$manifest$engine, "reference_marginal_cluster"))
    stop("Expected a result from cluster_model='joint_marginal'.", call. = FALSE)
  if (!is.numeric(level) || is.complex(level) || length(level) != 1L ||
      !is.finite(level) || level <= 0 || level >= 1)
    stop("level must be strictly between zero and one.", call. = FALSE)
  if (!is.null(formula) && (!is.null(response) || !is.null(covariates)))
    stop("Specify either formula or response/covariates, not both.", call. = FALSE)
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
  control <- catmech_marginal_cluster_control(spec$joint_control)
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
  if (!is.null(formula) || !is.null(response) || !is.null(covariates)) {
    value <- spec$value
    if (!is.character(value) || length(value) != 1L || is.na(value) || !nzchar(value) ||
        value == "(Intercept)" || value %in% predictors || !value %in% names(result$original) ||
        !is.numeric(result$original[[value]]) || is.complex(result$original[[value]]) ||
        !is.null(dim(result$original[[value]])))
      stop("The audited value column in spec$value is missing or malformed.", call. = FALSE)
    modeled <- c(value, predictors)
    target <- catmech_marginal_target(formula, response, covariates, value, modeled)
    if (!identical(target$response, value) || !identical(target$covariates, predictors)) {
      internal_map <- stats::setNames(c("T", if (length(predictors)) paste0("X", seq_along(predictors))), modeled)
      r_idx <- match(internal_map[[target$response]], internal)
      x_idx <- if (length(target$covariates)) match(unname(internal_map[target$covariates]), internal) else integer()
      terms <- c("(Intercept)", target$covariates)
      center <- unlist(transform$response_center[internal], use.names = FALSE)
      scale <- unlist(transform$response_scale[internal], use.names = FALSE)
      custom_flat <- matrix(NA_real_, nrow(flat), length(terms))
      if (!length(x_idx)) {
        custom_flat[, 1L] <- center[r_idx] + scale[r_idx] * flat[, paste0("mu:", internal[r_idx])]
      } else if (length(x_idx) == 1L) {
        r_lab <- internal[r_idx]
        x_lab <- internal[x_idx]
        pair_name <- if (r_idx >= x_idx) paste0(r_lab, ":", x_lab) else paste0(x_lab, ":", r_lab)
        custom_flat[, 2L] <- (flat[, paste0("sigma:", pair_name)] + flat[, paste0("psi:", pair_name)]) /
          (flat[, paste0("sigma:", x_lab, ":", x_lab)] + flat[, paste0("psi:", x_lab, ":", x_lab)]) *
          scale[r_idx] / scale[x_idx]
        custom_flat[, 1L] <- (center[r_idx] + scale[r_idx] * flat[, paste0("mu:", r_lab)]) -
          custom_flat[, 2L] * (center[x_idx] + scale[x_idx] * flat[, paste0("mu:", x_lab)])
      } else {
        total_columns <- paste0("sigma:", suffix)
        psi_columns <- paste0("psi:", suffix)
        mu_columns <- paste0("mu:", internal)
        for (i in seq_len(nrow(flat))) {
          z <- matrix(0, k, k)
          z[pairs] <- flat[i, total_columns] + flat[i, psi_columns]
          z <- z + t(z) - diag(diag(z), k)
          means_i <- center + scale * flat[i, mu_columns]
          beta_i <- as.vector(solve(z[x_idx, x_idx, drop = FALSE], z[x_idx, r_idx])) *
            scale[r_idx] / scale[x_idx]
          custom_flat[i, ] <- c(means_i[r_idx] - sum(beta_i * means_i[x_idx]), beta_i)
        }
      }
      coefficients <- array(custom_flat, c(dims[1L], dims[2L], length(terms)),
        dimnames = list(dimnames(tr)[[1L]], dimnames(tr)[[2L]], terms))
    }
  }
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

catmech_marginal_target <- function(formula, response, covariates, value, modeled) {
  if (!is.null(formula)) {
    if (!inherits(formula, "formula") || length(formula) != 3L)
      stop("formula must be a two-sided R formula.", call. = FALSE)
    lhs <- formula[[2L]]
    if (!is.symbol(lhs))
      stop("The formula response must be a single modeled variable name.", call. = FALSE)
    target_response <- as.character(lhs)
    if (!target_response %in% modeled)
      stop("The formula response must be the audited value or one of spec$predictors.", call. = FALSE)
    tt <- stats::terms(formula)
    if (!identical(attr(tt, "intercept"), 1L))
      stop("formula must include an intercept.", call. = FALSE)
    if (!is.null(attr(tt, "offset")))
      stop("formula must not include an offset.", call. = FALSE)
    vars <- as.list(attr(tt, "variables"))[-1L]
    if (!all(vapply(vars, is.symbol, logical(1L))))
      stop("formula must use plain variable names without transformations or calls.", call. = FALSE)
    var_names <- vapply(vars, as.character, character(1L))
    if (any(!var_names %in% modeled))
      stop("Every variable in formula must be the audited value or one of spec$predictors.", call. = FALSE)
    order <- attr(tt, "order")
    if (length(order) && any(order != 1L))
      stop("formula must be additive without interactions.", call. = FALSE)
    if (target_response %in% all.vars(formula[[3L]]))
      stop("The response variable cannot also appear as a covariate.", call. = FALSE)
    factors <- attr(tt, "factors")
    target_covariates <- if (!length(order)) character() else
      var_names[vapply(seq_len(ncol(factors)), function(j) which(factors[, j] != 0L), integer(1L))]
    if (any(!target_covariates %in% modeled) || anyDuplicated(target_covariates))
      stop("Every formula covariate must be a unique modeled variable.", call. = FALSE)
    return(list(response = target_response, covariates = target_covariates))
  }
  if (is.null(response)) {
    target_response <- value
  } else {
    if (!is.character(response) || length(response) != 1L || is.na(response) ||
        !response %in% modeled)
      stop("response must name the audited value or one of spec$predictors.", call. = FALSE)
    target_response <- response
  }
  allowed <- setdiff(modeled, target_response)
  if (is.null(covariates)) {
    target_covariates <- allowed
  } else {
    if (!is.character(covariates) || anyNA(covariates) || anyDuplicated(covariates) ||
        any(!covariates %in% allowed))
      stop("covariates must be unique names chosen from the remaining modeled variables.", call. = FALSE)
    target_covariates <- covariates
  }
  list(response = target_response, covariates = target_covariates)
}

#' Canonical prior sensitivity grid for joint marginal cluster models
#'
#' Automatically computes the prespecified 2-axis prior sensitivity protocol from
#' the Joint Marginal validation specification:
#' \enumerate{
#'   \item Covariance scale axis: \eqn{\text{prior\_scale} \in \{0.25, 1.0, 4.0\}}, holding regression scale fixed at 10.
#'   \item Regression scale axis: \eqn{\text{regression\_scale} \in \{5, 10, 20\}}, holding covariance scale fixed at 1.0.
#' }
#'
#' Evaluates the average absolute change in the posterior mean under the alternative
#' hyperparameters along each axis, and returns the primary sensitivity statistic
#' (the maximum of the two axis-level mean changes) along with Monte Carlo standard errors.
#'
#' @param fit A fitted \code{catmech_result} from \code{\link{catmech_reference}}
#'   with \code{cluster_model = "joint_marginal"}, or a named list of precomputed scenarios.
#' @param term Name of the target coefficient to assess. If \code{NULL},
#'   defaults to the first modeled predictor or formula covariate.
#' @param formula Optional formula passed to \code{\link{catmech_marginal_summary}} for custom outcome specifications.
#' @param prior_scales Numeric vector of covariance prior scales to evaluate (default \code{c(0.25, 1, 4)}).
#' @param regression_scales Numeric vector of regression prior scales to evaluate (default \code{c(5, 10, 20)}).
#' @param conf_level Credible level (default 0.95).
#' @param reference Reference value for direction checks (default 0).
#' @return An object of class \code{"catmech_marginal_sensitivity_grid"} containing:
#' \describe{
#'   \item{table}{Data frame of sensitivity results per scenario.}
#'   \item{axis_summary}{Summary comparing mean absolute shifts along each axis to MCSE bounds.}
#'   \item{statistic}{Primary sensitivity statistic: \eqn{\max(\text{mean\_abs\_change}_{\text{cov}}, \text{mean\_abs\_change}_{\text{reg}})}.}
#'   \item{mcse_mean}{Monte Carlo standard error of the posterior mean under the baseline prior.}
#'   \item{term}{The evaluated coefficient term.}
#'   \item{baseline}{Baseline scenario result row.}
#'   \item{status}{Sensitivity status from \code{\link{catmech_sensitivity}}.}
#'   \item{complete}{Logical indicating whether all evaluated scenarios succeeded.}
#'   \item{envelope}{Credible interval envelope across successful scenarios.}
#'   \item{raw_sensitivity}{Underlying output from \code{\link{catmech_sensitivity}}.}
#' }
#' @export
catmech_marginal_sensitivity_grid <- function(fit, term = NULL, formula = NULL,
                                              prior_scales = c(0.25, 1, 4),
                                              regression_scales = c(5, 10, 20),
                                              conf_level = 0.95, reference = 0) {
  if (!is.numeric(prior_scales) || !length(prior_scales) || anyNA(prior_scales) || any(prior_scales <= 0))
    stop("prior_scales must be a vector of positive numbers.", call. = FALSE)
  if (!is.numeric(regression_scales) || !length(regression_scales) || anyNA(regression_scales) || any(regression_scales <= 0))
    stop("regression_scales must be a vector of positive numbers.", call. = FALSE)
  if (!is.numeric(conf_level) || length(conf_level) != 1L || !is.finite(conf_level) ||
      conf_level <= 0 || conf_level >= 1)
    stop("conf_level must be between zero and one.", call. = FALSE)

  if (is.list(fit) && !inherits(fit, "catmech_result")) {
    results <- fit
    if (!length(results) || is.null(names(results)))
      stop("fit list must be a named list of catmech_result objects.", call. = FALSE)
    base_fit <- results[[1L]]
    if (!inherits(base_fit, "catmech_result"))
      stop("fit list elements must be catmech_result objects.", call. = FALSE)
    ctrl <- base_fit$spec$joint_control
    if (!is.list(ctrl)) ctrl <- list()
    base_ps <- if (!is.null(ctrl$prior_scale)) ctrl$prior_scale else 1
    base_rs <- if (!is.null(ctrl$regression_scale)) ctrl$regression_scale else 10
    cand_name <- sprintf("prior_%s_reg_%s", base_ps, base_rs)
    baseline_name <- if (cand_name %in% names(results)) cand_name else names(results)[1L]
  } else if (inherits(fit, "catmech_result")) {
    base_fit <- fit
    if (!identical(base_fit$spec$cluster_model, "joint_marginal"))
      stop("fit must be a result from cluster_model = 'joint_marginal'.", call. = FALSE)
    spec <- base_fit$spec
    orig <- base_fit$original
    m <- base_fit$manifest$m
    seed <- base_fit$manifest$seed
    ctrl <- spec$joint_control
    if (!is.list(ctrl)) ctrl <- list()
    base_ps <- if (!is.null(ctrl$prior_scale)) ctrl$prior_scale else 1
    base_rs <- if (!is.null(ctrl$regression_scale)) ctrl$regression_scale else 10
    baseline_name <- sprintf("prior_%s_reg_%s", base_ps, base_rs)

    results <- list()
    results[[baseline_name]] <- base_fit

    for (ps in prior_scales) {
      sc_name <- sprintf("prior_%s_reg_%s", ps, base_rs)
      if (identical(ps, base_ps)) {
        results[[sc_name]] <- base_fit
      } else if (!sc_name %in% names(results)) {
        sc_ctrl <- ctrl
        sc_ctrl$prior_scale <- ps
        sc_ctrl$regression_scale <- base_rs
        results[[sc_name]] <- tryCatch(
          catmech_reference(data = orig, value = spec$value, reference = spec$reference,
                            predictors = spec$predictors, method = if (!is.null(spec$method)) spec$method else "normal",
                            m = m, seed = seed, use_recorded = isTRUE(spec$use_recorded),
                            cluster = spec$cluster, cluster_model = "joint_marginal",
                            joint_control = sc_ctrl),
          error = identity
        )
      }
    }

    for (rs in regression_scales) {
      sc_name <- sprintf("prior_%s_reg_%s", base_ps, rs)
      if (identical(rs, base_rs)) {
        results[[sc_name]] <- base_fit
      } else if (!sc_name %in% names(results)) {
        sc_ctrl <- ctrl
        sc_ctrl$prior_scale <- base_ps
        sc_ctrl$regression_scale <- rs
        results[[sc_name]] <- tryCatch(
          catmech_reference(data = orig, value = spec$value, reference = spec$reference,
                            predictors = spec$predictors, method = if (!is.null(spec$method)) spec$method else "normal",
                            m = m, seed = seed, use_recorded = isTRUE(spec$use_recorded),
                            cluster = spec$cluster, cluster_model = "joint_marginal",
                            joint_control = sc_ctrl),
          error = identity
        )
      }
    }
  } else {
    stop("fit must be a catmech_result or a named list of catmech_results.", call. = FALSE)
  }

  target_term <- if (is.null(term)) {
    if (!is.null(formula)) {
      covs <- catmech_marginal_target(formula, NULL, NULL, base_fit$spec$value,
                                      c(base_fit$spec$value, base_fit$spec$predictors))$covariates
      if (length(covs)) covs[1L] else "(Intercept)"
    } else if (length(base_fit$spec$predictors)) {
      base_fit$spec$predictors[1L]
    } else {
      "(Intercept)"
    }
  } else {
    term
  }

  sens <- catmech_sensitivity(results = results, formula = formula, term = target_term,
                              reference = reference, conf_level = conf_level,
                              method = "marginal_summary")

  base_row <- sens$table[sens$table$scenario == baseline_name, ]
  base_est <- if (nrow(base_row) && base_row$status == "ok") base_row$estimate else NA_real_

  ctrl <- base_fit$spec$joint_control
  if (!is.list(ctrl)) ctrl <- list()
  base_ps <- if (!is.null(ctrl$prior_scale)) ctrl$prior_scale else 1
  base_rs <- if (!is.null(ctrl$regression_scale)) ctrl$regression_scale else 10

  cov_scenarios <- intersect(names(results), setdiff(sprintf("prior_%s_reg_%s", prior_scales, base_rs), baseline_name))
  cov_rows <- sens$table[sens$table$scenario %in% cov_scenarios & sens$table$status == "ok", ]
  cov_diffs <- if (is.finite(base_est) && nrow(cov_rows)) abs(cov_rows$estimate - base_est) else numeric(0)
  cov_mean_change <- if (length(cov_diffs)) mean(cov_diffs) else NA_real_
  cov_max_change <- if (length(cov_diffs)) max(cov_diffs) else NA_real_

  reg_scenarios <- intersect(names(results), setdiff(sprintf("prior_%s_reg_%s", base_ps, regression_scales), baseline_name))
  reg_rows <- sens$table[sens$table$scenario %in% reg_scenarios & sens$table$status == "ok", ]
  reg_diffs <- if (is.finite(base_est) && nrow(reg_rows)) abs(reg_rows$estimate - base_est) else numeric(0)
  reg_mean_change <- if (length(reg_diffs)) mean(reg_diffs) else NA_real_
  reg_max_change <- if (length(reg_diffs)) max(reg_diffs) else NA_real_

  statistic <- if (all(is.na(c(cov_mean_change, reg_mean_change)))) NA_real_ else
    max(c(cov_mean_change, reg_mean_change), na.rm = TRUE)

  base_summary <- tryCatch(
    catmech_marginal_summary(base_fit, level = conf_level, formula = formula),
    error = function(e) NULL
  )
  base_mcse <- if (!is.null(base_summary) && target_term %in% base_summary$term) {
    base_summary$mcse_mean[base_summary$term == target_term]
  } else {
    NA_real_
  }

  axis_summary <- data.frame(
    axis = c("covariance_scale", "regression_scale"),
    baseline = c(as.character(base_ps), as.character(base_rs)),
    alternatives = c(paste(setdiff(prior_scales, base_ps), collapse = ", "),
                     paste(setdiff(regression_scales, base_rs), collapse = ", ")),
    mean_abs_change = c(cov_mean_change, reg_mean_change),
    max_abs_change = c(cov_max_change, reg_max_change),
    ratio_to_mcse = if (is.finite(base_mcse) && base_mcse > 0) {
      c(cov_mean_change, reg_mean_change) / base_mcse
    } else {
      c(NA_real_, NA_real_)
    },
    stringsAsFactors = FALSE
  )

  structure(list(
    table = sens$table,
    axis_summary = axis_summary,
    statistic = statistic,
    mcse_mean = base_mcse,
    term = target_term,
    baseline = base_row,
    status = sens$status,
    complete = sens$complete,
    envelope = sens$successful_envelope,
    raw_sensitivity = sens
  ), class = c("catmech_marginal_sensitivity_grid", "list"))
}

#' Print a marginal sensitivity grid result
#'
#' @param x An object of class \code{catmech_marginal_sensitivity_grid}.
#' @param digits Number of decimal digits to print.
#' @param ... Unused.
#' @return The object \code{x}, invisibly.
#' @export
print.catmech_marginal_sensitivity_grid <- function(x, digits = 4L, ...) {
  cat("-- Catmech Joint Marginal Prior Sensitivity Grid --\n")
  b_est <- if (nrow(x$baseline) && is.finite(x$baseline$estimate)) round(x$baseline$estimate, digits) else "NA"
  b_se <- if (nrow(x$baseline) && is.finite(x$baseline$std.error)) round(x$baseline$std.error, digits) else "NA"
  b_mcse <- if (is.finite(x$mcse_mean)) round(x$mcse_mean, digits) else "NA"
  cat(sprintf("Target term: %s (baseline estimate: %s, posterior SD: %s, MCSE: %s)\n\n",
              x$term, b_est, b_se, b_mcse))

  cat("Axis summary:\n")
  df <- x$axis_summary
  df$mean_abs_change <- round(df$mean_abs_change, digits)
  df$max_abs_change <- round(df$max_abs_change, digits)
  df$ratio_to_mcse <- round(df$ratio_to_mcse, 2L)
  print(df, row.names = FALSE)
  cat("\n")

  stat_str <- if (is.finite(x$statistic)) round(x$statistic, digits) else "NA"
  mcse_ratio_str <- if (is.finite(x$statistic) && is.finite(x$mcse_mean) && x$mcse_mean > 0) {
    sprintf(" (%sx MCSE)", round(x$statistic / x$mcse_mean, 2L))
  } else ""
  cat(sprintf("Primary sensitivity statistic: %s%s\n", stat_str, mcse_ratio_str))

  env_str <- if (all(is.finite(x$envelope))) {
    sprintf("[%s, %s]", round(x$envelope[1L], digits), round(x$envelope[2L], digits))
  } else "[NA, NA]"
  cat(sprintf("Envelope: %s\n", env_str))
  cat(sprintf("Status: %s (complete: %s)\n", x$status, x$complete))

  if (!isTRUE(x$complete)) {
    cat("Note: One or more prior scenarios failed. Inspect $table for details.\n")
  }
  invisible(x)
}

#' Plot a marginal prior sensitivity grid
#'
#' Displays a 2-panel diagnostic plot showing posterior coefficient estimates and
#' credible intervals across covariance prior scales and regression prior scales,
#' compared against baseline estimates and Monte Carlo error bands.
#'
#' @param x A `catmech_marginal_sensitivity_grid` object.
#' @param y Unused.
#' @param col_points Color for scenario points (default `"navy"`).
#' @param col_base Color for baseline reference line (default `"firebrick"`).
#' @param ... Additional arguments passed to [graphics::par()].
#' @return Invisible `NULL`.
#' @export
plot.catmech_marginal_sensitivity_grid <- function(x, y, col_points = "navy",
                                                   col_base = "firebrick", ...) {
  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par), add = TRUE)
  graphics::par(mfrow = c(1, 2), mar = c(4.5, 4.5, 3.5, 1.5), ...)

  tbl <- x$table
  tbl_ok <- tbl[tbl$status == "ok", , drop = FALSE]
  b_est <- if (nrow(x$baseline) && is.finite(x$baseline$estimate)) x$baseline$estimate else NA_real_
  mcse <- if (is.finite(x$mcse_mean)) x$mcse_mean else 0
  bound <- 0.1 * mcse

  cov_names <- grep("^prior_[0-9.]+_reg_", tbl$scenario, value = TRUE)
  cov_scales <- as.numeric(gsub("^prior_([0-9.]+)_reg_.*", "\\1", cov_names))
  ord_cov <- order(cov_scales)
  cov_names <- cov_names[ord_cov]
  cov_scales <- cov_scales[ord_cov]
  cov_tbl <- tbl[match(cov_names, tbl$scenario), ]

  y_lim <- if (nrow(tbl_ok)) range(c(tbl_ok$lower, tbl_ok$upper, b_est + bound, b_est - bound), na.rm = TRUE) else c(-1, 1)

  graphics::plot(cov_scales, cov_tbl$estimate, type = "n", log = "x",
                 xlim = range(cov_scales, na.rm = TRUE), ylim = y_lim,
                 xlab = expression(paste("Covariance Prior Scale (", A^2, ")")),
                 ylab = paste("Estimate:", x$term),
                 main = "Covariance Prior Scale", xaxt = "n")
  graphics::axis(1, at = cov_scales, labels = cov_scales)
  if (is.finite(b_est)) {
    if (bound > 0) {
      graphics::rect(min(cov_scales) * 0.5, b_est - bound, max(cov_scales) * 2, b_est + bound,
                     col = grDevices::adjustcolor(col_base, alpha.f = 0.15), border = NA)
    }
    graphics::abline(h = b_est, col = col_base, lty = 2, lwd = 1.5)
  }
  for (i in seq_along(cov_scales)) {
    if (cov_tbl$status[i] == "ok") {
      graphics::arrows(cov_scales[i], cov_tbl$lower[i], cov_scales[i], cov_tbl$upper[i],
                       length = 0.04, angle = 90, code = 3, col = col_points, lwd = 1.2)
    }
  }
  graphics::points(cov_scales[cov_tbl$status == "ok"], cov_tbl$estimate[cov_tbl$status == "ok"],
                   pch = 19, col = col_points, cex = 1.2)

  reg_names <- grep("_reg_[0-9.]+$", tbl$scenario, value = TRUE)
  reg_scales <- as.numeric(gsub(".*_reg_([0-9.]+)$", "\\1", reg_names))
  ord_reg <- order(reg_scales)
  reg_names <- reg_names[ord_reg]
  reg_scales <- reg_scales[ord_reg]
  reg_tbl <- tbl[match(reg_names, tbl$scenario), ]

  graphics::plot(reg_scales, reg_tbl$estimate, type = "n", log = "x",
                 xlim = range(reg_scales, na.rm = TRUE), ylim = y_lim,
                 xlab = expression(paste("Regression Prior Scale (", lambda, ")")),
                 ylab = "", main = "Regression Prior Scale", xaxt = "n")
  graphics::axis(1, at = reg_scales, labels = reg_scales)
  if (is.finite(b_est)) {
    if (bound > 0) {
      graphics::rect(min(reg_scales) * 0.5, b_est - bound, max(reg_scales) * 2, b_est + bound,
                     col = grDevices::adjustcolor(col_base, alpha.f = 0.15), border = NA)
    }
    graphics::abline(h = b_est, col = col_base, lty = 2, lwd = 1.5)
  }
  for (i in seq_along(reg_scales)) {
    if (reg_tbl$status[i] == "ok") {
      graphics::arrows(reg_scales[i], reg_tbl$lower[i], reg_scales[i], reg_tbl$upper[i],
                       length = 0.04, angle = 90, code = 3, col = col_points, lwd = 1.2)
    }
  }
  graphics::points(reg_scales[reg_tbl$status == "ok"], reg_tbl$estimate[reg_tbl$status == "ok"],
                   pch = 19, col = col_points, cex = 1.2)

  invisible(NULL)
}
