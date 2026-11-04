# Multivariate Gaussian random-intercept model for (T,W,X), including the
# complete numeric predictors as responses. Unlike conditional regression,
# the marginal coefficient uses both between-site and within-site covariance.
catmech_marginal_cluster_prior <- function(control) {
  list(mean = "flat", covariance_family = "Huang_Wand",
    covariance_df = "dimension + 1", auxiliary_shape = .5,
    auxiliary_inverse_gamma_rate = 1 / control$prior_scale,
    covariance_conditional_scale = "4 * diag(1 / auxiliary)",
    marginal_sd_family = "half_t", marginal_sd_df = 2,
    marginal_sd_scale = sqrt(control$prior_scale),
    marginal_correlation = "uniform(-1,1)",
    covariance_components = "independent priors for residual and site covariance",
    scale_convention = "standardized response units; prior_scale multiplies squared half-t scale",
    scaling = "empirical observed-response means and standard deviations",
    dependence = "standard deviations and correlations are not jointly independent")
}

# Internal reproducible kernels; fixed covariance and inverse-Wishart options
# exist for independent numerical checks, not as user-selectable model controls.
catmech_marginal_gibbs <- function(input, seed, iterations, start = NULL,
                                  prior = list(family = "huang_wand", scale = 1), fixed = NULL) {
  check_count(iterations, "iterations", 1L)
  check_count(seed, "seed", 1L)
  if (!is.list(input) || !is.matrix(input$y) || !is.numeric(input$y) ||
      is.complex(input$y) || nrow(input$y) < 1L || ncol(input$y) < 1L ||
      any(!is.finite(input$y[!is.na(input$y)])) ||
      !is.numeric(input$subj) || length(input$subj) != nrow(input$y) ||
      anyNA(input$subj) || any(input$subj < 1 | input$subj != as.integer(input$subj)) ||
      !setequal(input$subj, seq_len(max(input$subj))))
    stop("Malformed marginal sampler input.", call. = FALSE)
  if (!is.list(prior) || !is.character(prior$family) || length(prior$family) != 1L ||
      !prior$family %in% c("huang_wand", "inverse_wishart") ||
      !is.numeric(prior$scale) || is.complex(prior$scale) || length(prior$scale) != 1L ||
      !is.finite(prior$scale) || prior$scale <= 0)
    stop("Malformed marginal covariance prior.", call. = FALSE)
  if (is.null(prior$interweave)) prior$interweave <- TRUE
  if (!is.logical(prior$interweave) || length(prior$interweave) != 1L || is.na(prior$interweave))
    stop("The internal interweave flag must be TRUE or FALSE.", call. = FALSE)
  if (!is.null(fixed)) {
    if (!is.list(fixed) || !all(c("sigma", "psi") %in% names(fixed)))
      stop("Both fixed covariance matrices are required.", call. = FALSE)
    for (name in c("sigma", "psi"))
      catmech_marginal_covariance(fixed[[name]], ncol(input$y))
  }
  restore <- preserve_rng()
  on.exit(restore(), add = TRUE)
  set.seed(seed)
  .Call("catmech_marginal_gibbs_cpp", list(y = input$y, subj = as.integer(input$subj)),
    as.integer(iterations), start, prior, fixed, PACKAGE = "catmech")
}

# Exact native scale-move ratio, exposed internally for a full-density/Jacobian
# oracle. The proposal acts on completed data, so no missing entries are allowed.
catmech_marginal_interweave_ratio <- function(y, subj, mu, u, sigma, psi, j, eta,
                                             prior = list(family = "huang_wand", scale = 1)) {
  if (!is.matrix(y) || !is.numeric(y) || is.complex(y) || any(!is.finite(y)) ||
      !is.matrix(u) || !is.numeric(u) || is.complex(u) || any(!is.finite(u)) ||
      ncol(y) != ncol(u) || !is.numeric(subj) || length(subj) != nrow(y) ||
      any(!is.finite(subj)) || any(subj != as.integer(subj)) ||
      !setequal(subj, seq_len(nrow(u))) || !is.numeric(mu) || is.complex(mu) ||
      length(mu) != ncol(y) || any(!is.finite(mu)) ||
      !is.numeric(j) || length(j) != 1L || !is.finite(j) || j != as.integer(j) ||
      j < 1L || j > ncol(y) || !is.numeric(eta) || is.complex(eta) ||
      length(eta) != 1L || !is.finite(eta) || !is.list(prior) ||
      !is.character(prior$family) || length(prior$family) != 1L ||
      !prior$family %in% c("huang_wand", "inverse_wishart") ||
      !is.numeric(prior$scale) || length(prior$scale) != 1L ||
      !is.finite(prior$scale) || prior$scale <= 0)
    stop("Malformed marginal scale-move oracle input.", call. = FALSE)
  catmech_marginal_covariance(sigma, ncol(y))
  catmech_marginal_covariance(psi, ncol(y))
  .Call("catmech_marginal_ratio_cpp", y, as.integer(subj), mu, u, sigma, psi,
    as.integer(j), eta, prior, PACKAGE = "catmech")
}

catmech_marginal_hw_draw <- function(n, dimension, scale = 1, seed = 1L) {
  check_count(n, "n", 1L)
  check_count(dimension, "dimension", 1L)
  check_count(seed, "seed", 1L)
  if (!is.numeric(scale) || is.complex(scale) || length(scale) != 1L ||
      !is.finite(scale) || scale <= 0)
    stop("scale must be a positive finite number.", call. = FALSE)
  restore <- preserve_rng()
  on.exit(restore(), add = TRUE)
  set.seed(seed)
  .Call("catmech_marginal_hw_cpp", as.integer(n), as.integer(dimension), scale,
    PACKAGE = "catmech")
}

catmech_marginal_covariance <- function(x, dimension) {
  if (!is.matrix(x) || !is.numeric(x) || is.complex(x) ||
      !identical(dim(x), rep(as.integer(dimension), 2L)) || any(!is.finite(x)) ||
      !isTRUE(all.equal(x, t(x), tolerance = 1e-12, check.attributes = FALSE)))
    stop("Marginal covariance must be a finite symmetric matrix of the correct size.", call. = FALSE)
  if (inherits(try(chol(x), silent = TRUE), "try-error"))
    stop("Marginal covariance must be positive definite.", call. = FALSE)
  invisible(x)
}

catmech_marginal_coefficients <- function(mu, sigma, psi, transform, predictors) {
  d <- length(predictors) + 2L
  catmech_marginal_covariance(sigma, d)
  catmech_marginal_covariance(psi, d)
  if (!is.numeric(mu) || is.complex(mu) || length(mu) != d || any(!is.finite(mu)))
    stop("Marginal means must be finite and match the response dimension.", call. = FALSE)
  labels <- c("T", "W", if (length(predictors)) paste0("X", seq_along(predictors)) else character())
  center <- unlist(transform$response_center[labels], use.names = FALSE)
  scale <- unlist(transform$response_scale[labels], use.names = FALSE)
  if (length(center) != d || length(scale) != d || any(!is.finite(c(center, scale))) || any(scale <= 0))
    stop("Marginal response transformations are malformed.", call. = FALSE)
  means <- center + scale * mu
  if (!length(predictors)) return(stats::setNames(means[1L], "(Intercept)"))
  total <- sigma + psi
  x <- seq.int(3L, d)
  beta <- as.vector(solve(total[x, x, drop = FALSE], total[x, 1L])) * scale[1L] / scale[x]
  stats::setNames(c(means[1L] - sum(beta * means[x]), beta), c("(Intercept)", predictors))
}

# Conditional target draws use observed components only, integrating and then
# sampling one shared site effect before drawing each missing T residual.
catmech_marginal_cluster_predict <- function(y, cluster, mu, sigma, psi, draw = TRUE) {
  if (!is.matrix(y) || !is.numeric(y) || is.complex(y) || ncol(y) < 2L || nrow(y) < 1L ||
      any(!is.finite(y[!is.na(y)])) || length(cluster) != nrow(y) || anyNA(cluster) ||
      !is.numeric(mu) || is.complex(mu) || length(mu) != ncol(y) || any(!is.finite(mu)) ||
      !is.logical(draw) || length(draw) != 1L || is.na(draw))
    stop("Malformed marginal prediction problem.", call. = FALSE)
  catmech_marginal_covariance(sigma, ncol(y))
  catmech_marginal_covariance(psi, ncol(y))
  group <- match(as.character(cluster), unique(as.character(cluster))) - 1L
  as.vector(.Call("catmech_marginal_predict_cpp", y, as.integer(group), mu, sigma, psi,
    draw, PACKAGE = "catmech"))
}

catmech_marginal_cluster_prepare <- function(data, value, gold, predictors, cluster) {
  if (!is.character(predictors) || anyNA(predictors) || anyDuplicated(predictors) ||
      any(!nzchar(predictors)) || any(predictors == "(Intercept)") ||
      any(!predictors %in% names(data)) ||
      any(!vapply(data[predictors], function(z) is.numeric(z) && is.null(dim(z)) && !is.complex(z) &&
        all(is.finite(z)), logical(1))))
    stop("joint_marginal requires complete finite numeric predictors with unique names other than '(Intercept)'.", call. = FALSE)
  x <- if (length(predictors)) as.matrix(data[predictors]) else matrix(numeric(), nrow(data), 0L)
  design <- cbind(`(Intercept)` = 1, x)
  catmech_model_design(design)
  y <- cbind(T = gold, W = data[[value]], x)
  labels <- c("T", "W", if (length(predictors)) paste0("X", seq_along(predictors)) else character())
  colnames(y) <- labels
  ids <- as.character(data[[cluster]])
  if (length(predictors)) {
    centered_x <- x
    for (j in seq_len(ncol(x)))
      centered_x[, j] <- x[, j] - ave(x[, j], ids, FUN = mean)
    if (qr(centered_x, tol = 1e-7)$rank < ncol(x))
      stop("joint_marginal requires independent within-site predictor variation; site-constant or within-site aliased predictors are unsupported.", call. = FALSE)
  }
  observed <- !is.na(y)
  paired <- observed[, 1L] & observed[, 2L]
  for (j in 1:2) {
    rows <- observed[, j]
    catmech_model_design(design[rows, , drop = FALSE])
    if (length(unique(ids[rows])) < 5L)
      stop("Marginal cluster fitting requires at least five observed sites for each measurement response.", call. = FALSE)
    site_design <- stats::model.matrix(~ factor(ids[rows]) - 1)
    rank <- qr(cbind(design[rows, , drop = FALSE], site_design), tol = 1e-7)$rank
    if (rank <= ncol(design) || sum(rows) <= rank)
      stop("Marginal cluster fitting needs within-site replication and between-site residual information for both measurements.", call. = FALSE)
  }
  if (sum(paired) < ncol(design) + 2L)
    stop("Marginal cluster fitting needs paired calibration rows with at least two residual degrees of freedom.", call. = FALSE)
  catmech_model_design(design[paired, , drop = FALSE])
  center <- colMeans(y, na.rm = TRUE)
  scale <- apply(y, 2L, stats::sd, na.rm = TRUE)
  if (any(!is.finite(scale)) || any(scale <= 0))
    stop("All marginal-model responses and predictors must have positive finite variation.", call. = FALSE)
  y <- sweep(sweep(y, 2L, center, "-"), 2L, scale, "/")
  rows <- which(rowSums(observed) > 0L)
  group <- match(ids[rows], unique(ids[rows]))
  list(y = y, cluster = ids, rows = rows, predictors = predictors, labels = labels,
    center = center, scale = scale,
    input = list(y = y[rows, , drop = FALSE], subj = as.integer(group)),
    prediction_input = list(y = y, subj = match(ids, unique(ids))),
    transform = list(response_center = as.list(center), response_scale = as.list(scale),
      predictor_names = as.list(predictors)))
}

catmech_marginal_cluster_trace <- function(problem, fit, keep) {
  d <- length(problem$labels)
  index <- which(lower.tri(matrix(0, d, d), diag = TRUE), arr.ind = TRUE)
  covariance_names <- paste(problem$labels[index[, 1L]], problem$labels[index[, 2L]], sep = ":")
  parts <- list(t(fit$mu[, keep, drop = FALSE]))
  names <- paste0("mu:", problem$labels)
  for (component in c("sigma", "psi")) {
    values <- vapply(seq_len(nrow(index)), function(k)
      fit[[component]][index[k, 1L], index[k, 2L], keep], numeric(length(keep)))
    parts[[length(parts) + 1L]] <- values
    names <- c(names, paste0(component, ":", covariance_names))
  }
  for (component in c("sigma", "psi")) {
    parts[[length(parts) + 1L]] <- t(fit[[paste0("aux_", component)]][, keep, drop = FALSE])
    names <- c(names, paste0("aux_", component, ":", problem$labels))
  }
  derived <- .Call("catmech_marginal_derived_cpp", problem$prediction_input, fit,
    as.integer(keep), unname(problem$center), unname(problem$scale), PACKAGE = "catmech")
  derived_names <- paste0("marginal:", c("(Intercept)", problem$predictors))
  missing <- is.na(problem$y[, 1L])
  if (any(missing & !is.na(problem$y[, 2L])))
    derived_names <- c(derived_names, "conditional_target_mean_recorded")
  if (any(missing & is.na(problem$y[, 2L])))
    derived_names <- c(derived_names, "conditional_target_mean_absent")
  parts[[length(parts) + 1L]] <- derived
  out <- do.call(cbind, parts)
  colnames(out) <- c(names, derived_names)
  if (any(!is.finite(out))) stop("Nonfinite marginal cluster parameter trace.", call. = FALSE)
  out
}

catmech_marginal_cluster_draws <- function(data, value, gold, predictors, m, cluster, control) {
  if (m > as.double(control$iterations) * control$chains)
    stop("m cannot exceed joint_control$iterations * joint_control$chains.", call. = FALSE)
  problem <- catmech_marginal_cluster_prepare(data, value, gold, predictors, cluster)
  d <- ncol(problem$y)
  seeds <- sample.int(.Machine$integer.max, control$chains)
  starts <- lapply(seq_len(control$chains), function(chain) {
    y <- problem$input$y
    y[is.na(y)] <- stats::rnorm(sum(is.na(y)))
    factor <- 2^seq(-1, 1, length.out = control$chains)[chain]
    list(mu = stats::rnorm(d, 0, .5), sigma = diag(factor, d),
      psi = diag(1 / factor, d), aux_sigma = rep(factor, d),
      aux_psi = rep(1 / factor, d), y = y)
  })
  prior <- list(family = "huang_wand", scale = control$prior_scale, interweave = TRUE)
  iterations <- control$iterations
  attempts <- list()
  repeat {
    chains <- lapply(seq_len(control$chains), function(chain)
      catmech_marginal_gibbs(problem$input, seeds[chain], control$burnin + iterations,
        starts[[chain]], prior))
    retained <- control$burnin + seq_len(iterations)
    parts <- lapply(chains, function(fit) catmech_marginal_cluster_trace(problem, fit, retained))
    traces <- array(NA_real_, c(iterations, control$chains, ncol(parts[[1L]])),
      dimnames = list(iteration = as.character(seq_len(iterations)),
        chain = as.character(seq_len(control$chains)), parameter = colnames(parts[[1L]])))
    for (chain in seq_len(control$chains)) traces[, chain, ] <- parts[[chain]]
    diagnostics <- catmech_model_chain_diagnostics(traces)
    diagnostics$engine <- "hw_marginal_cluster"
    attempts[[length(attempts) + 1L]] <- data.frame(iterations = iterations,
      convergence_pass = diagnostics$convergence_pass,
      max_rhat = max(diagnostics$summary$rhat), min_bulk_ess = min(diagnostics$summary$ess_bulk),
      min_tail_ess = min(diagnostics$summary$ess_tail))
    if (diagnostics$convergence_pass || iterations >= control$max_iterations) break
    iterations <- as.integer(min(2 * as.double(iterations), control$max_iterations))
  }
  diagnostics$attempts <- do.call(rbind, attempts)
  diagnostics$iterations_used <- iterations
  diagnostics$chain_seeds <- seeds
  diagnostics$chain_starts <- starts
  diagnostics$sampler_input <- problem$input
  diagnostics$sampler_prior <- prior
  diagnostics$input_rows <- problem$rows
  diagnostics$transform <- problem$transform
  diagnostics$backend <- "native_blocked_Gibbs"
  diagnostics$backend_isolation <- "R_seed_restored_per_chain_and_exact_restart"
  diagnostics$prior <- catmech_marginal_cluster_prior(control)
  diagnostics$scale_interweaving <- list(proposals_per_component = 2L, log_scale_sd = 1,
    auxiliary_update = "redrawn_after_collapsed_scale_moves",
    acceptance = do.call(rbind, lapply(seq_len(control$chains), function(chain)
      data.frame(chain = chain, component = problem$labels,
        proposed = as.vector(chains[[chain]]$scale_proposed),
        accepted = as.vector(chains[[chain]]$scale_accepted),
        acceptance_rate = as.vector(chains[[chain]]$scale_accepted / chains[[chain]]$scale_proposed)))))
  if (!diagnostics$convergence_pass) stop(structure(list(
    message = paste0("Marginal cluster model failed convergence checks (R-hat <= 1.05 and bulk/tail ESS >= ",
      100L * control$chains, " required). Inspect the condition's model_diagnostics."),
    call = NULL, model_diagnostics = diagnostics),
    class = c("catmech_convergence_error", "error", "condition")))
  allocation <- rep(m %/% control$chains, control$chains)
  if (m %% control$chains) allocation[seq_len(m %% control$chains)] <- allocation[seq_len(m %% control$chains)] + 1L
  draws <- matrix(NA_real_, nrow(data), m)
  selected <- list(); next_column <- 0L
  for (chain in seq_len(control$chains)) {
    keep <- if (!allocation[chain]) integer() else if (allocation[chain] == 1L)
      as.integer(ceiling(iterations / 2)) else as.integer(round(seq(1, iterations, length.out = allocation[chain])))
    for (post in keep) {
      next_column <- next_column + 1L
      k <- control$burnin + post
      target <- catmech_marginal_cluster_predict(problem$y, problem$cluster,
        chains[[chain]]$mu[, k], chains[[chain]]$sigma[, , k], chains[[chain]]$psi[, , k])
      draws[, next_column] <- target * problem$scale[1L] + problem$center[1L]
      draws[!is.na(gold), next_column] <- gold[!is.na(gold)]
      selected[[next_column]] <- data.frame(imputation = next_column, chain = chain, retained_iteration = post)
    }
  }
  diagnostics$selected_draws <- do.call(rbind, selected)
  list(draws = draws, model_diagnostics = diagnostics, transform = problem$transform)
}
