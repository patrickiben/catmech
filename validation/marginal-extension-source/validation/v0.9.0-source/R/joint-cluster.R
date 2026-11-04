# Explicit bivariate Gaussian random-intercept model, fitted by pan.
catmech_joint_cluster_control <- function(control = list()) {
  defaults <- list(chains = 4L, burnin = 1000L, iterations = 4000L,
                   max_iterations = 32000L, prior_scale = 1)
  if (!is.list(control) || (length(control) && (is.null(names(control)) ||
      anyNA(names(control)) || any(!nzchar(names(control))) || anyDuplicated(names(control)))))
    stop("joint_control must be a uniquely named list.", call. = FALSE)
  if (length(setdiff(names(control), names(defaults))))
    stop("Unknown joint_control fields: ", paste(setdiff(names(control), names(defaults)), collapse = ", "), call. = FALSE)
  defaults[names(control)] <- control
  for (name in c("chains", "burnin", "iterations", "max_iterations")) {
    minimum <- switch(name, chains = 4L, burnin = 0L, iterations = 100L, max_iterations = 100L)
    check_count(defaults[[name]], paste0("joint_control$", name), minimum)
    defaults[[name]] <- as.integer(defaults[[name]])
  }
  if (defaults$max_iterations < defaults$iterations ||
      as.double(defaults$burnin) + defaults$max_iterations > .Machine$integer.max)
    stop("joint_control requires max_iterations >= iterations and a supported total iteration count.", call. = FALSE)
  if (!is.numeric(defaults$prior_scale) || length(defaults$prior_scale) != 1L ||
      !is.finite(defaults$prior_scale) || defaults$prior_scale <= 0)
    stop("joint_control$prior_scale must be a positive finite number.", call. = FALSE)
  defaults
}

catmech_joint_cluster_prior <- function(control) {
  list(beta = "flat", covariance_family = "inverse_Wishart", covariance_df = 4L,
    covariance_scale = list(T = list(T = .5 * control$prior_scale, W = 0),
                            W = list(T = 0, W = .5 * control$prior_scale)),
    covariance_prior_mean = .5 * control$prior_scale,
    scale_convention = "standardized response units; same prior for residual and site covariance",
    scaling = "empirical observed-response means and standard deviations")
}

# Prepare the original observed components for repeated conditional calculations.
catmech_joint_cluster_prediction_problem <- function(y, x, cluster) {
  ids <- unique(as.character(cluster))
  group <- match(as.character(cluster), ids)
  g <- length(ids)
  present <- !is.na(y)
  both <- present[, 1L] & present[, 2L]
  only_t <- present[, 1L] & !present[, 2L]
  only_w <- !present[, 1L] & present[, 2L]
  grouped <- function(z, rows) {
    z <- as.matrix(z)
    out <- matrix(0, g, ncol(z))
    if (any(rows)) {
      sums <- rowsum(z[rows, , drop = FALSE], group[rows], reorder = TRUE)
      out[as.integer(rownames(sums)), ] <- sums
    }
    out
  }
  patterns <- lapply(list(both, only_t, only_w), function(rows) list(
    n = tabulate(group[rows], nbins = g), x = grouped(x, rows),
    y = grouped(replace(y, is.na(y), 0), rows)))
  list(y = y, x = x, group = group, ids = ids, patterns = patterns,
       missing = !present[, 1L], recorded = present[, 2L])
}

catmech_joint_cluster_condition <- function(problem, beta, sigma, psi) {
  ri <- chol2inv(chol(sigma))
  ui <- chol2inv(chol(psi))
  a <- problem$patterns[[1L]]
  t <- problem$patterns[[2L]]
  w <- problem$patterns[[3L]]
  h <- (a$y - a$x %*% beta) %*% ri
  h[, 1L] <- h[, 1L] + (t$y[, 1L] - as.vector(t$x %*% beta[, 1L])) / sigma[1L, 1L]
  h[, 2L] <- h[, 2L] + (w$y[, 2L] - as.vector(w$x %*% beta[, 2L])) / sigma[2L, 2L]
  p11 <- ui[1L, 1L] + a$n * ri[1L, 1L] + t$n / sigma[1L, 1L]
  p12 <- ui[1L, 2L] + a$n * ri[1L, 2L]
  p22 <- ui[2L, 2L] + a$n * ri[2L, 2L] + w$n / sigma[2L, 2L]
  # Invert through the Schur complement rather than subtracting two products.
  s22 <- 1 / (p22 - p12 * (p12 / p11))
  s12 <- -(p12 / p11) * s22
  s11 <- 1 / p11 + (p12 / p11)^2 * s22
  if (any(!is.finite(c(s11, s12, s22))) || any(s11 <= 0 | s22 <= 0))
    stop("Joint cluster conditional covariance is not numerically positive definite.", call. = FALSE)
  mean <- cbind(s11 * h[, 1L] + s12 * h[, 2L], s12 * h[, 1L] + s22 * h[, 2L])
  list(mean = mean, s11 = s11, s12 = s12, s22 = s22,
       row_variance = 1 / ri[1L, 1L], lambda = sigma[1L, 2L] / sigma[2L, 2L])
}

# Fixed-parameter predictive law. No imputed value enters the site posterior.
catmech_joint_cluster_predict <- function(y, x, cluster, beta, sigma, psi) {
  problem <- catmech_joint_cluster_prediction_problem(y, x, cluster)
  out <- as.numeric(y[, 1L])
  if (!any(problem$missing)) return(out)
  conditional <- catmech_joint_cluster_condition(problem, beta, sigma, psi)
  g <- length(problem$ids)
  z <- matrix(stats::rnorm(2L * g), g, 2L)
  u <- conditional$mean
  u[, 1L] <- u[, 1L] + sqrt(conditional$s11) * z[, 1L]
  remainder <- conditional$s22 - conditional$s12^2 / conditional$s11
  if (any(!is.finite(remainder)) || any(remainder <= 0))
    stop("Joint cluster conditional draw is numerically singular.", call. = FALSE)
  u[, 2L] <- u[, 2L] + conditional$s12 / sqrt(conditional$s11) * z[, 1L] + sqrt(remainder) * z[, 2L]
  location <- x %*% beta + u[problem$group, , drop = FALSE]
  recorded <- problem$missing & problem$recorded
  absent <- problem$missing & !problem$recorded
  out[recorded] <- stats::rnorm(sum(recorded),
    location[recorded, 1L] + conditional$lambda * (y[recorded, 2L] - location[recorded, 2L]),
    sqrt(conditional$row_variance))
  out[absent] <- stats::rnorm(sum(absent), location[absent, 1L], sqrt(sigma[1L, 1L]))
  if (any(!is.finite(out))) stop("Nonfinite joint cluster prediction.", call. = FALSE)
  out
}

catmech_joint_cluster_prepare <- function(data, value, gold, predictors, cluster) {
  y <- cbind(T = gold, W = data[[value]])
  x <- catmech_model_matrix(data, predictors)
  ids <- as.character(data[[cluster]])
  observed <- !is.na(y)
  paired <- observed[, 1L] & observed[, 2L]
  for (j in 1:2) {
    rows <- observed[, j]
    catmech_model_design(x[rows, , drop = FALSE])
    if (length(unique(ids[rows])) < 5L)
      stop("Joint cluster fitting requires at least five observed sites for each response, including five gold sites.", call. = FALSE)
    d <- stats::model.matrix(~ factor(ids[rows]) - 1)
    rank <- qr(cbind(x[rows, , drop = FALSE], d), tol = 1e-7)$rank
    if (rank <= ncol(x) || sum(rows) <= rank)
      stop("Joint cluster fitting needs within-site replication and between-site residual information for both observed responses.", call. = FALSE)
  }
  if (sum(paired) < ncol(x) + 2L)
    stop("Joint cluster fitting needs paired calibration rows with at least two residual degrees of freedom.", call. = FALSE)
  catmech_model_design(x[paired, , drop = FALSE])
  center <- colMeans(y, na.rm = TRUE)
  scale <- apply(y, 2L, stats::sd, na.rm = TRUE)
  if (any(!is.finite(scale)) || any(scale <= 0))
    stop("Both observed responses must have positive finite variation for joint scaling.", call. = FALSE)
  x_center <- rep(0, ncol(x)); x_scale <- rep(1, ncol(x))
  if (ncol(x) > 1L) {
    x_center[-1L] <- colMeans(x[, -1L, drop = FALSE])
    x_scale[-1L] <- apply(x[, -1L, drop = FALSE], 2L, stats::sd)
  }
  names(x_center) <- names(x_scale) <- colnames(x)
  y <- sweep(sweep(y, 2L, center, "-"), 2L, scale, "/")
  x <- sweep(sweep(x, 2L, x_center, "-"), 2L, x_scale, "/")
  active <- rowSums(observed) > 0L
  rows <- which(active)
  group <- match(ids[rows], sort(unique(ids[rows])))
  ordering <- order(group)
  rows <- rows[ordering]
  list(y = y, x = x, cluster = ids, rows = rows,
    input = list(y = y[rows, , drop = FALSE], pred = x[rows, , drop = FALSE], subj = group[ordering]),
    transform = list(response_center = as.list(center), response_scale = as.list(scale),
      predictor_center = as.list(x_center), predictor_scale = as.list(x_scale)),
    center = center, scale = scale,
    prediction = catmech_joint_cluster_prediction_problem(y, x, ids))
}

catmech_joint_cluster_trace <- function(problem, fit, keep) {
  p <- ncol(problem$x)
  beta_names <- c(paste0("beta_T:", colnames(problem$x)), paste0("beta_W:", colnames(problem$x)))
  names <- c(beta_names, "sigma_TT", "sigma_TW", "sigma_WW", "psi_TT", "psi_TW", "psi_WW")
  patterns <- list(recorded = problem$prediction$missing & problem$prediction$recorded,
                   absent = problem$prediction$missing & !problem$prediction$recorded)
  patterns <- patterns[vapply(patterns, any, logical(1))]
  names <- c(names, paste0("conditional_target_mean_", names(patterns)),
             "beta_T_original:(Intercept)")
  out <- matrix(NA_real_, length(keep), length(names), dimnames = list(NULL, names))
  out[, seq_len(2L * p)] <- matrix(aperm(fit$beta[, , keep, drop = FALSE], c(3L, 1L, 2L)),
                                  nrow = length(keep))
  out[, 2L * p + seq_len(6L)] <- cbind(fit$sigma[1L, 1L, keep], fit$sigma[1L, 2L, keep],
    fit$sigma[2L, 2L, keep], fit$psi[1L, 1L, keep], fit$psi[1L, 2L, keep], fit$psi[2L, 2L, keep])

  # A pattern mean depends on recipient design means and site counts, so it
  # does not require constructing a prediction for every row at every draw.
  g <- length(problem$prediction$ids)
  recipients <- lapply(patterns, function(rows) list(
    x = colMeans(problem$x[rows, , drop = FALSE]),
    weight = tabulate(problem$prediction$group[rows], nbins = g) / sum(rows),
    w = if (all(problem$prediction$recorded[rows])) mean(problem$y[rows, 2L]) else NA_real_))
  a <- problem$prediction$patterns[[1L]]
  t <- problem$prediction$patterns[[2L]]
  w <- problem$prediction$patterns[[3L]]
  # Bound temporary site-by-draw matrices even for many sites or long chains.
  batch_size <- max(1L, min(1024L, as.integer(floor(1e6 / max(g, 1L)))))
  for (first in seq.int(1L, length(keep), by = batch_size)) {
    rows <- seq.int(first, min(first + batch_size - 1L, length(keep)))
    beta_t <- t(out[rows, seq_len(p), drop = FALSE])
    beta_w <- t(out[rows, p + seq_len(p), drop = FALSE])
    sigma11 <- out[rows, 2L * p + 1L]; sigma12 <- out[rows, 2L * p + 2L]
    sigma22 <- out[rows, 2L * p + 3L]
    psi11 <- out[rows, 2L * p + 4L]; psi12 <- out[rows, 2L * p + 5L]
    psi22 <- out[rows, 2L * p + 6L]
    # The same two-dimensional precision algebra as the predictive helper,
    # evaluated across draws; Schur complements avoid determinant subtraction.
    ri11 <- 1 / (sigma11 - sigma12 * (sigma12 / sigma22))
    ri12 <- -(sigma12 / sigma22) * ri11
    ri22 <- 1 / sigma22 + (sigma12 / sigma22)^2 * ri11
    ui11 <- 1 / (psi11 - psi12 * (psi12 / psi22))
    ui12 <- -(psi12 / psi22) * ui11
    ui22 <- 1 / psi22 + (psi12 / psi22)^2 * ui11
    if (any(!is.finite(c(ri11, ri12, ri22, ui11, ui12, ui22))) ||
        any(c(sigma22, psi22, ri11, ui11) <= 0))
      stop("Joint cluster trace covariance is not numerically positive definite.", call. = FALSE)
    da_t <- a$y[, 1L] - a$x %*% beta_t
    da_w <- a$y[, 2L] - a$x %*% beta_w
    h_t <- sweep(da_t, 2L, ri11, "*") + sweep(da_w, 2L, ri12, "*") +
      sweep(t$y[, 1L] - t$x %*% beta_t, 2L, sigma11, "/")
    h_w <- sweep(da_t, 2L, ri12, "*") + sweep(da_w, 2L, ri22, "*") +
      sweep(w$y[, 2L] - w$x %*% beta_w, 2L, sigma22, "/")
    p11 <- matrix(ui11, g, length(rows), byrow = TRUE) + outer(a$n, ri11) + outer(t$n, 1 / sigma11)
    p12 <- matrix(ui12, g, length(rows), byrow = TRUE) + outer(a$n, ri12)
    p22 <- matrix(ui22, g, length(rows), byrow = TRUE) + outer(a$n, ri22) + outer(w$n, 1 / sigma22)
    schur <- p22 - p12 * (p12 / p11)
    if (any(!is.finite(schur)) || any(schur <= 0))
      stop("Joint cluster trace conditional covariance is not numerically positive definite.", call. = FALSE)
    u_w <- (h_w - p12 * (h_t / p11)) / schur
    u_t <- (h_t - p12 * u_w) / p11
    for (pattern in names(recipients)) {
      r <- recipients[[pattern]]
      means <- as.vector(crossprod(r$x, beta_t) + crossprod(r$weight, u_t))
      if (pattern == "recorded") means <- means + (sigma12 / sigma22) *
        as.vector(r$w - crossprod(r$x, beta_w) - crossprod(r$weight, u_w))
      out[rows, paste0("conditional_target_mean_", pattern)] <- means
    }
  }
  center <- unlist(problem$transform$predictor_center, use.names = FALSE)
  scale <- unlist(problem$transform$predictor_scale, use.names = FALSE)
  intercept <- out[, 1L]
  if (p > 1L) intercept <- intercept - as.vector(out[, 2L:p, drop = FALSE] %*% (center[-1L] / scale[-1L]))
  out[, "beta_T_original:(Intercept)"] <- problem$center[1L] + problem$scale[1L] * intercept
  if (any(!is.finite(out))) stop("Nonfinite joint cluster parameter trace.", call. = FALSE)
  out
}

# pan 2.0 keeps a cached normal variate outside its seeded uniform RNG. Each
# chain therefore runs in a new R process, including after an iteration retry.
catmech_joint_cluster_pan <- function(input, prior, seeds, starts, iter) {
  restore <- preserve_rng()
  workers <- NULL
  on.exit({
    if (!is.null(workers)) try(parallel::stopCluster(workers), silent = TRUE)
    restore()
  }, add = TRUE)
  run <- function(job, input, prior, iter, library_paths) {
    .libPaths(library_paths)
    pan::pan(input$y, input$subj, input$pred, xcol = seq_len(ncol(input$pred)),
      zcol = 1L, prior = prior, seed = job$seed, iter = iter, start = job$start)
  }
  environment(run) <- baseenv()
  jobs <- lapply(seq_along(seeds), function(j) list(seed = seeds[j], start = starts[[j]]))
  result <- vector("list", length(jobs))
  # Limit simultaneous workers to two, including under R CMD check --as-cran.
  # Every chain still gets a new process; a finished worker is never reused.
  for (first in seq.int(1L, length(jobs), by = 2L)) {
    batch <- seq.int(first, min(first + 1L, length(jobs)))
    # An explicit per-process port avoids inherited default-port collisions
    # when independent outer jobs were forked from the same R session.
    for (attempt in seq_len(20L)) {
      port <- as.integer(20000 + (as.double(Sys.getpid()) + sample.int(40000L, 1L)) %% 40000)
      candidate <- tryCatch(parallel::makePSOCKcluster(length(batch), port = port,
        setup_strategy = "parallel", setup_timeout = 60), error = identity)
      if (!inherits(candidate, "error")) {
        workers <- candidate
        break
      }
      if (!grepl("server socket", conditionMessage(candidate), fixed = TRUE)) stop(candidate)
    }
    if (is.null(workers)) stop("Could not allocate an isolated joint-sampler worker port.", call. = FALSE)
    result[batch] <- parallel::clusterApply(workers, jobs[batch], run, input = input,
      prior = prior, iter = iter, library_paths = .libPaths())
    parallel::stopCluster(workers)
    workers <- NULL
  }
  result
}

catmech_joint_cluster_draws <- function(data, value, gold, predictors, m, cluster, control) {
  if (m > as.double(control$iterations) * control$chains)
    stop("m cannot exceed joint_control$iterations * joint_control$chains.", call. = FALSE)
  problem <- catmech_joint_cluster_prepare(data, value, gold, predictors, cluster)
  p <- ncol(problem$x)
  prior <- list(a = 4, Binv = diag(.5 * control$prior_scale, 2L),
                c = 4, Dinv = diag(.5 * control$prior_scale, 2L))
  seeds <- sample.int(.Machine$integer.max, control$chains)
  starts <- lapply(seq_len(control$chains), function(chain) {
    y <- problem$input$y
    y[is.na(y)] <- stats::rnorm(sum(is.na(y)))
    factor <- 2^seq(-1, 1, length.out = control$chains)[chain]
    list(beta = matrix(stats::rnorm(2L * p, 0, .5), p, 2L),
         sigma = diag(factor, 2L), psi = diag(1 / factor, 2L), y = y)
  })
  iterations <- control$iterations
  attempts <- list()
  repeat {
    chains <- catmech_joint_cluster_pan(problem$input, prior, seeds, starts, control$burnin + iterations)
    retained <- control$burnin + seq_len(iterations)
    parts <- lapply(chains, function(fit) catmech_joint_cluster_trace(problem, fit, retained))
    traces <- array(NA_real_, c(iterations, control$chains, ncol(parts[[1L]])),
      dimnames = list(iteration = as.character(seq_len(iterations)),
        chain = as.character(seq_len(control$chains)), parameter = colnames(parts[[1L]])))
    for (chain in seq_len(control$chains)) traces[, chain, ] <- parts[[chain]]
    diagnostics <- catmech_model_chain_diagnostics(traces)
    diagnostics$engine <- "pan_joint_cluster"
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
  diagnostics$pan_input <- problem$input
  diagnostics$pan_prior <- prior
  diagnostics$input_rows <- problem$rows
  diagnostics$transform <- problem$transform
  diagnostics$pan_version <- as.character(utils::packageVersion("pan"))
  diagnostics$backend_isolation <- "fresh_R_worker_per_chain_and_attempt"
  if (!diagnostics$convergence_pass) stop(structure(list(
    message = paste0("Joint cluster model failed convergence checks (R-hat <= 1.05 and bulk/tail ESS >= ",
      100L * control$chains, " required). Inspect the condition's model_diagnostics."),
    call = NULL, model_diagnostics = diagnostics),
    class = c("catmech_convergence_error", "error", "condition")))
  allocation <- rep(m %/% control$chains, control$chains)
  if (m %% control$chains) allocation[seq_len(m %% control$chains)] <- allocation[seq_len(m %% control$chains)] + 1L
  draws <- matrix(NA_real_, nrow(data), m)
  selected <- list(); next_column <- 0L
  for (chain in seq_len(control$chains)) {
    keep <- if (!allocation[chain]) integer() else if (allocation[chain] == 1L) as.integer(ceiling(iterations / 2)) else
      as.integer(round(seq(1, iterations, length.out = allocation[chain])))
    for (post in keep) {
      next_column <- next_column + 1L
      k <- control$burnin + post
      beta <- matrix(chains[[chain]]$beta[, , k], p, 2L)
      target <- catmech_joint_cluster_predict(problem$y, problem$x, problem$cluster,
        beta, chains[[chain]]$sigma[, , k], chains[[chain]]$psi[, , k])
      draws[, next_column] <- target * problem$scale[1L] + problem$center[1L]
      draws[!is.na(gold), next_column] <- gold[!is.na(gold)]
      selected[[next_column]] <- data.frame(imputation = next_column, chain = chain, retained_iteration = post)
    }
  }
  diagnostics$selected_draws <- do.call(rbind, selected)
  list(draws = draws, model_diagnostics = diagnostics, transform = problem$transform)
}
