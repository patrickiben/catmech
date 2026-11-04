# Independent fixed-covariance Gaussian oracle. Integrating all site effects
# gives one observed-data covariance matrix. Flat-prior regression then has an
# exact normal posterior for (B_T, B_W). These are coefficient-summary tests,
# not simulation evidence about unknown-covariance fitting or interval coverage.
joint_analysis_oracle <- function(intercept_only = FALSE) {
  n <- 24L
  d <- data.frame(raw = .3 * seq_len(n) + sin(seq_len(n)),
    gold = .3 * seq_len(n) + sin(seq_len(n)) + .2 * cos(seq_len(n)),
    age = seq(8, 20, length.out = n),
    group = factor(rep(c("a", "b", "a"), length.out = n)),
    site = rep(LETTERS[1:6], each = 4L))
  d$gold[19:24] <- NA_real_
  d$raw[c(4, 8, 12, 16, 20, 24)] <- NA_real_
  predictors <- if (intercept_only) character() else c("age", "group")
  x <- if (intercept_only) matrix(1, n, 1L,
    dimnames = list(NULL, "(Intercept)")) else stats::model.matrix(~ age + group, d)
  p <- ncol(x)
  y <- cbind(T = d$gold, W = d$raw)
  sigma <- matrix(c(1.2, .45, .45, .8), 2L)
  psi <- matrix(c(.7, .3, .3, .5), 2L)
  k <- kronecker(diag(n), sigma) +
    kronecker(outer(d$site, d$site, `==`) * 1, psi)
  design <- matrix(0, 2L * n, 2L * p)
  design[seq(1L, 2L * n, 2L), seq_len(p)] <- x
  design[seq(2L, 2L * n, 2L), p + seq_len(p)] <- x
  observed <- which(!is.na(as.vector(t(y))))
  a <- design[observed, , drop = FALSE]
  ko <- k[observed, observed, drop = FALSE]
  covariance <- solve(crossprod(a, solve(ko, a)))
  location <- as.vector(covariance %*% crossprod(a, solve(ko, as.vector(t(y))[observed])))
  iterations <- 4000L
  chains <- 4L
  set.seed(if (intercept_only) 6102L else 6101L)
  original_draws <- sweep(matrix(stats::rnorm(iterations * chains * 2L * p),
    iterations * chains, 2L * p) %*% chol(covariance), 2L, location, "+")
  response_center <- colMeans(y, na.rm = TRUE)
  response_scale <- apply(y, 2L, stats::sd, na.rm = TRUE)
  predictor_center <- stats::setNames(rep(0, p), colnames(x))
  predictor_scale <- stats::setNames(rep(1, p), colnames(x))
  if (p > 1L) {
    predictor_center[-1L] <- colMeans(x[, -1L, drop = FALSE])
    predictor_scale[-1L] <- apply(x[, -1L, drop = FALSE], 2L, stats::sd)
  }
  standardized <- original_draws
  for (response in 1:2) {
    index <- (response - 1L) * p + seq_len(p)
    standardized[, index] <- sweep(original_draws[, index, drop = FALSE],
      2L, predictor_scale / response_scale[response], "*")
    standardized[, index[1L]] <- (original_draws[, index[1L]] -
      response_center[response] +
      as.vector(original_draws[, index, drop = FALSE] %*% predictor_center)) /
        response_scale[response]
  }
  sigma <- sigma / outer(response_scale, response_scale)
  psi <- psi / outer(response_scale, response_scale)
  names <- c(paste0("beta_T:", colnames(x)), paste0("beta_W:", colnames(x)),
    "sigma_TT", "sigma_TW", "sigma_WW", "psi_TT", "psi_TW", "psi_WW",
    "beta_T_original:(Intercept)")
  traces <- array(NA_real_, c(iterations, chains, length(names)),
    dimnames = list(iteration = as.character(seq_len(iterations)),
      chain = as.character(seq_len(chains)), parameter = names))
  for (j in seq_len(2L * p)) traces[, , j] <- standardized[, j]
  constants <- c(sigma[1, 1], sigma[1, 2], sigma[2, 2],
    psi[1, 1], psi[1, 2], psi[2, 2])
  for (j in seq_along(constants)) traces[, , 2L * p + j] <- constants[j]
  traces[, , "beta_T_original:(Intercept)"] <- original_draws[, 1L]
  transform <- list(response_center = as.list(response_center),
    response_scale = as.list(response_scale), predictor_center = as.list(predictor_center),
    predictor_scale = as.list(predictor_scale))
  fit <- structure(list(original = d, completed = list(d, d),
    manifest = list(engine = "reference_joint_cluster", m = 2L),
    spec = list(cluster_model = "joint", predictors = predictors, method = "normal",
      value = "raw", reference = "gold", cluster = "site",
      recorded_fit = "joint_bivariate_gaussian", joint_transform = transform,
      joint_control = list(chains = chains, burnin = 0L, iterations = iterations,
        max_iterations = iterations, prior_scale = 1)),
    model_diagnostics = list(engine = "pan_joint_cluster", convergence_pass = TRUE,
      traces = traces, iterations_used = iterations, transform = transform)),
    class = "catmech_result")
  list(fit = fit, location = location[seq_len(p)],
    covariance = covariance[seq_len(p), seq_len(p), drop = FALSE],
    original_draws = original_draws[, seq_len(p), drop = FALSE], x = x)
}

test_that("joint coefficient summaries agree with the observed-Gaussian oracle", {
  oracle <- joint_analysis_oracle()
  before <- .Random.seed
  out <- catmech_joint_summary(oracle$fit)
  expect_identical(.Random.seed, before)
  expect_identical(out$term, colnames(oracle$x))
  expect_equal(out$estimate, colMeans(oracle$original_draws), tolerance = 1e-12)
  expect_equal(out$posterior_sd, apply(oracle$original_draws, 2L, stats::sd), tolerance = 1e-12)
  expected_interval <- t(apply(oracle$original_draws, 2L, stats::quantile,
    probs = c(.025, .975), names = FALSE, type = 7L))
  expect_equal(unname(as.matrix(out[c("lower", "upper")])),
    unname(expected_interval), tolerance = 1e-12)
  repetitions <- nrow(oracle$original_draws)
  variance <- diag(oracle$covariance)
  # Six analytic MCSE bounds fixed before drawing: Gaussian means, variances,
  # and equal-tail quantiles from an independent observed-data GLS posterior.
  expect_lt(max(abs(out$estimate - oracle$location) /
    sqrt(variance / repetitions)), 6)
  expect_lt(max(abs(out$posterior_sd^2 - variance) /
    sqrt(2 * variance^2 / (repetitions - 1L))), 6)
  for (j in seq_along(variance)) {
    wanted <- oracle$location[j] + sqrt(variance[j]) * stats::qnorm(c(.025, .975))
    mcse <- sqrt(variance[j]) * sqrt(.025 * .975 / repetitions) /
      stats::dnorm(stats::qnorm(.975))
    expect_lt(max(abs(unlist(out[j, c("lower", "upper")]) - wanted)) / mcse, 6)
    samples <- matrix(oracle$original_draws[, j], nrow = 4000L)
    expect_equal(out$mcse_mean[j], posterior::mcse_mean(samples), tolerance = 1e-10)
    expect_equal(unlist(out[j, c("mcse_lower", "mcse_upper")], use.names = FALSE),
      posterior::mcse_quantile(samples, c(.025, .975), names = FALSE), tolerance = 1e-10)
    # Roundoff in the inverse transformation can break a tie between the two
    # central folded ranks, without materially changing this diagnostic.
    expect_equal(out$rhat[j], posterior::rhat(samples), tolerance = 1e-6)
    expect_equal(out$ess_bulk[j], posterior::ess_bulk(samples), tolerance = 1e-10)
    expect_equal(out$ess_tail[j], posterior::ess_tail(samples), tolerance = 1e-10)
  }
  expect_true(all(out$draws == 16000L & out$iterations == 4000L & out$chains == 4L))
  expect_true(all(out$interpretation == "conditional_joint_model_credible_interval"))
  expect_false(any(c("df", "std.error", "cluster_variance") %in% names(out)))
  # Omitting predictor centering from the intercept must fail this oracle.
  wrong_intercept <- oracle$fit$model_diagnostics$traces[, , "beta_T:(Intercept)"] *
    oracle$fit$spec$joint_transform$response_scale$T +
    oracle$fit$spec$joint_transform$response_center$T
  expect_gt(abs(mean(wrong_intercept) - oracle$location[1L]) /
    sqrt(variance[1L] / repetitions), 10)
})

test_that("joint summaries use all chains independently of completed datasets", {
  oracle <- joint_analysis_oracle()
  expected <- catmech_joint_summary(oracle$fit, level = .8)
  altered <- oracle$fit
  altered$completed <- NULL
  altered$manifest$m <- 100000L
  altered$model_diagnostics$selected_draws <- data.frame(chain = 1L, retained_iteration = 1L)
  expect_identical(catmech_joint_summary(altered, level = .8), expected)
  wanted <- t(apply(oracle$original_draws, 2L, stats::quantile,
    probs = c(.1, .9), names = FALSE, type = 7L))
  expect_equal(unname(as.matrix(expected[c("lower", "upper")])), unname(wanted), tolerance = 1e-12)
  expect_true(all(expected$level == .8))
  only <- joint_analysis_oracle(intercept_only = TRUE)
  out <- catmech_joint_summary(only$fit)
  expect_identical(out$term, "(Intercept)")
  expect_equal(out$estimate, mean(only$original_draws), tolerance = 1e-12)
  expect_equal(out$posterior_sd, stats::sd(as.vector(only$original_draws)), tolerance = 1e-12)
})

test_that("joint summaries refuse unsupported workflows and absent posterior fits", {
  fit <- joint_analysis_oracle()$fit
  expect_error(catmech_joint_summary(list()), "cluster_model='joint'", fixed = TRUE)
  other <- fit
  other$manifest$engine <- "reference_normal_recorded"
  expect_error(catmech_joint_summary(other), "cluster_model='joint'", fixed = TRUE)
  for (level in list(0, 1, NA_real_, Inf, c(.8, .9), "95%", .95 + 1i))
    expect_error(catmech_joint_summary(fit, level), "level")
  absent <- fit
  absent$spec$recorded_fit <- "not_fitted_all_reference_observed"
  absent$model_diagnostics$engine <- "reference_complete"
  absent$model_diagnostics$traces <- NULL
  expect_error(catmech_joint_summary(absent), "all-reference shortcut")
  fit$model_diagnostics$convergence_pass <- FALSE
  error <- tryCatch(catmech_joint_summary(fit), error = identity)
  expect_s3_class(error, "catmech_convergence_error")
  expect_identical(error$model_diagnostics, fit$model_diagnostics)
})

test_that("joint summaries reject malformed trace and transformation metadata", {
  fit <- joint_analysis_oracle()$fit
  change <- function(mutate, pattern) {
    bad <- mutate(fit)
    expect_error(catmech_joint_summary(bad), pattern)
  }
  change(function(z) { z$model_diagnostics$traces[1, 1, 1] <- NA_real_; z }, "finite")
  change(function(z) { z$model_diagnostics$traces <- matrix(1, 100, 4); z }, "array")
  change(function(z) { z$model_diagnostics$traces <- z$model_diagnostics$traces[, 1:3, , drop = FALSE]; z }, "four chains")
  change(function(z) { z$model_diagnostics$iterations_used <- 2000L; z }, "dimensions")
  change(function(z) { dimnames(z$model_diagnostics$traces)[[3]][2] <- "beta_T:wrong"; z }, "imputation design")
  change(function(z) { dimnames(z$model_diagnostics$traces)[[3]][2] <- "beta_T:(Intercept)"; z }, "unique")
  change(function(z) { z$model_diagnostics$traces[, , "sigma_TT"] <- -1; z }, "positive definite")
  change(function(z) { z$model_diagnostics$traces[, , "psi_TW"] <- 100; z }, "positive definite")
  change(function(z) { z$spec$joint_transform$response_scale$T <- 0; z }, "disagree")
  change(function(z) {
    z$spec$joint_transform$response_scale$T <- 0
    z$model_diagnostics$transform <- z$spec$joint_transform; z
  }, "response_scale")
  change(function(z) {
    z$spec$joint_transform$predictor_center[["(Intercept)"]] <- 1
    z$model_diagnostics$transform <- z$spec$joint_transform; z
  }, "unscaled intercept")
  change(function(z) { z$model_diagnostics$traces[, , "beta_T_original:(Intercept)"] <- 0; z }, "monitor conflicts")
})

test_that("joint summaries recheck the original-scale intercept combination", {
  fit <- joint_analysis_oracle()$fit
  trace <- fit$model_diagnostics$traces
  n <- dim(trace)[1L]
  centers <- unlist(fit$spec$joint_transform$predictor_center)
  scales <- unlist(fit$spec$joint_transform$predictor_scale)
  set.seed(7153)
  a <- matrix(stats::rnorm(4L * n), n, 4L)
  b <- matrix(stats::rnorm(4L * n), n, 4L)
  trace[, , "beta_T:age"] <- a
  trace[, , "beta_T:groupb"] <- b
  trace[, , "beta_T:(Intercept)"] <- centers[2] / scales[2] * a +
    centers[3] / scales[3] * b +
    matrix(rep(c(0, .01, .02, .03), each = n), n, 4L) +
    matrix(stats::rnorm(4L * n, sd = .0001), n, 4L)
  # Individual standardized coefficients pass, but their nearly cancelling
  # intercept combination has distinct chain means. Simulate an older result
  # without the newer explicit original-intercept monitor.
  kept <- !dimnames(trace)[[3L]] %in% "beta_T_original:(Intercept)"
  fit$model_diagnostics$traces <- trace[, , kept, drop = FALSE]
  for (name in paste0("beta_T:", c("(Intercept)", "age", "groupb"))) {
    expect_lt(posterior::rhat(trace[, , name]), 1.05)
    expect_gt(posterior::ess_bulk(trace[, , name]), 400)
    expect_gt(posterior::ess_tail(trace[, , name]), 400)
  }
  error <- tryCatch(catmech_joint_summary(fit), error = identity)
  expect_s3_class(error, "catmech_convergence_error")
  expect_false(error$model_diagnostics$convergence_pass)
  expect_identical(error$model_diagnostics$engine, "joint_population_fixed_effects")
  expect_gt(error$model_diagnostics$summary$rhat[1L], 1.05)
  expect_equal(dim(error$model_diagnostics$traces), c(n, 4L, 3L))
})
