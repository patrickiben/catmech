joint_trace_fixture <- function(predictors = c("x", "z"), pattern = "mixed") {
  index <- seq_len(48L)
  within <- rep(seq_len(6L), 8L)
  d <- data.frame(x = 20 + sin(index), z = -7 + cos(index / 3),
    site = rep(seq_len(8L), each = 6L),
    raw = 3 + sin(index / 2) + index / 50,
    gold = ifelse(within <= 3L, 10 + cos(index / 2) + index / 20, NA_real_))
  if (pattern == "mixed") d$raw[within %in% c(2L, 6L)] <- NA_real_
  if (pattern == "absent") d$raw[within > 3L] <- NA_real_
  catmech_joint_cluster_prepare(d, "raw", d$gold, predictors, "site")
}

joint_trace_parameters <- function(p, n) {
  beta <- array(sin(seq_len(p * 2L * n) / 7), c(p, 2L, n))
  sigma <- psi <- array(0, c(2L, 2L, n))
  for (j in seq_len(n)) {
    sigma[, , j] <- matrix(c(1 + j / (2 * n), .6 * sin(j), .6 * sin(j), 1.7), 2L)
    psi[, , j] <- matrix(c(.8, .4 * cos(j), .4 * cos(j), 1 + j / n), 2L)
  }
  list(beta = beta, sigma = sigma, psi = psi)
}

# Reference calculation uses the unchanged scalar conditional-prediction
# helper and averages individual predictions, independently of batch weights.
joint_trace_scalar_means <- function(problem, fit, keep) {
  patterns <- list(recorded = problem$prediction$missing & problem$prediction$recorded,
                  absent = problem$prediction$missing & !problem$prediction$recorded)
  patterns <- patterns[vapply(patterns, any, logical(1))]
  do.call(rbind, lapply(keep, function(k) {
    beta <- matrix(fit$beta[, , k], ncol(problem$x), 2L)
    conditional <- catmech_joint_cluster_condition(problem$prediction, beta,
      fit$sigma[, , k], fit$psi[, , k])
    mu <- problem$x %*% beta + conditional$mean[problem$prediction$group, , drop = FALSE]
    vapply(names(patterns), function(pattern) {
      rows <- patterns[[pattern]]
      if (pattern == "recorded") mean(mu[rows, 1L] + conditional$lambda *
        (problem$y[rows, 2L] - mu[rows, 2L])) else mean(mu[rows, 1L])
    }, numeric(1))
  }))
}

test_that("batched joint traces agree with individual conditional means", {
  problem <- joint_trace_fixture()
  fit <- joint_trace_parameters(ncol(problem$x), 2051L)
  keep <- c(2051L:1L, 1024L)
  set.seed(100)
  before <- .Random.seed
  actual <- catmech_joint_cluster_trace(problem, fit, keep)
  expect_identical(.Random.seed, before)
  expected <- joint_trace_scalar_means(problem, fit, keep)
  expect_equal(unname(actual[, paste0("conditional_target_mean_", colnames(expected)), drop = FALSE]),
    unname(expected), tolerance = 1e-12)
  expected_raw <- t(vapply(keep, function(k) c(as.vector(fit$beta[, , k]),
    fit$sigma[1L, 1L, k], fit$sigma[1L, 2L, k], fit$sigma[2L, 2L, k],
    fit$psi[1L, 1L, k], fit$psi[1L, 2L, k], fit$psi[2L, 2L, k]), numeric(12L)))
  expect_identical(unname(actual[, seq_len(12L)]), unname(expected_raw))
})

test_that("joint trace handles one draw, intercept-only design and one recipient pattern", {
  for (pattern in c("recorded", "absent", "mixed")) {
    problem <- joint_trace_fixture(character(), pattern)
    fit <- joint_trace_parameters(1L, 1L)
    actual <- catmech_joint_cluster_trace(problem, fit, 1L)
    expected <- joint_trace_scalar_means(problem, fit, 1L)
    expect_equal(unname(actual[, paste0("conditional_target_mean_", colnames(expected)), drop = FALSE]),
      unname(expected), tolerance = 1e-12)
    expect_equal(unname(actual[, "beta_T_original:(Intercept)"]),
      unname(problem$center[1L] + problem$scale[1L] * fit$beta[1L, 1L, 1L]), tolerance = 1e-14)
  }
})

test_that("joint original-scale intercept monitor evaluates the original design origin", {
  problem <- joint_trace_fixture()
  fit <- joint_trace_parameters(ncol(problem$x), 7L)
  actual <- catmech_joint_cluster_trace(problem, fit, seq_len(7L))
  center <- unlist(problem$transform$predictor_center)
  scale <- unlist(problem$transform$predictor_scale)
  original_origin <- c(1, (0 - center[-1L]) / scale[-1L])
  expected <- vapply(seq_len(7L), function(k) problem$center[1L] + problem$scale[1L] *
    sum(original_origin * fit$beta[, 1L, k]), numeric(1))
  expect_equal(unname(actual[, "beta_T_original:(Intercept)"]), unname(expected), tolerance = 1e-12)
  expect_gt(max(abs(expected - (problem$center[1L] + problem$scale[1L] * fit$beta[1L, 1L, ]))), 1)
})

test_that("joint trace rejects invalid covariance instead of returning misleading diagnostics", {
  problem <- joint_trace_fixture()
  fit <- joint_trace_parameters(ncol(problem$x), 1L)
  fit$sigma[, , 1L] <- matrix(c(1, 2, 2, 1), 2L)
  expect_error(catmech_joint_cluster_trace(problem, fit, 1L), "positive definite")
})

test_that("longer joint defaults retain explicit historical controls", {
  defaults <- catmech_joint_cluster_control()
  expect_identical(defaults$iterations, 4000L)
  expect_identical(defaults$max_iterations, 32000L)
  expect_identical(defaults[c("chains", "burnin", "prior_scale")],
    list(chains = 4L, burnin = 1000L, prior_scale = 1))
  old <- catmech_joint_cluster_control(list(iterations = 2000L, max_iterations = 8000L))
  expect_identical(old$iterations, 2000L)
  expect_identical(old$max_iterations, 8000L)
})
