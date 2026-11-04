# The predictive oracle uses one marginal covariance matrix for every observed
# and missing component. It deliberately does not use the implementation's
# per-site precision updates or conditional-row formulas.
joint_cluster_gaussian_oracle <- function(y, x, cluster, beta, sigma, psi) {
  n <- nrow(y)
  marginal <- kronecker(diag(n), sigma) +
    kronecker(outer(cluster, cluster, `==`) * 1, psi)
  mean <- as.vector(t(x %*% beta))
  observed <- which(!is.na(as.vector(t(y))))
  target <- 2L * which(is.na(y[, 1L])) - 1L
  cross <- marginal[target, observed, drop = FALSE]
  observed_covariance <- marginal[observed, observed, drop = FALSE]
  gain <- cross %*% solve(observed_covariance)
  list(rows = which(is.na(y[, 1L])), mean = as.vector(mean[target] +
    gain %*% (as.vector(t(y))[observed] - mean[observed])),
    covariance = marginal[target, target, drop = FALSE] -
      gain %*% t(cross))
}

joint_cluster_oracle_fixture <- function() {
  list(y = cbind(T = c(.8, NA, NA, 1.2, NA, NA, .7, NA, NA, NA, -.4),
                 W = c(1.1, 2.1, NA, NA, -1, NA, .2, 1.2, NA, NA, NA)),
    x = cbind(1, c(-.7, .2, .9, -.3, .5, -.8, .3, .6, -.1, .4, -.6)),
    cluster = c(rep("A", 4), rep("B", 2), rep("C", 2), rep("new", 2), "E"),
    beta = matrix(c(.3, .4, -.2, .6), 2L),
    sigma = matrix(c(1.2, .45, .45, .8), 2L),
    psi = matrix(c(.7, .3, .3, .5), 2L))
}

test_that("joint prediction agrees with a full Gaussian conditioning oracle", {
  f <- joint_cluster_oracle_fixture()
  oracle <- do.call(joint_cluster_gaussian_oracle, f)
  set.seed(171034)
  # Six analytic Monte Carlo standard errors were specified before drawing.
  # This checks every marginal mean and every pairwise covariance, including
  # recorded-present and recorded-missing targets in the same cluster.
  repetitions <- 10000L
  draws <- replicate(repetitions, do.call(catmech_joint_cluster_predict, f))
  known <- !is.na(f$y[, 1L])
  expect_true(all(is.finite(draws)))
  expect_equal(unname(draws[known, , drop = FALSE]),
    matrix(rep(f$y[known, 1L], repetitions), nrow = sum(known)), tolerance = 0)
  missing_draws <- t(draws[oracle$rows, , drop = FALSE])
  means <- colMeans(missing_draws)
  covariance <- stats::cov(missing_draws)
  mean_mcse <- sqrt(diag(oracle$covariance) / repetitions)
  covariance_mcse <- sqrt((outer(diag(oracle$covariance),
    diag(oracle$covariance)) + oracle$covariance^2) / (repetitions - 1L))
  expect_lt(max(abs(means - oracle$mean) / mean_mcse), 6)
  expect_lt(max(abs(covariance - oracle$covariance) / covariance_mcse), 6)

  # Prespecified negative controls establish that the oracle can detect
  # independent site draws for the two patterns and ignored W-only site data.
  a <- match(c(2L, 3L), oracle$rows)
  expect_gt(abs(oracle$covariance[a[1L], a[2L]]) /
    covariance_mcse[a[1L], a[2L]], 8)
  new <- match(c(9L, 10L), oracle$rows)
  expect_equal(oracle$covariance[new[1L], new[2L]], f$psi[1L, 1L],
    tolerance = 1e-12)
  expect_equal(diag(oracle$covariance)[new],
    rep(f$sigma[1L, 1L] + f$psi[1L, 1L], 2L), tolerance = 1e-12)
  w_only <- match(6L, oracle$rows)
  unconditional <- as.vector(f$x %*% f$beta)[6L]
  expect_gt(abs(oracle$mean[w_only] - unconditional) / mean_mcse[w_only], 8)
})

test_that("joint prediction conditions on gold-only anchors and restores row order", {
  f <- joint_cluster_oracle_fixture()
  oracle <- do.call(joint_cluster_gaussian_oracle, f)
  no_anchor <- f
  no_anchor$y[4L, 1L] <- NA_real_
  omitted <- do.call(joint_cluster_gaussian_oracle, no_anchor)
  a <- match(3L, oracle$rows)
  expect_gt(abs(oracle$mean[a] - omitted$mean[match(3L, omitted$rows)]), .01)
  permutation <- c(11L, 3L, 9L, 1L, 6L, 4L, 8L, 2L, 10L, 5L, 7L)
  shuffled <- f
  shuffled$y <- f$y[permutation, , drop = FALSE]
  shuffled$x <- f$x[permutation, , drop = FALSE]
  shuffled$cluster <- f$cluster[permutation]
  expected <- do.call(joint_cluster_gaussian_oracle, shuffled)
  order_in_original <- match(permutation[expected$rows], oracle$rows)
  expect_equal(expected$mean, oracle$mean[order_in_original], tolerance = 1e-12)
  expect_equal(expected$covariance,
    oracle$covariance[order_in_original, order_in_original], tolerance = 1e-12)
  set.seed(171035)
  draws <- replicate(2500L, do.call(catmech_joint_cluster_predict, shuffled))
  mean_mcse <- sqrt(diag(expected$covariance) / ncol(draws))
  expect_lt(max(abs(rowMeans(draws[expected$rows, , drop = FALSE]) -
    expected$mean) / mean_mcse), 6)
  known <- which(!is.na(shuffled$y[, 1L]))
  expect_equal(unname(draws[known, 1L]), unname(shuffled$y[known, 1L]), tolerance = 0)
})

joint_cluster_api_fixture <- function() {
  set.seed(30911)
  site <- rep(sprintf("site%02d", seq_len(24L)), each = 10L)
  x <- stats::rnorm(length(site))
  u <- matrix(stats::rnorm(48L), 24L, 2L) %*%
    chol(matrix(c(.7, .3, .3, .6), 2L))
  e <- matrix(stats::rnorm(480L), 240L, 2L) %*%
    chol(matrix(c(1, .35, .35, .8), 2L))
  truth <- .3 + .6 * x + u[rep(1:24, each = 10L), 1L] + e[, 1L]
  recorded <- -.2 + .8 * x + u[rep(1:24, each = 10L), 2L] + e[, 2L]
  within <- rep(seq_len(10L), 24L)
  gold <- ifelse(within <= 4L & site != "site24", truth, NA_real_)
  recorded[within %in% c(4L, 8L, 9L)] <- NA_real_
  # A completely unobserved site belongs to prediction, not covariance fitting.
  data.frame(raw = c(recorded, NA, NA), gold = c(gold, NA, NA),
    x = c(x, -.2, .4), site = c(site, "new", "new"), row = seq_len(242L))
}

test_that("joint mode preserves anchors, original data, mids and caller RNG", {
  d <- joint_cluster_api_fixture()
  set.seed(922)
  before <- .Random.seed
  fit <- catmech_reference(d, "raw", "gold", "x", use_recorded = TRUE,
    cluster = "site", cluster_model = "joint", m = 4L, seed = 24L)
  expect_identical(.Random.seed, before)
  expect_identical(fit$original, d)
  known <- !is.na(d$gold)
  for (j in seq_along(fit$completed)) {
    expect_identical(fit$completed[[j]]$raw[known], d$gold[known])
    expect_true(all(is.finite(fit$completed[[j]]$raw)))
    expect_identical(fit$completed[[j]][setdiff(names(d), "raw")],
      d[setdiff(names(d), "raw")])
    expect_equal(mice::complete(fit$mids, j)$raw, fit$completed[[j]]$raw)
  }
  values <- vapply(fit$completed, `[[`, numeric(nrow(d)), "raw")
  expect_true(all(apply(values[!known, , drop = FALSE], 1L, stats::var) > 0))
  expect_identical(fit$manifest$imputed_rows, which(!known))
  expect_true(fit$model_diagnostics$convergence_pass)
  expect_identical(fit$spec$cluster_model, "joint")

  # Independently reconstruct the standardized pan inputs and replay its first
  # chain. This is a wrapper contract check, distinct from the Gaussian oracle.
  metadata <- fit$model_diagnostics
  y <- cbind(T = d$gold, W = d$raw)
  y <- sweep(sweep(y, 2L, colMeans(y, na.rm = TRUE), "-"), 2L,
    apply(y, 2L, stats::sd, na.rm = TRUE), "/")
  x <- cbind(`(Intercept)` = 1, x = (d$x - mean(d$x)) / stats::sd(d$x))
  active <- which(rowSums(!is.na(y)) > 0L)
  groups <- match(d$site[active], sort(unique(d$site[active])))
  rows <- active[order(groups)]
  expect_equal(metadata$input_rows, rows, tolerance = 0)
  expect_equal(metadata$pan_input$y, y[rows, , drop = FALSE], tolerance = 0)
  expect_equal(unname(metadata$pan_input$pred), unname(x[rows, , drop = FALSE]),
    tolerance = 1e-14)
  expect_identical(metadata$pan_input$subj, groups[order(groups)])
  expect_equal(metadata$pan_prior, list(a = 4, Binv = diag(.5, 2L),
    c = 4, Dinv = diag(.5, 2L)), tolerance = 0)
  observed <- !is.na(y[rows, , drop = FALSE])
  for (start in metadata$chain_starts)
    expect_equal(start$y[observed], y[rows, , drop = FALSE][observed], tolerance = 0)
  # pan's saved Fortran normal cache is not reset by its seed argument. Use a
  # fresh process for the independent first-chain replay, as for a new fit.
  direct_args <- list(y = y[rows, , drop = FALSE], subj = groups[order(groups)],
    pred = x[rows, , drop = FALSE], xcol = 1:2, zcol = 1L,
    prior = list(a = 4, Binv = diag(.5, 2L), c = 4, Dinv = diag(.5, 2L)),
    seed = metadata$chain_seeds[1L],
    iter = 1000L + metadata$iterations_used, start = metadata$chain_starts[[1L]])
  replay <- function(args) {
    worker <- parallel::makePSOCKcluster(1L)
    on.exit(parallel::stopCluster(worker), add = TRUE)
    parallel::clusterCall(worker, function(paths) .libPaths(paths), .libPaths())
    parallel::clusterCall(worker, function(args) do.call(pan::pan, args), args)[[1L]]
  }
  direct <- replay(direct_args)
  kept <- 1000L + seq_len(metadata$iterations_used)
  expected_trace <- t(vapply(kept, function(k) c(as.vector(direct$beta[, , k]),
    direct$sigma[1L, 1L, k], direct$sigma[1L, 2L, k], direct$sigma[2L, 2L, k],
    direct$psi[1L, 1L, k], direct$psi[1L, 2L, k], direct$psi[2L, 2L, k]), numeric(10L)))
  expect_equal(unname(metadata$traces[, 1L, seq_len(10L)]),
    unname(expected_trace), tolerance = 1e-12)
  expect_true(all(metadata$selected_draws$retained_iteration >= 1L &
    metadata$selected_draws$retained_iteration <= metadata$iterations_used))
  expect_equal(tabulate(metadata$selected_draws$chain, nbins = 4L), rep(1L, 4L))

  # An unrelated call in the caller's process must not affect later fits.
  perturb <- direct_args
  perturb$iter <- 1L
  perturb$seed <- 839L
  invisible(do.call(pan::pan, perturb))
  repeated <- catmech_reference(d, "raw", "gold", "x", use_recorded = TRUE,
    cluster = "site", cluster_model = "joint", m = 4L, seed = 24L)
  expect_identical(repeated$completed, fit$completed)
  expect_identical(.Random.seed, before)
})

test_that("joint mode keeps caller row order with interleaved site labels", {
  d <- joint_cluster_api_fixture()
  permutation <- order(d$x)
  d <- d[permutation, , drop = FALSE]
  rownames(d) <- NULL
  fit <- catmech_reference(d, "raw", "gold", "x", use_recorded = TRUE,
    cluster = "site", cluster_model = "joint", m = 2L, seed = 25L)
  known <- !is.na(d$gold)
  for (out in fit$completed) {
    expect_identical(out$row, d$row)
    expect_identical(out$site, d$site)
    expect_identical(out$gold, d$gold)
    expect_identical(out$raw[known], d$gold[known])
    expect_true(all(is.finite(out$raw)))
  }
})

test_that("joint mode never treats calibration agreement as exact unaudited truth", {
  d <- joint_cluster_api_fixture()
  known <- !is.na(d$gold)
  d$raw[known] <- d$gold[known]
  fit <- catmech_reference(d, "raw", "gold", "x", use_recorded = TRUE,
    cluster = "site", cluster_model = "joint", m = 4L, seed = 26L)
  recorded <- !known & !is.na(d$raw)
  values <- vapply(fit$completed, `[[`, numeric(nrow(d)), "raw")
  expect_true(all(apply(values[recorded, , drop = FALSE], 1L, stats::var) > 0))
  expect_true(any(values[recorded, 1L] != d$raw[recorded]))
  expect_false(identical(fit$spec$recorded_fit,
    "recorded_kept_all_calibration_rows_agree"))
})

test_that("fully referenced joint requests need no stochastic model", {
  d <- data.frame(raw = c(100, NA, -100), gold = c(2, 2, 2),
    x = 1, site = c("one", "one", "two"))
  set.seed(101)
  before <- .Random.seed
  fit <- catmech_reference(d, "raw", "gold", "x", use_recorded = TRUE,
    cluster = "site", cluster_model = "joint", m = 2L, seed = 27L)
  expect_identical(.Random.seed, before)
  expect_length(fit$manifest$imputed_rows, 0L)
  for (out in fit$completed) expect_identical(out$raw, d$gold)
  expect_identical(fit$model_diagnostics$status, "not_applicable")
})

test_that("legacy complete-recorded draws are unchanged by explicit moment selection", {
  d <- joint_cluster_api_fixture()
  # Fill the recorded gaps solely for this compatibility check; no new model
  # claims are made about the moment estimator.
  d$raw[is.na(d$raw)] <- .2 + .5 * d$x[is.na(d$raw)]
  old <- catmech_reference(d, "raw", "gold", "x", use_recorded = TRUE,
    cluster = "site", m = 3L, seed = 28L)
  explicit <- catmech_reference(d, "raw", "gold", "x", use_recorded = TRUE,
    cluster = "site", cluster_model = "moment", m = 3L, seed = 28L)
  expect_identical(explicit$completed, old$completed)
  expect_identical(explicit$model_diagnostics, old$model_diagnostics)
})

test_that("joint controls and observed-data support are checked before sampling", {
  d <- joint_cluster_api_fixture()
  run <- function(data = d, m = 2L, ...) catmech_reference(data, "raw", "gold", "x",
    use_recorded = TRUE, cluster = "site", cluster_model = "joint", m = m, ...)
  for (control in list(list(prior_scale = 0), list(prior_scale = Inf),
    list(chains = 3L), list(iterations = 99L),
    list(iterations = 1000L, max_iterations = 999L), list(unknown = 1),
    list(1L), structure(list(1L, 2L), names = c("chains", "chains"))))
    expect_error(run(joint_control = control), "joint_control")
  expect_error(catmech_reference(d, "raw", "gold", "x", cluster_model = "joint"),
    "cluster|use_recorded")
  expect_error(run(m = 10000L), "m cannot exceed")
  e <- d
  e$raw[!is.na(e$gold)] <- NA_real_
  expect_error(run(e), "paired calibration")
  e <- d
  e$gold[!e$site %in% unique(e$site)[1:4]] <- NA_real_
  expect_error(run(e), "five.*sites")
  e <- d
  seen <- !is.na(e$gold)
  e$gold[seen & duplicated(e$site)] <- NA_real_
  expect_error(run(e), "within-site replication")
  e <- d
  e$x[!is.na(e$gold)] <- 1
  expect_error(run(e), "full column rank|constant predictor")
})

test_that("joint convergence refusal retains evidence and restores caller RNG", {
  d <- joint_cluster_api_fixture()
  set.seed(717)
  before <- .Random.seed
  error <- tryCatch(catmech_reference(d, "raw", "gold", "x", use_recorded = TRUE,
    cluster = "site", cluster_model = "joint", m = 2L, seed = 29L,
    joint_control = list(burnin = 0L, iterations = 100L,
      max_iterations = 100L, prior_scale = 2)), error = identity)
  expect_s3_class(error, "catmech_convergence_error")
  expect_identical(.Random.seed, before)
  diagnostic <- error$model_diagnostics
  expect_false(diagnostic$convergence_pass)
  expect_identical(diagnostic$iterations_used, 100L)
  expect_equal(dim(diagnostic$traces)[1:2], c(100L, 4L))
  expect_equal(nrow(diagnostic$attempts), 1L)
  expect_equal(diagnostic$pan_prior$Binv, diag(1, 2L), tolerance = 0)
  expect_equal(diagnostic$pan_prior$Dinv, diag(1, 2L), tolerance = 0)
  expect_identical(diagnostic$backend_isolation,
    "fresh_R_worker_per_chain_and_attempt")
})

test_that("joint holdout refits forward controls and hide each reference once", {
  d <- data.frame(raw = seq_len(18L) / 3, gold = c(seq_len(12L) / 2,
    rep(NA_real_, 6L)), x = seq(-1, 1, length.out = 18L),
    site = rep(c("A", "B", "C"), each = 6L))
  d$raw[c(3L, 9L, 15L)] <- NA_real_
  d$raw[c(2L, 4L, 6L)] <- d$gold[c(2L, 4L, 6L)]
  control <- list(chains = 4L, burnin = 200L, iterations = 400L,
    max_iterations = 800L, prior_scale = 4)
  fit <- structure(list(original = d,
    manifest = list(engine = "reference_joint_cluster"),
    spec = list(value = "raw", reference = "gold", predictors = "x",
      cluster = "site", cluster_model = "joint", joint_control = control,
      # Full-fit transforms must not be reused after a reference fold is hidden.
      joint_transform = list(response_center = list(T = 9999, W = 8888)))),
    class = "catmech_result")
  calls <- list()
  mock <- function(data, value, reference, predictors, method, m, seed,
                   use_recorded, cluster, cluster_model, joint_control, ...) {
    calls[[length(calls) + 1L]] <<- list(data = data, value = value,
      reference = reference, predictors = predictors, method = method,
      m = m, seed = seed, use_recorded = use_recorded, cluster = cluster,
      cluster_model = cluster_model, joint_control = joint_control,
      extra = list(...))
    # These deterministic dummy completions test forwarding only, not predictive
    # performance. They use no held-out reference value.
    list(completed = lapply(seq_len(m), function(j) {
      out <- data
      out[[value]] <- ifelse(is.na(data[[reference]]),
        1 + .4 * data$x + j / 10, data[[reference]])
      out
    }))
  }
  testthat::local_mocked_bindings(catmech_reference = mock, .package = "catmech")
  set.seed(318)
  before <- .Random.seed
  scored <- catmech_holdout_check(fit, folds = 3L, m = 20L, seed = 17L)
  expect_identical(.Random.seed, before)
  expect_length(calls, 3L)
  expect_length(scored$failures, 0L)
  known <- which(!is.na(d$gold))
  hidden <- vector("list", length(calls))
  for (j in seq_along(calls)) {
    call <- calls[[j]]
    expect_identical(call$value, "raw")
    expect_identical(call$reference, "gold")
    expect_identical(call$predictors, "x")
    expect_identical(call$method, "normal")
    expect_identical(call$use_recorded, TRUE)
    expect_identical(call$cluster, "site")
    expect_identical(call$cluster_model, "joint")
    expect_identical(call$joint_control, control)
    expect_identical(call$m, 20L)
    expect_identical(call$seed, 17L + j)
    expect_length(call$extra, 0L)
    expect_identical(call$data[setdiff(names(d), "gold")],
      d[setdiff(names(d), "gold")])
    hidden[[j]] <- known[is.na(call$data$gold[known])]
    expected <- d
    expected$gold[hidden[[j]]] <- NA_real_
    expect_identical(call$data, expected)
    expect_setequal(scored$rows$row[scored$rows$fold == j], hidden[[j]])
  }
  expect_identical(sort(unlist(hidden, use.names = FALSE)), known)
  expect_identical(sort(scored$rows$row), known)
  expect_identical(fit$original, d)
})
