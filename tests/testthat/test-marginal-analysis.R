marginal_public_fixture <- local({
  cache <- NULL
  function() {
    if (!is.null(cache)) return(cache)
    set.seed(822301)
    site <- rep(seq_len(24L), each = 10L)
    a <- stats::rnorm(24L)
    x <- a[site] + stats::rnorm(length(site))
    u <- matrix(stats::rnorm(48L), 24L, 2L) %*% chol(matrix(c(.5, .45, .45, .6), 2L))
    e <- matrix(stats::rnorm(480L), 240L, 2L) %*% chol(matrix(c(.5, .55, .55, .8), 2L))
    target <- .5 + .7 * x + .4 * a[site] + u[site, 1L] + e[, 1L]
    d <- data.frame(recorded = .8 + 1.05 * x + u[site, 2L] + e[, 2L],
      reference = target, x = x, site = as.character(site))
    d$reference[-sample(which(site <= 18L), 60L)] <- NA_real_
    d$recorded[seq(1L, 240L, 3L)] <- NA_real_
    fit <- catmech_reference(d, "recorded", "reference", "x", m = 8,
      seed = 822302, use_recorded = TRUE, cluster = "site", cluster_model = "joint_marginal")
    cache <<- list(fit = fit, data = d)
    cache
  }
})

test_that("the marginal engine integrates with completions and posterior summaries", {
  set.seed(822300)
  before <- .Random.seed
  fixture <- marginal_public_fixture()
  fit <- fixture$fit
  # Fixture construction itself uses RNG; the public call is checked separately.
  expect_identical(fit$original, fixture$data)
  expect_identical(fit$manifest$engine, "reference_marginal_cluster")
  expect_identical(fit$spec$cluster_model, "joint_marginal")
  known <- !is.na(fixture$data$reference)
  for (j in seq_along(fit$completed)) {
    z <- fit$completed[[j]]
    expect_equal(z$recorded[known], fixture$data$reference[known], tolerance = 0)
    expect_identical(z[-1L], fixture$data[-1L])
    expect_true(all(is.finite(z$recorded)))
    expect_equal(mice::complete(fit$mids, j)$recorded, z$recorded, tolerance = 0)
  }
  before <- .Random.seed
  out <- catmech_marginal_summary(fit)
  expect_identical(.Random.seed, before)
  expect_identical(out$term, c("(Intercept)", "x"))
  expect_true(all(out$interpretation == "conditional_joint_marginal_credible_interval"))
  expect_true(all(out$lower < out$estimate & out$upper > out$estimate))
  expect_true(all(out$rhat <= 1.05 & out$ess_bulk >= 400 & out$ess_tail >= 400))
  tr <- fit$model_diagnostics$traces
  expect_equal(out$estimate, vapply(paste0("marginal:", out$term),
    function(name) mean(tr[, , name]), numeric(1)), ignore_attr = TRUE)
  expect_equal(out$lower, vapply(paste0("marginal:", out$term),
    function(name) unname(stats::quantile(tr[, , name], .025)), numeric(1)), ignore_attr = TRUE)
  changed <- fit
  changed$completed <- NULL
  changed$manifest$m <- 900L
  changed$model_diagnostics$selected_draws <- NULL
  expect_identical(catmech_marginal_summary(changed), out)
  expect_error(catmech_joint_summary(fit), "cluster_model='joint'", fixed = TRUE)
  pooled <- catmech_pool_lm(fit, recorded ~ x)
  expect_true(all(is.finite(pooled$estimate)))
  expect_true(all(pooled$clusters == 24L))
})

test_that("marginal summaries refuse inconsistent metadata and wrong coefficients", {
  fit <- marginal_public_fixture()$fit
  expect_error(catmech_marginal_summary(list()), "Expected")
  for (level in list(0, 1, NA_real_, c(.8, .9), .9 + 0i))
    expect_error(catmech_marginal_summary(fit, level), "level")
  z <- fit; z$model_diagnostics$convergence_pass <- FALSE
  expect_error(catmech_marginal_summary(z), class = "catmech_convergence_error")
  z <- fit; z$model_diagnostics$traces[, , "marginal:x"] <-
    z$model_diagnostics$traces[, , "marginal:x"] + .1
  expect_error(catmech_marginal_summary(z), "conflict")
  z <- fit; z$model_diagnostics$traces[1, 1, "sigma:T:T"] <- -1
  expect_error(catmech_marginal_summary(z), "positive definite")
  z <- fit; z$spec$joint_transform$response_scale$T <- 2
  expect_error(catmech_marginal_summary(z), "inconsistent")
  z <- fit; z$model_diagnostics$iterations_used <- 200L
  expect_error(catmech_marginal_summary(z), "Malformed")
  z <- fit; z$model_diagnostics$traces[1, 1, 1] <- NA_real_
  expect_error(catmech_marginal_summary(z), "Malformed")
})

test_that("marginal public input restrictions and all-reference shortcut are explicit", {
  fixture <- marginal_public_fixture()
  d <- fixture$data
  call <- function(data, ...) catmech_reference(data, "recorded", "reference", "x",
    use_recorded = TRUE, cluster = "site", cluster_model = "joint_marginal", ...)
  z <- d; z$x <- factor(z$site)
  expect_error(call(z), "numeric|continuous")
  z <- d; z$x[1L] <- NA_real_
  expect_error(call(z), "complete|missing|finite")
  z <- d; z$x <- 1
  expect_error(call(z), "variation|constant|rank")
  expect_error(catmech_reference(d, "recorded", "reference", "x",
    cluster_model = "joint_marginal"), "requires|require")
  z <- d; z$reference <- seq_len(nrow(d)) / 10
  set.seed(822310); before <- .Random.seed
  all_reference <- call(z, m = 2, seed = 822311)
  expect_identical(.Random.seed, before)
  expect_equal(all_reference$completed[[1L]]$recorded, z$reference, tolerance = 0)
  expect_error(catmech_marginal_summary(all_reference), "all-reference")
})
