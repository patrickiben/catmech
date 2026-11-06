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

test_that("marginal refits preserve posterior traces across m and restore caller RNG", {
  fixture <- marginal_public_fixture()
  fit8 <- fixture$fit
  set.seed(822315L)
  before <- .Random.seed
  fit16 <- catmech_reference(fixture$data, "recorded", "reference", "x", m = 16,
    seed = 822302, use_recorded = TRUE, cluster = "site", cluster_model = "joint_marginal")
  expect_identical(.Random.seed, before)
  expect_identical(fit16$model_diagnostics$traces, fit8$model_diagnostics$traces)
  expect_identical(fit16$model_diagnostics$attempts, fit8$model_diagnostics$attempts)
  expect_identical(catmech_marginal_summary(fit16), catmech_marginal_summary(fit8))
  expect_length(fit16$completed, 16L)
  expect_identical(fit16$manifest$imputed_rows, which(is.na(fixture$data$reference)))
  expect_equal(mice::complete(fit16$mids, 16L)$recorded,
    fit16$completed[[16L]]$recorded, tolerance = 0)
  known <- !is.na(fixture$data$reference)
  expect_equal(fit16$completed[[16L]]$recorded[known], fixture$data$reference[known], tolerance = 0)
  expect_error(catmech_reference(fixture$data, "recorded", "reference", "x", m = 16001,
    seed = 822302, use_recorded = TRUE, cluster = "site", cluster_model = "joint_marginal"),
    "m cannot exceed")
  expect_identical(.Random.seed, before)
  expect_error(catmech_reference(fixture$data, "recorded", "reference", "x", m = 2,
    seed = 822302, use_recorded = TRUE, cluster = "site", cluster_model = "joint_marginal",
    joint_control = list(iterations = 200L, max_iterations = 100L)), "max_iterations")

  directory <- tempfile("catmech-marginal-result-")
  on.exit(unlink(directory, recursive = TRUE), add = TRUE)
  write_catmech(fit16, directory)
  restored <- readRDS(file.path(directory, "result.rds"))
  expect_identical(restored$model_diagnostics$traces, fit16$model_diagnostics$traces)
  expect_identical(catmech_marginal_summary(restored), catmech_marginal_summary(fit16))
  expect_identical(restored$original, fixture$data)
  manifest <- jsonlite::read_json(file.path(directory, "manifest.json"))
  expect_identical(manifest$engine, "reference_marginal_cluster")
  expect_identical(manifest$assumptions$cluster_model, "joint_marginal")
})

test_that("marginal holdout refits retain the engine without leaking hidden references", {
  fixture <- marginal_public_fixture()
  fit <- fixture$fit
  calls <- list()
  mock <- function(data, value, reference, predictors, method, m, seed,
    use_recorded, cluster, cluster_model, joint_control, ...) {
    calls[[length(calls) + 1L]] <<- list(data = data, value = value,
      reference = reference, predictors = predictors, method = method, m = m, seed = seed,
      use_recorded = use_recorded, cluster = cluster, cluster_model = cluster_model,
      joint_control = joint_control, extra = list(...))
    # Deterministic dummy predictions check forwarding and reference hiding,
    # not predictive accuracy. They do not access the cached original truth.
    list(completed = lapply(seq_len(m), function(j) {
      out <- data
      out[[value]] <- ifelse(is.na(data[[reference]]), 1 + .4 * data$x + j / 10, data[[reference]])
      out
    }))
  }
  testthat::local_mocked_bindings(catmech_reference = mock, .package = "catmech")
  set.seed(822319L); before <- .Random.seed
  scored <- catmech_holdout_check(fit, folds = 3L, m = 20L, seed = 822320L)
  expect_identical(.Random.seed, before)
  expect_length(calls, 3L)
  expect_length(scored$failures, 0L)
  known <- which(!is.na(fixture$data$reference))
  hidden <- vector("list", length(calls))
  for (j in seq_along(calls)) {
    call <- calls[[j]]
    expect_identical(call$cluster_model, "joint_marginal")
    expect_identical(call$joint_control, fit$spec$joint_control)
    expect_identical(call$use_recorded, TRUE)
    expect_identical(call$cluster, "site")
    expect_identical(call$predictors, "x")
    expect_identical(call$seed, 822320L + j)
    expect_length(call$extra, 0L)
    expect_identical(call$data[setdiff(names(fixture$data), "reference")],
      fixture$data[setdiff(names(fixture$data), "reference")])
    hidden[[j]] <- known[is.na(call$data$reference[known])]
    expect_setequal(scored$rows$row[scored$rows$fold == j], hidden[[j]])
  }
  expect_identical(sort(unlist(hidden, use.names = FALSE)), known)
  expect_identical(fit$original, fixture$data)
})

test_that("intercept-only marginal fits predict a site with no observed measurements", {
  d <- marginal_public_fixture()$data
  d <- rbind(d, data.frame(recorded = c(NA_real_, NA_real_),
    reference = c(NA_real_, NA_real_), x = c(0, 1), site = "unobserved"))
  set.seed(822321L); before <- .Random.seed
  fit <- catmech_reference(d, "recorded", "reference", character(), m = 4L,
    seed = 822322L, use_recorded = TRUE, cluster = "site", cluster_model = "joint_marginal")
  expect_identical(.Random.seed, before)
  out <- catmech_marginal_summary(fit)
  expect_identical(out$term, "(Intercept)")
  expect_true(all(is.finite(unlist(out[c("estimate", "posterior_sd", "lower", "upper")]))))
  expect_equal(ncol(fit$model_diagnostics$sampler_input$y), 2L)
  expect_equal(nrow(fit$model_diagnostics$sampler_input$y),
    sum(rowSums(!is.na(d[c("recorded", "reference")])) > 0L))
  new <- which(d$site == "unobserved")
  expect_false(any(new %in% fit$model_diagnostics$input_rows))
  draws <- vapply(fit$completed, function(z) z$recorded[new], numeric(length(new)))
  expect_true(all(is.finite(draws)))
  expect_true(all(apply(draws, 1L, stats::var) > 0))
  expect_equal(mice::complete(fit$mids, 4L)$recorded,
    fit$completed[[4L]]$recorded, tolerance = 0)
  known <- !is.na(d$reference)
  expect_equal(fit$completed[[4L]]$recorded[known], d$reference[known], tolerance = 0)
  expect_identical(fit$original, d)
})

test_that("conditional-regression prior controls are explicit and engine-specific", {
  fixture <- marginal_public_fixture()
  expect_identical(fixture$fit$spec$joint_control$regression_scale, 10)
  expect_identical(fixture$fit$spec$prior$covariance_family, "conditional_regression")
  expect_identical(fixture$fit$spec$prior$conditional_regression$proper, TRUE)
  call <- function(control, engine = "joint_marginal") catmech_reference(
    fixture$data, "recorded", "reference", "x", m = 2L, seed = 822323L,
    use_recorded = TRUE, cluster = "site", cluster_model = engine, joint_control = control)
  for (value in list(0, -1, NA_real_, Inf, c(5, 10), "wide", 1 + 0i, 1e200, 1e-200))
    expect_error(call(list(regression_scale = value)), "regression_scale")
  expect_error(call(list(regression_scale = 10), "joint"), "Unknown joint_control")
})

test_that("paired residual support excludes exact measurement degeneracy", {
  fixture <- marginal_public_fixture()
  d <- fixture$data
  call <- function(data) catmech_reference(data, "recorded", "reference", "x", m = 2L,
    seed = 822324L, use_recorded = TRUE, cluster = "site", cluster_model = "joint_marginal")
  paired <- !is.na(d$recorded) & !is.na(d$reference)
  equal <- d; equal$reference[paired] <- equal$recorded[paired]
  expect_error(call(equal), "paired|rank|variation|covariance")
  by_site <- d
  by_site$reference[paired] <- by_site$recorded[paired] +
    sin(as.integer(by_site$site[paired])) + .4 * by_site$x[paired]
  # Global residual scatter is full rank; the exact relation appears only
  # after accounting for site effects as well as X.
  expect_equal(qr(cbind(1, by_site$x[paired], by_site$reference[paired],
    by_site$recorded[paired]))$rank, 4L)
  expect_error(call(by_site), "paired|rank|variation|covariance")
  by_site$reference <- 1000 - 100 * by_site$reference
  by_site$recorded <- -20 + .01 * by_site$recorded
  by_site$x <- 20 - 3 * by_site$x
  expect_error(call(by_site), "paired|rank|variation|covariance")

  weak <- d; weak$reference[] <- NA_real_
  for (site in as.character(seq_len(18L))) {
    rows <- which(weak$site == site)
    first <- rows[which(!is.na(weak$recorded[rows]))[1L]]
    second <- rows[rows != first][1L]
    weak$recorded[second] <- NA_real_
    selected <- c(first, second)
    weak$reference[selected] <- .4 * weak$x[selected] + cos(selected) + sin(as.integer(site))
  }
  # The supported-data policy deliberately requires two paired within-site
  # residual degrees of freedom; this does not claim all weaker designs are
  # mathematically unidentified under every possible covariance prior.
  expect_error(call(weak), "paired|degrees|replication|support")
})

test_that("paired support checks are stable under measurement and predictor recoding", {
  d <- marginal_public_fixture()$data
  original <- catmech_marginal_cluster_prepare(d, "recorded", d$reference, "x", "site")
  recoded <- d
  recoded$reference <- 1000 - 100 * d$reference
  recoded$recorded <- -20 + .01 * d$recorded
  recoded$x <- 20 - 3 * d$x
  transformed <- catmech_marginal_cluster_prepare(recoded, "recorded", recoded$reference,
    "x", "site")
  expect_identical(transformed$rows, original$rows)
  expect_identical(transformed$input$subj, original$input$subj)
  expect_equal(transformed$input$y, sweep(original$input$y, 2L, c(-1, 1, -1), "*"),
    tolerance = 1e-10)
})

test_that("paired calibration at one site retains sufficient residual support", {
  set.seed(171100L)
  site <- rep(seq_len(9L), each = 8L)
  x <- stats::rnorm(length(site))
  d <- data.frame(recorded = 1.1 * x + stats::rnorm(length(x)),
    reference = .5 + .7 * x + stats::rnorm(length(x)), x = x, site = as.character(site))
  d$reference[site > 5L] <- NA_real_
  d$recorded[site %in% 2:5] <- NA_real_
  # Both responses have five observed sites, while the eight paired rows occur
  # at one site. Two independent residual measurement directions are available.
  prepared <- catmech_marginal_cluster_prepare(d, "recorded", d$reference, "x", "site")
  expect_identical(prepared$paired_support$paired_sites, 1L)
  expect_equal(prepared$paired_support$residual_df, 6L)
  expect_true(all(prepared$paired_support$standardized_residual_singular_values >
    prepared$paired_support$minimum_singular_value))
  expect_identical(prepared$rows, seq_len(nrow(d)))
})


test_that("marginal numerical budget extends without changing conditional defaults", {
  expect_identical(catmech:::catmech_marginal_cluster_control()$max_iterations, 64000L)
  expect_identical(catmech:::catmech_joint_cluster_control()$max_iterations, 32000L)
  expect_identical(catmech:::catmech_marginal_cluster_control(list(max_iterations=32000L))$max_iterations, 32000L)
})

test_that("marginal summary supports predictor-role formulas and rejects malformed targets", {
  fixture <- marginal_public_fixture()
  fit <- fixture$fit
  default_summary <- catmech_marginal_summary(fit)
  expect_identical(catmech_marginal_summary(fit, formula = recorded ~ x), default_summary)
  expect_identical(catmech_marginal_summary(fit, response = "recorded", covariates = "x"), default_summary)
  rev_summary <- catmech_marginal_summary(fit, formula = x ~ recorded)
  expect_identical(rev_summary, catmech_marginal_summary(fit, response = "x"))
  expect_identical(rev_summary, catmech_marginal_summary(fit, response = "x", covariates = "recorded"))
  expect_identical(rev_summary$term, c("(Intercept)", "recorded"))
  expect_true(all(is.finite(unlist(rev_summary[c("estimate", "posterior_sd", "lower", "upper")]))))
  expect_gt(rev_summary$estimate[2L], 0)
  int_summary <- catmech_marginal_summary(fit, formula = x ~ 1)
  expect_identical(int_summary$term, "(Intercept)")
  expect_warning(catmech_pool_lm(fit, x ~ recorded),
    "catmech_marginal_summary\\(result, formula = \\.\\.\\.\\)")

  expect_error(catmech_marginal_summary(fit, formula = x ~ recorded, response = "x"), "either formula or")
  expect_error(catmech_marginal_summary(fit, formula = ~ recorded), "two-sided")
  expect_error(catmech_marginal_summary(fit, formula = I(x + 1) ~ recorded), "single modeled variable")
  expect_error(catmech_marginal_summary(fit, formula = unknown ~ recorded), "audited value or one of spec\\$predictors")
  expect_error(catmech_marginal_summary(fit, formula = x ~ recorded - 1), "intercept")
  expect_error(catmech_marginal_summary(fit, formula = x ~ recorded + offset(recorded)), "offset")
  expect_error(catmech_marginal_summary(fit, formula = x ~ I(recorded^2)), "plain variable names")
  expect_error(catmech_marginal_summary(fit, formula = x ~ unknown), "audited value or one of spec\\$predictors")
  expect_error(catmech_marginal_summary(fit, formula = x ~ recorded:x), "additive")
  expect_error(catmech_marginal_summary(fit, formula = x ~ x + recorded), "cannot also appear as a covariate")
  expect_error(catmech_marginal_summary(fit, response = "unknown"), "audited value or one of spec\\$predictors")
  expect_error(catmech_marginal_summary(fit, response = "x", covariates = "x"), "remaining modeled variables")
  expect_error(catmech_marginal_summary(fit, response = "x", covariates = c("recorded", "recorded")), "unique")
})

test_that("adaptive chain continuation matches a single run of the same total length", {
  d <- marginal_public_fixture()$data
  doubled <- catmech_reference(d, "recorded", "reference", "x", m = 5L,
    seed = 822302L, use_recorded = TRUE, cluster = "site",
    cluster_model = "joint_marginal",
    joint_control = list(chains = 4L, burnin = 500L, iterations = 1000L, max_iterations = 8000L))
  direct <- catmech_reference(d, "recorded", "reference", "x", m = 5L,
    seed = 822302L, use_recorded = TRUE, cluster = "site",
    cluster_model = "joint_marginal",
    joint_control = list(chains = 4L, burnin = 500L, iterations = 8000L, max_iterations = 8000L))
  expect_identical(nrow(doubled$model_diagnostics$attempts), 4L)
  expect_identical(doubled$model_diagnostics$iterations_used, 8000L)
  expect_identical(doubled$model_diagnostics$traces, direct$model_diagnostics$traces)
  expect_identical(doubled$model_diagnostics$scale_interweaving, direct$model_diagnostics$scale_interweaving)
  expect_identical(doubled$completed, direct$completed)
  expect_identical(catmech_marginal_summary(doubled), catmech_marginal_summary(direct))
  expect_identical(catmech_marginal_summary(doubled, formula = x ~ recorded),
    catmech_marginal_summary(direct, formula = x ~ recorded))
})

test_that("catmech_sensitivity supports marginal_summary across prior scenarios", {
  fixture <- marginal_public_fixture()
  fit1 <- fixture$fit
  fit_err <- simpleError("Sampler failed convergence")
  res <- catmech_sensitivity(list(default = fit1, failed = fit_err),
    term = "x", method = "marginal_summary")
  expect_identical(nrow(res$table), 2L)
  expect_identical(res$table$status, c("ok", "failed"))
  expect_identical(res$status, "unresolved_incomplete_set")
  expect_true(is.finite(res$table$estimate[1L]))
  expect_true(is.na(res$table$estimate[2L]))
  expect_identical(res$table$direction[1L], "above")

  rev_res <- catmech_sensitivity(list(default = fit1), formula = x ~ recorded,
    term = "recorded", method = "marginal_summary")
  expect_identical(rev_res$table$status, "ok")
  expect_true(is.finite(rev_res$table$estimate[1L]))

  expect_error(catmech_sensitivity(list(a = fit1), formula = "not_formula", term = "x",
    method = "marginal_summary"), "formula must be an R formula or NULL")
})

test_that("catmech_marginal_sensitivity_grid computes 2-axis statistics and formats cleanly", {
  fixture <- marginal_public_fixture()
  fit1 <- fixture$fit

  # Error handling on invalid arguments
  expect_error(catmech_marginal_sensitivity_grid("not_a_fit"),
               "fit must be a catmech_result or a named list")
  expect_error(catmech_marginal_sensitivity_grid(fit1, prior_scales = -1),
               "prior_scales must be a vector of positive numbers")
  expect_error(catmech_marginal_sensitivity_grid(fit1, regression_scales = 0),
               "regression_scales must be a vector of positive numbers")
  expect_error(catmech_marginal_sensitivity_grid(fit1, conf_level = 1.5),
               "conf_level must be between zero and one")

  # Precomputed list of results
  list_grid <- catmech_marginal_sensitivity_grid(
    fit = list(prior_1_reg_10 = fit1, prior_2_reg_10 = fit1, prior_1_reg_20 = fit1),
    term = "x", prior_scales = c(1, 2), regression_scales = c(10, 20)
  )
  expect_s3_class(list_grid, "catmech_marginal_sensitivity_grid")
  expect_true(list_grid$complete)
  expect_identical(list_grid$term, "x")
  expect_equal(list_grid$statistic, 0, tolerance = 1e-10)
  expect_equal(list_grid$axis_summary$mean_abs_change, c(0, 0), tolerance = 1e-10)
  expect_true(is.finite(list_grid$mcse_mean))

  # Test print method
  printed <- capture.output(print(list_grid))
  expect_true(any(grepl("Catmech Joint Marginal Prior Sensitivity Grid", printed)))
  expect_true(any(grepl("Target term: x", printed)))
  expect_true(any(grepl("Primary sensitivity statistic", printed)))
  expect_true(any(grepl("complete: TRUE", printed)))

  # Grid with refitting on small grid
  grid_refit <- catmech_marginal_sensitivity_grid(
    fit = fit1, term = "x",
    prior_scales = c(1, 2), regression_scales = c(10, 20)
  )
  expect_true(grid_refit$complete)
  expect_true(is.finite(grid_refit$statistic))
  expect_true(grid_refit$statistic >= 0)
  expect_identical(nrow(grid_refit$axis_summary), 2L)
  expect_true(all(is.finite(grid_refit$axis_summary$mean_abs_change)))

  # Test plot method
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_silent(plot(grid_refit))
})
