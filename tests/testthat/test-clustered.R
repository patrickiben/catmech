cluster_fixture <- function(clusters = 30L, rows = 3L, audited = 12L, seed = 1L) {
  set.seed(seed)
  id <- rep(seq_len(clusters), each = rows)
  size <- rep(stats::rlnorm(clusters, 3, .5), each = rows)
  x <- size + stats::rnorm(clusters * rows, 0, 2)
  truth <- 1 + .5 * x + rep(stats::rnorm(clusters, 0, 2), each = rows) + stats::rnorm(clusters * rows)
  wrong <- stats::runif(clusters * rows) < .3
  recorded <- truth + wrong * stats::rnorm(clusters * rows, 1, 3)
  picked <- id %in% sample.int(clusters, audited)
  data.frame(id = id, group = rep(rep(c("a", "b"), length.out = clusters), each = rows), x = x,
             recorded = recorded, gold = ifelse(picked, truth, NA_real_), truth = truth)
}
manual_cr1 <- function(fit, ids) {
  x <- stats::model.matrix(fit); n <- nrow(x); k <- ncol(x); g <- length(unique(ids))
  bread <- solve(crossprod(x))
  meat <- crossprod(rowsum(x * stats::residuals(fit), ids))
  bread %*% meat %*% bread * (g / (g - 1)) * ((n - 1) / (n - k))
}

test_that("cluster needs use_recorded and valid columns", {
  d <- cluster_fixture()
  expect_error(catmech_reference(d, "recorded", "gold", "x", cluster = "id"), "only with use_recorded")
  expect_error(catmech_reference(d, "recorded", "gold", "x", use_recorded = TRUE, cluster = "nope"), "existing column")
  expect_error(catmech_reference(d, "recorded", "gold", "x", use_recorded = TRUE, cluster = "x"), "must differ")
  e <- d; e$id[1] <- NA
  expect_error(catmech_reference(e, "recorded", "gold", "x", use_recorded = TRUE, cluster = "id"), "missing values")
  expect_error(catmech_reference(d, "recorded", "gold", "x", use_recorded = TRUE, cluster = "id", cluster_strata = "nope"), "existing columns")
  e <- d; e$group[1] <- "b"; e$group[2] <- "a"
  expect_error(catmech_reference(e, "recorded", "gold", "x", use_recorded = TRUE, cluster = "id", cluster_strata = "group"), "one stratum")
  expect_error(catmech_reference(d, "recorded", "gold", "x", use_recorded = TRUE, cluster_strata = "group"), "cluster must name")
})

test_that("clustered imputation keeps references, fills the rest and is reproducible", {
  d <- cluster_fixture()
  fit <- catmech_reference(d, "recorded", "gold", "x", use_recorded = TRUE, cluster = "id", m = 5, seed = 7)
  known <- !is.na(d$gold)
  for (z in fit$completed) {
    expect_equal(z$recorded[known], d$gold[known])
    expect_false(anyNA(z$recorded))
    expect_identical(z$x, d$x)
  }
  values <- vapply(fit$completed, function(z) z$recorded[!known], numeric(sum(!known)))
  expect_gt(min(apply(values, 1, stats::sd)), 0)
  again <- catmech_reference(d, "recorded", "gold", "x", use_recorded = TRUE, cluster = "id", m = 5, seed = 7)
  expect_identical(lapply(fit$completed, `[[`, "recorded"), lapply(again$completed, `[[`, "recorded"))
  set.seed(11); before <- stats::runif(1)
  set.seed(11); catmech_reference(d, "recorded", "gold", "x", use_recorded = TRUE, cluster = "id", m = 2, seed = 3)
  expect_identical(stats::runif(1), before)
  expect_identical(fit$spec$cluster, "id")
  expect_identical(fit$spec$audited_clusters, 12L)
  expect_identical(fit$spec$cluster_uncertainty, "audited_cluster_bootstrap_each_imputation")
})

test_that("a fit without cluster records nothing about clusters", {
  d <- cluster_fixture()
  fit <- catmech_reference(d, "recorded", "gold", "x", use_recorded = TRUE, m = 2, seed = 7)
  expect_false(any(grepl("cluster", names(unclass(fit$spec)))))
  expect_false("clusters" %in% names(catmech_pool_lm(fit, recorded ~ x)))
})

test_that("clustered imputation refuses what it cannot support", {
  d <- cluster_fixture(audited = 4L)
  expect_error(catmech_reference(d, "recorded", "gold", "x", use_recorded = TRUE, cluster = "id", m = 2), "Fewer than five")
  d <- cluster_fixture()
  d$copy <- d$recorded
  expect_error(catmech_reference(d, "recorded", "gold", "copy", use_recorded = TRUE, cluster = "id", m = 2), "full column rank")
  d <- cluster_fixture(); d$recorded[!is.na(d$gold)] <- NA
  expect_error(catmech_reference(d, "recorded", "gold", "x", use_recorded = TRUE, cluster = "id", m = 2), "both recorded and reference")
})

test_that("clustered imputation handles missing recorded values and full agreement", {
  d <- cluster_fixture()
  gap <- which(is.na(d$gold))[1:5]
  d$recorded[gap] <- NA
  fit <- catmech_reference(d, "recorded", "gold", "x", use_recorded = TRUE, cluster = "id", m = 3, seed = 2)
  expect_false(any(vapply(fit$completed, function(z) anyNA(z$recorded), logical(1))))
  d <- cluster_fixture()
  d$recorded[!is.na(d$gold)] <- d$gold[!is.na(d$gold)]
  fit <- catmech_reference(d, "recorded", "gold", "x", use_recorded = TRUE, cluster = "id", m = 2, seed = 2)
  expect_identical(fit$spec$recorded_fit, "recorded_kept_all_calibration_rows_agree")
  expect_equal(fit$completed[[1]]$recorded, d$recorded)
})

test_that("strata keep every stratum in the resample, including a single audited cluster", {
  d <- cluster_fixture(clusters = 40L, audited = 40L)
  d$group <- ifelse(d$id <= 20L, "a", "b")
  d$gold[d$id > 8L & d$id <= 20L] <- NA
  d$gold[d$id > 21L] <- NA
  fit <- catmech_reference(d, "recorded", "gold", c("x", "group"), use_recorded = TRUE, cluster = "id",
                           cluster_strata = "group", m = 3, seed = 4)
  expect_length(fit$completed, 3L)
  expect_identical(fit$spec$cluster_strata, "group")
})

test_that("the pooled cluster-robust variance matches a manual CR1 calculation", {
  d <- cluster_fixture(audited = 30L)
  fit <- catmech_reference(d, "recorded", "gold", "x", use_recorded = TRUE, cluster = "id", m = 2, seed = 1)
  pooled <- catmech_pool_lm(fit, recorded ~ x)
  model <- stats::lm(gold ~ x, data = d)
  v <- manual_cr1(model, d$id)
  expect_equal(pooled$estimate, unname(stats::coef(model)))
  expect_equal(pooled$std.error, sqrt(diag(v)), ignore_attr = TRUE)
  expect_equal(pooled$df, rep(29, 2))
  expect_equal(pooled$clusters, rep(30L, 2))
  expect_true(all(endsWith(pooled$interpretation, "_cluster_robust")))
})

test_that("one row per cluster gives the unequal-variance form, and Rubin's rules use clusters minus one", {
  d <- cluster_fixture(clusters = 60L, rows = 1L, audited = 60L)
  fit <- catmech_reference(d, "recorded", "gold", "x", use_recorded = TRUE, m = 2, seed = 1)
  pooled <- catmech_pool_lm(fit, recorded ~ x, cluster = "id")
  model <- stats::lm(gold ~ x, data = d)
  x <- stats::model.matrix(model); e <- stats::residuals(model); bread <- solve(crossprod(x))
  v <- bread %*% crossprod(x * e) %*% bread * (60 / (60 - 2))
  expect_equal(pooled$std.error, sqrt(diag(v)), ignore_attr = TRUE)
  d <- cluster_fixture()
  fit <- catmech_reference(d, "recorded", "gold", "x", use_recorded = TRUE, cluster = "id", m = 4, seed = 9)
  pooled <- catmech_pool_lm(fit, recorded ~ x)
  fits <- lapply(fit$completed, function(z) stats::lm(recorded ~ x, data = z))
  q <- vapply(fits, function(f) unname(stats::coef(f)["x"]), numeric(1))
  u <- vapply(seq_along(fits), function(j) manual_cr1(fits[[j]], fit$completed[[j]]$id)["x", "x"], numeric(1))
  expected <- mice::pool.scalar(q, u, n = 30, k = 1)
  row <- pooled[pooled$term == "x", ]
  expect_equal(row$estimate, expected$qbar)
  expect_equal(row$std.error, sqrt(expected$t))
  expect_equal(row$df, expected$df)
})

test_that("pooling can use a cluster column for any fit, and refuses a bad one", {
  d <- cluster_fixture()
  fit <- catmech_reference(d, "recorded", "gold", "x", use_recorded = TRUE, m = 3, seed = 5)
  expect_true("clusters" %in% names(catmech_pool_lm(fit, recorded ~ x, cluster = "id")))
  expect_error(catmech_pool_lm(fit, recorded ~ x, cluster = "nope"), "must name a column")
  e <- d; e$one <- 1L
  fit <- catmech_reference(e, "recorded", "gold", "x", use_recorded = TRUE, m = 2, seed = 5)
  expect_error(catmech_pool_lm(fit, recorded ~ x, cluster = "one"), "at least two clusters")
  e$na <- ifelse(e$id == 1L, NA, e$id)
  fit <- catmech_reference(e, "recorded", "gold", "x", use_recorded = TRUE, m = 2, seed = 5)
  expect_error(catmech_pool_lm(fit, recorded ~ x, cluster = "na"), "missing values")
})

test_that("the hold-out check keeps the cluster option when it refits", {
  d <- cluster_fixture(clusters = 60L, audited = 40L)
  fit <- catmech_reference(d, "recorded", "gold", "x", use_recorded = TRUE, cluster = "id", m = 3, seed = 6)
  check <- suppressWarnings(catmech_holdout_check(fit, folds = 3L, m = 20L, seed = 6))
  expect_false(is.null(check$status))
  expect_length(check$failures, 0L)
})
