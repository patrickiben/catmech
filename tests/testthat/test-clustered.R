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
  expect_error(catmech_reference(d, "recorded", "gold", "x", use_recorded = TRUE, cluster = "id", cluster_strata = "group"), "unused argument")
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
  expect_identical(fit$spec$cluster_uncertainty, "random_intercept_drawn_each_imputation")
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

test_that("strata enter as predictors, including one with a single audited cluster", {
  d <- cluster_fixture(clusters = 40L, audited = 40L)
  d$group <- ifelse(d$id <= 20L, "a", "b")
  d$gold[d$id > 8L & d$id <= 20L] <- NA
  d$gold[d$id > 21L] <- NA
  fit <- catmech_reference(d, "recorded", "gold", c("x", "group"), use_recorded = TRUE, cluster = "id", m = 3, seed = 4)
  expect_length(fit$completed, 3L)
  expect_false("cluster_strata" %in% names(unclass(fit$spec)))
  expect_false(any(vapply(fit$completed, function(z) anyNA(z$recorded), logical(1))))
})

test_that("variance components match a one-way analysis of variance", {
  set.seed(31)
  groups <- 12L; each <- 4L; n <- groups * each
  cl <- factor(rep(seq_len(groups), each = each))
  y <- rep(stats::rnorm(groups, 0, 2), each = each) + stats::rnorm(n)
  moments <- catmech:::catmech_reference_cluster_moments(matrix(1, n, 1L), y, cl)
  means <- tapply(y, cl, mean)
  ms_between <- each * sum((means - mean(y))^2) / (groups - 1L)
  ms_within <- sum((y - means[cl])^2) / (n - groups)
  expect_equal(moments$df_within, n - groups)
  expect_equal(moments$df_between, groups - 1L)
  expect_equal(moments$sse_within / moments$df_within, ms_within)
  expect_equal(moments$ss_between / moments$df_between, ms_between)
  expect_equal(moments$trace, each * (groups - 1L))
  tau2 <- moments$df_between * (moments$ss_between / moments$df_between - moments$sse_within / moments$df_within) / moments$trace
  expect_equal(tau2, (ms_between - ms_within) / each)
})

test_that("the component trace holds with predictors and unequal clusters", {
  set.seed(32)
  cl <- factor(rep(1:8, times = c(2, 3, 4, 2, 5, 3, 2, 4)))
  n <- length(cl)
  X <- cbind(1, stats::rnorm(n))
  y <- as.vector(X %*% c(1, .5)) + stats::rnorm(8)[cl] + stats::rnorm(n)
  moments <- catmech:::catmech_reference_cluster_moments(X, y, cl)
  Z <- stats::model.matrix(~ cl - 1)
  M <- diag(n) - X %*% solve(crossprod(X), t(X))
  expect_equal(moments$trace, sum(diag(t(Z) %*% M %*% Z)))
  expect_equal(moments$sse_within, sum(stats::lm.fit(cbind(X, Z), y)$residuals^2))
})

test_that("generalized least squares matches explicit matrix algebra", {
  set.seed(33)
  cl <- factor(rep(1:7, times = c(1, 2, 3, 4, 2, 3, 5)))
  n <- length(cl)
  X <- cbind(1, stats::rnorm(n), stats::rexp(n))
  y <- as.vector(X %*% c(1, .5, -.2)) + stats::rnorm(7)[cl] + stats::rnorm(n)
  sigma2 <- 1.3; tau2 <- .7
  V <- sigma2 * diag(n) + tau2 * tcrossprod(stats::model.matrix(~ cl - 1))
  W <- solve(V)
  expected <- solve(t(X) %*% W %*% X)
  gls <- catmech:::catmech_reference_cluster_gls(X, y, cl, sigma2, tau2)
  expect_equal(unname(gls$coefficients), as.vector(expected %*% t(X) %*% W %*% y))
  R_inverse <- backsolve(qr.R(gls$qr), diag(3L))
  expect_equal(sigma2 * tcrossprod(R_inverse), expected, ignore_attr = TRUE)
  ols <- catmech:::catmech_reference_cluster_gls(X, y, cl, sigma2, 0)
  expect_equal(unname(ols$coefficients), unname(stats::lm.fit(X, y)$coefficients))
})

test_that("imputed rows carry the cluster correlation of the audited residuals", {
  set.seed(34)
  clusters <- 60L; rows <- 3L
  id <- rep(seq_len(clusters), each = rows)
  x <- stats::rnorm(clusters * rows)
  truth <- x + stats::rnorm(clusters * rows)
  shift <- rep(stats::rnorm(clusters, 0, 3), each = rows)
  d <- data.frame(id = id, x = x, recorded = truth + shift + stats::rnorm(clusters * rows, 0, .2))
  d$gold <- ifelse(id <= 30L, truth, NA_real_)
  fit <- catmech_reference(d, "recorded", "gold", "x", use_recorded = TRUE, cluster = "id", m = 200, seed = 5)
  values <- vapply(fit$completed, function(z) z$recorded[z$id == 40L], numeric(rows))
  noise <- values - rowMeans(values)
  correlation <- mean(stats::cor(t(noise))[upper.tri(diag(rows))])
  expect_gt(correlation, .5)
  independent <- d
  independent$recorded <- truth + stats::rnorm(clusters * rows, 0, 1)
  fit <- catmech_reference(independent, "recorded", "gold", "x", use_recorded = TRUE, cluster = "id", m = 200, seed = 5)
  values <- vapply(fit$completed, function(z) z$recorded[z$id == 40L], numeric(rows))
  noise <- values - rowMeans(values)
  expect_lt(abs(mean(stats::cor(t(noise))[upper.tri(diag(rows))])), .2)
})

test_that("an audited cluster with unaudited rows keeps its own intercept", {
  set.seed(35)
  clusters <- 40L; rows <- 4L
  id <- rep(seq_len(clusters), each = rows)
  x <- stats::rnorm(clusters * rows)
  effect <- stats::rnorm(clusters, 0, 3)
  effect[1L] <- 12
  truth <- x + effect[id] + stats::rnorm(clusters * rows, 0, .5)
  d <- data.frame(id = id, x = x, recorded = x + stats::rnorm(clusters * rows, 0, .5))
  d$gold <- ifelse(id <= 25L, truth, NA_real_)
  d$gold[d$id == 1L][4L] <- NA_real_
  fit <- catmech_reference(d, "recorded", "gold", "x", use_recorded = TRUE, cluster = "id", m = 100, seed = 6)
  row <- which(d$id == 1L)[4L]
  imputed <- vapply(fit$completed, function(z) z$recorded[row], numeric(1))
  expect_gt(mean(imputed), 9)
  other <- vapply(fit$completed, function(z) z$recorded[which(d$id == 30L)[1L]], numeric(1))
  expect_gt(stats::sd(other), 1.5)
})

test_that("with one row per cluster the intercept variance is zero and the draw still works", {
  d <- cluster_fixture(clusters = 60L, rows = 1L, audited = 30L)
  X <- cbind(1, d$x[!is.na(d$gold)], d$recorded[!is.na(d$gold)])
  moments <- catmech:::catmech_reference_cluster_moments(X, d$gold[!is.na(d$gold)], factor(d$id[!is.na(d$gold)]))
  expect_false(moments$separable)
  fit <- catmech_reference(d, "recorded", "gold", "x", use_recorded = TRUE, cluster = "id", m = 5, seed = 8)
  expect_false(any(vapply(fit$completed, function(z) anyNA(z$recorded), logical(1))))
  values <- vapply(fit$completed, function(z) z$recorded[is.na(d$gold)], numeric(sum(is.na(d$gold))))
  expect_gt(min(apply(values, 1, stats::sd)), 0)
})

test_that("audited rows with no residual variance are refused", {
  d <- cluster_fixture()
  d$gold[!is.na(d$gold)] <- 0
  expect_error(catmech_reference(d, "recorded", "gold", "x", use_recorded = TRUE, cluster = "id", m = 2, seed = 1), "no residual variance")
})

test_that("the pooled cluster-robust variance matches a manual CR1 calculation", {
  d <- cluster_fixture(audited = 30L)
  fit <- catmech_reference(d, "recorded", "gold", "x", use_recorded = TRUE, cluster = "id", m = 2, seed = 1)
  pooled <- catmech_pool_lm(fit, recorded ~ x, cluster_variance = "cr1s")
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
  pooled <- catmech_pool_lm(fit, recorded ~ x, cluster = "id", cluster_variance = "cr1s")
  model <- stats::lm(gold ~ x, data = d)
  x <- stats::model.matrix(model); e <- stats::residuals(model); bread <- solve(crossprod(x))
  v <- bread %*% crossprod(x * e) %*% bread * (60 / (60 - 2))
  expect_equal(pooled$std.error, sqrt(diag(v)), ignore_attr = TRUE)
  d <- cluster_fixture()
  fit <- catmech_reference(d, "recorded", "gold", "x", use_recorded = TRUE, cluster = "id", m = 4, seed = 9)
  pooled <- catmech_pool_lm(fit, recorded ~ x, cluster_variance = "cr1s")
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

test_that("the default cluster variance matches delete-one-cluster refits", {
  d <- cluster_fixture(audited = 30L)
  fit <- catmech_reference(d, "recorded", "gold", "x", use_recorded = TRUE, cluster = "id", m = 2, seed = 1)
  pooled <- catmech_pool_lm(fit, recorded ~ x)
  model <- stats::lm(gold ~ x, data = d)
  full <- stats::coef(model)
  units <- unique(d$id)
  deltas <- t(vapply(units, function(i) stats::coef(stats::lm(gold ~ x, data = d[d$id != i, ])) - full, numeric(2)))
  by_hand <- crossprod(deltas) * (length(units) - 1) / length(units)
  expect_equal(pooled$std.error, sqrt(diag(by_hand)), ignore_attr = TRUE)
  expect_equal(pooled$df, rep(29, 2))
  expect_true(all(pooled$cluster_variance == "cv3"))
  expect_error(catmech_pool_lm(fit, recorded ~ x, cluster_variance = "nope"), "should be one of")
})

test_that("the default cluster variance applies Rubin's rules with clusters minus one", {
  d <- cluster_fixture()
  fit <- catmech_reference(d, "recorded", "gold", "x", use_recorded = TRUE, cluster = "id", m = 4, seed = 9)
  pooled <- catmech_pool_lm(fit, recorded ~ x)
  row <- pooled[pooled$term == "x", ]
  fits <- lapply(fit$completed, function(z) stats::lm(recorded ~ x, data = z))
  q <- vapply(fits, function(f) unname(stats::coef(f)["x"]), numeric(1))
  u <- vapply(seq_along(fits), function(j) catmech:::catmech_cluster_jackknife(fits[[j]], fit$completed[[j]]$id)["x", "x"], numeric(1))
  expected <- mice::pool.scalar(q, u, n = 30, k = 1)
  expect_equal(row$std.error, sqrt(expected$t))
  expect_equal(row$df, expected$df)
})

test_that("a cluster whose rows alone determine a coefficient still gets a variance", {
  d <- data.frame(id = rep(1:6, each = 2), y = c(1, 2, 2, 4, 3, 5, 4, 7, 5, 8, 6, 9), x = c(0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6), z = c(1, 1, rep(0, 10)))
  model <- stats::lm(y ~ x + z, data = d)
  v <- catmech:::catmech_cluster_jackknife(model, d$id)
  expect_true(all(is.finite(v)))
})

test_that("pooling warns when a fit with cluster has the audited value as a predictor", {
  d <- cluster_fixture()
  fit <- catmech_reference(d, "recorded", "gold", "x", use_recorded = TRUE, cluster = "id", m = 2, seed = 3)
  expect_warning(catmech_pool_lm(fit, x ~ recorded), "biased slopes toward zero")
  expect_no_warning(catmech_pool_lm(fit, recorded ~ x))
  plain <- catmech_reference(d, "recorded", "gold", "x", use_recorded = TRUE, m = 2, seed = 3)
  expect_no_warning(catmech_pool_lm(plain, x ~ recorded, cluster = "id"))
})
