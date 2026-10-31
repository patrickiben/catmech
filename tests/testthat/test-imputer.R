imputer_fixture <- function() {
  set.seed(42)
  n <- 60L
  predictor <- stats::rnorm(n)
  value <- 5 + 0.15 * predictor + stats::rnorm(n, sd = 0.25)
  value[1L] <- NA_real_
  list(y = value, ry = !is.na(value), x = cbind(predictor = predictor),
       category = rep("Middle", n),
       bounds = make_interval_bounds(c(-Inf, 3, 7, Inf), c("Low", "Middle", "High")),
       flag = rep("No", n), corruption_probability = 0.8)
}

test_that("a single eligible donor cannot turn into an unrelated row index", {
  set.seed(7)
  index <- replicate(200L, catmech:::pick_pmm_donor(
    c(0, 0, 0, 0, 3), c(FALSE, FALSE, FALSE, FALSE, TRUE), 3, 5L
  ))
  expect_true(all(index == 5L))
  a <- imputer_fixture()
  a$y[5L] <- 9
  a$category[c(1L, 5L)] <- "High"
  set.seed(4)
  result <- do.call(mice.impute.catmech_threshold_corruption, a)
  expect_equal(as.numeric(result), 9)
  expect_equal(attr(result, "audit")$donor_row, 5L)
})

test_that("rank reduction uses independent columns and correct residual df", {
  a <- imputer_fixture()
  duplicated_x <- cbind(a$x, a$x[, 1L] * 2, constant = 1)
  set.seed(3)
  draw <- catmech:::catmech_bayes_lin_draw(a$y, a$ry, duplicated_x)
  expect_equal(draw$rank, 2L)
  expect_equal(draw$df, sum(a$ry) - 2L)
  expect_true(all(is.finite(draw$pred_draw)))
  expect_gt(draw$sigma, 0)
  set.seed(3)
  simple <- catmech:::catmech_bayes_lin_draw(a$y, a$ry, a$x)
  expect_equal(draw$pred_draw, simple$pred_draw, tolerance = 1e-12)
})

test_that("rank reduction cannot invent predictions for an unsupported level", {
  a <- imputer_fixture()
  a$x <- cbind(a$x, unsupported = as.numeric(seq_along(a$y) == 1L))
  expect_error(do.call(mice.impute.catmech_threshold_corruption, a),
               "outside the trusted regression column span")
  a$x[, "unsupported"] <- 0
  expect_silent(do.call(mice.impute.catmech_threshold_corruption, a))
})

test_that("matching compares fitted donor and drawn recipient predictions", {
  a <- imputer_fixture()
  a$n_donors <- 1L
  set.seed(19)
  draw <- catmech:::catmech_bayes_lin_draw(a$y, a$ry, a$x)
  donors <- which(a$ry)
  type_one <- donors[which.min(abs(draw$pred_hat[donors] - draw$pred_draw[1]))]
  type_two <- donors[which.min(abs(draw$pred_draw[donors] - draw$pred_draw[1]))]
  expect_false(type_one == type_two)
  set.seed(19)
  result <- do.call(mice.impute.catmech_threshold_corruption, a)
  expect_equal(attr(result, "audit")$donor_row, type_one)
})

test_that("tied donor predictions are not resolved by input row order", {
  set.seed(28)
  rows <- replicate(100, catmech:::pick_pmm_donor(rep(5, 6), rep(TRUE, 6), 5, 1))
  expect_gt(length(unique(rows)), 1L)
  expect_true(all(rows %in% 1:6))
})

test_that("intercept-only draws match the analytic normal model", {
  y <- c(1, 2, 4, 7, NA_real_)
  x <- matrix(numeric(), nrow = 5L, ncol = 0L)
  expected_mean <- mean(y[1:4])
  set.seed(101)
  expected_sigma <- sqrt(sum((y[1:4] - expected_mean)^2) / stats::rchisq(1, 3))
  expected_draw <- expected_mean + expected_sigma * stats::rnorm(1) / sqrt(4)
  set.seed(101)
  draw <- catmech:::catmech_bayes_lin_draw(y, !is.na(y), x)
  expect_equal(draw$pred_hat, rep(expected_mean, 5L))
  # The QR basis has a negative intercept diagonal. Its normal draw therefore
  # has the opposite sign under the same seed, with the same exact distribution.
  expect_equal(abs(draw$pred_draw - expected_mean),
               rep(abs(expected_draw - expected_mean), 5L))
  expect_equal(draw$sigma, expected_sigma)
  expect_equal(draw$df, 3L)
})

test_that("nonfinite predictors fail before incomplete or arbitrary output", {
  a <- imputer_fixture()
  a$x[1L, 1L] <- NA_real_
  expect_error(do.call(mice.impute.catmech_threshold_corruption, a), "finite predictors")
  a$x[1L, 1L] <- Inf
  expect_error(do.call(mice.impute.catmech_threshold_corruption, a), "finite predictors")
})

test_that("donor scarcity never re-admits flagged observations", {
  a <- imputer_fixture()
  a$category[1:2] <- "High"
  a$y[2] <- 9
  a$flag[2] <- "Yes"
  a$wy <- seq_along(a$y) %in% 1:2
  a$trusted <- rep(TRUE, length(a$y))
  expect_error(do.call(mice.impute.catmech_threshold_corruption, a),
               "Insufficient trusted within-category donors")
  a$fallback <- "normal"
  result <- do.call(mice.impute.catmech_threshold_corruption, a)
  expect_true(result[1] > 7)
  expect_true(attr(result, "audit")$used_fallback[1])
  expect_true(is.na(attr(result, "audit")$donor_row[1]))
})

test_that("min_donors is enforced instead of being an inert argument", {
  a <- imputer_fixture()
  a$category[c(1L, 5L)] <- "High"
  a$y[5L] <- 9
  a$min_donors <- 2L
  expect_error(do.call(mice.impute.catmech_threshold_corruption, a), "Insufficient trusted")
  a$fallback <- "normal"
  result <- do.call(mice.impute.catmech_threshold_corruption, a)
  expect_true(attr(result, "audit")$used_fallback)
  expect_equal(attr(result, "audit")$eligible_donors, 1L)
})

test_that("q endpoints preserve observed values or reconstruct with explicit audit", {
  a <- imputer_fixture()
  a$y[1L] <- 5.123456789
  a$ry[] <- TRUE
  a$wy <- seq_along(a$y) == 1L
  a$flag[1L] <- "Yes"
  a$corruption_probability <- 0
  kept <- do.call(mice.impute.catmech_threshold_corruption, a)
  expect_identical(as.numeric(kept), a$y[1L])
  expect_identical(attr(kept, "audit")$reason, "retained")
  a$corruption_probability <- 1
  replaced <- do.call(mice.impute.catmech_threshold_corruption, a)
  expect_identical(attr(replaced, "audit")$reason, "flag_draw")
  expect_true(attr(replaced, "audit")$reconstructed)
  expect_false(attr(replaced, "audit")$donor_row == 1L)
  a$flag[1L] <- NA_character_
  unknown <- do.call(mice.impute.catmech_threshold_corruption, a)
  expect_identical(as.numeric(unknown), a$y[1L])
})

test_that("a reconstruction decision need not produce a numerical change", {
  a <- imputer_fixture()
  a$y[c(1L, 5L)] <- 9
  a$category[c(1L, 5L)] <- "High"
  a$ry[] <- TRUE
  a$wy <- seq_along(a$y) == 1L
  a$flag[1L] <- "Yes"
  a$corruption_probability <- 1
  result <- do.call(mice.impute.catmech_threshold_corruption, a)
  expect_equal(as.numeric(result), 9)
  expect_true(attr(result, "audit")$reconstructed)
  expect_false(attr(result, "audit")$changed)
  expect_equal(attr(result, "audit")$donor_row, 5L)
})

test_that("exact category boundaries retain the package derivation convention", {
  a <- imputer_fixture()
  a$y[1:2] <- c(3, 7)
  a$category[1:2] <- as.character(catmech_derive_category(
    a$y[1:2], c(-Inf, 3, 7, Inf), c("Low", "Middle", "High")
  ))
  a$ry[] <- TRUE
  a$wy <- seq_along(a$y) %in% 1:2
  result <- do.call(mice.impute.catmech_threshold_corruption, a)
  expect_equal(as.numeric(result), c(3, 7))
  expect_true(all(attr(result, "audit")$reason == "retained"))
})

test_that("retention does not demand an unnecessary regression", {
  a <- imputer_fixture()
  a$y[] <- 5
  a$ry[] <- TRUE
  a$wy <- seq_along(a$y) == 1L
  a$flag[1L] <- "Yes"
  a$corruption_probability <- 0
  expect_equal(as.numeric(do.call(mice.impute.catmech_threshold_corruption, a)), 5)
  a$corruption_probability <- 1
  expect_error(do.call(mice.impute.catmech_threshold_corruption, a), "zero residual variance")
})

test_that("threshold conflicts are forced only if explicitly selected", {
  a <- imputer_fixture()
  a$y[1L] <- 20
  a$ry[] <- TRUE
  a$corruption_probability <- 0
  expect_error(do.call(mice.impute.catmech_threshold_corruption, a), "threshold conflicts")
  a$wy <- seq_along(a$y) == 1L
  result <- do.call(mice.impute.catmech_threshold_corruption, a)
  expect_true(value_in_bounds(as.numeric(result), a$bounds$Middle))
  expect_identical(attr(result, "audit")$reason, "threshold_conflict")
  expect_true(attr(result, "audit")$changed)
})

test_that("q must be valid without silent recycling", {
  a <- imputer_fixture()
  for (bad in list(-0.1, 1.1, NA_real_, Inf, c(0, 1), "0.8")) {
    a$corruption_probability <- bad
    expect_error(do.call(mice.impute.catmech_threshold_corruption, a),
                 "corruption_probability")
  }
})

test_that("unknown categories and unsupported fitting designs fail clearly", {
  a <- imputer_fixture()
  a$category[1] <- NA_character_
  expect_error(do.call(mice.impute.catmech_threshold_corruption, a), "known, nonmissing")
  a <- imputer_fixture()
  a$trusted <- seq_along(a$y) == 2L
  expect_error(do.call(mice.impute.catmech_threshold_corruption, a), "At least two trusted")
  a$trusted <- seq_along(a$y) %in% 2:3
  expect_error(do.call(mice.impute.catmech_threshold_corruption, a), "residual degrees")
})

test_that("narrow intervals are not widened by a fixed numerical epsilon", {
  set.seed(123)
  b <- list(lower = 1, upper = 1 + 1e-10,
            lower_inclusive = FALSE, upper_inclusive = TRUE)
  result <- catmech:::catmech_truncated_draw(1, 1, b)
  expect_gt(result, 1)
  expect_lte(result, b$upper)
})
