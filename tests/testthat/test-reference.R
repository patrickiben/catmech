reference_fixture <- function() {
  d <- simulate_catmech_lab(120L, seed = 37, corruption_probability = .5)
  d$gold <- NA_real_
  d$gold[seq(1L, nrow(d), by = 3L)] <- d$lab_true[seq(1L, nrow(d), by = 3L)]
  d[c("lab_true", "corrupted")] <- NULL
  d
}

test_that("trusted PMM delegates to the declared standard algorithm", {
  d <- reference_fixture()
  predictors <- c("age", "inflammation", "outcome", "lab_category")
  fit <- catmech_pmm(d, "lab_value", predictors, flag = "outlier_flag", m = 4, seed = 17)
  observed <- d$lab_value
  observed[d$outlier_flag == "Yes"] <- NA_real_
  ry <- !is.na(observed)
  x <- model.matrix(~ age + inflammation + outcome + lab_category, d)[, -1L, drop = FALSE]
  set.seed(17)
  for (j in 1:4) {
    expected <- observed
    expected[!ry] <- mice::mice.impute.pmm(observed, ry, x, wy = !ry,
      donors = 5, matchtype = 1, ridge = 1e-5)
    expect_equal(fit$completed[[j]]$lab_value, expected)
    expect_identical(fit$completed[[j]]$lab_category, d$lab_category)
    expect_equal(mice::complete(fit$mids, j)$lab_value, expected)
  }
  expect_identical(fit$original, d)
  expect_identical(fit$manifest$engine, "trusted_pmm")
  expect_identical(fit$manifest$imputed_rows, which(!ry))
})

test_that("matching provenance admits duplicate donor identity ambiguity", {
  d <- data.frame(y = c(1, 1, 2, 2, 3, 3, NA_real_), x = c(1, 3, 2, 4, 0, 5, 2.5))
  fit <- catmech_pmm(d, "y", "x", donors = 5, m = 3, seed = 21)
  expect_true(all(is.na(fit$audit$donor_row)))
  expect_true(all(fit$audit$donor_value_match_count == 2L))
  expect_false(any(fit$audit$donor_row_identified))
  expect_true(all(grepl(";", fit$audit$donor_value_rows, fixed = TRUE)))
  expect_true(all(fit$audit$eligible_donors == 6L))
  expect_true(all(fit$audit$candidate_pool_size == 5L))
})

test_that("PMM never silently reduces donor count or accepts unsupported design", {
  d <- data.frame(y = c(1, 2, 3, NA_real_), x = c(0, 1, 0, 1))
  expect_error(catmech_pmm(d, "y", "x", donors = 5), "Fewer trusted/reference")
  d$x <- c(0, 0, 0, 1)
  expect_error(catmech_pmm(d, "y", "x", donors = 2), "full column rank")
  d$x[2] <- NA_real_
  expect_error(catmech_pmm(d, "y", "x", donors = 2), "complete")
})

test_that("reference normal preserves anchors and ignores raw target at fixed seed", {
  d <- reference_fixture()
  pred <- c("age", "inflammation", "outcome")
  set.seed(21)
  before <- .Random.seed
  fit <- catmech_reference(d, "lab_value", "gold", pred, m = 4, seed = 24)
  expect_identical(.Random.seed, before)
  altered <- d
  altered$lab_value <- seq_len(nrow(d)) * 100
  other <- catmech_reference(altered, "lab_value", "gold", pred, m = 4, seed = 24)
  known <- !is.na(d$gold)
  for (j in 1:4) {
    expect_identical(fit$completed[[j]]$lab_value, other$completed[[j]]$lab_value)
    expect_equal(fit$completed[[j]]$lab_value[known], d$gold[known])
    expect_identical(fit$completed[[j]][setdiff(names(d), "lab_value")],
                     d[setdiff(names(d), "lab_value")])
    expect_equal(mice::complete(fit$mids, j)$lab_value, fit$completed[[j]]$lab_value)
  }
  expect_identical(fit$original, d)
  expect_identical(fit$mids$data$lab_value, d$gold)
  expect_identical(fit$manifest$engine, "reference_normal")
  expect_identical(fit$manifest$imputed_rows, which(!known))
  expect_true(all(fit$audit$reason[fit$audit$reference_observed] == "gold_reference"))
  expect_true(all(!fit$audit$reconstructed[fit$audit$reference_observed]))
})

test_that("reference PMM draws only from independently measured values", {
  d <- reference_fixture()
  fit <- catmech_reference(d, "lab_value", "gold", c("age", "outcome"),
                           method = "pmm", m = 4, seed = 39)
  known <- !is.na(d$gold)
  for (j in 1:4) {
    expect_true(all(fit$completed[[j]]$lab_value[!known] %in% d$gold[known]))
    expect_equal(fit$completed[[j]]$lab_value[known], d$gold[known])
  }
  expect_identical(fit$manifest$engine, "reference_pmm")
  expect_true(all(fit$audit$donor_row[!fit$audit$reference_observed] %in% which(known)))
})

test_that("full reference bypasses unnecessary stochastic fitting", {
  d <- data.frame(raw = c(100, -10, NA_real_), gold = rep(2, 3), x = 1:3,
                  constant = 0, single_level = factor(rep("one", 3)))
  for (method in c("normal", "pmm")) {
    fit <- catmech_reference(d, "raw", "gold", c("x", "constant", "single_level"),
                             method = method, m = 2, donors = 5)
    expect_equal(fit$completed[[1]]$raw, d$gold)
    expect_equal(mice::complete(fit$mids, 2)$raw, d$gold)
    expect_length(fit$manifest$imputed_rows, 0)
    expect_true(all(fit$diagnostics$requested == 0))
    expect_true(all(fit$audit$reason == "gold_reference"))
  }
  d$constant[1] <- NA_real_
  expect_error(catmech_reference(d, "raw", "gold", "constant", m = 2), "complete")
  d$constant[1] <- Inf
  expect_error(catmech_reference(d, "raw", "gold", "constant", m = 2), "finite")
})

test_that("observed-only PMM has no stochastic fit or donor requirement", {
  d <- data.frame(y = c(1, 2, 3), x = 1:3, constant = 1,
                  single_level = factor(rep("one", 3)))
  fit <- catmech_pmm(d, "y", c("x", "constant", "single_level"), donors = 5, m = 2)
  expect_equal(fit$completed[[1]], d)
  expect_equal(mice::complete(fit$mids, 2)$y, d$y)
  expect_length(fit$manifest$imputed_rows, 0)
  expect_equal(nrow(fit$audit), 0L)
  d$y[3] <- NA_real_
  expect_error(catmech_pmm(d, "y", c("constant", "single_level"), donors = 2), "constant predictor")
})

test_that("both PMM APIs restore RNG and retain public call provenance", {
  d <- reference_fixture()
  set.seed(941)
  before <- .Random.seed
  fit <- catmech_pmm(d, "lab_value", c("age", "outcome"), flag = "outlier_flag", m = 2)
  expect_identical(.Random.seed, before)
  expect_identical(fit$mids$call[[1]], as.name("catmech_pmm"))
  fit <- catmech_reference(d, "lab_value", "gold", c("age", "outcome"), method = "pmm", m = 2)
  expect_identical(.Random.seed, before)
  expect_identical(fit$mids$call[[1]], as.name("catmech_reference"))
  expect_identical(fit$manifest$claim_level, "sensitivity_only")
  expect_identical(fit$spec$predictors, c("age", "outcome"))
})

test_that("partial reference PMM requires explicitly labeled sensitivity pooling", {
  d <- reference_fixture()
  fit <- catmech_reference(d, "lab_value", "gold", c("age", "outcome"),
                           method = "pmm", m = 4, seed = 39)
  expect_identical(fit$manifest$claim_level, "sensitivity_only")
  expect_error(catmech_pool_lm(fit, outcome ~ lab_value + age), "allow_sensitivity")
  pooled <- catmech_pool_lm(fit, outcome ~ lab_value + age, allow_sensitivity = TRUE)
  expect_true(all(pooled$interpretation == "sensitivity_only"))
  report <- catmech_sensitivity(list(reference_pmm = fit), outcome ~ lab_value + age,
                                "lab_value", expected_scenarios = "reference_pmm")
  expect_true(report$complete)
  expect_identical(report$table$interpretation, "sensitivity_only")
  normal <- catmech_reference(d, "lab_value", "gold", c("age", "outcome"), m = 4)
  expect_true(all(catmech_pool_lm(normal, outcome ~ lab_value + age)$interpretation ==
                    "conditional_reference_model"))
})

test_that("complete reference PMM permits ordinary complete-data pooling", {
  d <- simulate_catmech_lab(80, seed = 48)
  d$gold <- d$lab_true
  d[c("lab_true", "corrupted")] <- NULL
  fit <- catmech_reference(d, "lab_value", "gold", c("age", "outcome"),
                           method = "pmm", m = 3)
  expect_length(fit$manifest$imputed_rows, 0L)
  expect_identical(fit$manifest$claim_level, "conditional_reference_matching")
  pooled <- catmech_pool_lm(fit, outcome ~ lab_value + age)
  expect_true(all(pooled$interpretation == "conditional_reference_matching"))
  expected <- lm(outcome ~ lab_value + age, data = fit$completed[[1]])
  expect_equal(pooled$df, rep(df.residual(expected), nrow(pooled)))
  expect_equal(unname(as.matrix(pooled[c("lower", "upper")])), unname(confint(expected)))
})

test_that("reference input and support failures remain explicit", {
  d <- reference_fixture()
  expect_error(catmech_reference(d, "lab_value", "absent", "age"), "existing column")
  expect_error(catmech_reference(d, "lab_value", "gold", c("age", "gold")), "exclude")
  d$gold[] <- NA_real_
  expect_error(catmech_reference(d, "lab_value", "gold", "age"), "independently observed")
  d$gold[1:3] <- 1:3
  d$only_missing <- 1
  d$only_missing[1:3] <- 0
  expect_error(catmech_reference(d, "lab_value", "gold", "only_missing"), "full column rank")
  d$proxy <- d$lab_value
  expect_error(catmech_reference(d, "lab_value", "gold", "proxy"), "complete")
  d$gold[1] <- Inf
  expect_error(catmech_reference(d, "lab_value", "gold", "age"), "finite numeric")
})

test_that("reference audit distinguishes known corrections from imputations", {
  d <- reference_fixture()
  known <- which(!is.na(d$gold))
  d$lab_value[known[1:2]] <- c(-100, NA_real_)
  fit <- catmech_reference(d, "lab_value", "gold", c("age", "outcome"), m = 2)
  a <- fit$audit[fit$audit$imputation == 1, ]
  expect_true(a$changed[a$row == known[1]])
  expect_true(is.na(a$changed[a$row == known[2]]))
  expect_true(all(a$reference_observed[match(known, a$row)]))
  expect_equal(fit$diagnostics$reference_observed, rep(length(known), 2))
  expect_equal(fit$manifest$data_hash,
    digest::digest(serialize(d, NULL, version = 2), algo = "sha256", serialize = FALSE))
  expect_equal(fit$manifest$assumption_hash, catmech_hash(fit$spec))
})

test_that("use_recorded imputes from recorded values and handles every pattern", {
  d <- reference_fixture()
  predictors <- c("age", "outcome")
  set.seed(12)
  before <- .Random.seed
  fit <- catmech_reference(d, "lab_value", "gold", predictors, use_recorded = TRUE, m = 4, seed = 31)
  expect_identical(.Random.seed, before)
  known <- !is.na(d$gold)
  for (j in 1:4) {
    expect_equal(fit$completed[[j]]$lab_value[known], d$gold[known])
    expect_true(all(is.finite(fit$completed[[j]]$lab_value)))
    expect_equal(mice::complete(fit$mids, j)$lab_value, fit$completed[[j]]$lab_value)
  }
  expect_identical(fit$manifest$engine, "reference_normal_recorded")
  expect_identical(fit$spec$recorded_fit, "linear_given_recorded")
  expect_identical(fit$manifest$imputed_rows, which(!known))
  expect_true(all(catmech_pool_lm(fit, outcome ~ lab_value + age)$interpretation == "conditional_reference_model"))
  # The recorded value informs the imputation, so changing it changes the result.
  altered <- d
  altered$lab_value[!known & !is.na(d$lab_value)] <- altered$lab_value[!known & !is.na(d$lab_value)] + 5
  other <- catmech_reference(altered, "lab_value", "gold", predictors, use_recorded = TRUE, m = 4, seed = 31)
  expect_false(identical(fit$completed[[1]]$lab_value, other$completed[[1]]$lab_value))
  # The plain workflow is unchanged by the new argument's existence.
  plain <- catmech_reference(d, "lab_value", "gold", predictors, m = 4, seed = 31)
  expect_null(plain$spec$use_recorded)
  expect_identical(plain$manifest$engine, "reference_normal")
})

test_that("use_recorded keeps recorded values when calibration rows all agree", {
  d <- simulate_catmech_lab(120L, seed = 44, flag_rate = 0, corruption_probability = 0)
  d$gold <- NA_real_
  d$gold[seq(1L, 120L, 3L)] <- d$lab_true[seq(1L, 120L, 3L)]
  d[c("lab_true", "corrupted")] <- NULL
  fit <- catmech_reference(d, "lab_value", "gold", c("age", "outcome"), use_recorded = TRUE, m = 3, seed = 2)
  open <- is.na(d$gold) & !is.na(d$lab_value)
  expect_identical(fit$spec$recorded_fit, "recorded_kept_all_calibration_rows_agree")
  expect_equal(fit$completed[[2]]$lab_value[open], d$lab_value[open])
  # Records with no recorded value are still imputed.
  expect_true(all(is.finite(fit$completed[[2]]$lab_value)))
  expect_false(anyNA(fit$completed[[2]]$lab_value))
})

test_that("use_recorded refuses unsupported requests", {
  d <- reference_fixture()
  expect_error(catmech_reference(d, "lab_value", "gold", "age", method = "pmm", use_recorded = TRUE), "only to method")
  expect_error(catmech_reference(d, "lab_value", "gold", "age", use_recorded = NA), "TRUE or FALSE")
  apart <- d
  apart$gold <- ifelse(is.na(d$lab_value), 1 + d$age / 50, NA_real_)
  expect_error(catmech_reference(apart, "lab_value", "gold", "age", use_recorded = TRUE), "both recorded and reference")
})
