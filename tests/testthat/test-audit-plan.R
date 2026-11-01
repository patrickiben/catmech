audit_fixture <- function(n = 400L, audited = 120L, seed = 5L, noise = .5) {
  restore <- preserve_rng()
  on.exit(restore(), add = TRUE)
  set.seed(seed)
  x <- stats::rnorm(n)
  truth <- .5 * x + stats::rnorm(n)
  gold <- rep(NA_real_, n)
  rows <- sample.int(n, audited)
  gold[rows] <- truth[rows]
  data.frame(y = .8 * truth + .3 * x + stats::rnorm(n), x = x, recorded = truth + stats::rnorm(n, 0, noise),
             gold = gold, truth = truth)
}
plan_quietly <- function(...) suppressWarnings(catmech_audit_plan(...))
plan_args <- list(value = "recorded", reference = "gold", predictors = c("x", "y"), formula = y ~ recorded + x,
                  term = "recorded", bootstrap = 100L, m = 5L)
run_plan <- function(data, ...) do.call(plan_quietly, utils::modifyList(c(list(data = data), plan_args), list(...)))

test_that("k has the known closed form under a normal model, with and without the recorded value", {
  restore <- preserve_rng()
  on.exit(restore(), add = TRUE)
  set.seed(61)
  n <- 60000L
  S <- 1; s2 <- 1; b <- .8; d2 <- .25
  truth <- stats::rnorm(n, 0, sqrt(S))
  frame <- data.frame(y = b * truth + stats::rnorm(n, 0, sqrt(s2)), recorded = truth + stats::rnorm(n, 0, sqrt(d2)), gold = truth)
  expected <- function(tau2) tau2 / S * (1 + b^2 * S / s2 - 2 * b^2 * tau2 / s2)
  values <- catmech_audit_pool_k(frame, "recorded", "gold", "y", y ~ recorded, "recorded")
  expect_equal(values[["k"]], expected(1 / (1 / S + 1 / d2 + b^2 / s2)), tolerance = .02)
  expect_equal(values[["discard_k"]], expected(1 / (1 / S + b^2 / s2)), tolerance = .02)
  # With homoskedastic data, the empirical version estimates the same quantity.
  expect_equal(values[["classical"]], values[["k"]], tolerance = .05)
  expect_equal(values[["discard_classical"]], values[["discard_k"]], tolerance = .05)
})

test_that("k has the known closed form when the audited value is the outcome", {
  restore <- preserve_rng()
  on.exit(restore(), add = TRUE)
  set.seed(62)
  n <- 60000L
  x <- stats::rnorm(n)
  truth <- 1 + .5 * x + stats::rnorm(n)
  frame <- data.frame(x = x, recorded = truth + stats::rnorm(n, 0, .5), gold = truth)
  values <- catmech_audit_pool_k(frame, "recorded", "gold", "x", recorded ~ x, "x")
  expect_equal(values[["k"]], 1 / (1 + 1 / .25), tolerance = .02)
  expect_equal(values[["discard_k"]], 1, tolerance = .02)
})

test_that("the closed-form imputation variance agrees with brute-force simulation", {
  frame <- audit_fixture()
  known <- !is.na(frame$gold)
  a <- catmech_audit_analysis(frame, "recorded", y ~ recorded + x, "recorded", ifelse(known, frame$gold, frame$recorded))
  G0 <- catmech_model_matrix(frame, c("x", "y"))
  G <- cbind(G0, .catmech_recorded = frame$recorded)
  p <- catmech_audit_imputation(which(known), G, G0, frame$gold, rep(TRUE, nrow(frame)))
  parts <- catmech_audit_k(a, p$mu, p$tau2, known, frame$gold, detail = TRUE)
  restore <- preserve_rng()
  on.exit(restore(), add = TRUE)
  set.seed(63)
  draws <- vapply(seq_len(20000L), function(i) {
    y <- p$mu + sqrt(p$tau2) * stats::rnorm(length(p$mu))
    (parts$A0 + parts$A1 * y) * (parts$E0 + parts$E1 * y)
  }, numeric(length(p$mu)))
  # The denominator is s2 times weight[j], and weight[j] is the mean squared weight of the expected design.
  center <- parts$A0 + parts$A1 * ifelse(known, frame$gold, p$mu)
  simulated <- mean(apply(draws, 1L, stats::var)) / (parts$s2 * mean(center^2 + parts$A1^2 * ifelse(known, 0, p$tau2)))
  expect_equal(parts$k, simulated, tolerance = .03)
})

test_that("predicted widths track actual linear-MI widths on normal simulated data", {
  restore <- preserve_rng()
  on.exit(restore(), add = TRUE)
  set.seed(64)
  n <- 3000L
  x <- stats::rnorm(n)
  truth <- .5 * x + stats::rnorm(n)
  frame <- data.frame(y = .8 * truth + .3 * x + stats::rnorm(n), x = x, recorded = truth + stats::rnorm(n, 0, .5))
  audit_order <- sample.int(n)
  with_audit <- function(size) {
    d <- frame
    d$gold <- NA_real_
    d$gold[audit_order[seq_len(size)]] <- truth[audit_order[seq_len(size)]]
    d
  }
  plan <- catmech_audit_plan(with_audit(600L), "recorded", "gold", c("x", "y"), y ~ recorded + x, "recorded",
                             fractions = c(.1, .5), bootstrap = 200L, seed = 3L)
  expect_identical(plan$flags, character())
  expect_identical(plan$status, "no_flags")
  actual <- vapply(c(300L, 1500L), function(size) {
    fit <- catmech_reference(with_audit(size), "recorded", "gold", c("x", "y"), use_recorded = TRUE, m = 40L, seed = 8L)
    row <- catmech_pool_lm(fit, y ~ recorded + x)
    row <- row[row$term == "recorded", ]
    row$upper - row$lower
  }, numeric(1))
  expect_equal(plan$table$width, actual, tolerance = .12)
  oracle <- stats::lm(y ~ recorded + x, data = transform(frame, recorded = truth))
  full <- unname(diff(stats::confint(oracle)["recorded", ]))
  expect_equal(plan$pilot$full_width, full, tolerance = .05)
  expect_true(all(plan$table$relative_lower <= plan$table$relative_width))
  expect_true(all(plan$table$relative_width <= plan$table$relative_upper))
})

test_that("the plan table and answer follow the formula", {
  plan <- run_plan(audit_fixture(), fractions = c(.5, .1, .25, 1), target_ratio = 1.1)
  tb <- plan$table
  expect_identical(tb$fraction, c(.1, .25, .5, 1))
  expect_identical(tb$records, c(40, 100, 200, 400))
  expect_identical(tb$additional_records, c(0, 0, 80, 280))
  k <- plan$k$estimate[plan$k$approach == "recorded"]
  k0 <- plan$k$estimate[plan$k$approach == "discard"]
  expect_equal(tb$relative_width, sqrt(1 + k * (1 - tb$fraction) / tb$fraction))
  expect_equal(tb$discard_relative_width, sqrt(1 + k0 * (1 - tb$fraction) / tb$fraction))
  expect_equal(tb$width, plan$pilot$full_width * tb$relative_width)
  expect_equal(tb$discard_width, plan$pilot$discard_full_width * tb$discard_relative_width)
  expect_equal(tb$relative_width[4], 1)
  expect_true(all(diff(tb$relative_width) < 0))
  expect_true(all(tb$discard_relative_width >= tb$relative_width))
  expect_equal(tb$predictive_lower[4], 1)
  expect_true(k0 > k)
  answer <- plan$answer
  expect_equal(sqrt(1 + k * (1 - answer$fraction) / answer$fraction), 1.1)
  expect_identical(answer$records, as.integer(ceiling(answer$fraction * 400)))
  expect_true(answer$records_lower <= answer$records && answer$records <= answer$records_upper)
  expect_identical(answer$smallest_grid_fraction, min(tb$fraction[tb$relative_width <= 1.1]))
  expect_identical(tb$small_audit, c(FALSE, FALSE, FALSE, FALSE))
  expect_true(run_plan(audit_fixture(), fractions = .05)$table$small_audit)
  expect_true(is.na(run_plan(audit_fixture(), fractions = .1, target_ratio = 1.0001)$answer$smallest_grid_fraction))
})

test_that("the pilot fit is the one catmech_reference and catmech_pool_lm return", {
  frame <- audit_fixture()
  plan <- run_plan(frame, seed = 19L)
  fit <- catmech_reference(frame, "recorded", "gold", c("x", "y"), use_recorded = TRUE, m = 5L, seed = 19L)
  row <- catmech_pool_lm(fit, y ~ recorded + x)
  row <- row[row$term == "recorded", ]
  expect_equal(plan$pilot$width, row$upper - row$lower)
  expect_equal(plan$pilot$within_variance, row$within_variance)
  expect_equal(plan$pilot$full_width, 2 * stats::qt(.975, 397) * sqrt(row$within_variance))
})

test_that("the bootstrap is seeded and restores the caller's RNG", {
  frame <- audit_fixture()
  set.seed(71)
  before <- get(".Random.seed", envir = .GlobalEnv)
  first <- run_plan(frame, seed = 2L)
  expect_identical(before, get(".Random.seed", envir = .GlobalEnv))
  second <- run_plan(frame, seed = 2L)
  third <- run_plan(frame, seed = 9L)
  expect_identical(first, second)
  expect_identical(first$k$estimate, third$k$estimate)
  expect_false(identical(first$k$lower, third$k$lower))
  if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) rm(".Random.seed", envir = .GlobalEnv)
  run_plan(frame, seed = 2L)
  expect_false(exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
})

test_that("records without a recorded value use the model without it", {
  frame <- audit_fixture()
  absent <- seq(3L, nrow(frame), by = 10L)
  frame$recorded[absent] <- NA_real_
  plan <- run_plan(frame)
  complete <- run_plan(audit_fixture())
  expect_gt(plan$k$estimate[1], complete$k$estimate[1])
  expect_lt(plan$k$estimate[1], plan$k$estimate[2])
  expect_equal(plan$pilot$calibration_rows, sum(!is.na(frame$gold) & !is.na(frame$recorded)))
  expect_true(all(is.finite(as.matrix(plan$table))))
})

test_that("the analysis can be another coefficient, a factor or an interaction", {
  frame <- audit_fixture()
  frame$arm <- factor(rep(c("a", "b"), length.out = nrow(frame)))
  plan <- run_plan(frame, predictors = c("x", "y", "arm"), formula = y ~ recorded * arm + x, term = "recorded:armb")
  expect_true(is.finite(plan$k$estimate[1]))
  other <- run_plan(frame, predictors = c("x", "y", "arm"), formula = y ~ recorded + arm + x, term = "armb")
  expect_true(other$k$estimate[1] < plan$k$estimate[1])
  outcome <- run_plan(frame, predictors = c("x", "arm"), formula = recorded ~ x + arm, term = "armb")
  expect_true(outcome$k$estimate[1] > 0 && !"outcome_not_in_predictors" %in% outcome$flags)
})

test_that("unusable pilots and analyses are refused", {
  frame <- audit_fixture()
  expect_error(catmech_audit_plan(transform(frame, gold = replace(gold, which(!is.na(gold))[-(1:29)], NA)), "recorded", "gold", c("x", "y"), y ~ recorded + x, "recorded"), "needs at least 30")
  many <- transform(frame, gold = replace(gold, which(!is.na(gold))[-(1:35)], NA))
  for (i in 1:15) many[[paste0("z", i)]] <- sin(i * seq_len(nrow(many)))
  expect_error(catmech_audit_plan(many, "recorded", "gold", c("x", "y", paste0("z", 1:15)), y ~ recorded + x, "recorded"), "needs at least 39")
  expect_error(catmech_audit_plan(transform(frame, gold = NA_real_), "recorded", "gold", c("x", "y"), y ~ recorded + x, "recorded"), "At least one reference")
  expect_error(catmech_audit_plan(transform(frame, gold = truth), "recorded", "gold", c("x", "y"), y ~ recorded + x, "recorded"), "nothing to plan")
  expect_error(catmech_audit_plan(transform(frame, recorded = NA_real_), "recorded", "gold", c("x", "y"), y ~ recorded + x, "recorded"), "No record")
  expect_error(catmech_audit_plan(frame, "recorded", "gold", c("x", "y"), y ~ log(recorded) + x, "log(recorded)"), "linear in")
  expect_error(catmech_audit_plan(frame, "recorded", "gold", c("x", "y"), y ~ I(recorded^2) + x, "I(recorded^2)"), "linear")
  expect_error(catmech_audit_plan(frame, "recorded", "gold", c("x", "y"), y ~ x, "x"), "does not use")
  expect_error(catmech_audit_plan(frame, "recorded", "gold", c("x", "y"), y ~ recorded + x, "absent"), "not a coefficient")
  expect_error(catmech_audit_plan(transform(frame, w = replace(x, 1L, NA)), "recorded", "gold", c("x", "y"), y ~ recorded + w, "recorded"), "must be complete")
  expect_error(catmech_audit_plan(frame, "recorded", "gold", c("x", "y"), "y ~ recorded", "recorded"), "two-sided")
  expect_error(catmech_audit_plan(frame, "recorded", "gold", c("x", "y"), y ~ recorded, c("a", "b")), "one coefficient")
  expect_error(catmech_audit_plan(frame, "recorded", "recorded", c("x", "y"), y ~ recorded, "recorded"), "excluded|exclude|different")
  expect_error(catmech_audit_plan(frame, "recorded", "gold", c("x", "y"), y ~ recorded, "recorded", fractions = c(0, .5)), "fractions")
  expect_error(catmech_audit_plan(frame, "recorded", "gold", c("x", "y"), y ~ recorded, "recorded", target_ratio = 1), "target_ratio")
  expect_error(catmech_audit_plan(frame, "recorded", "gold", c("x", "y"), y ~ recorded, "recorded", conf_level = 1), "conf_level")
  expect_error(catmech_audit_plan(frame, "recorded", "gold", c("x", "y"), y ~ recorded, "recorded", bootstrap = 10L), "bootstrap")
  rare <- frame
  rare$site <- factor(ifelse(is.na(rare$gold) & seq_len(nrow(rare)) %% 2L == 0L, "b", "a"))
  expect_error(suppressWarnings(catmech_audit_plan(rare, "recorded", "gold", c("x", "y", "site"), y ~ recorded + x, "recorded", bootstrap = 100L)), "full column rank")
})

test_that("each risk to the plan raises a clear warning and a flag", {
  frame <- audit_fixture()
  expect_warning(small <- catmech_audit_plan(transform(frame, gold = replace(gold, which(!is.na(gold))[-(1:50)], NA)), "recorded", "gold", c("x", "y"), y ~ recorded + x, "recorded", bootstrap = 100L, m = 5L), "Small pilot")
  expect_identical(small$flags, "small_pilot")
  expect_identical(small$status, "review_flags")
  expect_warning(catmech_audit_plan(frame, "recorded", "gold", "x", y ~ recorded + x, "recorded", bootstrap = 100L, m = 5L), "not in predictors")
  # Exact-or-different records.
  exact <- frame
  exact$recorded <- exact$truth
  expect_warning(none <- catmech_audit_plan(exact, "recorded", "gold", c("x", "y"), y ~ recorded + x, "recorded", bootstrap = 100L, m = 5L), "No pilot row disagrees")
  expect_identical(none$k$estimate[1], 0)
  expect_identical(none$k$upper[1], 0)
  expect_true(all(none$table$relative_width == 1))
  expect_identical(none$pilot$recorded_fit, "recorded_kept_all_calibration_rows_agree")
  few <- exact
  wrong <- which(!is.na(few$gold))[1:6]
  few$recorded[wrong] <- few$truth[wrong] + .3
  expect_warning(some <- catmech_audit_plan(few, "recorded", "gold", c("x", "y"), y ~ recorded + x, "recorded", bootstrap = 100L, m = 5L), "Only 6 pilot rows disagree")
  expect_true(some$k$estimate[1] > 0)
  expect_identical(some$flags, "few_disagreements")
  # Gross errors on unaudited rows, then on pilot rows.
  gross <- frame
  gross$recorded[which(is.na(gross$gold))[1:5]] <- 50
  expect_warning(g1 <- catmech_audit_plan(gross, "recorded", "gold", c("x", "y"), y ~ recorded + x, "recorded", bootstrap = 100L, m = 5L), "gross errors")
  expect_true("gross_errors" %in% g1$flags)
  gross <- frame
  gross$recorded[which(!is.na(gross$gold))[1:2]] <- 12
  expect_warning(g2 <- catmech_audit_plan(gross, "recorded", "gold", c("x", "y"), y ~ recorded + x, "recorded", bootstrap = 100L, m = 5L), "gross errors")
  expect_true(g2$pilot$largest_discrepancy_ratio > 5)
  # A pilot chosen on the recorded value.
  chosen <- frame
  chosen$gold <- ifelse(abs(chosen$recorded - mean(chosen$recorded)) > 1, chosen$truth, NA_real_)
  expect_warning(s <- catmech_audit_plan(chosen, "recorded", "gold", c("x", "y"), y ~ recorded + x, "recorded", bootstrap = 100L, m = 5L), "simple random pilot")
  expect_true("pilot_looks_selective" %in% s$flags)
})
