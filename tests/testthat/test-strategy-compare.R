test_that("catmech_compare_strategies works across all 4 strategies and plots", {
  set.seed(101)
  n <- 60
  d <- data.frame(
    outcome = rnorm(n),
    y_bin = rbinom(n, 1, 0.5),
    lab_value = rnorm(n, mean = 5, sd = 1),
    lab_ref = rnorm(n, mean = 5, sd = 1),
    age = rnorm(n, mean = 50, sd = 10)
  )
  # 20% audit (12 records)
  audited <- sample.int(n, 12)
  d$lab_ref[-audited] <- NA_real_

  # Strategy comparison for linear model
  comp <- catmech_compare_strategies(
    data = d,
    value = "lab_value",
    reference = "lab_ref",
    formula = outcome ~ lab_value + age,
    m = 3L,
    seed = 123L
  )

  expect_s3_class(comp, "catmech_strategy_comparison")
  expect_true(is.data.frame(comp))
  expected_strats <- c("validation_subsample", "naive_trust_recorded", "linear_mi_recorded", "contamination_model")
  expect_true(all(expected_strats %in% comp$strategy))
  expect_true(all(c("term", "estimate", "std.error", "ci_width", "rel_efficiency") %in% names(comp)))
  expect_equal(comp$n_total[1L], n)
  expect_equal(comp$n_audited[1L], 12L)

  # Check print method
  expect_output(print(comp), "Catmech Audit Strategy Comparison")
  expect_output(print(comp), "validation_subsample")

  # Check plot method
  tmp <- tempfile(fileext = ".png")
  grDevices::png(tmp)
  expect_no_error(plot(comp, term = "lab_value"))
  grDevices::dev.off()
  unlink(tmp)

  # Single-predictor formula (outcome ~ lab_value) with missing un-audited lab_value
  d_single <- d
  un_aud <- which(is.na(d_single$lab_ref))[1L]
  d_single$lab_value[un_aud] <- NA_real_
  comp_single <- catmech_compare_strategies(
    data = d_single,
    value = "lab_value",
    reference = "lab_ref",
    formula = outcome ~ lab_value,
    m = 2L,
    seed = 456L
  )
  expect_s3_class(comp_single, "catmech_strategy_comparison")
  expect_true(all(expected_strats %in% comp_single$strategy))

  # Error handling: insufficient audits
  d_tiny <- d
  d_tiny$lab_ref[!is.na(d_tiny$lab_ref)][1:10] <- NA_real_ # only 2 audits left
  expect_error(
    catmech_compare_strategies(d_tiny, "lab_value", "lab_ref", outcome ~ lab_value),
    "At least 5 audited records"
  )
})
