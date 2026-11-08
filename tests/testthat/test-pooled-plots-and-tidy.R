test_that("plot, tidy and glance work on catmech_pooled_lm, glm and coxph", {
  set.seed(42)
  n <- 80
  d <- data.frame(
    y = rnorm(n),
    y_bin = rbinom(n, 1, 0.4),
    time = rexp(n, rate = 0.1) + 0.1,
    status = rbinom(n, 1, 0.6),
    recorded = rnorm(n),
    reference = rnorm(n),
    age = rnorm(n)
  )
  # 25% audit
  audited <- sample.int(n, 20)
  d$reference[-audited] <- NA_real_

  # recorded is a predictor in every analysis below, so each outcome joins the imputation
  # model: y, y_bin, and for the Cox model the event indicator and Nelson-Aalen hazard.
  d$hazard <- mice::nelsonaalen(d, time, status)
  res <- catmech_reference(d, value = "recorded", reference = "reference",
                           predictors = c("age", "y", "y_bin", "status", "hazard"),
                           use_recorded = TRUE, m = 3L, seed = 123L)

  # Pooled LM
  pooled_lm <- catmech_pool_lm(res, y ~ recorded + age)
  expect_s3_class(pooled_lm, "catmech_pooled_lm")
  expect_output(print(pooled_lm), "Catmech Pooled Linear Model")

  # Tidy and Glance LM
  t_lm <- tidy(pooled_lm)
  expect_true(is.data.frame(t_lm))
  expect_true(all(c("term", "estimate", "std.error", "statistic", "p.value", "conf.low", "conf.high") %in% names(t_lm)))
  expect_equal(nrow(t_lm), 3L)

  g_lm <- glance(pooled_lm)
  expect_true(is.data.frame(g_lm))
  expect_equal(g_lm$nobs, n)
  expect_equal(g_lm$m, 3L)

  # Plot LM (verify graphics call succeeds without error)
  tmp <- tempfile(fileext = ".png")
  grDevices::png(tmp)
  expect_no_error(plot(pooled_lm))
  expect_no_error(plot(pooled_lm, intercept = TRUE))
  grDevices::dev.off()
  unlink(tmp)

  # Pooled GLM (logistic)
  pooled_glm <- catmech_pool_glm(res, y_bin ~ recorded + age, family = stats::binomial(), exponentiate = TRUE)
  expect_s3_class(pooled_glm, "catmech_pooled_glm")

  t_glm <- tidy(pooled_glm, exponentiate = TRUE)
  expect_true(is.data.frame(t_glm))
  expect_true(all(c("term", "estimate", "conf.low", "conf.high") %in% names(t_glm)))

  g_glm <- glance(pooled_glm)
  expect_equal(g_glm$family, "binomial")

  tmp2 <- tempfile(fileext = ".png")
  grDevices::png(tmp2)
  expect_no_error(plot(pooled_glm))
  expect_no_error(plot(pooled_glm, exponentiate = FALSE))
  grDevices::dev.off()
  unlink(tmp2)

  # Pooled Cox PH
  pooled_cox <- catmech_pool_coxph(res, survival::Surv(time, status) ~ recorded + age, exponentiate = TRUE)
  expect_s3_class(pooled_cox, "catmech_pooled_coxph")

  t_cox <- tidy(pooled_cox, exponentiate = TRUE)
  expect_true(is.data.frame(t_cox))
  g_cox <- glance(pooled_cox)
  expect_equal(g_cox$nobs, n)
  expect_equal(g_cox$n_events, sum(d$status))

  tmp3 <- tempfile(fileext = ".png")
  grDevices::png(tmp3)
  expect_no_error(plot(pooled_cox))
  grDevices::dev.off()
  unlink(tmp3)
})
