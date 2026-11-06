analysis_fixture <- function() {
  d <- simulate_catmech_lab(180,seed=21)
  d <- d[setdiff(names(d),c("lab_true","corrupted"))]
  catmech_model(d,"lab_value",c("age","sex","inflammation","outcome"),
    flag="outlier_flag",m=8,seed=18)
}

test_that("pooled model agrees with an independent scalar calculation", {
  fit <- analysis_fixture()
  f <- outcome ~ lab_value + age + sex + inflammation
  tab <- catmech_pool_lm(fit,f)
  models <- lapply(fit$completed,function(d) lm(f,data=d))
  q <- vapply(models,function(z) coef(z)["lab_value"],numeric(1))
  u <- vapply(models,function(z) vcov(z)["lab_value","lab_value"],numeric(1))
  wanted_var <- mean(u)+(1+1/length(q))*var(q)
  row <- tab[tab$term=="lab_value",]
  expect_equal(row$estimate,mean(q))
  expect_equal(row$std.error^2,wanted_var)
  expect_equal(row$interpretation,"conditional_model_based")
  expect_lt(row$lower,row$upper)
})

test_that("pooling refuses unconverged and legacy results by default", {
  fit <- analysis_fixture()
  fit$model_diagnostics$convergence_pass <- FALSE
  expect_error(catmech_pool_lm(fit,outcome~lab_value),"convergence")
  fit <- analysis_fixture()
  fit$manifest$engine <- NULL
  expect_error(catmech_pool_lm(fit,outcome~lab_value),"sensitivity")
  tab <- catmech_pool_lm(fit,outcome~lab_value,allow_sensitivity=TRUE)
  expect_true(all(tab$interpretation=="sensitivity_only"))
})

test_that("assumption audit never hides failed scenarios", {
  fit <- analysis_fixture()
  result <- catmech_sensitivity(list(ignored=fit,exact=simpleError("category conflict")),
                                outcome~lab_value+age+sex+inflammation,"lab_value")
  expect_false(result$complete)
  expect_equal(result$status,"unresolved_incomplete_set")
  expect_equal(result$table$status,c("ok","failed"))
  expect_match(result$table$message[2],"category conflict")
  expect_match(result$interpretation,"not a calibrated")
  expect_true(all(is.finite(result$successful_envelope)))
  none <- catmech_sensitivity(list(a=simpleError("failure")),outcome~lab_value,"lab_value")
  expect_true(all(is.na(none$successful_envelope)))
})

test_that("scenario envelope is descriptive and can show a changed conclusion", {
  fit <- analysis_fixture()
  flipped <- fit
  flipped$completed <- lapply(flipped$completed,function(d) {d$outcome <- -d$outcome;d})
  result <- catmech_sensitivity(list(positive=fit,negative=flipped),
    outcome~lab_value+age+sex+inflammation,"lab_value")
  expect_true(result$complete)
  expect_equal(result$status,"assumption_sensitive")
  expect_lt(result$successful_envelope[1],0)
  expect_gt(result$successful_envelope[2],0)
  stable <- catmech_sensitivity(list(a=fit,b=fit),outcome~lab_value+age+sex+inflammation,"lab_value")
  expect_equal(stable$status,"direction_stable_over_supplied_scenarios")
})

test_that("analysis missingness and singular coefficients cannot disappear", {
  fit <- analysis_fixture()
  fit$completed[[1]]$outcome[1] <- NA
  expect_error(catmech_pool_lm(fit,outcome~lab_value),"missing")
  fit <- analysis_fixture()
  fit$completed <- lapply(fit$completed,function(d) {d$duplicate <- d$lab_value;d})
  expect_error(catmech_pool_lm(fit,outcome~lab_value+duplicate),"estimable")
  expect_error(catmech_sensitivity(list(a=fit),"not a formula","lab_value"),"formula")
  expect_error(catmech_pool_lm(fit,outcome~lab_value,conf_level=1),"strictly")
  fit <- analysis_fixture()
  fit$completed <- lapply(fit$completed,function(d) {d$outcome <- 0;d})
  expect_error(catmech_pool_lm(fit,outcome~lab_value),"degenerate")
})

test_that("prespecified scenarios and common data cannot be silently omitted", {
  fit <- analysis_fixture()
  report <- catmech_sensitivity(list(a=fit),outcome~lab_value,"lab_value",
                                expected_scenarios=c("a","missing"))
  expect_equal(report$status,"unresolved_incomplete_set")
  expect_equal(report$table$scenario,c("a","missing"))
  expect_match(report$table$message[2],"not supplied")
  changed <- fit
  changed$original$age[1] <- changed$original$age[1]+1
  expect_error(catmech_sensitivity(list(a=fit,b=changed),outcome~lab_value,"lab_value"),"same original data")
  expect_error(catmech_sensitivity(list(a=fit),outcome~lab_value,"lab_value",
                                  expected_scenarios="b"),"unplanned")
})

test_that("complete-data pooling uses the ordinary regression interval", {
  d <- simulate_catmech_lab(90, seed=43, missing_rate=0, flag_rate=0)
  result <- catmech_model(d,"lab_value",c("age","sex","inflammation","outcome"),m=4)
  formula <- outcome ~ lab_value + age + sex + inflammation
  expected <- lm(formula,data=d)
  pooled <- catmech_pool_lm(result,formula)
  expect_equal(pooled$df,rep(df.residual(expected),nrow(pooled)))
  expect_equal(unname(as.matrix(pooled[c("lower","upper")])),unname(confint(expected)))
  result$completed[[2]]$lab_value[1] <- result$completed[[2]]$lab_value[1]+1
  expect_error(catmech_pool_lm(result,formula),"No-imputation metadata")
})

test_that("zero empirical between-imputation variance is not complete data", {
  result <- analysis_fixture()
  result$completed <- rep(result$completed[1],length(result$completed))
  formula <- outcome ~ lab_value + age + sex + inflammation
  pooled <- catmech_pool_lm(result,formula)
  fit <- lm(formula,data=result$completed[[1]])
  expect_true(all(pooled$df < df.residual(fit)))
  result$manifest$engine <- "unregistered_future_engine"
  expect_error(catmech_pool_lm(result,formula,allow_sensitivity=TRUE),"Unrecognized")
})

test_that("unsupported analysis shapes fail explicitly instead of returning NULL", {
  result <- analysis_fixture()
  expect_error(catmech_pool_lm(result,cbind(outcome,age)~lab_value),"univariate linear models")
  expect_error(catmech_pool_lm(result,outcome~0),"at least one named coefficient")
  expect_equal(catmech_pool_lm(result,outcome~1)$term,"(Intercept)")
  expect_equal(catmech_pool_lm(result,cbind(outcome)~lab_value),
               catmech_pool_lm(result,outcome~lab_value))
})

test_that("S3 print and summary methods provide clean, informative console output", {
  fit <- analysis_fixture()
  txt <- utils::capture.output(val <- print(fit))
  expect_identical(val, fit)
  expect_true(any(grepl("<catmech_result: gaussian_model>", txt)))
  expect_true(any(grepl("Target value:      lab_value", txt)))
  expect_true(any(grepl("Predictors:        age, sex, inflammation, outcome", txt)))
  expect_true(any(grepl("Completions \\(m\\):   8", txt)))

  s <- summary(fit)
  expect_s3_class(s, "summary.catmech_result")
  s_txt <- utils::capture.output(s_val <- print(s))
  expect_identical(s_val, s)
  expect_true(any(grepl("<summary.catmech_result: gaussian_model>", s_txt)))
  expect_true(any(grepl("Observations:    180", s_txt)))

  diag <- fit$model_diagnostics
  d_txt <- utils::capture.output(d_val <- print(diag))
  expect_identical(d_val, diag)
  expect_true(any(grepl("<catmech_model_diagnostics:", d_txt)))

  # catmech_sensitivity S3 print method
  sens <- catmech_sensitivity(list(default = fit), formula = outcome ~ lab_value + age, term = "lab_value")
  expect_s3_class(sens, "catmech_sensitivity")
  sens_txt <- utils::capture.output(sens_val <- print(sens))
  expect_identical(sens_val, sens)
  expect_true(any(grepl("Catmech Sensitivity Analysis", sens_txt)))
  expect_true(any(grepl("Successful envelope:", sens_txt)))
  expect_true(any(grepl("default", sens_txt)))
})

test_that("catmech_pool_glm pools logistic and Poisson models with Rubin's rules", {
  fit <- analysis_fixture()
  for (j in seq_along(fit$completed)) {
    fit$completed[[j]]$y_bin <- as.numeric(fit$completed[[j]]$outcome > stats::median(fit$completed[[j]]$outcome))
    fit$completed[[j]]$y_cnt <- stats::rpois(nrow(fit$completed[[j]]), lambda = pmax(0.1, abs(fit$completed[[j]]$lab_value)))
    fit$completed[[j]]$site <- rep(seq_len(10L), length.out = nrow(fit$completed[[j]]))
  }
  # Logistic regression
  res_bin <- catmech_pool_glm(fit, y_bin ~ lab_value + age, family = stats::binomial(), exponentiate = TRUE)
  expect_s3_class(res_bin, "catmech_pooled_glm")
  expect_true("estimate_exp" %in% names(res_bin))
  expect_equal(res_bin$estimate_exp, exp(res_bin$estimate), tolerance = 1e-10)
  expect_equal(res_bin$lower_exp, exp(res_bin$lower), tolerance = 1e-10)
  expect_true(all(res_bin$df > 0))
  expect_true(all(res_bin$std.error > 0))

  # Test print method
  printed <- utils::capture.output(val <- print(res_bin))
  expect_identical(val, res_bin)
  expect_true(any(grepl("Catmech Pooled GLM: binomial", printed)))
  expect_true(any(grepl("estimate_exp", printed)))

  # Poisson regression
  res_cnt <- catmech_pool_glm(fit, y_cnt ~ lab_value + age, family = stats::poisson())
  expect_identical(res_cnt$family[1L], "poisson")

  # Cluster-robust GLM covariance
  res_clust <- catmech_pool_glm(fit, y_bin ~ lab_value + age, family = stats::binomial(), cluster = "site")
  expect_true("clusters" %in% names(res_clust))
  expect_identical(res_clust$clusters[1L], 10L)

  # Sensitivity analysis with pool_glm
  sens_glm <- catmech_sensitivity(list(default = fit), formula = y_bin ~ lab_value + age,
                                  term = "lab_value", method = "pool_glm", family = stats::binomial())
  expect_s3_class(sens_glm, "catmech_sensitivity")
  expect_identical(sens_glm$table$status[1L], "ok")
  expect_true(is.finite(sens_glm$table$estimate[1L]))
})

test_that("catmech_pool_coxph pools Cox proportional hazards models with Rubin's rules", {
  fit <- analysis_fixture()
  set.seed(123)
  n <- nrow(fit$completed[[1L]])
  time_fixed <- stats::rexp(n, rate = 0.1)
  status_fixed <- stats::rbinom(n, 1, prob = 0.6)
  site_fixed <- rep(seq_len(10L), length.out = n)

  for (j in seq_along(fit$completed)) {
    fit$completed[[j]]$surv_time <- time_fixed
    fit$completed[[j]]$surv_event <- status_fixed
    fit$completed[[j]]$site <- site_fixed
  }

  f_cox <- survival::Surv(surv_time, surv_event) ~ lab_value + age
  res_cox <- catmech_pool_coxph(fit, f_cox, exponentiate = TRUE)

  expect_s3_class(res_cox, "catmech_pooled_coxph")
  expect_true(all(c("estimate", "std.error", "df", "lower", "upper",
                    "estimate_exp", "lower_exp", "upper_exp",
                    "within_variance", "between_variance", "nevent") %in% names(res_cox)))
  expect_equal(res_cox$estimate_exp, exp(res_cox$estimate), tolerance = 1e-10)
  expect_equal(res_cox$lower_exp, exp(res_cox$lower), tolerance = 1e-10)
  expect_equal(res_cox$upper_exp, exp(res_cox$upper), tolerance = 1e-10)
  expect_true(all(res_cox$df > 0))
  expect_true(all(res_cox$std.error > 0))
  expect_equal(res_cox$nevent[1L], sum(status_fixed))

  # Test print method
  printed <- utils::capture.output(val <- print(res_cox))
  expect_identical(val, res_cox)
  expect_true(any(grepl("Catmech Pooled Cox PH Model", printed)))
  expect_true(any(grepl("Events:", printed)))
  expect_true(any(grepl("estimate_exp", printed)))

  # Cluster-robust Cox covariance
  res_clust <- catmech_pool_coxph(fit, f_cox, cluster = "site")
  expect_true("clusters" %in% names(res_clust))
  expect_identical(res_clust$clusters[1L], 10L)
  expect_true(grepl("_cluster_robust", res_clust$interpretation[1L]))

  # Sensitivity analysis with pool_coxph
  sens_cox <- catmech_sensitivity(list(default = fit), formula = f_cox,
                                  term = "lab_value", method = "pool_coxph")
  expect_s3_class(sens_cox, "catmech_sensitivity")
  expect_identical(sens_cox$table$status[1L], "ok")
  expect_true(is.finite(sens_cox$table$estimate[1L]))

  # Negative control: missing cluster column
  expect_error(catmech_pool_coxph(fit, f_cox, cluster = "nonexistent_col"),
               "cluster must name a column")
})
