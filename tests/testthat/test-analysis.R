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
