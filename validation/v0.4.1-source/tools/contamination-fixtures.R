# Scenario generators for the v0.4 contamination study. Not exported package API.
# Each case starts from the unchanged cohort in tools/lab-fixture.R. Mechanisms not
# in that fixture are applied here with their own seeded streams, so the imputers
# only ever see `input`; simulated truth stays in the evaluator's fields.

contamination_stream_bases <- function(phase) {
  if (phase == "heldout") {
    list(selection = 180100000, noise = 186100000, mechanism = 188100000, imputation = 190000000)
  } else {
    list(selection = 170100000, noise = 176100000, mechanism = 178100000, imputation = 175000000)
  }
}

contamination_scenarios <- function(primary_reps, stress_reps) {
  cell <- function(scenario, family, fixture = "missed_flags", fraction = .25, selection = "random",
                   noise = 0, mechanism = "none", corruption_strata = "flag", error_strata = NULL,
                   agreement_tol = 0, reference_sd = 0, methods = NULL, holdout = FALSE) {
    list(scenario = scenario, family = family, fixture = fixture, fraction = fraction,
         selection = selection, noise = noise, mechanism = mechanism,
         corruption_strata = corruption_strata, error_strata = error_strata,
         agreement_tol = agreement_tol, reference_sd = reference_sd,
         methods = methods, holdout = holdout,
         reps = if (family == "primary") primary_reps else stress_reps)
  }
  list(
    cell("clean", "primary", fixture = "clean", holdout = TRUE),
    cell("flagged_errors", "primary", fixture = "flagged_errors"),
    cell("missed_flags", "primary", holdout = TRUE),
    cell("site_shift_declared", "primary", fixture = "site_shift",
         corruption_strata = "SITEID", error_strata = "SITEID"),
    cell("rounded_records", "support", fixture = "rounded_records", agreement_tol = .005),
    cell("reference_05pct", "support", fraction = .05),
    cell("reference_10pct", "support", fraction = .10),
    cell("reference_50pct", "support", fraction = .50),
    cell("site_shift_undeclared", "support", fixture = "site_shift", holdout = TRUE),
    cell("targeted_audit", "credibility", selection = "targeted"),
    cell("nonignorable_clean", "credibility", fixture = "clean", selection = "mnar"),
    cell("nonignorable_missed", "credibility", selection = "mnar"),
    cell("noisy_declared", "credibility", noise = .08, reference_sd = .08),
    cell("noisy_wrong_sd", "credibility", noise = .08, reference_sd = .04, methods = "contamination"),
    cell("noisy_undeclared", "credibility", noise = .08, methods = "contamination"),
    cell("heavy_tail_target", "wrong_model", mechanism = "heavy_tail", holdout = TRUE),
    cell("nonlinear_target", "wrong_model", mechanism = "nonlinear", holdout = TRUE),
    cell("heteroskedastic_target", "wrong_model", mechanism = "heteroskedastic"),
    cell("error_bimodal", "wrong_model", mechanism = "error_bimodal", holdout = TRUE),
    cell("error_gross_and_small", "wrong_model", mechanism = "error_gross_and_small"),
    cell("error_multiplicative", "wrong_model", mechanism = "error_multiplicative"),
    cell("differential_outcome", "wrong_model", mechanism = "differential_outcome", holdout = TRUE),
    cell("value_dependent_error", "wrong_model", mechanism = "value_dependent_error"),
    cell("flag_from_recorded", "wrong_model", mechanism = "flag_from_recorded", holdout = TRUE),
    cell("assay_noise_tolerance", "wrong_model", mechanism = "assay_noise", agreement_tol = .06),
    cell("assay_noise_no_tolerance", "wrong_model", mechanism = "assay_noise", methods = "contamination"),
    cell("recorded_missing_mnar", "wrong_model", mechanism = "recorded_missing_mnar"))
}

contamination_logit_intercept <- function(score, target) {
  stats::uniroot(function(a) mean(stats::plogis(a + score)) - target, c(-30, 30))$root
}

# Redraws who is missing, who is wrong and who is flagged, with the same error
# size and flag operating characteristics as the fixture's missed_flags scenario.
contamination_redraw_recording <- function(target, missing_probability, corrupt_probability) {
  n <- length(target)
  hidden <- stats::runif(n) < missing_probability
  wrong <- !hidden & stats::runif(n) < corrupt_probability
  shift <- stats::rnorm(n, .35, .07)
  flagged <- !hidden & stats::runif(n) < ifelse(wrong, .60, .05)
  recorded <- target + ifelse(wrong, shift, 0)
  recorded[hidden] <- NA_real_
  list(recorded = recorded, flagged = flagged)
}

make_contamination_case <- function(s, r, phase, master, make_fixture) {
  bases <- contamination_stream_bases(phase)
  f <- make_fixture(n = 600L, replicate = r, phase = phase, seed = master, scenario = s$fixture)
  d <- f$data
  truth <- f$truth$LAB_TRUE
  recorded <- d$LBSTRESN
  outcome <- d$outcome
  flagged <- as.character(d$QCFLAG) == "Yes"
  site <- as.integer(d$SITEID)
  centre <- 0.95 + 0.001 * (d$age - 50) + 0.025 * (d$sex == "M") + 0.025 * (d$inflammation - 2) +
    c(-0.025, 0, 0.025)[site]
  observed <- !is.na(recorded)
  set.seed(as.integer(bases$mechanism + r))
  replace_truth <- function(new_truth) {
    # Same recording shifts and the same outcome noise, so the slope stays at 0.8.
    change <- new_truth - truth
    recorded <<- recorded + change
    outcome <<- outcome + 0.8 * change
    truth <<- new_truth
  }
  redraw_positive <- function(draw_residual) {
    value <- centre + draw_residual(length(centre))
    while (any(value <= 0)) {
      low <- which(value <= 0)
      value[low] <- centre[low] + draw_residual(length(centre))[low]
    }
    value
  }
  shift_rows <- which(observed & f$truth$RANDOM_ERROR)
  if (s$mechanism == "heavy_tail") {
    replace_truth(redraw_positive(function(k) .2 * stats::rt(k, 3) / sqrt(3)))
  } else if (s$mechanism == "nonlinear") {
    replace_truth(redraw_positive(function(k) stats::rnorm(k, 0, .2)) + .06 * ((d$inflammation - 2)^2 - 1))
  } else if (s$mechanism == "heteroskedastic") {
    scale <- exp(.5 * (d$inflammation - 2))
    scale <- scale / sqrt(mean(scale^2))
    replace_truth(redraw_positive(function(k) .2 * scale * stats::rnorm(k)))
  } else if (s$mechanism == "error_bimodal") {
    sign <- ifelse(stats::runif(length(shift_rows)) < .5, -1, 1)
    recorded[shift_rows] <- truth[shift_rows] + sign * f$truth$VALUE_SHIFT[shift_rows]
  } else if (s$mechanism == "error_gross_and_small") {
    gross <- stats::runif(length(shift_rows)) < .25
    recorded[shift_rows] <- ifelse(gross, 10 * truth[shift_rows],
      truth[shift_rows] + stats::rnorm(length(shift_rows), .10, .05))
  } else if (s$mechanism == "error_multiplicative") {
    recorded[shift_rows] <- 1.4 * truth[shift_rows]
  } else if (s$mechanism == "differential_outcome") {
    score <- 1.5 * as.numeric(scale(outcome))
    made <- contamination_redraw_recording(truth, .2, stats::plogis(contamination_logit_intercept(score, .08) + score))
    recorded <- made$recorded
    flagged <- made$flagged
  } else if (s$mechanism == "value_dependent_error") {
    score <- 1.5 * as.numeric(scale(truth))
    made <- contamination_redraw_recording(truth, .2, stats::plogis(contamination_logit_intercept(score, .08) + score))
    recorded <- made$recorded
    flagged <- made$flagged
  } else if (s$mechanism == "flag_from_recorded") {
    flagged <- !is.na(recorded) & (recorded < .55 | recorded > 1.35)
  } else if (s$mechanism == "assay_noise") {
    clean <- which(observed & !f$truth$RANDOM_ERROR)
    recorded[clean] <- truth[clean] + stats::rnorm(length(clean), 0, .02)
  } else if (s$mechanism == "recorded_missing_mnar") {
    score <- 2 * as.numeric(scale(truth))
    made <- contamination_redraw_recording(truth, stats::plogis(contamination_logit_intercept(score, .2) + score), .08)
    recorded <- made$recorded
    flagged <- made$flagged
  } else if (s$mechanism != "none") stop("Unknown mechanism: ", s$mechanism)
  predictors <- c("age", "sex", "inflammation", "SITEID", "outcome")
  input <- d[, setdiff(predictors, "outcome")]
  input$outcome <- outcome
  input$lab_value <- recorded
  input$flag <- factor(ifelse(flagged, "Yes", "No"), levels = c("No", "Yes"))
  set.seed(as.integer(bases$selection + r))
  u <- stats::runif(nrow(input))
  probability <- rep(s$fraction, nrow(input))
  analyst_probability <- probability
  if (s$selection == "targeted") {
    z <- (recorded - mean(recorded, na.rm = TRUE)) / stats::sd(recorded, na.rm = TRUE)
    probability <- ifelse(is.na(recorded), .15, ifelse(flagged, .80, .15 + .50 * (abs(z) > 1.5)))
    analyst_probability <- probability
  } else if (s$selection == "mnar") {
    mx <- stats::model.matrix(stats::reformulate(predictors), data = input)
    residual <- as.numeric(scale(stats::lm.fit(mx, truth)$residuals))
    probability <- stats::plogis(contamination_logit_intercept(2 * residual, .25) + 2 * residual)
    analyst_probability <- rep(.25, nrow(input))
  }
  gold <- truth
  if (s$noise > 0) {
    set.seed(as.integer(bases$noise + r))
    gold <- gold + stats::rnorm(nrow(input), 0, s$noise)
  }
  gold[u >= probability] <- NA_real_
  input$gold_value <- gold
  analysis <- setdiff(predictors, "outcome")
  list(input = input, truth = truth, slope = f$manifest$known_slope, predictors = predictors,
       corrupted = !is.na(recorded) & abs(recorded - truth) > s$agreement_tol + 1e-8,
       sampling_probability = analyst_probability,
       analysis_x = stats::model.matrix(stats::reformulate(analysis), data = input),
       model_x = stats::model.matrix(stats::reformulate(predictors), data = input)[, -1L, drop = FALSE],
       selection_seed = as.integer(bases$selection + r),
       imputation_seed = as.integer(bases$imputation + r),
       stream = paste0("master", master, ":", phase, ":", r))
}
