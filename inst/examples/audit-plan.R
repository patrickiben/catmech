# Run after installing catmech: Rscript inst/examples/audit-plan.R
# All records and references in this example are synthetic. A random 25% pilot estimates
# how much precision each larger audit would buy for one regression coefficient.
set.seed(181)
audit_data <- data.frame(age = stats::rnorm(600, 50, 10), site = rep(c("a", "b", "c"), 200))
true_value <- 1 + .01 * (audit_data$age - 50) + stats::rnorm(600, 0, .3)
audit_data$outcome <- 2 + .8 * true_value + stats::rnorm(600, 0, .5)
# About 15% of recorded values are wrong, by a moderate amount.
wrong <- stats::runif(600) < .15
audit_data$recorded <- true_value + wrong * stats::rnorm(600, .1, .25)

# Reference values exist only on a simple random pilot.
pilot_rows <- sample.int(600, 150)
audit_data$reference <- NA_real_
audit_data$reference[pilot_rows] <- true_value[pilot_rows]

audit_plan <- catmech::catmech_audit_plan(
  audit_data, value = "recorded", reference = "reference",
  predictors = c("age", "site", "outcome"),
  formula = outcome ~ recorded + age + site, term = "recorded",
  fractions = c(.10, .25, .50, .75), target_ratio = 1.25, seed = 183
)

# k, the share of the analysis variance the auxiliaries cannot recover, with its interval.
print(audit_plan$k[c("approach", "estimate", "lower", "upper")], digits = 3)
# Records to audit against expected width relative to a full audit, keeping and discarding
# the recorded value. The predictive columns show what one audit of that size may realize.
print(audit_plan$table[c("fraction", "records", "relative_width", "relative_lower", "relative_upper",
                         "predictive_lower", "predictive_upper", "discard_relative_width")], digits = 3)
# Records needed to come within 25% of the full-audit width.
print(audit_plan$answer, digits = 3)
print(audit_plan$assumptions)
