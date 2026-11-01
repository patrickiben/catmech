# Run after installing catmech: Rscript inst/examples/contamination.R
# All records and references are synthetic. Real use requires independently obtained
# reference measurements of the same quantity, time point and units.
contamination_data <- catmech::simulate_catmech_lab(n = 400, seed = 181, corruption_probability = 0.5)

# A simple random 25% audit: expose truth only on measured rows, then remove the
# simulator's hidden truth and corruption indicator.
set.seed(182)
audited <- sort(sample.int(nrow(contamination_data), nrow(contamination_data) / 4))
contamination_data$reference_value <- NA_real_
contamination_data$reference_value[audited] <- contamination_data$lab_true[audited]
contamination_data[c("lab_true", "corrupted")] <- NULL
contamination_predictors <- c("age", "sex", "inflammation", "outcome")
contamination_formula <- outcome ~ lab_value + age + sex + inflammation

# Step 1: read what the references say before choosing a model.
contamination_audit <- catmech::catmech_recording_audit(
  contamination_data, "lab_value", "reference_value", by = "outlier_flag"
)
print(contamination_audit$overall)
print(contamination_audit$by)
print(contamination_audit$notes)

# Step 2: use recorded values without trusting them. Flagged and unflagged records
# get separate corruption rates, because the audit shows they differ.
contamination_fit <- catmech::catmech_contamination(
  contamination_data, "lab_value", "reference_value", contamination_predictors,
  corruption_strata = "outlier_flag", m = 20, seed = 183
)
print(contamination_fit$contamination$parameters)
print(utils::head(contamination_fit$contamination$corruption[
  order(-contamination_fit$contamination$corruption$posterior_corruption), ]))

# Step 3: score the model on references it did not see. A clean result is not proof.
contamination_check <- catmech::catmech_holdout_check(contamination_fit, folds = 5, seed = 184)
print(contamination_check$status)
print(contamination_check$summary)

# Step 4: compare with the workflow that ignores recorded values, and with declared
# departures from ignorable reference selection. These are assumptions to report.
contamination_scenarios <- list(
  ignores_recorded = catmech::catmech_reference(
    contamination_data, "lab_value", "reference_value", contamination_predictors,
    m = 20, seed = 185
  ),
  contamination = contamination_fit,
  unreferenced_lower = catmech::catmech_contamination(
    contamination_data, "lab_value", "reference_value", contamination_predictors,
    corruption_strata = "outlier_flag", selection_delta = -0.25, m = 20, seed = 186
  ),
  unreferenced_higher = catmech::catmech_contamination(
    contamination_data, "lab_value", "reference_value", contamination_predictors,
    corruption_strata = "outlier_flag", selection_delta = 0.25, m = 20, seed = 187
  )
)
contamination_comparison <- catmech::catmech_sensitivity(
  contamination_scenarios, contamination_formula, "lab_value",
  expected_scenarios = names(contamination_scenarios)
)
print(contamination_comparison$table)
cat("Scenario-set status:", contamination_comparison$status, "\n")
