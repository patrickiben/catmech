# Run after installing catmech: Rscript inst/examples/reference.R
# All records and gold measurements in this example are synthetic. Real use
# requires independently obtained reference measurements of the same quantity,
# time point and units. Do not turn flags, categories or imputed values into gold.
reference_data <- catmech::simulate_catmech_lab(
  n = 300, seed = 171, corruption_probability = 0.5,
  category_error_rate = 0.1
)

# A simple random 25% validation subset: expose truth only at measured anchors,
# then immediately remove the simulator's hidden truth and corruption indicator.
set.seed(172)
reference_rows <- sort(sample.int(nrow(reference_data), nrow(reference_data) / 4))
reference_data$gold_value <- NA_real_
reference_data$gold_value[reference_rows] <- reference_data$lab_true[reference_rows]
reference_data[c("lab_true", "corrupted")] <- NULL
reference_predictors <- c("age", "sex", "inflammation", "outcome")
reference_formula <- outcome ~ lab_value + age + sex + inflammation
reference_expected <- c("flag_normal", "flag_pmm", "flag_pmm_category",
                        "reference_normal", "reference_pmm")

# Every workflow receives the same original data. These are different declared
# assumptions and information sets, not a competition to choose a favorable fit.
# The flag workflows trust unflagged observed lab values. Neither uses gold.
reference_scenarios <- list(
  flag_normal = tryCatch(catmech::catmech_model(
    reference_data, "lab_value", reference_predictors,
    flag = "outlier_flag", category_mode = "ignore", m = 20, seed = 173
  ), error = identity),
  flag_pmm = tryCatch(catmech::catmech_pmm(
    reference_data, "lab_value", reference_predictors,
    flag = "outlier_flag", donors = 5, m = 20, seed = 174
  ), error = identity),
  # A category is an ordinary predictor here; no interval is forced. Comparing
  # this with flag_pmm changes information, so any difference is not solely PMM.
  flag_pmm_category = tryCatch(catmech::catmech_pmm(
    reference_data, "lab_value", c(reference_predictors, "lab_category"),
    flag = "outlier_flag", donors = 5, m = 20, seed = 175
  ), error = identity),
  # Only independently measured gold values anchor these reference workflows.
  # All non-reference targets are unavailable even when an unflagged raw value
  # exists. The raw measurement is not automatically used as a predictor.
  reference_normal = tryCatch(catmech::catmech_reference(
    reference_data, "lab_value", "gold_value", reference_predictors,
    method = "normal", m = 20, seed = 176
  ), error = identity),
  # Experimental sensitivity scenario: held-out package studies found poor
  # interval coverage for partial-reference PMM, even with random accurate gold.
  # The algorithm remains standard PMM; this is a conservative package policy,
  # not a claim that every PMM use is invalid.
  reference_pmm = tryCatch(catmech::catmech_reference(
    reference_data, "lab_value", "gold_value", reference_predictors,
    method = "pmm", donors = 5, m = 20, seed = 177
  ), error = identity)
)
# This function explicitly computes sensitivity output. Direct pooled inference
# for partial-reference PMM is refused unless allow_sensitivity = TRUE is given
# to catmech_pool_lm(); its interpretation remains "sensitivity_only".
reference_comparison <- catmech::catmech_sensitivity(
  reference_scenarios, formula = reference_formula, term = "lab_value",
  expected_scenarios = reference_expected
)
print(reference_comparison$table)
cat("Scenario-set status:", reference_comparison$status, "\n")
cat(reference_comparison$interpretation, "\n")

# Additional measured information changes the identification assumptions;
# normality, reference selection and donor support still matter. One synthetic
# result is not evidence of general coverage or that reference PMM is superior.
# Raw lab values remain available in every successful fit's $original. In its
# $mids adapter, observed target values are gold anchors, with all others missing.
# To export a successful reference result to a new or empty directory:
# catmech::write_catmech(reference_scenarios$reference_normal, "output/reference-example")
