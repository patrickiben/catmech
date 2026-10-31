# Run after installing catmech: Rscript inst/examples/model.R
# All records are synthetic. This example demonstrates an analysis conditional
# on stated assumptions; it does not validate clinical use or identify truth.
example_data <- catmech::simulate_catmech_lab(
  n = 160, seed = 7, corruption_probability = 0.5,
  category_error_rate = 0.1
)
example_data <- example_data[setdiff(names(example_data), c("lab_true", "corrupted"))]
example_predictors <- c("age", "sex", "inflammation", "outcome")
example_formula <- outcome ~ lab_value + age + sex + inflammation

# Main analysis: discard every flagged measurement and ignore the category.
# The downstream outcome is included among the imputation predictors.
primary <- catmech::catmech_model(
  example_data, value = "lab_value", predictors = example_predictors,
  flag = "outlier_flag", category_mode = "ignore", m = 20, seed = 203
)
primary_analysis <- catmech::catmech_pool_lm(primary, example_formula)
print(primary_analysis[primary_analysis$term == "lab_value", ])

# Prespecify all scenarios before looking at their estimates. The noisy matrix
# is known only because this is a synthetic example with a specified generator.
# In an application it requires external validation or a sensitivity declaration.
expected <- c("ignore_categories", "known_noisy_categories", "exact_categories")
category_levels <- c("Low", "Medium", "High")
confusion <- matrix(0.05, 3, 3,
  dimnames = list(true = category_levels, observed = category_levels))
diag(confusion) <- 0.9
category_arguments <- list(
  data = example_data, value = "lab_value", predictors = example_predictors,
  flag = "outlier_flag", category = "lab_category", levels = category_levels,
  thresholds = c(-Inf, 5.5, 7, Inf), m = 20, seed = 204
)
scenarios <- list(
  ignore_categories = primary,
  known_noisy_categories = tryCatch(
    do.call(catmech::catmech_model,
      c(category_arguments, list(category_mode = "noisy", confusion = confusion))),
    error = identity
  ),
  exact_categories = tryCatch(
    do.call(catmech::catmech_model,
      c(category_arguments, list(category_mode = "exact"))),
    error = identity
  )
)

# The exact-category assumption contradicts some trusted synthetic measurements.
# Its caught error remains in the table; a failed scenario is not dropped.
comparison <- catmech::catmech_sensitivity(
  scenarios, formula = example_formula, term = "lab_value",
  expected_scenarios = expected
)
print(comparison$table)
cat("Scenario-set status:", comparison$status, "\n")
cat(comparison$interpretation, "\n")

# To export the primary analysis, choose a new or empty directory:
# catmech::write_catmech(primary, "output/model-example")
