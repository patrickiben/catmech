# Runs the real catmech_audit_plan() on the same CSV and test vectors the page uses, for comparison with the browser.
# Run from the repository root, after installing the package to a scratch library:
#   CATMECH_LIB=/path/to/scratch/lib Rscript docs/audit-planner/qa/r-cross-check.R
library(catmech, lib.loc = Sys.getenv("CATMECH_LIB"))
d <- read.csv("docs/audit-planner/qa/example-pilot-nhanes.csv", na.strings = c("NA", "."), stringsAsFactors = FALSE)
plan <- suppressWarnings(catmech_audit_plan(d, value = "recorded_kg", reference = "reference_kg",
  predictors = c("age", "sex", "race", "height", "sbp"), formula = sbp ~ recorded_kg + age + sex + race + height,
  term = "recorded_kg", target_ratio = 1.25, seed = 1))
out <- list(k = plan$k, table = plan$table[, c("fraction", "records", "additional_records", "relative_width", "relative_lower", "relative_upper", "predictive_lower", "predictive_upper", "discard_relative_width")],
  answer = plan$answer, flags = plan$flags, pilot = plan$pilot[, c("records", "reference_rows", "calibration_rows", "disagreements", "largest_discrepancy_ratio", "largest_recorded_z", "selection_z", "bootstrap_usable")])
print(out)
con <- file("docs/audit-planner/qa/r-cross-check.txt", "w"); sink(con); cat("catmech_audit_plan() on qa/example-pilot-nhanes.csv\n"); print(out); sink(); close(con)
vals <- c(k = plan$k$estimate[1], discard_k = plan$k$estimate[2], classical = plan$k$classical[1], discard_classical = plan$k$classical[2],
  k_lower = plan$k$lower[1], k_upper = plan$k$upper[1], discard_lower = plan$k$lower[2], discard_upper = plan$k$upper[2],
  records = plan$answer$records, records_lower = plan$answer$records_lower, records_upper = plan$answer$records_upper,
  disagreements = plan$pilot$disagreements, selection_z = plan$pilot$selection_z, largest_recorded_z = plan$pilot$largest_recorded_z, bootstrap_usable = plan$pilot$bootstrap_usable)
vals <- c(vals, setNames(plan$table$relative_width, paste0("relative_", plan$table$fraction)), setNames(plan$table$discard_relative_width, paste0("discard_relative_", plan$table$fraction)),
  setNames(plan$table$predictive_lower, paste0("predictive_lower_", plan$table$fraction)), setNames(plan$table$predictive_upper, paste0("predictive_upper_", plan$table$fraction)))
writeLines(sprintf("%s %.15g", names(vals), vals), "docs/audit-planner/qa/r-cross-check-values.txt")
