# Run in a fresh process after installing the built source archive:
# Rscript tools/smoke-installed.R /path/to/isolated/R-library
args <- commandArgs(trailingOnly = TRUE)
if (length(args)) .libPaths(c(args[1], .libPaths()))
stopifnot(!"package:catmech" %in% search())
stopifnot(utils::packageVersion("catmech") >= "0.3.1")
stopifnot(all(c("catmech_pmm", "catmech_reference") %in%
                getNamespaceExports("catmech")))

d <- catmech::simulate_catmech_lab(200, seed = 7, corruption_probability = 0.5)
d <- d[setdiff(names(d), c("lab_true", "corrupted"))]
predictors <- c("age", "sex", "inflammation", "outcome")
analysis_formula <- outcome ~ lab_value + age + sex + inflammation
flagged_observed <- !is.na(d$lab_value) & as.character(d$outlier_flag) == "Yes"
stopifnot(any(flagged_observed))

check_alignment <- function(fit, original) {
  stopifnot(identical(fit$original, original))
  for (j in seq_along(fit$completed)) {
    adapter <- mice::complete(fit$mids, j)
    stopifnot(identical(adapter$lab_value, fit$completed[[j]]$lab_value))
    stopifnot(identical(adapter$lab_value[flagged_observed],
                        fit$completed[[j]]$lab_value[flagged_observed]))
    stopifnot(identical(fit$completed[[j]]$lab_category, original$lab_category))
  }
  invisible(TRUE)
}

check_export <- function(fit, label) {
  destination <- tempfile(paste0("catmech-installed-", label, "-"))
  catmech::write_catmech(fit, destination)
  restored <- readRDS(file.path(destination, "result.rds"))
  stopifnot(identical(restored$original, fit$original),
            identical(restored$completed, fit$completed),
            identical(restored$model_diagnostics, fit$model_diagnostics))
  manifest <- jsonlite::read_json(file.path(destination, "manifest.json"),
                                  simplifyVector = FALSE)
  stopifnot(identical(manifest$assumption_hash, catmech::catmech_hash(fit$spec)))
  canonical <- paste(readLines(file.path(destination, "assumptions.canonical.json")),
                     collapse = "\n")
  stopifnot(identical(digest::digest(canonical, algo = "sha256", serialize = FALSE),
                     manifest$assumption_hash))
  invisible(manifest)
}

# Principal model, namespace-only, including observed flagged cells in mids.
primary <- catmech::catmech_model(d, "lab_value", predictors,
  flag = "outlier_flag", m = 5, seed = 22)
check_alignment(primary, d)
stopifnot(primary$model_diagnostics$convergence_pass,
          all(primary$audit$reconstructed))
pooled <- catmech::catmech_pool_lm(primary, analysis_formula)
adapter_pooled <- summary(mice::pool(with(primary$mids,
  lm(outcome ~ lab_value + age + sex + inflammation))))
stopifnot(all(is.finite(pooled$estimate)),
  isTRUE(all.equal(pooled$estimate, adapter_pooled$estimate, tolerance = 1e-10)))
primary_manifest <- check_export(primary, "ignore")
stopifnot(identical(primary_manifest$engine, "gaussian_model"),
          "thresholds" %in% names(primary_manifest$assumptions),
          "confusion" %in% names(primary_manifest$assumptions),
          is.null(primary_manifest$assumptions$thresholds),
          is.null(primary_manifest$assumptions$confusion))

# Exercise both Gibbs likelihood modes on known, uncorrupted category labels.
category_levels <- c("Low", "Medium", "High")
category_arguments <- list(data = d, value = "lab_value", predictors = predictors,
  flag = "outlier_flag", category = "lab_category", levels = category_levels,
  thresholds = c(-Inf, 5.5, 7, Inf), m = 4, seed = 24,
  burnin = 200, iterations = 400, chains = 4)
exact <- do.call(catmech::catmech_model,
                 c(category_arguments, list(category_mode = "exact")))
check_alignment(exact, d)
stopifnot(exact$model_diagnostics$convergence_pass,
          all(exact$diagnostics$threshold_violations == 0))
confusion <- matrix(.05, 3, 3, dimnames = list(category_levels, category_levels))
diag(confusion) <- .9
noisy <- do.call(catmech::catmech_model,
  c(category_arguments, list(category_mode = "noisy", confusion = confusion)))
check_alignment(noisy, d)
stopifnot(noisy$model_diagnostics$convergence_pass)
noisy_manifest <- check_export(noisy, "noisy")
stopifnot(identical(noisy_manifest$assumptions$thresholds[[1]], "-Inf"),
          identical(noisy_manifest$assumptions$thresholds[[4]], "Inf"),
          identical(names(noisy_manifest$assumptions$confusion), category_levels))
for (level in category_levels) {
  row <- noisy_manifest$assumptions$confusion[[level]]
  stopifnot(identical(names(row), category_levels),
    isTRUE(all.equal(unname(unlist(row)), unname(confusion[level, ]))))
}

# A captured API error remains an unresolved member of the finite scenario set.
unsupported <- tryCatch(catmech::catmech_model(d, "lab_value", predictors,
  flag = "outlier_flag", category_mode = "ignore", confusion = confusion), error = identity)
comparison <- catmech::catmech_sensitivity(
  list(ignore = primary, exact = exact, noisy = noisy, invalid = unsupported),
  formula = analysis_formula, term = "lab_value",
  expected_scenarios = c("ignore", "exact", "noisy", "invalid")
)
stopifnot(inherits(unsupported, "error"), !comparison$complete,
          identical(comparison$status, "unresolved_incomplete_set"),
          identical(comparison$table$status[comparison$table$scenario == "invalid"], "failed"))
stopifnot(!"package:catmech" %in% search())

# Keep the legacy path tested and require explicit pooling opt-in.
spec <- catmech::read_catmech_spec(system.file("specs", "lab.yaml", package = "catmech"))
legacy <- catmech::catmech_impute(d, spec, predictors, m = 3, seed = 22)
check_alignment(legacy, d)
stopifnot(all(legacy$diagnostics$threshold_violations == 0))
refused <- tryCatch(catmech::catmech_pool_lm(legacy, analysis_formula), error = identity)
stopifnot(inherits(refused, "error"))
legacy_pooled <- catmech::catmech_pool_lm(legacy, analysis_formula, allow_sensitivity = TRUE)
stopifnot(all(legacy_pooled$interpretation == "sensitivity_only"))
check_export(legacy, "legacy")
cat("Installed namespace-only model, legacy, pooling and JSON/RDS checks: PASS\n")

# Standard PMM, including category as an explicitly chosen ordinary predictor.
trusted_pmm <- catmech::catmech_pmm(d, "lab_value", c(predictors, "lab_category"),
  flag = "outlier_flag", donors = 5, m = 4, seed = 35)
check_alignment(trusted_pmm, d)
trusted_rows <- which(!is.na(d$lab_value) & as.character(d$outlier_flag) != "Yes")
stopifnot(all(trusted_pmm$audit$eligible_donors == length(trusted_rows)),
          all(trusted_pmm$audit$candidate_pool_size == 5L),
          identical(trusted_pmm$manifest$engine, "trusted_pmm"))
for (j in seq_along(trusted_pmm$completed)) {
  stopifnot(all(trusted_pmm$completed[[j]]$lab_value[trusted_pmm$manifest$imputed_rows] %in%
                  d$lab_value[trusted_rows]))
}
check_export(trusted_pmm, "trusted-pmm")

# Reference measurements are exposed only in a randomly selected 25% subset.
reference_raw <- catmech::simulate_catmech_lab(200, seed = 7, corruption_probability = 0.5)
set.seed(36)
gold_rows <- sort(sample.int(nrow(reference_raw), nrow(reference_raw) / 4))
reference_data <- reference_raw[setdiff(names(reference_raw), c("lab_true", "corrupted"))]
reference_data$gold_value <- NA_real_
reference_data$gold_value[gold_rows] <- reference_raw$lab_true[gold_rows]
rm(reference_raw)
known <- !is.na(reference_data$gold_value)
reference_fits <- lapply(c("normal", "pmm"), function(method) {
  fit <- catmech::catmech_reference(reference_data, "lab_value", "gold_value", predictors,
    method = method, m = 4, seed = 37)
  check_alignment(fit, reference_data)
  stopifnot(identical(fit$mids$data$lab_value, reference_data$gold_value),
            identical(fit$manifest$imputed_rows, which(!known)),
            identical(fit$manifest$reference_rows, which(known)),
            all(fit$audit$reason[fit$audit$reference_observed] == "gold_reference"),
            all(!fit$audit$reconstructed[fit$audit$reference_observed]),
            all(fit$audit$reconstructed[!fit$audit$reference_observed]))
  for (j in seq_along(fit$completed)) {
    stopifnot(identical(fit$completed[[j]]$lab_value[known], reference_data$gold_value[known]),
              identical(fit$completed[[j]]$gold_value, reference_data$gold_value))
    if (method == "pmm")
      stopifnot(all(fit$completed[[j]]$lab_value[!known] %in% reference_data$gold_value[known]))
  }
  if (method == "pmm") {
    refused <- tryCatch(catmech::catmech_pool_lm(fit, analysis_formula), error = identity)
    stopifnot(inherits(refused, "error"),
              grepl("allow_sensitivity", conditionMessage(refused), fixed = TRUE))
    analysis <- catmech::catmech_pool_lm(fit, analysis_formula, allow_sensitivity = TRUE)
    stopifnot(all(analysis$interpretation == "sensitivity_only"),
              identical(fit$manifest$claim_level, "sensitivity_only"))
  } else analysis <- catmech::catmech_pool_lm(fit, analysis_formula)
  stopifnot(all(is.finite(analysis$estimate)))
  manifest <- check_export(fit, paste0("reference-", method))
  stopifnot(identical(manifest$engine, paste0("reference_", method)),
            "donors" %in% names(manifest$assumptions))
  if (method == "normal") stopifnot(is.null(manifest$assumptions$donors))
  else stopifnot(identical(manifest$claim_level, "sensitivity_only"))
  fit
})
names(reference_fits) <- c("normal", "pmm")

# A completely measured reference endpoint has no stochastic target cells.
# The pooled coefficient, SE, df and t interval equal a direct gold-data lm.
complete_reference <- catmech::simulate_catmech_lab(80, seed = 38)
complete_reference$gold_value <- complete_reference$lab_true
complete_reference[c("lab_true", "corrupted")] <- NULL
gold_analysis <- complete_reference
gold_analysis$lab_value <- gold_analysis$gold_value
ordinary <- stats::lm(analysis_formula, data = gold_analysis)
ordinary_ci <- stats::confint(ordinary)
for (method in c("normal", "pmm")) {
  fit <- catmech::catmech_reference(complete_reference, "lab_value", "gold_value", predictors,
    method = method, m = 3, seed = 39)
  pooled <- catmech::catmech_pool_lm(fit, analysis_formula)
  ord <- match(pooled$term, names(stats::coef(ordinary)))
  stopifnot(length(fit$manifest$imputed_rows) == 0L,
            !any(pooled$interpretation == "sensitivity_only"),
            identical(mice::complete(fit$mids, 2)$lab_value, complete_reference$gold_value),
            isTRUE(all.equal(pooled$estimate, unname(stats::coef(ordinary)[ord]))),
            isTRUE(all.equal(pooled$std.error, unname(sqrt(diag(stats::vcov(ordinary)))[ord]))),
            all(pooled$df == stats::df.residual(ordinary)),
            isTRUE(all.equal(pooled$lower, unname(ordinary_ci[ord, 1]))),
            isTRUE(all.equal(pooled$upper, unname(ordinary_ci[ord, 2]))))
  manifest <- check_export(fit, paste0("complete-reference-", method))
  stopifnot("imputed_rows" %in% names(manifest),
            is.list(manifest$imputed_rows), length(manifest$imputed_rows) == 0L)
}
stopifnot(!"package:catmech" %in% search())
cat("Installed namespace-only PMM, reference, adapter, audit, export and complete-reference checks: PASS\n")

# Run all bundled examples. The model/reference examples use namespace calls;
# the older lab example intentionally attaches the package.
model_example <- new.env(parent = globalenv())
source(system.file("examples", "model.R", package = "catmech"), local = model_example)
stopifnot(!"package:catmech" %in% search(),
          inherits(model_example$scenarios$exact_categories, "error"),
          !model_example$comparison$complete,
          identical(model_example$comparison$status, "unresolved_incomplete_set"))
reference_example <- new.env(parent = globalenv())
source(system.file("examples", "reference.R", package = "catmech"), local = reference_example)
stopifnot(!"package:catmech" %in% search(), reference_example$reference_comparison$complete,
          !any(c("lab_true", "corrupted") %in% names(reference_example$reference_data)),
          sum(!is.na(reference_example$reference_data$gold_value)) ==
            nrow(reference_example$reference_data) / 4,
          identical(names(reference_example$reference_scenarios), reference_example$reference_expected))
stopifnot(identical(reference_example$reference_comparison$table$interpretation[
  reference_example$reference_comparison$table$scenario == "reference_pmm"], "sensitivity_only"))
source(system.file("examples", "lab.R", package = "catmech"), local = new.env(parent = globalenv()))
contamination_example <- new.env(parent = baseenv())
source(system.file("examples", "contamination.R", package = "catmech"), local = contamination_example)
stopifnot(inherits(contamination_example$contamination_fit, "catmech_result"),
          identical(contamination_example$contamination_fit$manifest$engine, "reference_contamination"))
audit_plan_example <- new.env(parent = baseenv())
source(system.file("examples", "audit-plan.R", package = "catmech"), local = audit_plan_example)
stopifnot(identical(audit_plan_example$audit_plan$status, "no_flags"),
          identical(audit_plan_example$audit_plan$k$approach, c("recorded", "discard")),
          all(diff(audit_plan_example$audit_plan$table$relative_width) < 0),
          audit_plan_example$audit_plan$answer$records <= nrow(audit_plan_example$audit_data))
cat("Bundled main-model, reference, legacy, contamination and audit-plan examples: PASS\n")
