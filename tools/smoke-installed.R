# Run in a fresh process after installing the built source archive:
# Rscript tools/smoke-installed.R /path/to/isolated/R-library
args <- commandArgs(trailingOnly = TRUE)
if (length(args)) .libPaths(c(args[1], .libPaths()))
stopifnot(!"package:catmech" %in% search())
stopifnot(utils::packageVersion("catmech") >= "0.2.0")

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

# Run both bundled examples. The main model example uses namespace calls;
# the older lab example intentionally attaches the package.
model_example <- new.env(parent = globalenv())
source(system.file("examples", "model.R", package = "catmech"), local = model_example)
stopifnot(!"package:catmech" %in% search(),
          inherits(model_example$scenarios$exact_categories, "error"),
          !model_example$comparison$complete,
          identical(model_example$comparison$status, "unresolved_incomplete_set"))
source(system.file("examples", "lab.R", package = "catmech"), local = new.env(parent = globalenv()))
cat("Bundled main-model and legacy examples: PASS\n")
