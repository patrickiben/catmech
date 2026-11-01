#!/usr/bin/env Rscript
# Compare supported scalar pooling with frozen v0.3.0 and exercise v0.3.1 guards.
# The frozen function is sourced into an isolated environment; completions are
# shared between old/new pooling calls, so this isolates the inference changes.
args <- commandArgs(trailingOnly = TRUE)
out <- if (length(args)) args[[1L]] else "validation/pooling-guards"
if (dir.exists(out) && length(list.files(out, all.files = TRUE, no.. = TRUE)))
  stop("Refusing to overwrite a nonempty evidence directory: ", out)
frozen_file <- "validation/v0.3.0-source/R/analysis.R"
stopifnot(file.exists("DESCRIPTION"), file.exists(frozen_file))
pkgload::load_all(".", quiet = TRUE)
source_files <- c(sort(list.files("R", pattern = "\\.R$", full.names = TRUE)),
                  "DESCRIPTION", frozen_file,
                  "validation/v0.3.0-source/SOURCE.txt", "tools/check-pooling-guards.R")
hash_files <- function(paths) vapply(paths, function(p)
  digest::digest(file = p, algo = "sha256"), character(1))
source_before <- hash_files(source_files)
frozen <- new.env(parent = baseenv())
sys.source(frozen_file, envir = frozen)
old_pool <- frozen$catmech_pool_lm

generated <- simulate_catmech_lab(150L, seed = 44031L,
  missing_rate = .2, flag_rate = .1, corruption_probability = 1)
truth <- generated$lab_true
d <- generated[setdiff(names(generated), c("lab_true", "corrupted"))]
d$gold_value <- NA_real_
anchor <- seq(1L, nrow(d), by = 3L)
d$gold_value[anchor] <- truth[anchor]
predictors <- c("age", "sex", "inflammation", "outcome")
full <- d
full$lab_value <- truth
full$gold_value <- truth
full$outlier_flag[] <- "No"
spec <- catmech_spec("lab_value", "lab_category", c("Low", "Medium", "High"),
  c(-Inf, 5.5, 7, Inf), flag = "outlier_flag", corruption_probability = 1,
  conflict_policy = "reconstruct")
fits <- list(
  gaussian_partial = catmech_model(d, "lab_value", predictors,
    flag = "outlier_flag", m = 5L, seed = 70L),
  trusted_pmm_partial = catmech_pmm(d, "lab_value", predictors,
    flag = "outlier_flag", m = 5L, seed = 71L),
  reference_normal_partial = catmech_reference(d, "lab_value", "gold_value",
    predictors, method = "normal", m = 5L, seed = 72L),
  reference_pmm_partial = catmech_reference(d, "lab_value", "gold_value",
    predictors, method = "pmm", m = 5L, seed = 73L),
  legacy_optin = catmech_impute(d, spec, predictors, m = 5L, seed = 74L),
  gaussian_complete = catmech_model(full, "lab_value", predictors,
    flag = "outlier_flag", m = 3L, seed = 75L),
  trusted_pmm_complete = catmech_pmm(full, "lab_value", predictors,
    flag = "outlier_flag", m = 3L, seed = 76L),
  reference_normal_complete = catmech_reference(full, "lab_value", "gold_value",
    predictors, method = "normal", m = 3L, seed = 77L),
  reference_pmm_complete = catmech_reference(full, "lab_value", "gold_value",
    predictors, method = "pmm", m = 3L, seed = 78L))
formulas <- list(adjusted = outcome ~ lab_value + age + sex + inflammation,
                 intercept_only = outcome ~ 1,
                 single_column_response = cbind(outcome) ~ lab_value)
numeric_fields <- c("estimate", "std.error", "df", "lower", "upper",
                    "within_variance", "between_variance", "m", "n")
comparisons <- list()
for (name in names(fits)) for (form in names(formulas)) for (level in c(.9, .95)) {
  optin <- name %in% c("legacy_optin", "reference_pmm_partial")
  before <- old_pool(fits[[name]], formulas[[form]], conf_level = level,
                     allow_sensitivity = optin)
  after <- catmech_pool_lm(fits[[name]], formulas[[form]], conf_level = level,
                          allow_sensitivity = optin)
  exact <- identical(before[numeric_fields], after[numeric_fields]) &&
    identical(before$term, after$term)
  delta <- max(abs(unlist(before[numeric_fields], use.names = FALSE) -
                     unlist(after[numeric_fields], use.names = FALSE)))
  comparisons[[length(comparisons) + 1L]] <- data.frame(
    workflow = name, analysis = form, confidence_level = level,
    explicit_sensitivity_optin = optin, coefficient_count = nrow(after),
    max_absolute_difference = delta, bit_identical = exact,
    previous_interpretation = paste(unique(before$interpretation), collapse = ";"),
    current_interpretation = paste(unique(after$interpretation), collapse = ";"))
}
comparison_table <- do.call(rbind, comparisons)

checks <- list()
check <- function(name, pass, detail) {
  checks[[length(checks) + 1L]] <<- data.frame(
    check = name, pass = isTRUE(pass), detail = detail)
}
fails <- function(expr) {
  result <- tryCatch(force(expr), error = identity)
  list(error = inherits(result, "error"),
       message = if (inherits(result, "error")) conditionMessage(result) else "Returned without error")
}
base <- fits$gaussian_partial
invalid <- list(multivariate_response = cbind(outcome, age) ~ lab_value,
                no_coefficients = outcome ~ 0,
                offset_only = outcome ~ 0 + offset(age),
                no_response = ~ lab_value,
                aliased_coefficient = outcome ~ lab_value + I(lab_value))
for (name in names(invalid)) {
  result <- fails(catmech_pool_lm(base, invalid[[name]]))
  check(paste0("reject_", name), result$error, result$message)
}
check("old_multivariate_silent_null_reproduced",
      is.null(old_pool(base, invalid$multivariate_response)),
      "Frozen v0.3.0 returned NULL for a multivariate lm response; v0.3.1 must refuse it explicitly.")
check("old_zero_coefficient_silent_null_reproduced",
      is.null(old_pool(base, invalid$no_coefficients)),
      "Frozen v0.3.0 returned NULL for a zero-coefficient model; v0.3.1 must refuse it explicitly.")
single <- catmech_pool_lm(base, cbind(outcome) ~ lab_value)
ordinary <- catmech_pool_lm(base, outcome ~ lab_value)
check("one_column_response_is_supported_scalar_analysis",
      identical(single, ordinary),
      "R fits a one-column cbind response as a scalar lm; its pooled result equals the ordinary scalar response.")
check("intercept_only_is_supported",
      identical(catmech_pool_lm(base, outcome ~ 1)$term, "(Intercept)"),
      "An intercept is an estimable coefficient and must not be confused with a zero-coefficient model.")
partial_gate <- fails(catmech_pool_lm(fits$reference_pmm_partial, formulas$adjusted))
check("partial_reference_pmm_default_refused", partial_gate$error, partial_gate$message)
partial <- catmech_pool_lm(fits$reference_pmm_partial, formulas$adjusted,
                          allow_sensitivity = TRUE)
check("partial_reference_pmm_optin_marked_sensitivity",
      all(grepl("sensitivity", partial$interpretation, fixed = TRUE)),
      paste(unique(partial$interpretation), collapse = ";"))
legacy_gate <- fails(catmech_pool_lm(fits$legacy_optin, formulas$adjusted))
check("legacy_default_refused", legacy_gate$error, legacy_gate$message)
full_pmm <- catmech_pool_lm(fits$reference_pmm_complete, formulas$adjusted)
gold_lm <- stats::lm(formulas$adjusted, data = full)
gold_summary <- summary(gold_lm)$coefficients
gold_ci <- stats::confint(gold_lm)
position <- match(full_pmm$term, rownames(gold_summary))
check("complete_reference_pmm_unrestricted_complete_data_endpoint",
      isTRUE(all.equal(full_pmm$estimate, unname(gold_summary[position, "Estimate"]), tolerance = 1e-12)) &&
      isTRUE(all.equal(full_pmm$std.error, unname(gold_summary[position, "Std. Error"]), tolerance = 1e-12)) &&
      identical(full_pmm$df, rep(as.numeric(stats::df.residual(gold_lm)), nrow(full_pmm))) &&
      isTRUE(all.equal(full_pmm$lower, unname(gold_ci[position, 1]), tolerance = 1e-12)) &&
      isTRUE(all.equal(full_pmm$upper, unname(gold_ci[position, 2]), tolerance = 1e-12)),
      "All-reference PMM does no imputation and remains available without opt-in, with ordinary gold-data lm uncertainty.")
equal_draws <- fits$reference_pmm_partial
equal_draws$completed <- rep(equal_draws$completed[1L], length(equal_draws$completed))
equal_gate <- fails(catmech_pool_lm(equal_draws, formulas$adjusted))
check("zero_empirical_between_variance_does_not_bypass_gate",
      equal_gate$error, equal_gate$message)
bad_meta <- fits$reference_pmm_partial
bad_meta$manifest$imputed_rows <- integer()
bad_meta$completed[[2L]]$lab_value[2L] <- bad_meta$completed[[2L]]$lab_value[2L] + 1
meta_gate <- fails(catmech_pool_lm(bad_meta, formulas$adjusted))
check("inconsistent_no_imputation_metadata_refused", meta_gate$error, meta_gate$message)
guard_table <- do.call(rbind, checks)

source_after <- hash_files(source_files)
if (!identical(source_before, source_after))
  stop("Sources changed during comparison; no evidence written.")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
utils::write.csv(comparison_table, file.path(out, "scalar-numerical-equivalence.csv"), row.names = FALSE)
utils::write.csv(guard_table, file.path(out, "guard-controls.csv"), row.names = FALSE)
utils::write.csv(data.frame(file = source_files, sha256 = unname(source_before)),
                 file.path(out, "source-hashes.csv"), row.names = FALSE)
capture.output(sessionInfo(), file = file.path(out, "sessionInfo.txt"))
passed <- all(comparison_table$bit_identical) && all(guard_table$pass)
writeLines(c(
  "# Pooling guards and frozen-source numerical equivalence", "",
  paste0("Overall result: **", if (passed) "PASS" else "FAIL", "**."), "",
  paste0(nrow(comparison_table), " supported scalar comparisons use identical completion objects in the frozen v0.3.0 pool and current pool."),
  "They cover all five workflow engines, four complete-data endpoints, three supported scalar analysis forms, and 90%/95% confidence levels.",
  "The comparison requires identical coefficient names and bit-identical estimate, standard error, degrees of freedom, limits, within/between variance, m, and n.",
  "Partial reference PMM and the legacy engine receive explicit sensitivity opt-in; changed interpretation labels are recorded separately.", "",
  "Guard checks reproduce the old silent-NULL multiresponse/zero-coefficient bugs and verify explicit rejection.",
  "A single-column response and an intercept-only analysis remain supported. Partial-reference PMM is refused by default,",
  "even with zero empirical between-imputation variance. Fully referenced PMM remains ordinary complete-data analysis.",
  "The frozen source is from commit ea5ec6ae4688c5ff8a2f2086cc6f45e69cb5f21c; all current/frozen source hashes were unchanged during the check.", "",
  "These are engineering and numerical-equivalence checks. They do not establish coverage or rehabilitate the undercoverage of reference PMM."
), file.path(out, "README.md"))
print(comparison_table[, c("workflow", "analysis", "confidence_level", "bit_identical")])
print(guard_table[, c("check", "pass")])
if (!passed) stop("One or more pooling guard checks failed; evidence retained.")
