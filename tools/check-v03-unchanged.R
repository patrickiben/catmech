#!/usr/bin/env Rscript
# Confirms the v0.4 additions leave existing workflows numerically unchanged. Runs the
# same calls on the frozen v0.3.1 sources and the current sources in separate processes,
# then compares completions, audits and pooled results exactly.
# Usage: Rscript tools/check-v03-unchanged.R [OUTPUT_FILE]
args <- commandArgs(trailingOnly = TRUE)
file_arg <- grep("^--file=", commandArgs(), value = TRUE)
root <- normalizePath(file.path(dirname(sub("^--file=", "", file_arg[1])), ".."))
output <- if (length(args)) args[1] else file.path(root, "validation", "v0.3-unchanged-check.txt")
worker <- tempfile(fileext = ".R")
writeLines(c(
  "args <- commandArgs(TRUE)",
  "pkgload::load_all(args[1], quiet = TRUE)",
  "predictors <- c('age', 'inflammation', 'outcome')",
  "out <- list()",
  "for (seed in c(11L, 12L, 13L)) {",
  "  d <- simulate_catmech_lab(200, seed = seed, corruption_probability = .5, category_error_rate = .1)",
  "  d$gold <- NA_real_",
  "  d$gold[seq(1L, 200L, 3L)] <- d$lab_true[seq(1L, 200L, 3L)]",
  "  d[c('lab_true', 'corrupted')] <- NULL",
  "  labels <- c('Low', 'Medium', 'High')",
  "  H <- matrix(.05, 3, 3, dimnames = list(labels, labels)); diag(H) <- .9",
  "  fits <- list(",
  "    model = catmech_model(d, 'lab_value', predictors, flag = 'outlier_flag', m = 5, seed = seed),",
  "    noisy = catmech_model(d, 'lab_value', predictors, flag = 'outlier_flag', category = 'lab_category', thresholds = c(-Inf, 5.5, 7, Inf),",
  "      levels = labels, category_mode = 'noisy', confusion = H, m = 5, seed = seed),",
  "    pmm = catmech_pmm(d, 'lab_value', predictors, flag = 'outlier_flag', m = 5, seed = seed),",
  "    reference_normal = catmech_reference(d, 'lab_value', 'gold', predictors, m = 5, seed = seed),",
  "    reference_pmm = catmech_reference(d, 'lab_value', 'gold', predictors, method = 'pmm', m = 5, seed = seed),",
  "    legacy = catmech_impute(d, catmech_spec('lab_value', 'lab_category', labels, c(-Inf, 5.5, 7, Inf), flag = 'outlier_flag',",
  "      corruption_probability = .5, conflict_policy = 'reconstruct'), predictors, m = 5, seed = seed))",
  "  for (name in names(fits)) {",
  "    fit <- fits[[name]]",
  "    out[[paste(name, seed)]] <- list(completed = fit$completed, audit = fit$audit, diagnostics = fit$diagnostics,",
  "      assumption_hash = fit$manifest$assumption_hash, claim_level = fit$manifest$claim_level,",
  "      pooled = catmech_pool_lm(fit, outcome ~ lab_value + age, allow_sensitivity = TRUE))",
  "  }",
  "}",
  "saveRDS(out, args[2], version = 2)"), worker)
run <- function(source_root) {
  path <- tempfile(fileext = ".rds")
  status <- system2(file.path(R.home("bin"), "Rscript"), c(shQuote(worker), shQuote(source_root), shQuote(path)))
  if (status != 0L) stop("Worker failed for ", source_root)
  readRDS(path)
}
old <- run(file.path(root, "validation", "v0.3.1-source"))
new <- run(root)
stopifnot(identical(names(old), names(new)))
same <- vapply(names(old), function(key) identical(old[[key]], new[[key]]), logical(1))
lines <- c(paste("Comparison of v0.3.1 frozen sources with current sources:", sum(same), "of", length(same), "workflow runs identical in completions, audits, diagnostics, assumption hashes, claim levels and pooled results."),
           paste(names(same), ifelse(same, "identical", "DIFFERENT")))
writeLines(lines, output)
cat(lines[1], "\n")
if (!all(same)) stop("A previously released workflow changed numerically.")
