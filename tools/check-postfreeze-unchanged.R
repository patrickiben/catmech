#!/usr/bin/env Rscript
# The linear-MI option and hold-out advice came after the frozen v0.4 runs. This confirms
# they leave those runs' numerics alone, by running the same calls on the frozen and the
# current sources in separate processes and comparing them exactly.
# Usage: Rscript tools/check-postfreeze-unchanged.R [OUTPUT_FILE]
args <- commandArgs(trailingOnly = TRUE)
file_arg <- grep("^--file=", commandArgs(), value = TRUE)
root <- normalizePath(file.path(dirname(sub("^--file=", "", file_arg[1])), ".."))
output <- if (length(args)) args[1] else file.path(root, "validation", "v0.4.0-postfreeze-check.txt")
worker <- tempfile(fileext = ".R")
writeLines(c(
  "args <- commandArgs(TRUE)",
  "pkgload::load_all(args[1], quiet = TRUE)",
  "predictors <- c('age', 'inflammation', 'outcome')",
  "out <- list()",
  "for (seed in c(21L, 22L, 23L)) {",
  "  d <- simulate_catmech_lab(240, seed = seed, corruption_probability = .5)",
  "  d$gold <- NA_real_",
  "  d$gold[seq(1L, 240L, 3L)] <- d$lab_true[seq(1L, 240L, 3L)]",
  "  d[c('lab_true', 'corrupted')] <- NULL",
  "  noisy <- d",
  "  noisy$gold <- noisy$gold + .1 * sin(seq_len(240))",
  "  fits <- list(",
  "    contamination = catmech_contamination(d, 'lab_value', 'gold', predictors, corruption_strata = 'outlier_flag', m = 5, seed = seed),",
  "    shifted = catmech_contamination(d, 'lab_value', 'gold', predictors, corruption_strata = 'outlier_flag', selection_delta = .3, m = 5, seed = seed),",
  "    noisy = catmech_contamination(noisy, 'lab_value', 'gold', predictors, reference_sd = .1, m = 5, seed = seed, max_iterations = 2400),",
  "    bootstrap = catmech_reference(d, 'lab_value', 'gold', predictors, method = 'pmm_abb', m = 5, seed = seed))",
  "  for (name in names(fits)) {",
  "    fit <- fits[[name]]",
  "    out[[paste(name, seed)]] <- list(completed = fit$completed, audit = fit$audit, diagnostics = fit$diagnostics,",
  "      contamination = fit$contamination, assumption_hash = fit$manifest$assumption_hash,",
  "      pooled = catmech_pool_lm(fit, outcome ~ lab_value + age))",
  "  }",
  "}",
  "saveRDS(out, args[2], version = 2)"), worker)
run <- function(source_root) {
  path <- tempfile(fileext = ".rds")
  status <- system2(file.path(R.home("bin"), "Rscript"), c(shQuote(worker), shQuote(source_root), shQuote(path)))
  if (status != 0L) stop("Worker failed for ", source_root)
  readRDS(path)
}
old <- run(file.path(root, "validation", "v0.4.0-source"))
new <- run(root)
stopifnot(identical(names(old), names(new)))
same <- vapply(names(old), function(key) identical(old[[key]], new[[key]]), logical(1))
lines <- c(paste("Comparison of frozen v0.4.0 study sources with current sources:", sum(same), "of", length(same), "workflow runs identical."),
           paste(names(same), ifelse(same, "identical", "DIFFERENT")))
writeLines(lines, output)
cat(lines[1], "\n")
if (!all(same)) stop("A frozen study workflow changed numerically.")
