# Checks that the 1.0.0 rewrite of the marginal Gibbs sampler (src/marginal.cpp and the
# chain continuation in R/marginal-sampler.R) reproduces the 0.10.1 draws exactly.
#
# Usage: Rscript tools/review-sampler-equivalence-2026-11-07.R <old library> <new library>
#   Install 0.10.1 (commit 662252a) into the old library and the build to test into the
#   new one. Each library is loaded in its own R process.
#
# The fit is the public fixture of tests/testthat/test-marginal-analysis.R with a short
# first run, so the sampler must double: from 1,000 iterations it converges after
# doubling, and from 100 it is refused at the 800-iteration cap. Completions, posterior
# summaries and every attempt's diagnostics must be identical.

args <- commandArgs(TRUE)
if (length(args) != 2L) stop("Give the old and the new library paths.", call. = FALSE)

worker <- function(lib, iterations, file) {
  code <- sprintf('
    suppressPackageStartupMessages(library(catmech, lib.loc = %s))
    set.seed(822301)
    site <- rep(seq_len(24L), each = 10L)
    a <- stats::rnorm(24L)
    x <- a[site] + stats::rnorm(length(site))
    u <- matrix(stats::rnorm(48L), 24L, 2L) %%*%% chol(matrix(c(.5, .45, .45, .6), 2L))
    e <- matrix(stats::rnorm(480L), 240L, 2L) %%*%% chol(matrix(c(.5, .55, .55, .8), 2L))
    target <- .5 + .7 * x + .4 * a[site] + u[site, 1L] + e[, 1L]
    d <- data.frame(recorded = .8 + 1.05 * x + u[site, 2L] + e[, 2L],
      reference = target, x = x, site = as.character(site))
    d$reference[-sample(which(site <= 18L), 60L)] <- NA_real_
    d$recorded[seq(1L, 240L, 3L)] <- NA_real_
    out <- tryCatch({
      fit <- suppressWarnings(catmech_reference(d, "recorded", "reference", "x", m = 8,
        seed = 822302, use_recorded = TRUE, cluster = "site", cluster_model = "joint_marginal",
        joint_control = list(iterations = %dL, burnin = 100L, max_iterations = %dL)))
      list(completed = fit$completed, summary = catmech_marginal_summary(fit),
           attempts = fit$model_diagnostics$attempts, traces = fit$model_diagnostics$traces)
    }, error = function(e) list(error = conditionMessage(e),
      attempts = e$model_diagnostics$attempts, traces = e$model_diagnostics$traces))
    saveRDS(out, %s)', deparse(lib), iterations, 8L * iterations, deparse(file))
  script <- tempfile(fileext = ".R")
  writeLines(code, script)
  status <- system2(file.path(R.home("bin"), "Rscript"), c("--vanilla", script))
  if (status != 0L) stop("Worker failed for ", lib, call. = FALSE)
  readRDS(file)
}

for (iterations in c(1000L, 100L)) {
  old <- worker(args[1L], iterations, tempfile(fileext = ".rds"))
  new <- worker(args[2L], iterations, tempfile(fileext = ".rds"))
  same <- vapply(names(old), function(k) identical(old[[k]], new[[k]]), logical(1L))
  cat(sprintf("start %d iterations: %s; attempts %s\n", iterations,
              if (all(same) && setequal(names(old), names(new))) "IDENTICAL" else
                paste("DIFFERS in", paste(names(same)[!same], collapse = ", ")),
              paste(old$attempts$iterations, collapse = " -> ")))
}
