#!/usr/bin/env Rscript
# Exact numerical/assumption compatibility for supported v0.7.2 workflows.
# Usage: Rscript tools/check-v072-unchanged.R OUTPUT_FILE
args <- commandArgs(TRUE)
file_arg <- grep("^--file=", commandArgs(), value = TRUE)
root <- normalizePath(file.path(dirname(sub("^--file=", "", file_arg[1])), ".."))
output <- if (length(args)) args[1] else tempfile(fileext = ".txt")
worker <- tempfile(fileext = ".R")
writeLines(c(
  "args <- commandArgs(TRUE)",
  "pkgload::load_all(args[1], quiet=TRUE)",
  "out <- list()",
  "for (seed in 881:883) {",
  "  set.seed(seed)",
  "  d <- data.frame(site=rep(1:30, each=5), x=rnorm(150))",
  "  truth <- 1+.6*d$x+rnorm(30)[d$site]+rnorm(150)",
  "  d$recorded <- .8*truth+rnorm(30)[d$site]+rnorm(150)",
  "  d$gold <- ifelse(d$site<=15, truth, NA_real_)",
  "  for (scenario in c('plain','complete_W','separate_sites','no_recipient_W','agreement','full_gold')) {",
  "    z <- d",
  "    if (scenario=='separate_sites') z$recorded[z$site>=23] <- NA_real_",
  "    if (scenario=='no_recipient_W') z$recorded[z$site>15] <- NA_real_",
  "    if (scenario=='agreement') z$recorded <- truth",
  "    if (scenario=='full_gold') z$gold <- truth",
  "    fit <- suppressWarnings(catmech_reference(z, 'recorded', 'gold', 'x',",
  "      use_recorded=TRUE, cluster=if(scenario=='plain') NULL else 'site', m=10, seed=seed))",
  "    out[[paste(scenario, seed)]] <- list(completed=fit$completed, audit=fit$audit,",
  "      diagnostics=fit$diagnostics, model_diagnostics=fit$model_diagnostics, spec=fit$spec,",
  "      assumption_hash=fit$manifest$assumption_hash, claim_level=fit$manifest$claim_level,",
  "      pooled=suppressWarnings(catmech_pool_lm(fit, recorded~x, cluster='site')))",
  "  }",
  "}",
  "saveRDS(out, args[2], version=2)"), worker)
run <- function(source_root) {
  path <- tempfile(fileext = ".rds")
  on.exit(unlink(path))
  status <- system2(file.path(R.home("bin"), "Rscript"),
    c(shQuote(worker), shQuote(source_root), shQuote(path)))
  if (status != 0L) stop("Worker failed for ", source_root)
  readRDS(path)
}
old <- run(file.path(root, "validation", "v0.7.2-source"))
new <- run(root)
unlink(worker)
stopifnot(identical(names(old), names(new)))
same <- vapply(names(old), function(key) identical(old[[key]], new[[key]]), logical(1))
lines <- c(paste("Frozen v0.7.2 versus current:", sum(same), "of", length(same),
                 "supported workflow runs identical."),
           paste(names(same), ifelse(same, "identical", "DIFFERENT")))
writeLines(lines, output)
cat(lines[1], "\n")
if (!all(same)) stop("A supported default reference workflow changed.")
