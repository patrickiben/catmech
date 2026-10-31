# Run from a fresh R process after installing the built source archive.
args <- commandArgs(trailingOnly=TRUE)
if (length(args)) .libPaths(c(args[1], .libPaths()))
stopifnot(!"package:catmech" %in% search())
d <- catmech::simulate_catmech_lab(200, seed=7)
d <- d[setdiff(names(d),c("lab_true","corrupted"))]
s <- catmech::read_catmech_spec(system.file("specs","lab.yaml",package="catmech"))
fit <- catmech::catmech_impute(d,s,c("age","sex","inflammation","outcome"),m=3,seed=22)
stopifnot(identical(fit$original,d), all(fit$diagnostics$threshold_violations==0))
for (j in seq_along(fit$completed)) {
  stopifnot(identical(mice::complete(fit$mids,j)$lab_value,fit$completed[[j]]$lab_value))
}
pooled <- summary(mice::pool(with(fit$mids,
  lm(outcome ~ lab_value + age + sex + inflammation))))
stopifnot("lab_value" %in% pooled$term, all(is.finite(pooled$estimate)))
destination <- tempfile("catmech-installed-output-")
catmech::write_catmech(fit,destination)
stopifnot(identical(readRDS(file.path(destination,"result.rds"))$original,d))
cat("Installed, namespace-only API and analysis/export smoke checks: PASS\n")
source(system.file("examples","lab.R",package="catmech"))
