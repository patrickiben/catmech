# Runs every R snippet the guide's method picker can show, on a small simulated dataset with whole-site audits.
# Usage, from the repository root:
#   node docs/learn/qa/picker-snippets.js docs/learn/index.html SNIPPETS.json
#   Rscript tools/learn-picker-check.R SNIPPETS.json > docs/learn/qa/picker-check.txt
suppressPackageStartupMessages({
  if (requireNamespace("pkgload", quietly = TRUE) && file.exists("DESCRIPTION")) pkgload::load_all(".", quiet = TRUE) else library(catmech)
})
snips <- jsonlite::fromJSON(commandArgs(TRUE)[1], simplifyDataFrame = FALSE)
set.seed(7)
n <- 600; site <- rep(sprintf("s%02d", 1:30), each = 20)
age <- round(runif(n, 20, 80)); sex <- factor(sample(c("F", "M"), n, TRUE))
truth <- 5 + 0.02 * age + 0.3 * (sex == "M") + rnorm(n, 0, 0.8)
outcome <- 1 + 0.5 * truth + 0.01 * age + rnorm(n)
wrong <- runif(n) < 0.08
lab_value <- truth + ifelse(wrong, rnorm(n, 0.6, 0.2), 0)
audited <- site %in% sprintf("s%02d", seq(1, 30, 3))
dat <- data.frame(site, age, sex, outcome, lab_value, reference_value = ifelse(audited, truth, NA_real_))
sd_checked <- 0.2; sd_guess <- 0.2; delta_guess <- 0.2
ok <- 0
for (s in snips) {
  label <- sprintf("%-13s %-9s sites=%s", s$method, s$role, s$sites)
  res <- tryCatch({ withCallingHandlers(eval(parse(text = s$code), envir = new.env(parent = globalenv())),
    warning = function(w) { message("  warning: ", conditionMessage(w)); invokeRestart("muffleWarning") }); "OK" },
    error = function(e) paste("ERROR:", conditionMessage(e)))
  if (res == "OK") ok <- ok + 1
  cat(sprintf("%s  %s\n", label, res))
}
cat(ok, "of", length(snips), "snippets ran\n")
