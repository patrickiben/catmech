# Numbers for the worked example in docs/learn/index.html: the README's NHANES example, run as written,
# plus the other methods on the same audit. An illustration, not part of any frozen study.
# Usage, from the repository root: Rscript tools/learn-guide-example.R
suppressPackageStartupMessages({
  if (requireNamespace("pkgload", quietly = TRUE) && file.exists("DESCRIPTION")) pkgload::load_all(".", quiet = TRUE) else library(catmech)
})
d <- read.csv("validation/real-data/nhanes-pool.csv")
d <- d[d$split == "evaluation", ]
set.seed(1)
audited <- sample(nrow(d), nrow(d) / 4)
d$audit_kg <- NA_real_
d$audit_kg[audited] <- d$reference_kg[audited]
predictors <- c("age", "sex", "race", "height", "sbp")
model <- sbp ~ recorded_kg + age + sex + race + height
show <- function(label, r) sprintf("%s: %.3f, 95%% interval %.3f to %.3f, width %.3f", label, r[1], r[2], r[3], r[3] - r[2])
pooled <- function(fit) { p <- catmech_pool_lm(fit, model); p <- p[p$term == "recorded_kg", ]; c(p$estimate, p$lower, p$upper) }
ols <- function(data, value) {
  f <- stats::lm(stats::as.formula(sub("recorded_kg", value, deparse(model))), data = data)
  ci <- stats::confint(f)[value, ]
  c(unname(stats::coef(f)[value]), unname(ci[1]), unname(ci[2]))
}
gap <- abs(d$recorded_kg[audited] - d$reference_kg[audited])
both <- !is.na(gap)
measured <- transform(d, measured_kg = reference_kg)
lines <- c(
  sprintf("Evaluation pool: %d people, %d with a self-reported weight; audited: %d", nrow(d), sum(!is.na(d$recorded_kg)), length(audited)),
  sprintf("Audited people with both weights whose self-report differs from the measurement by more than 1 kg: %d of %d (%.1f%%)", sum(gap[both] > 1), sum(both), 100 * mean(gap[both] > 1)),
  sprintf("Median absolute difference on those rows: %.1f kg", stats::median(gap[both])),
  show("Fill in from the recorded value (catmech_reference, use_recorded = TRUE)", pooled(catmech_reference(d, "recorded_kg", "audit_kg", predictors = predictors, use_recorded = TRUE, m = 20, seed = 123))),
  show("Fill in without the recorded value (catmech_reference default)", pooled(catmech_reference(d, "recorded_kg", "audit_kg", predictors = predictors, m = 20, seed = 123))),
  show("Audited rows only", ols(measured[audited, ], "measured_kg")),
  show("Trust the recorded value, rows that have one", ols(d, "recorded_kg")),
  show("Every measured weight", ols(measured, "measured_kg"))
)
writeLines(lines, "docs/learn/qa/guide-example.txt")
writeLines(lines)
