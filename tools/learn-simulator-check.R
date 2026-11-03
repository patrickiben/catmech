# Runs the audit simulator's studies through the R package, so docs/learn/qa/compare.js can check
# the browser engine against it. Usage, from the repository root:
#   node docs/learn/qa/make-datasets.js STUDIES.json
#   Rscript tools/learn-simulator-check.R STUDIES.json
# Writes docs/learn/qa/r-results.csv, r-exact.json, r-draw.json and r-session.txt.
suppressPackageStartupMessages({
  if (requireNamespace("pkgload", quietly = TRUE) && file.exists("DESCRIPTION")) pkgload::load_all(".", quiet = TRUE) else library(catmech)
})
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: Rscript tools/learn-simulator-check.R STUDIES.json", call. = FALSE)
studies <- jsonlite::fromJSON(args[[1]], simplifyVector = TRUE, simplifyDataFrame = FALSE)
out <- "docs/learn/qa"

frame <- function(s) data.frame(age = s$age, sbp = s$sbp, recorded = s$recorded, truth = s$truth,
  checked = ifelse(s$audited, s$truth, NA_real_))
ols <- function(model, d) {
  fit <- stats::lm(model, data = d)
  ci <- stats::confint(fit)[2, ]
  c(estimate = unname(stats::coef(fit)[2]), se = unname(sqrt(diag(stats::vcov(fit)))[2]),
    df = fit$df.residual, lower = unname(ci[1]), upper = unname(ci[2]))
}
imputed <- function(d, recorded, m, seed) {
  fit <- suppressWarnings(catmech_reference(d, "recorded", "checked", predictors = c("age", "sbp"),
    use_recorded = recorded, m = m, seed = seed))
  pooled <- suppressWarnings(catmech_pool_lm(fit, sbp ~ recorded + age))
  list(fit = fit, row = pooled[pooled$term == "recorded", ])
}

rows <- list()
for (s in studies) {
  d <- frame(s)
  audited <- d[s$audited, ]
  add <- function(method, v) rows[[length(rows) + 1L]] <<- data.frame(id = s$id, error = s$error,
    selection = s$selection, method = method, estimate = v[["estimate"]], se = v[["se"]], df = v[["df"]],
    lower = v[["lower"]], upper = v[["upper"]])
  add("trust", ols(sbp ~ recorded + age, d))
  add("audited", ols(sbp ~ truth + age, audited))
  add("oracle", ols(sbp ~ truth + age, d))
  for (method in c("default", "linear")) {
    r <- imputed(d, method == "linear", 20L, 1000L + s$id)$row
    add(method, c(estimate = r$estimate, se = r$std.error, df = r$df, lower = r$lower, upper = r$upper))
  }
}
results <- do.call(rbind, rows)
utils::write.csv(results, file.path(out, "r-results.csv"), row.names = FALSE)

# Exact case: the package's own completed copies, so the browser's analysis and pooling can be checked
# to rounding error.
exact <- lapply(c(0L, 40L), function(id) {
  s <- studies[[id + 1L]]
  d <- frame(s)
  lapply(stats::setNames(c("default", "linear"), c("default", "linear")), function(method) {
    r <- imputed(d, method == "linear", 5L, 77L + id)
    list(id = id, copies = lapply(r$fit$completed, function(z) z$recorded),
      pooled = list(estimate = r$row$estimate, se = r$row$std.error, df = r$row$df, lower = r$row$lower, upper = r$row$upper))
  })
})
writeLines(jsonlite::toJSON(exact, digits = NA, auto_unbox = TRUE), file.path(out, "r-exact.json"))

# One deterministic draw through the package's own internal steps (catmech_model_design, catmech_model_fit and
# the draw in catmech_model_parameter_draw), with fixed chi-square, coefficient and row deviates, so the
# browser's fillIn can be checked exactly rather than on average.
ns <- asNamespace("catmech")
s <- studies[[1L]]
d <- frame(s)
aud <- as.logical(s$audited)
X <- cbind(1, d$age, d$sbp, d$recorded)
design <- get("catmech_model_design", ns)(X[aud, , drop = FALSE])
fit <- get("catmech_model_fit", ns)(design, d$truth[aud])
chi <- 0.9 * design$df
z <- c(0.3, -1.1, 0.7, 0.2)
e <- sin(seq_len(sum(!aud)))
sigma <- sqrt(fit$sse / chi)
noise <- numeric(design$p)
noise[design$qr$pivot] <- backsolve(design$R, z)
beta <- fit$beta + sigma * noise
draw <- list(id = 0L, beta_hat = fit$beta, sse = fit$sse, df = design$df, R = design$R, pivot = design$qr$pivot,
  chisq = chi, z = z, e = e, expected = as.vector(X[!aud, , drop = FALSE] %*% beta) + sigma * e)
writeLines(jsonlite::toJSON(draw, digits = NA, auto_unbox = TRUE), file.path(out, "r-draw.json"))
writeLines(utils::capture.output(utils::sessionInfo()), file.path(out, "r-session.txt"))
cat("Wrote", nrow(results), "rows and two exact cases to", out, "\n")
