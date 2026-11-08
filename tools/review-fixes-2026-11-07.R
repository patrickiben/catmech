# Reproduces the numbers in validation/REVIEW_FIXES_2026-11-07.md.
#
# Usage: Rscript tools/review-fixes-2026-11-07.R [library] [output.json]
#   library      optional library path holding the catmech build to test
#   output.json  defaults to validation/review-fixes-2026-11-07/results.json
#
# Four checks, each on simulated data with a known answer:
#   1. A logistic slope when the audited value is a predictor, imputed with and
#      without the outcome (200 data sets).
#   2. catmech_compare_strategies() with the predictors its 1.0.0 default chose
#      (outcome left out) and with its new default (outcome kept) (100 data sets).
#   3. The cluster-robust GLM standard error under an estimated dispersion, with
#      the 1.0.0 bread (vcov) and the corrected bread (cov.unscaled): overdispersed
#      counts (dispersion above 1) and a Gamma outcome (dispersion below 1).
#   4. catmech_audit_power() defaults against the double-sampling regression
#      estimator for linear imputation.
#   5. Unclustered catmech_pool_glm() and catmech_pool_coxph() against mice::pool().
# Runtime is several minutes; the contamination model in check 2 dominates.

args <- commandArgs(TRUE)
lib <- if (length(args) >= 1L && nzchar(args[1L])) args[1L] else NULL
out <- if (length(args) >= 2L) args[2L] else file.path("validation", "review-fixes-2026-11-07", "results.json")
suppressPackageStartupMessages(library(catmech, lib.loc = lib))
quiet <- function(expr) suppressWarnings(suppressMessages(expr))
mcse <- function(x) stats::sd(x) / sqrt(length(x))

# Shared design: true value t, recorded w = t + noise, auxiliary x, 40% random audit.
simulate <- function(n, seed, outcome = c("binary", "gaussian"), slope = 0.8) {
  outcome <- match.arg(outcome)
  set.seed(seed)
  t <- stats::rnorm(n)
  x <- 0.3 * t + stats::rnorm(n)
  w <- t + stats::rnorm(n, sd = 0.7)
  y <- if (outcome == "binary") stats::rbinom(n, 1, stats::plogis(-0.3 + slope * t + 0.3 * x)) else
    1 + slope * t + 0.3 * x + stats::rnorm(n)
  list(t = t, data = data.frame(lab = w, lab_ref = ifelse(seq_len(n) %in% sample(n, 0.4 * n), t, NA_real_),
                                x = x, y = y))
}

# 1. Logistic regression with the audited value as a predictor.
glm_rows <- lapply(seq_len(200L), function(r) {
  s <- simulate(600L, 5000L + r, "binary")
  d <- s$data
  gold <- stats::coef(stats::glm(y ~ t + x, family = stats::binomial(), data = cbind(d, t = s$t)))[["t"]]
  naive <- stats::coef(stats::glm(y ~ lab + x, family = stats::binomial(), data = d))[["lab"]]
  fit <- function(predictors) {
    res <- quiet(catmech_reference(d, "lab", "lab_ref", predictors, m = 10, seed = r, use_recorded = TRUE))
    p <- quiet(catmech_pool_glm(res, y ~ lab + x, family = stats::binomial()))
    p[p$term == "lab", c("estimate", "lower", "upper")]
  }
  without <- fit("x")
  with <- fit(c("x", "y"))
  c(gold = gold, naive = naive, without_outcome = without$estimate, with_outcome = with$estimate,
    covered_without = without$lower <= 0.8 && 0.8 <= without$upper,
    covered_with = with$lower <= 0.8 && 0.8 <= with$upper)
})
glm_rows <- do.call(rbind, glm_rows)
check1 <- list(true_log_odds_ratio = 0.8, data_sets = nrow(glm_rows), n = 600L, audit = 0.4, m = 10L,
  mean = as.list(colMeans(glm_rows)), mcse = as.list(apply(glm_rows, 2, mcse)))

# 2. Strategy comparison with a continuous outcome.
cmp_rows <- lapply(seq_len(100L), function(r) {
  d <- simulate(600L, 9000L + r, "gaussian")$data
  slopes <- function(predictors) {
    cmp <- as.data.frame(quiet(catmech_compare_strategies(d, "lab", "lab_ref", y ~ lab + x,
                                                          predictors = predictors, m = 10, seed = r)))
    s <- cmp[cmp$term == "lab", ]
    stats::setNames(s$estimate, s$strategy)
  }
  old <- slopes("x")    # what the 1.0.0 default detected for y ~ lab + x
  new <- slopes(NULL)   # the new default keeps the outcome
  c(validation_subsample = new[["validation_subsample"]],
    naive_trust_recorded = new[["naive_trust_recorded"]],
    linear_mi_old_default = old[["linear_mi_recorded"]],
    linear_mi_new_default = new[["linear_mi_recorded"]],
    contamination_old_default = old[["contamination_model"]],
    contamination_new_default = new[["contamination_model"]])
})
cmp_rows <- do.call(rbind, cmp_rows)
check2 <- list(true_slope = 0.8, data_sets = nrow(cmp_rows), n = 600L, audit = 0.4, m = 10L,
  mean = as.list(colMeans(cmp_rows)), mcse = as.list(apply(cmp_rows, 2, mcse)),
  note = paste("The contamination model assumes some records are exactly right;",
               "here every record has error, so its remaining bias is expected misspecification."))

# 3. Cluster-robust GLM standard error under an estimated dispersion.
set.seed(1)
g <- rep(1:30, each = 20)
xx <- stats::rnorm(600)
mu <- exp(0.5 + 0.3 * xx + stats::rnorm(30, sd = 0.3)[g])
yy <- stats::rnbinom(600, mu = mu, size = 1.5)
old_bread_se <- function(f) {
  X <- stats::model.matrix(f)
  meat <- crossprod(rowsum(X * (f$weights * stats::residuals(f, "working")), g))
  v <- stats::vcov(f) %*% meat %*% stats::vcov(f) * (30 / 29) * (599 / 598)
  sqrt(v[2L, 2L])
}
mu_gamma <- exp(0.5 + 0.3 * xx + stats::rnorm(30, sd = 0.3)[g])
y_gamma <- stats::rgamma(600, shape = 3.7, rate = 3.7 / mu_gamma)
fits3 <- list(poisson = stats::glm(yy ~ xx, family = stats::poisson()),
              quasipoisson = stats::glm(yy ~ xx, family = stats::quasipoisson()),
              gamma_log = stats::glm(y_gamma ~ xx, family = stats::Gamma(link = "log")))
check3 <- lapply(fits3, function(f) {
  list(dispersion = summary(f)$dispersion, se_1_0_0 = old_bread_se(f),
       se_fixed = sqrt(catmech:::catmech_glm_cluster_vcov(f, g)[2L, 2L]))
})

# 4. Audit-size defaults against two-phase regression theory.
pwr <- catmech_audit_power(N = 1000, seed = 1)
k_linear <- (1 - 0.85^2) * (1 - 0.30)
f <- pwr$recommended$fraction
check4 <- list(recommended_fraction = f, expected_ratio = pwr$recommended$expected_ratio,
  k_default = pwr$k_summary$estimate, k_linear = k_linear,
  linear_ratio_at_recommended = sqrt(1 + k_linear * (1 - f) / f),
  linear_fraction_for_1_25 = k_linear / (1.25^2 - 1 + k_linear))

# 5. Unclustered GLM and Cox pooling against mice::pool().
set.seed(3)
n5 <- 400
t5 <- stats::rnorm(n5)
d5 <- data.frame(lab = t5 + stats::rnorm(n5, sd = 0.5), ref = ifelse(seq_len(n5) <= 150, t5, NA_real_),
                 x = stats::rnorm(n5), y = stats::rbinom(n5, 1, stats::plogis(0.5 * t5)),
                 tm = stats::rexp(n5, exp(0.4 * t5)), st = stats::rbinom(n5, 1, 0.7))
d5$hz <- mice::nelsonaalen(d5, tm, st)
res5 <- quiet(catmech_reference(d5, "lab", "ref", c("x", "y", "st", "hz"), m = 10, seed = 4, use_recorded = TRUE))
pg <- catmech_pool_glm(res5, y ~ lab + x)
mg <- summary(mice::pool(mice::as.mira(lapply(res5$completed, function(z)
  stats::glm(y ~ lab + x, family = stats::binomial(), data = z)))))
pc <- catmech_pool_coxph(res5, survival::Surv(tm, st) ~ lab + x)
mc <- summary(mice::pool(mice::as.mira(lapply(res5$completed, function(z)
  survival::coxph(survival::Surv(tm, st) ~ lab + x, data = z)))))
check5 <- list(glm_max_abs_diff = max(abs(c(pg$estimate - mg$estimate, pg$std.error - mg$std.error, pg$df - mg$df))),
               coxph_max_abs_diff = max(abs(c(pc$estimate - mc$estimate, pc$std.error - mc$std.error, pc$df - mc$df))),
               mice_version = as.character(utils::packageVersion("mice")))

git <- tryCatch(system2("git", c("describe", "--always", "--dirty", "--abbrev=12"), stdout = TRUE, stderr = FALSE),
                error = function(e) NA_character_, warning = function(w) NA_character_)
result <- list(script = "tools/review-fixes-2026-11-07.R", source_tree = if (length(git)) git[1L] else NA_character_,
  catmech_version = as.character(utils::packageVersion("catmech", lib.loc = lib)),
  r_version = paste(R.version$major, R.version$minor, sep = "."),
  glm_outcome_in_imputation = check1, compare_strategies_default = check2,
  glm_cluster_dispersion = check3, audit_power_defaults = check4, pooling_vs_mice = check5)
dir.create(dirname(out), recursive = TRUE, showWarnings = FALSE)
jsonlite::write_json(result, out, auto_unbox = TRUE, digits = 6, pretty = TRUE)
cat("Wrote", out, "\n")
str(result[-(1:4)], digits.d = 4)
