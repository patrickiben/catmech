#!/usr/bin/env Rscript
# Independent posterior quadrature check for the Gaussian category Gibbs model.
# No MCMC sampling or package likelihood routines are used by the quadrature.
# Run from the repository root; optional argument is a new output directory.

args <- commandArgs(trailingOnly = TRUE)
directory <- if (length(args)) args[[1L]] else "validation/model-oracle"
if (dir.exists(directory) && length(list.files(directory, all.files = TRUE, no.. = TRUE)))
  stop("Oracle output directory must be new or empty.")
if (!dir.exists(directory) && !dir.create(directory, recursive = TRUE))
  stop("Could not create oracle output directory.")
suppressPackageStartupMessages(pkgload::load_all(".", quiet = TRUE))
sources <- c("R/model.R", "R/spec.R", "R/impute.R", "tools/check-model-oracle.R")
source_hash <- vapply(sources, function(p) digest::digest(file = p, algo = "sha256"), character(1))

# Standard-normal and normalized Gamma(shape, scale=1) quadrature via their
# orthogonal-polynomial Jacobi matrices. Squared first eigenvector entries are
# weights summing to one. This uses only base R eigendecomposition.
normal_rule <- function(order) {
  J <- matrix(0, order, order)
  off <- sqrt(seq_len(order - 1L))
  J[cbind(seq_len(order - 1L), 2:order)] <- off
  J[cbind(2:order, seq_len(order - 1L))] <- off
  e <- eigen(J, symmetric = TRUE)
  list(node = e$values, weight = e$vectors[1L, ]^2)
}
gamma_rule <- function(order, shape) {
  alpha <- shape - 1
  J <- diag(2 * (seq_len(order) - 1) + alpha + 1)
  k <- seq_len(order - 1L)
  off <- sqrt(k * (k + alpha))
  J[cbind(k, k + 1L)] <- off
  J[cbind(k + 1L, k)] <- off
  e <- eigen(J, symmetric = TRUE)
  stopifnot(all(e$values > 0))
  list(node = e$values, weight = e$vectors[1L, ]^2)
}

observed <- c(-1.4, -1.0, -0.7, -0.5, -0.25, -0.1, 0.1, 0.3, 0.55, 0.8, 1.1, 1.5)
n_high <- 6L
n_low <- 2L
observed_mean <- mean(observed)
sse <- sum((observed - observed_mean)^2)
df <- length(observed) - 1L

quadrature <- function(order, use_categories = TRUE) {
  zn <- normal_rule(order)
  gn <- gamma_rule(order, df / 2)
  grid <- expand.grid(i = seq_len(order), j = seq_len(order))
  # Under the trusted-data reference posterior, W=SSE/sigma² is chi-square(df)
  # and Z=sqrt(n)*(mu-ybar)/sigma is independent standard normal.
  sigma <- sqrt(sse / (2 * gn$node[grid$j]))
  mu <- observed_mean + sigma * zn$node[grid$i] / sqrt(length(observed))
  ratio <- mu / sigma
  log_weight <- log(zn$weight[grid$i]) + log(gn$weight[grid$j])
  if (use_categories) log_weight <- log_weight +
    n_high * stats::pnorm(ratio, log.p = TRUE) +
    n_low * stats::pnorm(-ratio, log.p = TRUE)
  weight <- exp(log_weight - max(log_weight))
  weight <- weight / sum(weight)
  mu_mean <- sum(weight * mu)
  # Mean of the eight interval-only latent measurements conditional on mu,sigma.
  positive_mills <- exp(stats::dnorm(ratio, log = TRUE) - stats::pnorm(ratio, log.p = TRUE))
  negative_mills <- exp(stats::dnorm(ratio, log = TRUE) - stats::pnorm(-ratio, log.p = TRUE))
  latent_mean <- mu + sigma * (n_high * positive_mills - n_low * negative_mills) /
    (n_high + n_low)
  c(mu_mean = mu_mean, mu_variance = sum(weight * (mu - mu_mean)^2),
    sigma_squared_mean = sum(weight * sigma^2),
    latent_mean = sum(weight * latent_mean))
}

orders <- c(64L, 128L, 256L)
values <- lapply(orders, quadrature)
oracle <- values[[length(values)]]
quadrature_change <- abs(oracle - values[[length(values) - 1L]])
convergence <- data.frame(quantity = names(oracle),
  order_64 = values[[1L]], order_128 = values[[2L]], order_256 = oracle,
  last_change = quadrature_change, row.names = NULL)
utils::write.csv(convergence, file.path(directory, "quadrature-convergence.csv"), row.names = FALSE)
if (any(quadrature_change > 2e-6)) stop("Independent quadrature did not stabilize.")

without_categories <- quadrature(256L, FALSE)
trusted_exact <- c(mu_mean = observed_mean,
  mu_variance = sse / (length(observed) * (df - 2)),
  sigma_squared_mean = sse / (df - 2))
known_answer_error <- abs(without_categories[names(trusted_exact)] - trusted_exact)
if (any(known_answer_error > 2e-6)) stop("Quadrature fails the conjugate no-category oracle.")

data <- data.frame(value = c(observed, rep(NA_real_, n_high + n_low)),
  band = factor(c(ifelse(observed <= 0, "Low", "High"), rep("High", n_high), rep("Low", n_low)),
                levels = c("Low", "High")))
fit <- catmech_model(data, "value", character(0), category = "band",
  thresholds = c(-Inf, 0, Inf), levels = c("Low", "High"), category_mode = "exact",
  m = 20L, burnin = 1500L, iterations = 5000L, chains = 4L, seed = 20261101L)
traces <- fit$model_diagnostics$traces
mu <- traces[, , "beta:(Intercept)"]
sigma_squared <- traces[, , "sigma"]^2
quantities <- list(mu_mean = mu, mu_variance = (mu - oracle[["mu_mean"]])^2,
                   sigma_squared_mean = sigma_squared, latent_mean = traces[, , "latent_mean"])
estimate <- vapply(quantities, mean, numeric(1))
mcse <- vapply(quantities, posterior::mcse_mean, numeric(1))
checks <- data.frame(quantity = names(oracle), oracle = unname(oracle),
  mcmc = unname(estimate), mcse = unname(mcse),
  quadrature_error_bound = unname(2 * quadrature_change),
  absolute_difference = unname(abs(estimate - oracle)), row.names = NULL)
checks$tolerance <- 5 * checks$mcse + checks$quadrature_error_bound
checks$pass <- is.finite(checks$mcse) & checks$absolute_difference <= checks$tolerance
# Deliberately wrong posterior ignores all interval-only records. It must be
# distinguishable on the mean; otherwise this fixture would not test their use.
negative_control <- data.frame(
  control = "omit_category_likelihood", correct_mu_mean = oracle[["mu_mean"]],
  wrong_mu_mean = without_categories[["mu_mean"]],
  separation_in_mcse = abs(oracle[["mu_mean"]] - without_categories[["mu_mean"]]) / mcse[["mu_mean"]])
negative_control$detected <- negative_control$separation_in_mcse > 10
after_hash <- vapply(sources, function(p) digest::digest(file = p, algo = "sha256"), character(1))
if (!identical(source_hash, after_hash)) stop("Source changed during the oracle run.")
utils::write.csv(checks, file.path(directory, "posterior-moment-checks.csv"), row.names = FALSE)
utils::write.csv(negative_control, file.path(directory, "negative-control.csv"), row.names = FALSE)
utils::write.csv(data.frame(quantity = names(known_answer_error), absolute_error = known_answer_error),
                 file.path(directory, "conjugate-oracle-check.csv"), row.names = FALSE)
utils::write.csv(fit$model_diagnostics$summary, file.path(directory, "chain-diagnostics.csv"), row.names = FALSE)
saveRDS(fit$model_diagnostics, file.path(directory, "chain-diagnostics.rds"), version = 2)
utils::write.csv(data.frame(path = sources, sha256 = source_hash),
                 file.path(directory, "source-hashes.csv"), row.names = FALSE)
writeLines(utils::capture.output(utils::sessionInfo()), file.path(directory, "sessionInfo.txt"))
writeLines(c(
  "# Independent Gaussian category-model posterior check", "",
  "Twelve trusted observations and eight interval-only observations were used in an intercept-only model.",
  "The reference prior was p(mu,sigma²) proportional to 1/sigma². Six interval-only values were positive and two nonpositive.", "",
  "The oracle integrates the trusted conjugate posterior over independent standard-normal and chi-square coordinates, weighted by the exact category likelihood. Base R Gaussian quadrature uses no package likelihood or Gibbs-update code.",
  "Quadrature orders 64, 128 and 256 were compared, and the zero-category case was checked against closed-form posterior moments.", "",
  sprintf("Four Gibbs chains used 1,500 burn-in and 5,000 retained iterations each. All four posterior moment checks passed: %s.", all(checks$pass)),
  sprintf("The known-wrong posterior omitting the category likelihood was detected: %s.", all(negative_control$detected)),
  "The acceptance tolerance is five Monte Carlo standard errors plus twice the last quadrature-order change, fixed before observing the chain results.", "",
  "This checks the sampler's posterior moments for one known example. It does not establish frequentist confidence-interval coverage, robustness to misspecification, or posterior calibration across applications.", "",
  "Reproduce from the repository root: `Rscript tools/check-model-oracle.R validation/new-model-oracle`.",
  "See posterior-moment-checks.csv, quadrature-convergence.csv, conjugate-oracle-check.csv, negative-control.csv and source-hashes.csv."
), file.path(directory, "README.md"))
print(checks, row.names = FALSE)
if (!all(checks$pass) || !all(negative_control$detected)) stop("Posterior oracle check failed.")
