# Independent dense observed-data Gaussian calculations. These deliberately do
# not use the implementation's site-precision or missing-pattern calculations.
marginal_oracle_fixture <- function() {
  list(y = cbind(T = c(.8, NA, NA, 1.2, NA, NA, .7, NA, NA, NA, -.4),
                 W = c(1.1, 2.1, NA, NA, -1, NA, .2, 1.2, NA, NA, NA),
                 x = c(-.7, .2, .9, -.3, .5, -.8, .3, .6, -.1, .4, -.6)),
    cluster = c(rep("A", 4), rep("B", 2), rep("C", 2), rep("new", 2), "E"),
    mu = c(.3, -.2, .1),
    sigma = matrix(c(.99, 1.285, .7, 1.285, 1.9025, 1.05, .7, 1.05, 1), 3L),
    psi = matrix(c(1.71, 1.605, 1.1, 1.605, 1.7025, 1.05, 1.1, 1.05, 1), 3L))
}

marginal_dense_oracle <- function(y, cluster, mu, sigma, psi) {
  n <- nrow(y); d <- ncol(y)
  covariance <- kronecker(diag(n), sigma) +
    kronecker(outer(cluster, cluster, `==`) * 1, psi)
  design <- kronecker(matrix(1, n, 1L), diag(d))
  mean <- as.vector(design %*% mu)
  values <- as.vector(t(y))
  observed <- which(!is.na(values))
  target <- d * (which(is.na(y[, 1L])) - 1L) + 1L
  observed_covariance <- covariance[observed, observed, drop = FALSE]
  gain <- covariance[target, observed, drop = FALSE] %*% solve(observed_covariance)
  a <- design[observed, , drop = FALSE]
  mu_covariance <- solve(crossprod(a, solve(observed_covariance, a)))
  list(rows = which(is.na(y[, 1L])),
    mean = as.vector(mean[target] + gain %*% (values[observed] - mean[observed])),
    covariance = covariance[target, target, drop = FALSE] -
      gain %*% covariance[observed, target, drop = FALSE],
    mu_mean = as.vector(mu_covariance %*%
      crossprod(a, solve(observed_covariance, values[observed]))),
    mu_covariance = mu_covariance)
}

test_that("marginal predictions agree with independent dense Gaussian conditioning", {
  f <- marginal_oracle_fixture()
  oracle <- do.call(marginal_dense_oracle, f)
  set.seed(170801L)
  repetitions <- 10000L
  samples <- replicate(repetitions, do.call(catmech_marginal_cluster_predict, f))
  known <- !is.na(f$y[, 1L])
  expect_true(all(is.finite(samples)))
  expect_equal(unname(samples[known, , drop = FALSE]),
    matrix(rep(f$y[known, 1L], repetitions), nrow = sum(known)), tolerance = 0)
  missing_samples <- t(samples[oracle$rows, , drop = FALSE])
  mean_mcse <- sqrt(diag(oracle$covariance) / repetitions)
  covariance_mcse <- sqrt((outer(diag(oracle$covariance),
    diag(oracle$covariance)) + oracle$covariance^2) / (repetitions - 1L))
  expect_lt(max(abs(colMeans(missing_samples) - oracle$mean) / mean_mcse), 6)
  expect_lt(max(abs(stats::cov(missing_samples) - oracle$covariance) / covariance_mcse), 6)

  # Mixed recording patterns must share the same site draw. A site with no T
  # or W still has observed X: treating it as wholly unobserved is incorrect.
  mixed <- match(c(2L, 3L), oracle$rows)
  expect_gt(abs(oracle$covariance[mixed[1L], mixed[2L]]) /
    covariance_mcse[mixed[1L], mixed[2L]], 8)
  unaudited <- match(c(9L, 10L), oracle$rows)
  expect_lt(oracle$covariance[unaudited[1L], unaudited[2L]], f$psi[1L, 1L])
  expect_lt(max(diag(oracle$covariance)[unaudited]), f$sigma[1L, 1L] + f$psi[1L, 1L])
  expect_gt(max(abs(oracle$mean[unaudited] - f$mu[1L]) / mean_mcse[unaudited]), 8)
})

test_that("marginal prediction handles row order and several continuous predictors", {
  f <- marginal_oracle_fixture()
  # Independent linear construction gives a valid four-response covariance;
  # the fourth variable contributes information beyond the first predictor.
  add <- c(.2, -.3, .4)
  extend <- function(v) rbind(cbind(v, v %*% add),
    c(as.vector(crossprod(add, v)), as.numeric(crossprod(add, v %*% add)) + .7))
  f$sigma <- extend(f$sigma)
  f$psi <- extend(f$psi)
  f$mu <- c(f$mu, .4)
  f$y <- cbind(f$y, x2 = seq(-.9, .9, length.out = nrow(f$y)))
  permutation <- c(11L, 3L, 9L, 1L, 6L, 4L, 8L, 2L, 10L, 5L, 7L)
  original <- do.call(marginal_dense_oracle, f)
  f$y <- f$y[permutation, , drop = FALSE]
  f$cluster <- f$cluster[permutation]
  oracle <- do.call(marginal_dense_oracle, f)
  locations <- match(permutation[oracle$rows], original$rows)
  expect_equal(oracle$mean, original$mean[locations], tolerance = 1e-12)
  expect_equal(oracle$covariance, original$covariance[locations, locations], tolerance = 1e-12)
  set.seed(170802L)
  repetitions <- 4000L
  samples <- replicate(repetitions, do.call(catmech_marginal_cluster_predict, f))
  mean_mcse <- sqrt(diag(oracle$covariance) / repetitions)
  expect_lt(max(abs(rowMeans(samples[oracle$rows, , drop = FALSE]) -
    oracle$mean) / mean_mcse), 6)
})

test_that("fixed-covariance Gibbs mean draws agree with the observed likelihood", {
  f <- marginal_oracle_fixture()
  oracle <- do.call(marginal_dense_oracle, f)
  input <- list(y = f$y, subj = match(f$cluster, unique(f$cluster)))
  retained <- 6000L; burnin <- 1000L; chains <- 4L; d <- ncol(f$y)
  samples <- array(NA_real_, c(retained, chains, d))
  for (chain in seq_len(chains)) {
    fit <- catmech_marginal_gibbs(input, seed = 170810L + chain,
      iterations = retained + burnin, fixed = list(sigma = f$sigma, psi = f$psi))
    samples[, chain, ] <- t(fit$mu[, burnin + seq_len(retained), drop = FALSE])
    expect_true(all(fit$scale_proposed == 0 & fit$scale_accepted == 0))
    expect_equal(fit$sigma[, , c(1L, burnin + retained)],
      array(rep(f$sigma, 2L), c(d, d, 2L)), tolerance = 0)
    expect_equal(fit$psi[, , c(1L, burnin + retained)],
      array(rep(f$psi, 2L), c(d, d, 2L)), tolerance = 0)
  }
  centered <- samples
  for (j in seq_len(d)) {
    centered[, , j] <- samples[, , j] - oracle$mu_mean[j]
    expect_lt(abs(mean(centered[, , j])) / posterior::mcse_mean(samples[, , j]), 6)
    expect_lt(posterior::rhat(samples[, , j]), 1.02)
    expect_gt(posterior::ess_bulk(samples[, , j]), 1000)
    for (k in seq_len(j)) {
      product <- centered[, , j] * (samples[, , k] - oracle$mu_mean[k])
      expect_lt(abs(mean(product) - oracle$mu_covariance[j, k]) /
        posterior::mcse_mean(product), 6)
    }
  }
  # A fixed-covariance marginal intercept is a known linear contrast of mu.
  covariance <- f$sigma + f$psi
  slopes <- solve(covariance[-c(1L, 2L), -c(1L, 2L), drop = FALSE],
    covariance[-c(1L, 2L), 1L])
  contrast <- c(1, 0, -slopes)
  intercept <- samples[, , 1L] - slopes * samples[, , 3L]
  location <- sum(contrast * oracle$mu_mean)
  variance <- as.numeric(crossprod(contrast, oracle$mu_covariance %*% contrast))
  expect_lt(abs(mean(intercept) - location) / posterior::mcse_mean(intercept), 6)
  square <- (intercept - location)^2
  expect_lt(abs(mean(square) - variance) / posterior::mcse_mean(square), 6)
  expect_equal(unname(slopes), .9, tolerance = 1e-12)
  expect_equal(f$sigma[3L, 1L] / f$sigma[3L, 3L], .7, tolerance = 1e-12)
  expect_equal(f$psi[3L, 1L] / f$psi[3L, 3L], 1.1, tolerance = 1e-12)
})

# Huang and Wand (2013), doi:10.1214/13-BA815, equation (2):
# Q|a ~ IW(nu+d-1, 2*nu*diag(1/a)); a_k ~ IG(1/2, 1/A^2).
# At nu=2 each SD/A is half-t2 and each pairwise correlation is Uniform(-1,1).
# In the conditional-regression prior these properties apply to the S and
# Omega blocks, not to the assembled full covariance. These independent CDF
# checks avoid variance moments, which do not exist for the half-t2 prior.
test_that("Huang-Wand helper has its declared block marginal distributions", {
  n <- 12000L; d <- 3L; scale <- 4
  samples <- catmech_marginal_hw_draw(n, d, scale = scale, seed = 170830L)
  covariance <- samples$covariance
  expect_equal(dim(covariance), c(d, d, n))
  expect_true(all(is.finite(covariance)))
  expect_equal(covariance, aperm(covariance, c(2L, 1L, 3L)), tolerance = 1e-10)
  # Distributional rather than moment checks also remain stable in this prior's
  # heavy tails. All bounds are six exact binomial Monte Carlo standard errors.
  probabilities <- c(.25, .5, .75, .9)
  thresholds <- stats::qt((probabilities + 1) / 2, df = 2)
  for (j in seq_len(d)) {
    standard_deviations <- sqrt(covariance[j, j, ]) / sqrt(scale)
    expect_true(all(standard_deviations > 0))
    for (k in seq_along(probabilities)) {
      actual <- mean(standard_deviations <= thresholds[k])
      mcse <- sqrt(probabilities[k] * (1 - probabilities[k]) / n)
      expect_lt(abs(actual - probabilities[k]) / mcse, 6)
    }
    if (j > 1L) for (k in seq_len(j - 1L)) {
      rho <- covariance[j, k, ] / sqrt(covariance[j, j, ]) / sqrt(covariance[k, k, ])
      expect_true(all(abs(rho) < 1))
      for (probability in c(.25, .5, .75)) {
        actual <- mean(rho <= 2 * probability - 1)
        expect_lt(abs(actual - probability) /
          sqrt(probability * (1 - probability) / n), 6)
      }
    }
  }
})

# Independent implementation comparison for the unknown-covariance posterior.
# At dimension two, IW(d+2, .5 I) matches pan's documented IW(4, .5 I)
# priors. The mean is intercept-only and flat in both implementations. This
# does not claim that pan implements the hierarchical prior; it tests the
# common missing-data Gaussian likelihood and all covariance Gibbs updates.
test_that("unknown-covariance Gibbs posterior agrees with independent pan sampling", {
  set.seed(170850L)
  group <- rep(seq_len(24L), each = 6L)
  within <- rep(seq_len(6L), 24L)
  sigma <- matrix(c(1.1, .4, .4, .8), 2L)
  psi <- matrix(c(.6, .25, .25, .5), 2L)
  site <- matrix(stats::rnorm(48L), 24L, 2L) %*% chol(psi)
  y <- matrix(stats::rnorm(288L), 144L, 2L) %*% chol(sigma) + site[group, ]
  y <- sweep(y, 2L, c(.4, -.3), "+")
  y[within > 2L | group > 18L, 1L] <- NA_real_
  y[within %in% c(2L, 5L), 2L] <- NA_real_
  active <- rowSums(!is.na(y)) > 0
  y <- y[active, , drop = FALSE]; group <- group[active]
  retained <- 10000L; burnin <- 2000L; chains <- 4L
  keep <- burnin + seq_len(retained)
  prior <- list(a = 4, Binv = diag(.5, 2L), c = 4, Dinv = diag(.5, 2L))
  starts <- lapply(seq_len(chains), function(chain) {
    initial <- y; initial[is.na(initial)] <- 0
    list(beta = matrix(c(-.2, .2) * chain, 1L),
      sigma = diag(.5 + .3 * chain, 2L), psi = diag(1.8 - .3 * chain, 2L), y = initial)
  })
  reference <- catmech_joint_cluster_pan(
    list(y = y, subj = group, pred = matrix(1, nrow(y), 1L)),
    prior, seeds = 170860L + seq_len(chains), starts = starts,
    iter = burnin + retained)
  candidate <- lapply(seq_len(chains), function(chain) catmech_marginal_gibbs(
    list(y = y, subj = group), seed = 170880L + chain, iterations = burnin + retained,
    prior = list(family = "inverse_wishart", scale = 1)))
  traces <- function(fits, reference = FALSE) {
    output <- array(NA_real_, c(retained, chains, 8L))
    for (chain in seq_len(chains)) {
      fit <- fits[[chain]]
      means <- if (reference) fit$beta[1L, , keep] else fit$mu[, keep]
      output[, chain, ] <- cbind(t(means), fit$sigma[1L, 1L, keep],
        fit$sigma[1L, 2L, keep], fit$sigma[2L, 2L, keep],
        fit$psi[1L, 1L, keep], fit$psi[1L, 2L, keep], fit$psi[2L, 2L, keep])
    }
    output
  }
  a <- traces(candidate); b <- traces(reference, reference = TRUE)
  for (j in seq_len(dim(a)[3L])) {
    expect_lt(posterior::rhat(a[, , j]), 1.05)
    expect_lt(posterior::rhat(b[, , j]), 1.05)
    mcse <- sqrt(posterior::mcse_mean(a[, , j])^2 + posterior::mcse_mean(b[, , j])^2)
    expect_lt(abs(mean(a[, , j]) - mean(b[, , j])) / mcse, 6)
    # Center both second moments at one fixed constant. The difference retains
    # the same zero expectation without separately estimated centering terms.
    aa <- a[, , j]^2; bb <- b[, , j]^2
    mcse <- sqrt(posterior::mcse_mean(aa)^2 + posterior::mcse_mean(bb)^2)
    expect_lt(abs(mean(aa) - mean(bb)) / mcse, 6)
    probabilities <- c(.025, .975)
    qa <- stats::quantile(a[, , j], probabilities, names = FALSE)
    qb <- stats::quantile(b[, , j], probabilities, names = FALSE)
    mcse <- sqrt(posterior::mcse_quantile(a[, , j], probabilities, names = FALSE)^2 +
      posterior::mcse_quantile(b[, , j], probabilities, names = FALSE)^2)
    expect_lt(max(abs(qa - qb) / mcse), 6)
  }
})

# Construct parameter draws with known marginal coefficients. Their within-site
# and between-site coefficients deliberately differ. This tests a distribution's
# transformation and summary, not affine invariance of priors under refitting.
marginal_affine_oracle <- function() {
  set.seed(170910L)
  n <- 96L; groups <- rep(seq_len(12L), each = 8L)
  x <- cbind(5 + stats::rnorm(n), -2 + stats::rnorm(n))
  d <- data.frame(check.names = FALSE, "recorded value" = .8 + x[, 1L] + stats::rnorm(n),
    "true value" = 1.2 + .6 * x[, 1L] - .4 * x[, 2L] + stats::rnorm(n),
    T = x[, 1L], W = x[, 2L], "site label" = as.character(groups))
  d[["true value"]][rep(1:8, 12L) > 4L | groups > 10L] <- NA_real_
  d[["recorded value"]][rep(1:8, 12L) %in% c(3L, 7L) | groups == 12L] <- NA_real_
  size <- 4000L
  beta <- cbind(stats::rnorm(size, .6, .12), stats::rnorm(size, -.4, .10))
  alpha <- stats::rnorm(size, 1.2, .15)
  predictor_mu <- cbind(stats::rnorm(size, 5, .15), stats::rnorm(size, -2, .12))
  mu <- rbind(alpha + rowSums(beta * predictor_mu),
    .8 + .9 * predictor_mu[, 1L] + .2 * predictor_mu[, 2L], t(predictor_mu))
  sigma <- psi <- array(NA_real_, c(4L, 4L, size))
  s <- matrix(c(1.3, .45, .45, .9), 2L)
  residual <- matrix(c(.8, .2, .2, .6), 2L)
  site <- matrix(c(.4, .15, .15, .5), 2L)
  block <- function(x, e) {
    v <- matrix(0, 4L, 4L)
    v[1:2, 1:2] <- x; v[3:4, 3:4] <- e
    v
  }
  loading <- function(b, w) rbind(c(b, 1, 0), c(w, 0, 1), c(1, 0, 0, 0), c(0, 1, 0, 0))
  within <- block(s, residual); between <- block(.5 * s, site)
  for (i in seq_len(size)) {
    a <- loading(beta[i, ] + c(.2, -.125), c(.9, .2))
    b <- loading(beta[i, ] - c(.4, -.25), c(1.1, -.3))
    sigma[, , i] <- a %*% within %*% t(a)
    psi[, , i] <- b %*% between %*% t(b)
  }
  list(data = d, mu = mu, sigma = sigma, psi = psi,
    coefficients = cbind(alpha, beta), predictors = c("T", "W"))
}

marginal_affine_result <- function(f) {
  d <- f$data
  problem <- catmech_marginal_cluster_prepare(d, "recorded value", d[["true value"]],
    f$predictors, "site label")
  means <- sweep(sweep(f$mu, 1L, problem$center, "-"), 1L, problem$scale, "/")
  size <- ncol(means)
  sigma <- psi <- array(NA_real_, c(4L, 4L, size))
  scaling <- outer(problem$scale, problem$scale)
  for (i in seq_len(size)) {
    sigma[, , i] <- f$sigma[, , i] / scaling
    psi[, , i] <- f$psi[, , i] / scaling
  }
  fit <- list(mu = means, sigma = sigma, psi = psi,
    aux_sigma = matrix(1, 4L, size), aux_psi = matrix(1, 4L, size))
  trace <- catmech_marginal_cluster_trace(problem, fit, seq_len(size))
  traces <- array(trace, c(size / 4L, 4L, ncol(trace)),
    dimnames = list(NULL, NULL, colnames(trace)))
  result <- structure(list(original = d,
    manifest = list(engine = "reference_marginal_cluster", m = 2L),
    spec = list(value = "recorded value", cluster_model = "joint_marginal", predictors = f$predictors,
      joint_control = list(chains = 4L, burnin = 0L, iterations = size / 4L,
        max_iterations = size / 4L), joint_transform = problem$transform,
      recorded_fit = "joint_marginal_gaussian"),
    model_diagnostics = list(engine = "hw_marginal_cluster", convergence_pass = TRUE,
      traces = traces, iterations_used = as.integer(size / 4L), transform = problem$transform)),
    class = "catmech_result")
  list(result = result, problem = problem)
}

test_that("multiple-predictor marginal summaries recover known affine coefficients", {
  f <- marginal_affine_oracle()
  prepared <- marginal_affine_result(f)
  out <- catmech_marginal_summary(prepared$result)
  expect_identical(out$term, c("(Intercept)", "T", "W"))
  expect_equal(out$estimate, colMeans(f$coefficients), tolerance = 1e-12, ignore_attr = TRUE)
  expect_equal(out$posterior_sd, apply(f$coefficients, 2L, stats::sd),
    tolerance = 1e-12, ignore_attr = TRUE)
  expected <- t(apply(f$coefficients, 2L, stats::quantile, c(.025, .975), names = FALSE))
  expect_equal(unname(as.matrix(out[c("lower", "upper")])), unname(expected), tolerance = 1e-12)
  # The final site has only predictors observed; its rows belong in the
  # likelihood for the predictor distribution and in the prediction problem.
  expect_identical(prepared$problem$rows, seq_len(nrow(f$data)))
  rows <- f$data[["site label"]] == "12"
  expect_true(all(is.na(prepared$problem$input$y[rows, 1:2])))
  expect_true(all(is.finite(prepared$problem$input$y[rows, 3:4])))

  # Apply X' = a + B X and T' = c + h T. A non-diagonal B forces the
  # multivariable covariance solve; h<0 also exercises negative recoding.
  b <- matrix(c(2, -.3, .6, 1.4), 2L)
  a <- c(4, -7); h <- -2; shift <- 3
  q <- diag(4L); q[1:2, 1:2] <- diag(h, 2L); q[3:4, 3:4] <- b
  recoded <- f
  recoded$predictors <- c("first predictor", "second: predictor")
  recoded$data <- f$data
  recoded$data[["recorded value"]] <- shift + h * f$data[["recorded value"]]
  recoded$data[["true value"]] <- shift + h * f$data[["true value"]]
  recoded$data[c("T", "W")] <- NULL
  xx <- sweep(as.matrix(f$data[c("T", "W")]) %*% t(b), 2L, a, "+")
  recoded$data[recoded$predictors] <- as.data.frame(xx)
  recoded$mu <- sweep(q %*% f$mu, 1L, c(shift, shift, a), "+")
  for (i in seq_len(ncol(f$mu))) {
    recoded$sigma[, , i] <- q %*% f$sigma[, , i] %*% t(q)
    recoded$psi[, , i] <- q %*% f$psi[, , i] %*% t(q)
  }
  beta <- h * f$coefficients[, 2:3] %*% solve(b)
  recoded$coefficients <- cbind(shift + h * f$coefficients[, 1L] - as.vector(beta %*% a), beta)
  re <- marginal_affine_result(recoded)
  actual <- catmech_marginal_summary(re$result)
  expect_identical(actual$term, c("(Intercept)", recoded$predictors))
  expect_equal(actual$estimate, colMeans(recoded$coefficients), tolerance = 1e-11, ignore_attr = TRUE)
  expect_equal(actual$posterior_sd, apply(recoded$coefficients, 2L, stats::sd),
    tolerance = 1e-11, ignore_attr = TRUE)
  expected <- t(apply(recoded$coefficients, 2L, stats::quantile, c(.025, .975), names = FALSE))
  expect_equal(unname(as.matrix(actual[c("lower", "upper")])), unname(expected), tolerance = 1e-11)

  bad <- f$data; bad$T <- factor(bad$T > median(bad$T))
  expect_error(catmech_reference(bad, "recorded value", "true value", c("T", "W"),
    cluster_model = "joint_marginal", use_recorded = TRUE, cluster = "site label"), "numeric")
  bad <- f$data; names(bad)[names(bad) == "T"] <- "(Intercept)"
  expect_error(catmech_reference(bad, "recorded value", "true value", c("(Intercept)", "W"),
    cluster_model = "joint_marginal", use_recorded = TRUE, cluster = "site label"), "Intercept")
})

test_that("marginal summaries with the audited value as a predictor match exact covariance solves", {
  skip_on_cran() # slow; GitHub Actions runs it with NOT_CRAN=true
  f <- marginal_affine_oracle()
  prepared <- marginal_affine_result(f)
  size <- ncol(f$mu)
  oracle_draws <- function(r_idx, x_idx) {
    out <- matrix(NA_real_, size, 1L + length(x_idx))
    for (i in seq_len(size)) {
      v <- f$sigma[, , i] + f$psi[, , i]
      m <- f$mu[, i]
      if (!length(x_idx)) {
        out[i, 1L] <- m[r_idx]
      } else {
        b <- as.vector(solve(v[x_idx, x_idx, drop = FALSE], v[x_idx, r_idx]))
        out[i, ] <- c(m[r_idx] - sum(b * m[x_idx]), b)
      }
    }
    out
  }
  full_expected <- oracle_draws(3L, c(1L, 4L))
  by_formula <- catmech_marginal_summary(prepared$result, formula = T ~ `recorded value` + W)
  by_args <- catmech_marginal_summary(prepared$result, response = "T",
    covariates = c("recorded value", "W"))
  by_default_covs <- catmech_marginal_summary(prepared$result, response = "T")
  expect_identical(by_formula, by_args)
  expect_identical(by_formula, by_default_covs)
  expect_identical(by_formula$term, c("(Intercept)", "recorded value", "W"))
  expect_equal(by_formula$estimate, colMeans(full_expected), tolerance = 1e-11, ignore_attr = TRUE)
  expect_equal(by_formula$posterior_sd, apply(full_expected, 2L, stats::sd),
    tolerance = 1e-11, ignore_attr = TRUE)
  q_full <- t(apply(full_expected, 2L, stats::quantile, c(.025, .975), names = FALSE))
  expect_equal(unname(as.matrix(by_formula[c("lower", "upper")])), unname(q_full), tolerance = 1e-11)

  uni_expected <- oracle_draws(3L, 1L)
  uni_summary <- catmech_marginal_summary(prepared$result, formula = T ~ `recorded value`)
  expect_identical(uni_summary$term, c("(Intercept)", "recorded value"))
  expect_equal(uni_summary$estimate, colMeans(uni_expected), tolerance = 1e-11, ignore_attr = TRUE)
  expect_equal(uni_summary$posterior_sd, apply(uni_expected, 2L, stats::sd),
    tolerance = 1e-11, ignore_attr = TRUE)

  int_expected <- oracle_draws(3L, integer())
  int_summary <- catmech_marginal_summary(prepared$result, formula = T ~ 1)
  expect_identical(int_summary,
    catmech_marginal_summary(prepared$result, response = "T", covariates = character()))
  expect_identical(int_summary$term, "(Intercept)")
  expect_equal(int_summary$estimate, mean(int_expected), tolerance = 1e-12)
  expect_equal(int_summary$posterior_sd, stats::sd(int_expected), tolerance = 1e-12)
})

# Evaluate full probability densities independently of the native optimized
# scale-move ratio. Constants cancel because the model and dimension are fixed.
marginal_oracle_log_normal <- function(residual, covariance) {
  residual <- as.matrix(residual)
  -.5 * (nrow(residual) * ncol(residual) * log(2 * pi) +
    nrow(residual) * as.numeric(determinant(covariance, logarithm = TRUE)$modulus) +
    sum((residual %*% solve(covariance)) * residual))
}
marginal_oracle_log_covariance_prior <- function(covariance, family, scale) {
  d <- nrow(covariance)
  logdet <- as.numeric(determinant(covariance, logarithm = TRUE)$modulus)
  precision <- solve(covariance)
  if (family == "huang_wand")
    -(2 + 2 * d) / 2 * logdet - (2 + d) / 2 * sum(log(2 * diag(precision) + 1 / scale))
  else -(d + 2 + d + 1) / 2 * logdet - .25 * scale * sum(diag(precision))
}

test_that("interweaving ratios match full densities and the complete Jacobian", {
  group <- c(1L, 2L, 1L, 3L, 3L, 3L, 2L, 3L, 1L)
  y <- matrix(c(.2, -.4, .9, 1.1, -.7, .3, .8, .6, -.2,
    -.1, .3, .5, -.8, .2, 1.2, .4, -.4, .7), 9L, 2L)
  mu <- c(.3, -.2)
  u <- matrix(c(.2, -.3, .5, -.1, .4, .2), 3L, 2L)
  sigma <- matrix(c(.8, -.25, -.25, .6), 2L)
  psi <- matrix(c(.6, .2, .2, .9), 2L)
  extend <- function(v) {
    loading <- c(.2, -.3)
    rbind(cbind(v, v %*% loading),
      c(as.vector(crossprod(loading, v)), as.numeric(crossprod(loading, v %*% loading)) + .7))
  }
  sigma <- extend(sigma); psi <- extend(psi)
  y <- cbind(y, c(.3, -.5, .1, .7, .2, -.8, .6, -.2, .4))
  u <- cbind(u, c(.1, -.2, .4)); mu <- c(mu, .1)
  density <- function(u, psi, family, scale) {
    residual <- sweep(y - u[group, ], 2L, mu, "-")
    marginal_oracle_log_normal(residual, sigma) +
      marginal_oracle_log_normal(u, psi) +
      marginal_oracle_log_covariance_prior(psi, family, scale)
  }
  for (j in seq_along(mu)) for (eta in c(-.8, .55)) {
    transform <- diag(length(mu)); transform[j, j] <- exp(eta)
    proposed_u <- u %*% transform
    proposed_psi <- transform %*% psi %*% transform
    jacobian <- (nrow(u) + ncol(u) + 1) * eta
    for (family in c("huang_wand", "inverse_wishart")) for (scale in c(.25, 1, 4)) {
      expected <- density(proposed_u, proposed_psi, family, scale) -
        density(u, psi, family, scale) + jacobian
      actual <- catmech_marginal_interweave_ratio(y, group, mu, u, sigma, psi,
        j, eta, prior = list(family = family, scale = scale))
      expect_equal(actual, expected, tolerance = 1e-11)
    }
    # Independently differentiate the full deterministic state transformation:
    # all U entries and the nonredundant symmetric Psi entries participate.
    indices <- which(lower.tri(psi, diag = TRUE), arr.ind = TRUE)
    state <- c(u, psi[indices])
    move <- function(state) {
      uu <- matrix(state[seq_along(u)], nrow(u), ncol(u))
      pp <- matrix(0, length(mu), length(mu))
      pp[indices] <- state[length(u) + seq_len(nrow(indices))]
      pp <- pp + t(pp) - diag(diag(pp))
      pp <- transform %*% pp %*% transform
      c(uu %*% transform, pp[indices])
    }
    step <- 1e-5
    derivative <- vapply(seq_along(state), function(k) {
      upper <- lower <- state; upper[k] <- upper[k] + step; lower[k] <- lower[k] - step
      (move(upper) - move(lower)) / (2 * step)
    }, numeric(length(state)))
    expect_equal(as.numeric(determinant(derivative, logarithm = TRUE)$modulus),
      jacobian, tolerance = 1e-9)
  }
})

test_that("unfactorized Huang-Wand posterior is unchanged by scale interweaving", {
  set.seed(170940L)
  group <- rep(seq_len(24L), each = 10L)
  sigma <- matrix(c(.6, .15, .15, .5), 2L)
  psi <- matrix(c(.8, .3, .3, .7), 2L)
  site <- matrix(stats::rnorm(48L), 24L, 2L) %*% chol(psi)
  y <- sweep(matrix(stats::rnorm(480L), 240L, 2L) %*% chol(sigma) + site[group, ],
    2L, c(.4, -.3), "+")
  y[seq(7L, 240L, 7L), 1L] <- NA_real_
  y[seq(10L, 240L, 10L), 2L] <- NA_real_
  retained <- 8000L; burnin <- 2000L; chains <- 4L
  keep <- burnin + seq_len(retained)
  run <- function(interweave, seed) {
    output <- array(NA_real_, c(retained, chains, 12L))
    for (chain in seq_len(chains)) {
      f <- catmech_marginal_gibbs(list(y = y, subj = group), seed = seed + chain,
        iterations = burnin + retained,
        prior = list(family = "huang_wand", scale = 1, interweave = interweave))
      expect_equal(as.vector(f$scale_proposed), rep(if (interweave) 2 * (burnin + retained) else 0, 2L))
      expect_true(all(f$scale_accepted >= 0 & f$scale_accepted <= f$scale_proposed))
      output[, chain, ] <- cbind(t(f$mu[, keep]), f$sigma[1L, 1L, keep],
        f$sigma[1L, 2L, keep], f$sigma[2L, 2L, keep], f$psi[1L, 1L, keep],
        f$psi[1L, 2L, keep], f$psi[2L, 2L, keep],
        t(log(f$aux_sigma[, keep])), t(log(f$aux_psi[, keep])))
    }
    output
  }
  a <- run(TRUE, 170950L); b <- run(FALSE, 170970L)
  for (j in seq_len(dim(a)[3L])) {
    expect_lt(posterior::rhat(a[, , j]), 1.05)
    expect_lt(posterior::rhat(b[, , j]), 1.05)
    mcse <- sqrt(posterior::mcse_mean(a[, , j])^2 + posterior::mcse_mean(b[, , j])^2)
    expect_lt(abs(mean(a[, , j]) - mean(b[, , j])) / mcse, 6)
    probabilities <- c(.025, .975)
    qa <- stats::quantile(a[, , j], probabilities, names = FALSE)
    qb <- stats::quantile(b[, , j], probabilities, names = FALSE)
    mcse <- sqrt(posterior::mcse_quantile(a[, , j], probabilities, names = FALSE)^2 +
      posterior::mcse_quantile(b[, , j], probabilities, names = FALSE)^2)
    expect_lt(max(abs(qa - qb) / mcse), 6)
  }
})
