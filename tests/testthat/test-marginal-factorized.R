# Independent conditional-distribution oracles for the proper regression
# factorization of a covariance matrix. Values are centered row/site vectors,
# not cohorts used in repeated-sampling validation.
factorized_oracle_block <- function(k) {
  set.seed(171010L + k)
  x <- matrix(stats::rnorm(60L * k), 60L, k)
  if (k == 2L) x[, 2L] <- .55 * x[, 1L] + .8 * x[, 2L]
  beta <- matrix(c(.8, -.4, -.3, .6)[seq_len(2L * k)], k, 2L)
  y <- x %*% beta + matrix(stats::rnorm(120L), 60L, 2L) %*%
    chol(matrix(c(.7, -.2, -.2, .9), 2L))
  list(values = cbind(y, x), auxiliary = seq(.8, 1.4, length.out = k + 2L),
    regression_scale = 2.5)
}

factorized_inverse_wishart_variance <- function(scale, df) {
  p <- nrow(scale)
  # Marginal variance of each covariance entry under IW(df, scale).
  ((df - p + 1) * scale^2 + (df - p - 1) * outer(diag(scale), diag(scale))) /
    ((df - p) * (df - p - 1)^2 * (df - p - 3))
}

test_that("factorized covariance blocks match exact IW and matrix-normal conditionals", {
  for (k in 1:2) {
    f <- factorized_oracle_block(k)
    y <- f$values[, 1:2, drop = FALSE]
    x <- f$values[, 2L + seq_len(k), drop = FALSE]
    nobs <- nrow(x); samples <- 10000L
    precision <- crossprod(x) + diag(1 / f$regression_scale^2, k)
    location <- solve(precision, crossprod(x, y))
    s_scale <- crossprod(x) + diag(4 / f$auxiliary[2L + seq_len(k)], k)
    omega_scale <- crossprod(y) - crossprod(y, x) %*% location +
      diag(4 / f$auxiliary[1:2], 2L)
    expected_s <- s_scale / nobs
    expected_omega <- omega_scale / nobs
    out <- catmech_marginal_factorized_draw(f$values, f$auxiliary,
      regression_scale = f$regression_scale, n = samples, seed = 171020L + k)
    expect_equal(dim(out$covariance), c(k + 2L, k + 2L, samples))
    s <- array(NA_real_, c(k, k, samples))
    omega <- array(NA_real_, c(2L, 2L, samples))
    coefficients <- innovations <- matrix(NA_real_, samples, 2L * k)
    for (i in seq_len(samples)) {
      q <- out$covariance[, , i]
      s[, , i] <- q[2L + seq_len(k), 2L + seq_len(k), drop = FALSE]
      b <- q[1:2, 2L + seq_len(k), drop = FALSE] %*%
        solve(matrix(s[, , i], k, k))
      omega[, , i] <- q[1:2, 1:2] - b %*% matrix(s[, , i], k, k) %*% t(b)
      coefficients[i, ] <- as.vector(t(b))
      innovations[i, ] <- as.vector(chol(precision) %*% (t(b) - location) %*%
        solve(chol(omega[, , i])))
    }
    mcse_s <- sqrt(factorized_inverse_wishart_variance(s_scale, nobs + k + 1) / samples)
    mcse_omega <- sqrt(factorized_inverse_wishart_variance(omega_scale, nobs + 3) / samples)
    expect_lt(max(abs(apply(s, c(1L, 2L), mean) - expected_s) / mcse_s), 6)
    expect_lt(max(abs(apply(omega, c(1L, 2L), mean) - expected_omega) / mcse_omega), 6)
    variance <- kronecker(expected_omega, solve(precision))
    expect_lt(max(abs(colMeans(coefficients) - as.vector(location)) /
      sqrt(diag(variance) / samples)), 6)
    expect_lt(max(abs(colMeans(innovations))) * sqrt(samples), 6)
    covariance_mcse <- sqrt((outer(rep(1, 2L * k), rep(1, 2L * k)) + diag(2L * k)) /
      (samples - 1L))
    expect_lt(max(abs(stats::cov(innovations) - diag(2L * k)) / covariance_mcse), 6)
    # S and Omega are independent given the completed vectors and fixed
    # auxiliaries. A cross-block dependence in this sampler would be spurious.
    cross <- (s[1L, 1L, ] - expected_s[1L, 1L]) *
      (omega[1L, 1L, ] - expected_omega[1L, 1L])
    expect_lt(abs(mean(cross)) / (stats::sd(cross) / sqrt(samples)), 6)
  }
})

factorized_oracle_log_hw <- function(q, scale) {
  d <- nrow(q)
  -(2 + 2 * d) / 2 * as.numeric(determinant(q, logarithm = TRUE)$modulus) -
    (2 + d) / 2 * sum(log(2 * diag(solve(q)) + 1 / scale))
}
factorized_oracle_log_prior <- function(q, scale, regression_scale) {
  d <- nrow(q); k <- d - 2L
  if (!k) return(factorized_oracle_log_hw(q, scale))
  xx <- 2L + seq_len(k)
  s <- q[xx, xx, drop = FALSE]
  b <- q[1:2, xx, drop = FALSE] %*% solve(s)
  omega <- q[1:2, 1:2] - b %*% s %*% t(b)
  factorized_oracle_log_hw(s, scale) + factorized_oracle_log_hw(omega, scale) -
    k / 2 * as.numeric(determinant(omega, logarithm = TRUE)$modulus) -
    2 * k * log(regression_scale) -
    sum((solve(omega, b)) * b) / (2 * regression_scale^2) -
    2 * as.numeric(determinant(s, logarithm = TRUE)$modulus)
}
factorized_oracle_log_normal <- function(residual, q) {
  -.5 * (nrow(residual) * ncol(residual) * log(2 * pi) +
    nrow(residual) * as.numeric(determinant(q, logarithm = TRUE)$modulus) +
    sum((residual %*% solve(q)) * residual))
}

test_that("factorized scale moves use the density on full covariance coordinates", {
  set.seed(171030L)
  group <- c(1L, 2L, 1L, 3L, 3L, 3L, 2L, 3L, 1L)
  for (d in 3:4) {
    z <- matrix(stats::rnorm(d * d), d, d)
    sigma <- crossprod(z) / d + diag(.5, d)
    z <- matrix(stats::rnorm(d * d), d, d)
    psi <- crossprod(z) / d + diag(.5, d)
    y <- matrix(stats::rnorm(length(group) * d), length(group), d)
    u <- matrix(stats::rnorm(3L * d, sd = .3), 3L, d)
    mu <- seq(-.3, .4, length.out = d)
    for (scale in c(.25, 4)) for (regression_scale in c(5, 20)) {
      density <- function(u, psi) {
        residual <- sweep(y - u[group, ], 2L, mu, "-")
        factorized_oracle_log_normal(residual, sigma) +
          factorized_oracle_log_normal(u, psi) +
          factorized_oracle_log_prior(psi, scale, regression_scale)
      }
      for (j in seq_len(d)) for (eta in c(-.65, .8)) {
        transform <- diag(d); transform[j, j] <- exp(eta)
        expected <- density(u %*% transform, transform %*% psi %*% transform) -
          density(u, psi) + (nrow(u) + d + 1) * eta
        actual <- catmech_marginal_interweave_ratio(y, group, mu, u, sigma, psi,
          j, eta, prior = list(family = "conditional_regression", scale = scale,
            regression_scale = regression_scale))
        expect_equal(actual, expected, tolerance = 1e-10)
      }
    }
  }
})

test_that("the covariance factorization Jacobian equals determinant(S) squared", {
  for (k in 1:2) {
    s <- if (k == 1L) matrix(.7, 1L) else matrix(c(.7, .2, .2, 1.2), 2L)
    b <- matrix(seq(-.3, .8, length.out = 2L * k), 2L, k)
    omega <- matrix(c(.6, -.15, -.15, .8), 2L)
    is <- which(lower.tri(s, diag = TRUE), arr.ind = TRUE)
    io <- which(lower.tri(omega, diag = TRUE), arr.ind = TRUE)
    iq <- which(lower.tri(matrix(0, k + 2L, k + 2L), diag = TRUE), arr.ind = TRUE)
    state <- c(s[is], b, omega[io])
    rebuild <- function(state) {
      ss <- matrix(0, k, k); ss[is] <- state[seq_len(nrow(is))]
      ss <- ss + t(ss) - diag(diag(ss), k)
      bb <- matrix(state[nrow(is) + seq_along(b)], 2L, k)
      oo <- matrix(0, 2L, 2L)
      oo[io] <- state[nrow(is) + length(b) + seq_len(nrow(io))]
      oo <- oo + t(oo) - diag(diag(oo), 2L)
      q <- rbind(cbind(oo + bb %*% ss %*% t(bb), bb %*% ss), cbind(ss %*% t(bb), ss))
      q[iq]
    }
    step <- 1e-5
    derivative <- vapply(seq_along(state), function(j) {
      upper <- lower <- state; upper[j] <- upper[j] + step; lower[j] <- lower[j] - step
      (rebuild(upper) - rebuild(lower)) / (2 * step)
    }, numeric(length(state)))
    expect_equal(as.numeric(determinant(derivative, logarithm = TRUE)$modulus),
      2 * as.numeric(determinant(s, logarithm = TRUE)$modulus), tolerance = 1e-9)
  }
})

test_that("factorized-prior posterior is unchanged by interweaving", {
  set.seed(171040L)
  group <- rep(seq_len(24L), each = 8L)
  sigma <- matrix(c(.99, 1.285, .7, 1.285, 1.9025, 1.05, .7, 1.05, 1), 3L)
  psi <- matrix(c(1.71, 1.605, 1.1, 1.605, 1.7025, 1.05, 1.1, 1.05, 1), 3L)
  site <- matrix(stats::rnorm(72L), 24L, 3L) %*% chol(psi)
  y <- sweep(matrix(stats::rnorm(length(group) * 3L), length(group), 3L) %*%
    chol(sigma) + site[group, ], 2L, c(.4, -.3, .1), "+")
  y[seq(7L, nrow(y), 7L), 1L] <- NA_real_
  y[seq(5L, nrow(y), 5L), 2L] <- NA_real_
  retained <- 8000L; burnin <- 2000L; chains <- 4L
  keep <- burnin + seq_len(retained)
  indices <- which(lower.tri(sigma, diag = TRUE), arr.ind = TRUE)
  run <- function(interweave, seed) {
    output <- array(NA_real_, c(retained, chains, 23L))
    for (chain in seq_len(chains)) {
      f <- catmech_marginal_gibbs(list(y = y, subj = group), seed = seed + chain,
        iterations = burnin + retained, prior = list(family = "conditional_regression",
          scale = 1, regression_scale = 10, interweave = interweave))
      covariance_traces <- function(name) vapply(seq_len(nrow(indices)), function(j)
        f[[name]][indices[j, 1L], indices[j, 2L], keep], numeric(retained))
      slope <- (f$sigma[3L, 1L, keep] + f$psi[3L, 1L, keep]) /
        (f$sigma[3L, 3L, keep] + f$psi[3L, 3L, keep])
      intercept <- f$mu[1L, keep] - slope * f$mu[3L, keep]
      output[, chain, ] <- cbind(t(f$mu[, keep]), covariance_traces("sigma"),
        covariance_traces("psi"), t(log(f$aux_sigma[, keep])),
        t(log(f$aux_psi[, keep])), intercept, slope)
    }
    output
  }
  a <- run(TRUE, 171050L); b <- run(FALSE, 171070L)
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
