# Candidate replacements for catmech_reference_cluster_fill(), with its arguments, and the rescaled-bootstrap comparator.
# Locked in validation/CLUSTERED_FIX_PLAN.md. Only base R and stats, no added dependency.

fix_check_audit <- function(x, fit_rows, ids) {
  if (qr(x[fit_rows, , drop = FALSE], tol = 1e-7)$rank != ncol(x) || length(fit_rows) <= ncol(x))
    stop("The trusted regression design must have full column rank; remove aliased or unsupported predictors.", call. = FALSE)
  if (length(unique(ids[fit_rows])) < 5L)
    stop("Fewer than five audited clusters: cluster uncertainty is refused.", call. = FALSE)
}

# Rows of D are scaled delete-one-cluster coefficient differences, so crossprod(D) is the stratified jackknife covariance.
fix_jackknife_differences <- function(x, y, fit_rows, ids, strata) {
  by_cluster <- split(fit_rows, as.character(ids[fit_rows]))
  groups <- if (is.null(strata)) list(names(by_cluster)) else
    split(names(by_cluster), as.character(strata[vapply(by_cluster, `[`, integer(1), 1L)]))
  p <- ncol(x)
  full <- stats::lm.fit(x[fit_rows, , drop = FALSE], y[fit_rows])$coefficients
  differences <- list()
  lone <- 0L
  for (members in groups) {
    kept <- list()
    if (length(members) > 1L) for (cluster in members) {
      rows <- setdiff(fit_rows, by_cluster[[cluster]])
      refit <- stats::lm.fit(x[rows, , drop = FALSE], y[rows])
      if (refit$rank == p && length(rows) > p) kept[[length(kept) + 1L]] <- refit$coefficients - full
    }
    if (!length(kept)) {
      lone <- lone + 1L
      next
    }
    differences <- c(differences, lapply(kept, function(d) d * sqrt((length(members) - 1) / length(kept))))
  }
  if (!length(differences))
    stop("No stratum has an audited cluster that can be left out: cluster uncertainty cannot be estimated.", call. = FALSE)
  do.call(rbind, differences) * sqrt(length(groups) / (length(groups) - lone))
}

jackknife_fill <- function(x, y, fit_rows, fill_rows, ids, strata, m) {
  fix_check_audit(x, fit_rows, ids)
  n <- length(fit_rows); p <- ncol(x)
  fit <- stats::lm.fit(x[fit_rows, , drop = FALSE], y[fit_rows])
  s2 <- sum(fit$residuals^2) / (n - p)
  if (!is.finite(s2) || s2 <= 0) stop("The audited rows leave no residual variance.", call. = FALSE)
  D <- fix_jackknife_differences(x, y, fit_rows, ids, strata)
  out <- matrix(NA_real_, length(fill_rows), m)
  for (j in seq_len(m)) {
    beta <- fit$coefficients + as.vector(crossprod(D, stats::rnorm(nrow(D))))
    sigma2 <- s2 * (n - p) / stats::rchisq(1L, n - p)
    out[, j] <- as.vector(x[fill_rows, , drop = FALSE] %*% beta) + stats::rnorm(length(fill_rows), 0, sqrt(sigma2))
  }
  out
}

# Moment (Henderson III) variance components, drawn in each imputation, then generalized least squares given them.
random_intercept_fill <- function(x, y, fit_rows, fill_rows, ids, strata, m) {
  fix_check_audit(x, fit_rows, ids)
  X <- x[fit_rows, , drop = FALSE]; Y <- y[fit_rows]; n <- nrow(X); p <- ncol(X)
  cl <- factor(as.character(ids[fit_rows])); ni <- tabulate(cl)
  ols <- stats::lm.fit(X, Y)
  sse_ols <- sum(ols$residuals^2)
  if (!is.finite(sse_ols) || sse_ols <= 0) stop("The audited rows leave no residual variance.", call. = FALSE)
  D <- stats::model.matrix(~ cl - 1)
  within <- stats::lm.fit(cbind(X, D), Y)
  df_w <- n - within$rank; df_b <- within$rank - p
  sse_w <- sum(within$residuals^2); ss_b <- max(0, sse_ols - sse_w)
  separable <- df_w > 0L && df_b > 0L && sse_w > 0
  c_value <- if (separable) n - sum(crossprod(qr.Q(ols$qr), D)^2) else NA_real_
  xbar <- rowsum(X, cl, reorder = TRUE) / ni; ybar <- as.vector(rowsum(Y, cl, reorder = TRUE)) / ni
  fill_clusters <- as.character(ids[fill_rows]); fill_units <- unique(fill_clusters)
  out <- matrix(NA_real_, length(fill_rows), m)
  for (j in seq_len(m)) {
    draw <- NULL
    for (attempt in seq_len(100L)) {
      if (separable) {
        sigma2 <- sse_w / stats::rchisq(1L, df_w)
        tau2 <- max(0, df_b * (ss_b / stats::rchisq(1L, df_b) - sigma2) / c_value)
      } else {
        sigma2 <- sse_ols / stats::rchisq(1L, n - p)
        tau2 <- 0
      }
      lambda <- 1 - sqrt(sigma2 / (sigma2 + ni * tau2))
      gls <- stats::lm.fit(X - (lambda * xbar)[as.integer(cl), , drop = FALSE], Y - (lambda * ybar)[as.integer(cl)])
      if (gls$rank != p) next
      beta <- gls$coefficients + as.vector(backsolve(qr.R(gls$qr), stats::rnorm(p))) * sqrt(sigma2)
      effect <- stats::rnorm(length(fill_units), 0, sqrt(tau2)); names(effect) <- fill_units
      seen <- fill_units[fill_units %in% levels(cl)]
      if (length(seen) && tau2 > 0) {
        k <- ni[match(seen, levels(cl))]
        residual_mean <- as.vector(rowsum(Y - as.vector(X %*% beta), cl, reorder = TRUE)) / ni
        effect[seen] <- k * tau2 / (sigma2 + k * tau2) * residual_mean[match(seen, levels(cl))] +
          stats::rnorm(length(seen), 0, sqrt(tau2 * sigma2 / (sigma2 + k * tau2)))
      }
      draw <- as.vector(x[fill_rows, , drop = FALSE] %*% beta) + unname(effect[fill_clusters]) + stats::rnorm(length(fill_rows), 0, sqrt(sigma2))
      break
    }
    if (is.null(draw)) stop("Audited clusters lost design support in imputation; audit more clusters.", call. = FALSE)
    out[, j] <- draw
  }
  out
}

# Comparator from the review of the first test: n_h - 1 clusters drawn in each stratum, weighted n_h / (n_h - 1).
rao_wu_fill <- function(x, y, fit_rows, fill_rows, ids, strata, m) {
  fix_check_audit(x, fit_rows, ids)
  by_cluster <- split(fit_rows, as.character(ids[fit_rows]))
  groups <- if (is.null(strata)) list(names(by_cluster)) else
    split(names(by_cluster), as.character(strata[vapply(by_cluster, `[`, integer(1), 1L)]))
  out <- matrix(NA_real_, length(fill_rows), m)
  for (j in seq_len(m)) {
    draw <- NULL
    for (attempt in seq_len(200L)) {
      rows <- integer(); w <- numeric()
      for (u in groups) {
        nh <- length(u)
        pick <- if (nh == 1L) u else u[sample.int(nh, nh - 1L, replace = TRUE)]
        got <- by_cluster[pick]
        rows <- c(rows, unlist(got, use.names = FALSE)); w <- c(w, rep(if (nh == 1L) 1 else nh / (nh - 1), sum(lengths(got))))
      }
      fit <- stats::lm.wfit(x[rows, , drop = FALSE], y[rows], w)
      if (fit$rank != ncol(x) || sum(w) <= ncol(x)) next
      variance <- sum(w * fit$residuals^2) / (sum(w) - ncol(x))
      if (!is.finite(variance) || variance <= 0) next
      draw <- as.vector(x[fill_rows, , drop = FALSE] %*% fit$coefficients) + stats::rnorm(length(fill_rows), 0, sqrt(variance))
      break
    }
    if (is.null(draw)) stop("Audited clusters lost design support in resampling; supply cluster_strata or audit more clusters.", call. = FALSE)
    out[, j] <- draw
  }
  out
}
