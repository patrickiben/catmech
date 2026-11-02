# Two-phase design-based estimator for a regression coefficient on audited values.
#
# Phase 1 = all clusters in the cohort (a sample from the superpopulation of clusters).
# Phase 2 = audited clusters, stratified SRS, design weights N_h / n_h from the counts.
#
# Nuisance terms (strata, and the covariate when the audited value is a predictor) are
# partialled out on phase-1 data instead of being re-estimated from one or two audited
# clusters per stratum. The unknown normal-equation totals are then estimated with survey
# weights calibrated (linear/GREG, constant within cluster) to phase-1 totals of auxiliaries
# built from the recorded value:
#   audited outcome:    beta = sum(xt * truth) / sum(xt^2); aux = xt*lab, xt^2, xt
#   audited predictor:  beta = sum(ts * yt) / sum(ts^2);    aux = ths*yt, ths^2
# where ths is a working imputation of the truth given lab, covariate, outcome and a shrunken
# stratum offset, residualised on phase-1 nuisance terms. If truth equals the working value the
# estimate reduces exactly to the phase-1 fit.
#
# Variance: V = V1 + V2. For a regression estimator the superpopulation variance is
# lambda^2 S_a^2 / G + S_r^2 / n, so the fpc on phase 2 cancels against the residual part of
# phase 1 and no fpc is used.
#   V1 = delete-cluster (leverage-corrected) variance of the auxiliary part on all clusters.
#        For count-like values it is floored by a smooth mean-variance working model.
#   V2 = delete-one-audited-cluster jackknife with re-imputation and re-calibration, all audited
#        clusters in one group (collapsed strata), so singleton strata need no special rule.
# df = Satterthwaite of G - H - 1 and n - 1.
suppressPackageStartupMessages(library(survey))

.sv_csum <- function(v, cl, G) {
  r <- numeric(G); s <- rowsum(v, cl, reorder = TRUE); r[as.integer(rownames(s))] <- s; r
}
.sv_clamp <- function(v, c0) pmax(pmin(v, c0), -c0)

# Closed form of survey::calibrate(calfun = "linear", aggregate.stage = 1) on cluster-level
# auxiliaries, used for jackknife replicates (speed) and as fallback.
.sv_cal <- function(w, A, tot) {
  keep <- w > 0; Ak <- A[keep, , drop = FALSE]; wk <- w[keep]
  M <- crossprod(Ak * sqrt(wk)); rhs <- tot - colSums(Ak * wk)
  sc <- sqrt(pmax(diag(M), 1e-300)); M2 <- M / outer(sc, sc)
  e <- eigen(M2, symmetric = TRUE); ok <- e$values > 1e-8 * max(e$values)
  inv <- e$vectors[, ok, drop = FALSE] %*% (t(e$vectors[, ok, drop = FALSE]) / e$values[ok])
  lam <- (inv %*% (rhs / sc)) / sc
  g <- numeric(length(w)); g[keep] <- wk * (1 + drop(Ak %*% lam)); g
}

# Point-estimate weights through the survey package; cluster-level result of length G.
.sv_cal_survey <- function(w, A, tot, h) {
  id <- which(w > 0); k <- ncol(A)
  dd <- data.frame(id = id, h = factor(h[id]), w = w[id], A[id, , drop = FALSE])
  nm <- paste0("a", seq_len(k)); names(dd)[4:(3 + k)] <- nm
  names(tot) <- nm
  g <- tryCatch({
    des <- if (nlevels(dd$h) > 1L) svydesign(ids = ~id, strata = ~h, weights = ~w, data = dd)
           else svydesign(ids = ~id, weights = ~w, data = dd)
    cal <- calibrate(des, as.formula(paste("~ 0 +", paste(nm, collapse = "+"))), population = tot,
                     calfun = "linear", aggregate.stage = 1)
    1 / cal$prob
  }, error = function(e) NULL, warning = function(w) NULL)
  if (is.null(g) || anyNA(g) || !all(is.finite(g))) return(.sv_cal(w, A, tot))
  out <- numeric(length(w)); out[id] <- g; out
}

.sv_jk1 <- function(idx, w, est) {
  n <- length(idx)
  t <- vapply(idx, function(c0) {
    w1 <- w; w1[idx] <- w1[idx] * n / (n - 1); w1[c0] <- 0
    tryCatch(est(w1, c0), error = function(e) NA_real_)
  }, 0)
  t <- t[is.finite(t)]
  if (length(t) < max(3L, n - 2L)) stop("jackknife failed")
  (n - 1) / n * sum((t - mean(t))^2)
}

.sv_jk_p1 <- function(pc, lev, h) {
  V <- 0
  for (k in unique(h)) {
    i <- h == k; m <- sum(i)
    if (m > 1) { x <- pc[i] / pmax(1 - lev[i], .1); V <- V + (m - 1) / m * sum((x - mean(x))^2) }
  }
  V
}

# Count-like recorded values: Var(cluster residual total) = c1 m^2 + c2 m over clusters.
.sv_mv <- function(xt, fit, el, b1, T0, cl, G, lev) {
  nc <- .sv_csum(rep(1, length(xt)), cl, G); m <- pmax(.sv_csum(fit, cl, G), 1e-3)
  ec <- .sv_csum(el, cl, G) / pmax(1 - lev, .1); sc <- .sv_csum(xt, cl, G)
  cf <- qr.coef(qr(cbind(m^2, m) / (m + 1)), ec^2 / (m + 1))
  if (anyNA(cf)) return(0)
  cf <- pmax(cf, 0)
  sum((sc / nc)^2 * (cf[1] * m^2 + cf[2] * m)) * b1^2 / T0^2
}

fit_method <- function(p, seed) {   # deterministic: seed is not needed
  d <- p$input; a <- as.logical(p$audited); f <- p$formula
  cl <- as.integer(factor(d$cluster_id)); G <- max(cl); N <- nrow(d)
  st <- if ("stratum" %in% names(d)) droplevels(as.factor(d$stratum)) else factor(rep("all", N))
  first <- !duplicated(cl); h <- integer(G); h[cl[first]] <- as.integer(st[first])
  aud <- logical(G); aud[cl[a]] <- TRUE
  Nh <- tabulate(h); nh <- tabulate(h[aud], nbins = length(Nh)); H <- length(Nh); n <- sum(aud)
  if (any(nh == 0L)) stop("a stratum has no audited cluster")
  if (n < 4L) stop("too few audited clusters")
  if (anyNA(d$gold_value[a]) || anyNA(d$lab_value) || anyNA(d$covariate)) stop("missing values")
  w0 <- ifelse(aud, Nh[h] / nh[h], 0); lab <- d$lab_value
  gold <- d$gold_value; gold[!a] <- 0
  exact <- all(gold[a] == lab[a]); idx <- which(aud)
  countlike <- all(lab >= 0) && all(lab == round(lab))
  isout <- identical(all.vars(f[[2]]), "lab_value")
  X <- model.matrix(f[-2], d); j <- match(p$term, colnames(X))
  if (is.na(j)) stop("term not in model")
  Z <- X[, -j, drop = FALSE]; QZ <- qr(Z)
  Z <- Z[, QZ$pivot[seq_len(QZ$rank)], drop = FALSE]; QZ <- qr(Z)
  df1 <- max(2, G - H - 1); V2 <- 0

  if (isout) {
    x <- X[, j]; xt <- qr.resid(QZ, x); T0 <- sum(xt^2)
    if (T0 < 1e-10) stop("covariate has no variation after adjustment")
    vc <- .sv_csum(xt^2, cl, G); lev <- vc / T0 + 1 / Nh[h]
    el <- qr.resid(qr(cbind(Z, xt)), lab)
    if (exact) { th <- sum(xt * lab) / T0; b1 <- 1 }
    else {
      A <- cbind(.sv_csum(xt * lab, cl, G), vc, .sv_csum(xt, cl, G)); tot <- colSums(A)
      u <- .sv_csum(xt * gold, cl, G)
      est <- function(w, c0 = 0) sum(.sv_cal(w, A, tot) * u) / T0
      th <- sum(.sv_cal_survey(w0, A, tot, h) * u) / T0
      k <- w0 > 0
      B <- qr.coef(qr(A[k, , drop = FALSE] * sqrt(w0[k])), u[k] * sqrt(w0[k])); B[is.na(B)] <- 0
      b1 <- B[1]
      V2 <- .sv_jk1(idx, w0, est)
    }
    pc <- .sv_csum(xt * b1 * el, cl, G) / T0
    V1 <- .sv_jk_p1(pc, lev, h)
    if (countlike) V1 <- max(V1, .sv_mv(xt, lab - el, el, b1, T0, cl, G, lev))
  } else {
    resp <- all.vars(f[[2]]); yt <- qr.resid(QZ, d[[resp]])
    Pr <- cbind(1, lab, d$covariate, d[[resp]])
    if (sum(qr.resid(QZ, lab)^2) < 1e-10) stop("recorded value has no variation after adjustment")
    build <- function(keep, wr) {
      b <- lm.wfit(Pr[keep, , drop = FALSE], gold[keep], wr[keep])$coefficients; b[is.na(b)] <- 0
      t0 <- drop(Pr %*% b)
      if (H > 1L && !exact) {        # empirical-Bayes stratum offset of the imputation residual
        rc <- tapply(gold[keep] - t0[keep], cl[keep], mean); hc <- h[as.integer(names(rc))]
        nk <- tabulate(hc, nbins = H); hm <- numeric(H); s <- tapply(rc, hc, mean)
        hm[as.integer(names(s))] <- s; nn <- sum(nk); nz <- nk > 0; Hk <- sum(nz)
        if (nn > Hk && Hk > 1L) {
          sw <- sum((rc - hm[hc])^2) / (nn - Hk); mg <- sum(nk * hm) / nn
          bs <- sum(nk * (hm - mg)^2) / (Hk - 1); nt <- (nn - sum(nk^2) / nn) / (Hk - 1)
          tau <- max(0, (bs - sw) / nt)
          if (tau > 0) t0 <- t0 + (ifelse(nz, tau / (tau + sw / pmax(nk, 1)), 0) * (hm - mg))[h[cl]]
        }
      }
      fitZ <- t0 - qr.resid(QZ, t0); ths <- t0 - fitZ; ts <- gold - fitZ
      if (!exact && !countlike) {    # no extrapolation of cluster shifts beyond the audited range
        cm <- .sv_csum(ths, cl, G) / .sv_csum(rep(1, N), cl, G)
        ths <- ths - (cm - .sv_clamp(cm, max(abs(cm[unique(cl[keep])]))))[cl]
      }
      list(A = cbind(.sv_csum(ths * yt, cl, G), .sv_csum(ths^2, cl, G)),
           u1 = .sv_csum(ts * yt, cl, G), u2 = .sv_csum(ts^2, cl, G), ths = ths)
    }
    b0 <- build(a, w0[cl])
    if (exact) { th <- sum(b0$ths * yt) / sum(b0$ths^2); D <- sum(b0$ths^2) }
    else {
      g <- .sv_cal_survey(w0, b0$A, colSums(b0$A), h); D <- sum(g * b0$u2); th <- sum(g * b0$u1) / D
      est <- function(w, c0 = 0) {
        b <- build(a & cl != c0, w[cl]); g1 <- .sv_cal(w, b$A, colSums(b$A)); D1 <- sum(g1 * b$u2)
        if (D1 > 0) sum(g1 * b$u1) / D1 else NA_real_
      }
      V2 <- .sv_jk1(idx, w0, est)
    }
    if (!is.finite(D) || D <= 0) stop("degenerate audited predictor")
    vc <- .sv_csum(b0$ths^2, cl, G); vz <- .sv_csum(rowSums(qr.Q(QZ)^2), cl, G)
    pc <- .sv_csum(b0$ths * (yt - th * b0$ths), cl, G) / D
    V1 <- .sv_jk_p1(pc, vc / sum(vc) + vz, h)
  }
  V <- V1 + V2
  if (!is.finite(th) || !is.finite(V) || V <= 0) stop("no usable variance")
  df <- if (V2 > 0 && V1 > 0) V^2 / (V1^2 / df1 + V2^2 / (n - 1L)) else if (V2 > 0) n - 1L else df1
  list(estimate = unname(th), se = sqrt(V), df = df)
}
