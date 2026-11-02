# Multiple imputation with a two-level (Mundlak) imputation model, run through mice.
#
# Imputation model for the error d = gold - recorded on audited clusters:
#   d_ij = Zb_g gamma_b + Zw_ij gamma_w + u_g + e_ij
# Zb = cluster means of (recorded, covariate, outcome if allowed), plus strata with
# empirical-Bayes partial pooling. Zw = within-cluster deviations of the same variables.
# Separate within and between slopes keep second moments right at both levels, which a
# common-slope random intercept does not when errors shift whole clusters.
# The single incomplete variable means a Gibbs chain is not needed: the posterior factors
# into a within and a between regression, so mice runs one pass and each draw is exact
# (balanced clusters; approximate when cluster sizes differ).
#
# Analysis: OLS on completed data with CR2 variance and Bell-McCaffrey df.
# Pooling: Rubin's total variance. Reference df by Satterthwaite, combining the analysis df
# with the df of the between-cluster regression on audited clusters, because the between-
# imputation variance is estimated from those few clusters and not from m.

.ind2l_prepare <- function(y, ry, W, Zr, strat, cl, deflate) {
  all_v <- cbind(W = W, Zr)
  cidx <- match(cl, unique(cl)); G <- max(cidx); ns <- tabulate(cidx, G)
  means <- rowsum(all_v, cidx, reorder = TRUE) / ns
  within <- all_v - means[cidx, , drop = FALSE]
  d <- y - W
  gaud <- as.vector(rowsum(as.numeric(ry), cidx, reorder = TRUE)) > 0
  if (any(ry != gaud[cidx])) stop("audit must cover whole clusters")
  dmean <- numeric(G); dmean[gaud] <- as.vector(rowsum(d[ry], cidx[ry], reorder = TRUE)) / ns[gaud]
  n_a <- sum(ry); G_a <- sum(gaud)
  scale_y <- 1 + sum(y[ry]^2)
  # Within-cluster regression
  dw <- d[ry] - dmean[cidx[ry]]
  Zw <- within[ry, , drop = FALSE]
  cs <- sqrt(colSums(Zw^2)); keep_w <- which(cs > 1e-10 * (1 + sqrt(colSums(all_v^2))))
  rw <- 0L; gw <- numeric(0); Rw <- NULL; kw <- integer(0); RSSw <- sum(dw^2)
  if (length(keep_w)) {
    qw <- qr(Zw[, keep_w, drop = FALSE], tol = 1e-9); rw <- qw$rank
    kw <- keep_w[qw$pivot[seq_len(rw)]]
    qw <- qr(Zw[, kw, drop = FALSE]); Rw <- qr.R(qw)
    gw <- backsolve(Rw, qr.qty(qw, dw)[seq_len(rw)])
    RSSw <- sum((dw - drop(Zw[, kw, drop = FALSE] %*% gw))^2)
  }
  df_w <- n_a - G_a - rw
  # Between-cluster regression, row-weighted to match the OLS analysis
  kbar <- 1 / mean(1 / ns[gaud])
  wt <- sqrt(ns[gaud] / kbar)
  Zb <- cbind(`(I)` = 1, means)
  Zba <- Zb[gaud, , drop = FALSE] * wt; db <- dmean[gaud] * wt
  qb0 <- qr(Zba, tol = 1e-9); k0 <- qb0$pivot[seq_len(qb0$rank)]
  Zb <- Zb[, k0, drop = FALSE]; Zba <- Zba[, k0, drop = FALSE]
  p0 <- ncol(Zba); df_b <- G_a - p0; edf <- 0
  if (df_b < 3) stop("too few audited clusters for the imputation model")
  if (!is.null(strat)) {
    sc <- droplevels(strat[match(seq_len(G), cidx)])
    fs <- droplevels(sc[gaud]); Sh <- nlevels(fs)
    if (Sh >= 2L && G_a - Sh >= 2L) {
      r0 <- qr.resid(qr(Zba), db)
      s0 <- sum(r0^2) / df_b
      nh <- tabulate(fs); mh <- as.vector(rowsum(r0, fs, reorder = TRUE)) / nh
      ssb <- sum(nh * mh^2); n0 <- (G_a - sum(nh^2) / G_a) / (Sh - 1)
      ssw <- (sum(r0^2) - ssb) / (G_a - Sh)
      om2 <- max(0, (ssb / (Sh - 1) - ssw) / n0)
      if (om2 > 1e-12 * (s0 + 1e-300)) {
        S <- stats::model.matrix(~ sc - 1)
        S <- S[, colSums(S[gaud, , drop = FALSE]) > 0, drop = FALSE]
        Sa <- S[gaud, , drop = FALSE] * wt
        lam <- max(s0, 1e-12) / om2
        full <- cbind(Zba, Sa)
        edf <- sum(diag(solve(crossprod(full) + diag(c(rep(0, p0), rep(lam, ncol(S)))), crossprod(full)))) - p0
        Zba <- rbind(full, cbind(matrix(0, ncol(S), p0), sqrt(lam) * diag(ncol(S))))
        db <- c(db, rep(0, ncol(S)))
        Zb <- cbind(Zb, S)
        df_b <- G_a - p0 - edf
      }
    }
  }
  qb <- qr(Zba); Rb <- qr.R(qb); rb <- ncol(Zba)
  if (qb$rank < rb) stop("rank-deficient between-cluster imputation model")
  gb <- backsolve(Rb, qr.qty(qb, db)[seq_len(rb)])
  RSSb <- sum((db - drop(Zba %*% gb))^2)
  exact <- RSSw + RSSb < 1e-14 * scale_y
  # Deflation for an imputed predictor: draws of gamma and of the variances inflate the
  # variance of the imputed value by about 2 p / G plus the inverse chi-square excess, and
  # that excess attenuates a slope on the imputed value. No such effect when it is the outcome.
  f_b <- f_w <- 1
  if (deflate) {
    f_b <- max(0.1, (df_b - 2) / df_b - 2 * (p0 - 1 + edf) / G_a)
    if (df_w > 2) f_w <- max(0.1, (df_w - 2) / df_w - 2 * rw / (n_a - G_a))
  }
  mis <- which(!ry)
  list(mis = mis, W = W, cidx = cidx, G = G, exact = exact, df_b = df_b, df_w = df_w, rw = rw, rb = rb,
       RSSw = RSSw, RSSb = RSSb, gw = gw, gb = gb, Rw = Rw, Rb = Rb, kbar = kbar, f_b = f_b, f_w = f_w,
       Zbm = Zb[cidx[mis], , drop = FALSE], Zwm = within[mis, kw, drop = FALSE])
}

.ind2l_draw <- function(st) {
  s2 <- if (st$df_w > 0) st$RSSw / stats::rchisq(1, st$df_w) else 0
  gws <- if (st$rw > 0) st$gw + sqrt(s2) * backsolve(st$Rw, stats::rnorm(st$rw)) else numeric(0)
  # sigma_b^2 = tau^2 + sigma^2 / k, so its posterior is truncated below at sigma^2 / k
  lower <- s2 / st$kbar
  cmax <- if (lower > 0) st$RSSb / lower else Inf
  pmax_chi <- stats::pchisq(cmax, st$df_b)
  sb2 <- if (pmax_chi <= 0) lower else st$RSSb / stats::qchisq(stats::runif(1, 0, pmax_chi), st$df_b)
  sb2 <- max(sb2, lower)
  gbs <- st$gb + sqrt(sb2) * backsolve(st$Rb, stats::rnorm(st$rb))
  t2 <- max(0, st$f_b * sb2 - st$f_w * lower)
  u <- stats::rnorm(st$G, 0, sqrt(t2))
  fit <- drop(st$Zbm %*% gbs)
  if (st$rw > 0) fit <- fit + drop(st$Zwm %*% gws)
  st$W[st$mis] + fit + u[st$cidx[st$mis]] + stats::rnorm(length(st$mis), 0, sqrt(st$f_w * s2))
}

# Elementary mice method. The state is fixed: the only incomplete variable depends on complete data.
mice.impute.ind2l <- function(y, ry, x, wy = NULL, ind2l_state, ...) .ind2l_draw(ind2l_state)

# CR2 pieces for one design matrix.
.cr2_setup <- function(X, term, gidx, want_df = TRUE, vtype = "cr3") {
  qx <- qr(X); rk <- qx$rank; X <- X[, qx$pivot[seq_len(rk)], drop = FALSE]
  j <- match(term, colnames(X)); if (is.na(j)) stop("term aliased in the analysis model")
  XtXi <- chol2inv(qr.R(qr(X))); a <- XtXi[j, ]
  L <- X %*% XtXi; xa <- drop(X %*% a); w <- numeric(nrow(X)); G <- max(gidx)
  rows <- split(seq_len(nrow(X)), factor(gidx, levels = seq_len(G)))
  for (idx in rows) {
    n_g <- length(idx)
    H <- L[idx, , drop = FALSE] %*% t(X[idx, , drop = FALSE])
    if (n_g == 1L) {
      w[idx] <- xa[idx] / (if (vtype == "cr3") max(1 - H, 1e-8) else sqrt(max(1 - H, 1e-8)))
    } else {
      ev <- eigen(diag(n_g) - H, symmetric = TRUE)
      pw <- if (vtype == "cr3") 1 else 0.5
      inv <- ifelse(ev$values > 1e-8, ev$values^(-pw), 0)
      w[idx] <- ev$vectors %*% (inv * crossprod(ev$vectors, xa[idx]))
    }
  }
  df <- NA_real_
  if (want_df) {
    Wf <- matrix(0, nrow(X), G); Wf[cbind(seq_len(nrow(X)), gidx)] <- w
    Mw <- Wf - X %*% (XtXi %*% crossprod(X, Wf))
    S <- crossprod(Mw); df <- sum(diag(S))^2 / sum(S^2)
  }
  list(X = X, XtXi = XtXi, a = a, w = w, df = df, j = j)
}

.cr2_apply <- function(cs, Y, gidx) {
  Y <- as.matrix(Y)
  b <- cs$XtXi %*% crossprod(cs$X, Y)
  E <- Y - cs$X %*% b
  u <- colSums(rowsum(cs$w * E, gidx, reorder = FALSE)^2)
  list(est = b[cs$j, ], u = u)
}

.ind2l_analyse <- function(d, formula, term, Ymat, gidx, vtype = "cr3") {
  rhs <- "lab_value" %in% all.vars(formula[[3]])
  m <- ncol(Ymat)
  if (!rhs) {
    d$lab_value <- Ymat[, 1]
    mf <- stats::model.frame(formula, d, na.action = stats::na.fail); X <- stats::model.matrix(formula, mf)
    cs <- .cr2_setup(X, term, gidx, TRUE, vtype); r <- .cr2_apply(cs, Ymat, gidx)
    return(list(est = r$est, u = r$u, df_com = cs$df))
  }
  est <- u <- numeric(m); df_com <- NA_real_
  d$lab_value <- Ymat[, 1]
  mf <- stats::model.frame(formula, d, na.action = stats::na.fail); X <- stats::model.matrix(formula, mf); yy <- stats::model.response(mf)
  lv <- which(colnames(X) == "lab_value")
  for (i in seq_len(m)) {
    X[, lv] <- Ymat[, i]
    cs <- .cr2_setup(X, term, gidx, i == 1L, vtype)
    if (i == 1L) df_com <- cs$df
    r <- .cr2_apply(cs, yy, gidx); est[i] <- r$est; u[i] <- r$u
  }
  list(est = est, u = u, df_com = df_com)
}

.ind2l_pool <- function(a, df_b) {
  m <- length(a$est); qbar <- mean(a$est); ubar <- mean(a$u); B <- stats::var(a$est)
  if (m < 2L || B <= 0) return(list(estimate = qbar, se = sqrt(ubar), df = a$df_com))
  Bp <- (1 + 1 / m) * B
  nu <- (ubar + Bp)^2 / (ubar^2 / a$df_com + Bp^2 * (1 / df_b + 1 / (m - 1)))
  list(estimate = qbar, se = sqrt(ubar + Bp), df = nu)
}

fit_method <- function(p, seed, m = 30L, vtype = "cr3") {
  set.seed(seed)
  d <- p$input
  if (!requireNamespace("mice", quietly = TRUE)) stop("mice not available")
  ry <- as.logical(p$audited)
  if (sum(ry) < 1L) stop("no audited rows")
  gidx <- match(d$cluster_id, unique(d$cluster_id))
  rhs <- "lab_value" %in% all.vars(p$formula[[3]])
  preds <- intersect(setdiff(p$predictors, c("stratum", "lab_value", "gold_value", "cluster_id")), names(d))
  Zr <- as.matrix(d[, preds, drop = FALSE])
  strat <- if ("stratum" %in% names(d)) d$stratum else NULL
  if (all(ry)) {
    a <- .ind2l_analyse(d, p$formula, p$term, matrix(d$gold_value, ncol = 1L), gidx, vtype)
    return(list(estimate = a$est, se = sqrt(a$u), df = a$df_com))
  }
  st <- .ind2l_prepare(d$gold_value, ry, d$lab_value, Zr, strat, d$cluster_id, deflate = rhs)
  if (st$exact) {
    # Recorded equals truth on every audited row: nothing to model, analyse the recorded values.
    Y <- matrix(ifelse(ry, d$gold_value, d$lab_value), ncol = 1L)
    a <- .ind2l_analyse(d, p$formula, p$term, Y, gidx, vtype)
    return(list(estimate = a$est, se = sqrt(a$u), df = a$df_com))
  }
  # mice looks elementary methods up from its own namespace, which reaches only the global search path.
  assign("mice.impute.ind2l", mice.impute.ind2l, envir = globalenv())
  md <- data.frame(gold_value = d$gold_value, lab_value = d$lab_value)
  pm <- matrix(0L, 2L, 2L, dimnames = list(names(md), names(md))); pm["gold_value", "lab_value"] <- 1L
  imp <- suppressWarnings(mice::mice(md, m = m, maxit = 1L, method = c("ind2l", ""), predictorMatrix = pm,
    printFlag = FALSE, remove.collinear = FALSE, remove.constant = FALSE,
    eps = 0, seed = NA, blots = list(gold_value = list(ind2l_state = st))))
  Y <- matrix(d$gold_value, nrow(d), m); Y[!ry, ] <- as.matrix(imp$imp$gold_value)
  a <- .ind2l_analyse(d, p$formula, p$term, Y, gidx, vtype)
  .ind2l_pool(a, st$df_b)
}
