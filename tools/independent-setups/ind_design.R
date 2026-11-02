# ind_design: design-based correction of the recorded analysis, written directly.
#
# Estimator (a coefficient-level difference/GREG estimator with one tuned scalar):
#   estimate = T + kappa * (F - R)
#     T = design-weighted OLS on audited rows with true values    (weights N_h / n_h)
#     R = the same regression on audited rows with recorded values
#     F = OLS on all rows with recorded values
#   kappa = regression through the origin of T's cluster scores on R's, audited clusters only,
#   clipped to [0, 2]. If recorded equals truth on every audited row, kappa = 1 and estimate = F.
#
# Variance (superpopulation, no finite-population correction, since the target is the coefficient
# in the population the cohort represents): sum over all clusters of (g_c + kappa * h_c)^2, where
#   g_c = T score - kappa * R score (audited clusters), h_c = F score (all clusters).
#   Scores use CR2 leverage adjustment, pseudo-inverse where a stratum effect absorbs a cluster,
#   or CR3 when the response is a count, because variance then rises with leverage.
#   g_c is also divided by (1 - leverage in the kappa regression), a leave-one-out correction for
#   tuning; leverage is scaled by n_h / (n_h - 1) since stratum effects use up contrasts.
# Degrees of freedom: Bell-McCaffrey Satterthwaite for the combined statistic, each part scaled to
#   its observed scores, capped by an empirical value 3 v^2 / sum(s^4) that falls with one dominant
#   cluster.

.idz_fit <- function(X, y, w, term) {
  sw <- sqrt(w); qx <- qr(X * sw, tol = 1e-9); r <- qx$rank; piv <- qx$pivot[seq_len(r)]
  Q <- qr.Q(qx)[, seq_len(r), drop = FALSE]; R <- qr.R(qx)[seq_len(r), seq_len(r), drop = FALSE]
  yw <- y * sw; qty <- crossprod(Q, yw); coefs <- backsolve(R, qty)
  pos <- match(term, colnames(X)[piv])
  if (is.na(pos)) return(list(ok = FALSE))
  u <- numeric(r); u[pos] <- 1
  list(ok = TRUE, Q = Q, e = as.numeric(yw - Q %*% qty),
       L = as.numeric(Q %*% backsolve(R, u, transpose = TRUE)), est = coefs[pos])
}

.idz_adjust <- function(fit, rows, power) lapply(rows, function(i) {
  ee <- eigen(diag(length(i)) - tcrossprod(fit$Q[i, , drop = FALSE]), symmetric = TRUE)
  d <- ifelse(ee$values > 1e-8, pmax(ee$values, 1e-8)^(-power), 0)
  ee$vectors %*% (d * t(ee$vectors))
})

.idz_scores <- function(fit, rows, adj) vapply(seq_along(rows), function(j) {
  i <- rows[[j]]; sum(fit$L[i] * (adj[[j]] %*% fit$e[i]))
}, 0)

# q_c' M q_d under a working iid model, for the Bell-McCaffrey degrees of freedom.
.idz_cross <- function(fit, rows, adj) {
  G <- length(rows); B <- matrix(0, ncol(fit$Q), G); dg <- numeric(G)
  for (j in seq_len(G)) {
    i <- rows[[j]]; qc <- as.numeric(adj[[j]] %*% fit$L[i])
    dg[j] <- sum(qc^2); B[, j] <- crossprod(fit$Q[i, , drop = FALSE], qc)
  }
  out <- -crossprod(B); diag(out) <- diag(out) + dg; out
}

.idz_core <- function(p) {
  d <- p$input; f <- p$formula; term <- p$term; n <- nrow(d)
  if (anyNA(d$lab_value)) stop("recorded value missing")
  au <- as.logical(p$audited) & !is.na(d$gold_value)
  if (sum(au) < 3L) stop("too few audited rows")
  cid <- as.character(d$cluster_id); cl_all <- unique(cid); aud_cl <- unique(cid[au])
  y_name <- all.vars(f)[1L]; outcome_case <- identical(y_name, "lab_value")
  has_s <- "stratum" %in% names(d)
  strat <- if (has_s) as.character(d$stratum) else rep("all", n)
  str_cl <- strat[match(cl_all, cid)]
  Nh <- table(str_cl); nh <- table(factor(str_cl[cl_all %in% aud_cl], levels = names(Nh)))
  if (any(nh == 0L)) stop("a stratum has no audited cluster")
  w_row <- as.numeric((Nh / nh)[strat])

  rhs <- stats::delete.response(stats::terms(f))
  mk <- function(v, idx) {
    dd <- d[idx, , drop = FALSE]; dd$lab_value <- v[idx]
    stats::model.matrix(rhs, dd, xlev = if (has_s) list(stratum = levels(factor(d$stratum))))
  }
  X_all <- mk(d$lab_value, seq_len(n)); X_R <- X_all[au, , drop = FALSE]
  if (outcome_case) {
    X_T <- X_R; y_T <- d$gold_value[au]; y_R <- d$lab_value[au]; y_F <- d$lab_value
  } else {
    X_T <- mk(d$gold_value, which(au)); y_F <- d[[y_name]]; y_T <- y_R <- y_F[au]
    if (anyNA(y_F)) stop("outcome missing")
  }
  cid_a <- cid[au]
  rows_all <- unname(split(seq_len(n), factor(cid, levels = cl_all)))
  rows_a <- unname(split(seq_along(cid_a), factor(cid_a, levels = aud_cl)))
  map_a <- match(aud_cl, cl_all); G <- length(cl_all); Ga <- length(aud_cl)

  yy <- if (outcome_case) d$lab_value else y_F
  power <- if (all(yy >= 0) && all(abs(yy - round(yy)) < 1e-8)) 1 else 0.5

  fT <- .idz_fit(X_T, y_T, w_row[au], term)
  if (!fT$ok) stop("term not estimable in the audited clusters")
  fR <- .idz_fit(X_R, y_R, w_row[au], term)
  fF <- .idz_fit(X_all, y_F, rep(1, n), term)
  aT <- .idz_adjust(fT, rows_a, power)
  pT <- .idz_scores(fT, rows_a, aT)
  kap <- 0; ctrl <- 0; pR <- numeric(Ga); h <- numeric(G)
  if (fR$ok && fF$ok) {
    aR <- if (outcome_case) aT else .idz_adjust(fR, rows_a, power)
    aF <- .idz_adjust(fF, rows_all, power)
    pR <- .idz_scores(fR, rows_a, aR); pF <- .idz_scores(fF, rows_all, aF)
    ctrl <- fF$est - fR$est
    kap <- if (sum(pR^2) > 0) max(0, min(2, sum(pT * pR) / sum(pR^2))) else 1
    h <- kap * pF
  }
  est <- fT$est + kap * ctrl
  g <- pT - kap * pR
  if (kap > 0 && sum(pR^2) > 0) {
    ng <- as.numeric(nh[str_cl[map_a]])
    lev <- pmin(ifelse(ng > 1, ng / (ng - 1), 1) * pR^2 / sum(pR^2), 0.6)
    g <- g / (1 - lev)
  }
  sa <- numeric(G); sa[map_a] <- g
  s <- sa + h
  v <- sum(s^2)
  if (!is.finite(est) || !is.finite(v) || v <= 0) stop("no estimable variance")

  Ca <- .idz_cross(fT, rows_a, aT)
  C <- matrix(0, G, G); C[map_a, map_a] <- sum(g^2) / max(sum(diag(Ca)), 1e-300) * Ca
  if (kap > 0) {
    Cf <- .idz_cross(fF, rows_all, aF)
    C <- C + sum(h^2) / max(sum(diag(Cf)), 1e-300) * Cf
  }
  df <- max(1, min(sum(diag(C))^2 / sum(C^2), 3 * v^2 / sum(s^4)))
  list(estimate = est, se = sqrt(v), df = df, kappa = kap)
}

fit_method <- function(p, seed) {
  set.seed(seed)  # the method has no randomness; kept for the interface
  r <- .idz_core(p)
  list(estimate = r$estimate, se = r$se, df = r$df)
}
