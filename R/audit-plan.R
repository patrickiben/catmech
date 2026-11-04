#' Plan how many records to audit from a random pilot
#'
#' Answers "how many records should we audit?" for linear multiple imputation with
#' the recorded value, [catmech_reference()] with `use_recorded = TRUE`. It predicts
#' how much wider the interval for one regression coefficient will be at each audit
#' fraction than with a full audit, using the classical two-phase form
#' `width(f) = width_full * sqrt(1 + k * (1 - f) / f)`. It predicts precision only. It says
#' nothing about bias or coverage, and it does not replace [catmech_recording_audit()].
#'
#' @param data A data frame holding every record in the cohort, audited or not.
#' @param value Name of the numeric recorded value, with `NA` where nothing was recorded.
#' @param reference Name of the numeric reference value, observed on the pilot audit
#'   only. The pilot must be a random sample of the records.
#' @param predictors Complete auxiliary columns available on every record, as in
#'   [catmech_reference()]. Include the analysis outcome.
#' @param formula The planned linear analysis, as passed to [catmech_pool_lm()]. It must
#'   be linear in `value`, so transform the column first rather than writing `log(value)`.
#' @param term Name of the coefficient the audit is meant to protect.
#' @param fractions Audit fractions to report, each above zero and at most one.
#' @param target_ratio Width relative to a full audit that counts as good enough, above
#'   one. Default 1.25, meaning within 25\% of the full-audit width.
#' @param conf_level Level of the analysis interval and of the plan's intervals.
#' @param bootstrap Bootstrap replicates of the pilot, at least 100.
#' @param m Imputations for the pilot's own fits, at least two.
#' @param seed Nonnegative integer seed. The caller's RNG state is restored.
#' @return A list with:
#'   \describe{
#'     \item{`status`}{`"no_flags"` or `"review_flags"`.}
#'     \item{`flags`}{Short names of every warning raised.}
#'     \item{`k`}{One row each for `recorded` (keep the recorded value, linear MI) and
#'       `discard` (the [catmech_reference()] default): `estimate`, `lower`, `upper`,
#'       `classical`, the empirical version described below, and `rubin_check`, the value
#'       implied by the pilot fit's own between-imputation variance.}
#'     \item{`table`}{One row per fraction: `fraction`, `records`, `additional_records`
#'       beyond the pilot, `small_audit`, which marks a size the function would refuse as a
#'       pilot, `relative_width` with `relative_lower` and `relative_upper` from the
#'       bootstrap interval for `k`, `predictive_lower` and `predictive_upper`
#'       for the width a single audit of that size may realize, `width` with `width_lower`
#'       and `width_upper` anchored on the pilot fit, and the same for discarding recorded
#'       values: `discard_relative_width`, `discard_relative_lower`,
#'       `discard_relative_upper` and `discard_width`.}
#'     \item{`answer`}{Records needed to reach `target_ratio`: `fraction` and `records` from
#'       the estimate, `records_lower` and `records_upper` from the interval for `k`, and
#'       `smallest_grid_fraction`, the smallest of `fractions` predicted to reach it, or `NA`.}
#'     \item{`pilot`}{The pilot's size, disagreements, fits and diagnostics.}
#'     \item{`assumptions`}{The assumptions the plan rests on, as text.}
#'   }
#' @details
#' Classical theory for an efficient two-phase estimator gives
#' `Var(f) / Var(1) = 1 + k * (1 - f) / f` with `k = Var(psi - E[psi | aux]) / Var(psi)`,
#' where `psi` is the influence function of the coefficient using true values and `aux` is
#' what is known on every record. This is not a new result.
#'
#' The planner evaluates that ratio under the working models, because the interval from
#' Rubin's rules comes from them. The numerator is `Var(psi | aux)` under the normal linear
#' imputation model fitted on the pilot's reference rows, averaged over all records. The
#' denominator is `Var(psi)` under the homoskedastic analysis model. The result has a closed
#' form and needs no simulation. The `classical` column instead uses empirical variances on
#' the pilot rows. On unequal-variance data the two differ: the first describes the reported
#' interval and the second the true sampling variance. Use `estimate` to plan interval width.
#'
#' `width` is anchored on the pilot's own fit: its within-imputation variance estimates the
#' full-audit variance, and `relative_width` scales the result. The pilot's total width is not
#' the anchor, because its between-imputation variance is too noisy with few imputations.
#'
#' Assumptions: the pilot is a simple random sample of these records; the analysis is linear
#' MI with Rubin's rules through [catmech_pool_lm()]; the reference is exact and measures the
#' same quantity; the later audit is random and is imputed from the same model. The formula
#' is a large-sample one. It understates small audits, most of all when discarding recorded
#' values, because little data then drives a large share of the variance.
#'
#' The intervals answer different questions. `relative_lower` and `relative_upper` come from a
#' bias-corrected and accelerated bootstrap for `k`, so they describe the expected width. A
#' single small audit can land far from that, because its own residual variance is noisy;
#' `predictive_lower` and `predictive_upper` try to include this. A pilot can miss rare errors
#' altogether, and then both ranges are too low and too narrow, so neither is a guarantee.
#'
#' The function refuses a pilot with fewer than 30 rows holding both values or fewer than 20
#' residual degrees of freedom in the imputation model, and a pilot design that is rank
#' deficient. It warns below 100 rows or 10 disagreements, when no pilot row disagrees, when
#' the pilot looks selective on the recorded value (either z above 3.29, for its mean or spread),
#' and when the outcome is missing from `predictors`. It also warns when gross errors could
#' distort a linear fit: a pilot discrepancy above 5 robust SDs of the reference, or a recorded
#' value above 10 robust SDs from the median. Linear MI failed in the synthetic gross-error cell
#' of the contamination study, so do not rely on a plan carrying that warning.
#' @export
catmech_audit_plan <- function(data, value, reference, predictors, formula, term,
                               fractions = c(.05, .10, .25, .50, .75), target_ratio = 1.25,
                               conf_level = .95, bootstrap = 500L, m = 20L, seed = 1L) {
  catmech_reference_inputs(data, value, predictors, m, seed, 1L, excluded = c(value, reference))
  if (!is.character(reference) || length(reference) != 1L || is.na(reference) ||
      !nzchar(reference) || !reference %in% names(data) || reference == value)
    stop("reference must name an existing column different from the recorded value.", call. = FALSE)
  gold <- data[[reference]]
  if (!is.numeric(gold) || any(!is.na(gold) & !is.finite(gold)))
    stop("reference must contain finite numeric values or NA.", call. = FALSE)
  if (!inherits(formula, "formula") || length(formula) != 3L)
    stop("formula must be a two-sided R formula.", call. = FALSE)
  if (!is.character(term) || length(term) != 1L || is.na(term) || !nzchar(term))
    stop("term must name one coefficient.", call. = FALSE)
  if (!is.numeric(fractions) || !length(fractions) || anyNA(fractions) ||
      any(fractions <= 0 | fractions > 1) || anyDuplicated(fractions))
    stop("fractions must be distinct numbers above zero and at most one.", call. = FALSE)
  if (!is.numeric(target_ratio) || length(target_ratio) != 1L || !is.finite(target_ratio) || target_ratio <= 1)
    stop("target_ratio must be one number above one.", call. = FALSE)
  if (!is.numeric(conf_level) || length(conf_level) != 1L || !is.finite(conf_level) ||
      conf_level <= 0 || conf_level >= 1)
    stop("conf_level must be strictly between zero and one.", call. = FALSE)
  check_count(bootstrap, "bootstrap", 100L)
  if (seed > .Machine$integer.max - 3L) stop("seed is too large.", call. = FALSE)
  fractions <- sort(fractions)
  W <- data[[value]]
  known <- !is.na(gold)
  recorded <- !is.na(W)
  n <- nrow(data)
  if (!any(known)) stop("At least one reference value is required.", call. = FALSE)
  if (all(known)) stop("Every record already has a reference value, so there is nothing to plan.", call. = FALSE)
  if (!any(recorded)) stop("No record has a recorded value.", call. = FALSE)
  calibration <- known & recorded
  fill <- ifelse(known, gold, ifelse(recorded, W, mean(gold[known])))
  analysis <- catmech_audit_analysis(data, value, formula, term, fill)
  G0 <- catmech_model_matrix(data, predictors)
  G <- cbind(G0, .catmech_recorded = ifelse(recorded, W, 0))
  needed <- max(30L, ncol(G) + 20L)
  if (sum(calibration) < needed)
    stop("The pilot has ", sum(calibration), " rows with both values and needs at least ", needed,
         " to estimate k: 30, and 20 more than the ", ncol(G), " imputation coefficients. Audit more random records.",
         call. = FALSE)

  flags <- character()
  flag <- function(name, message) {
    flags <<- c(flags, name)
    warning(message, call. = FALSE)
  }
  agree <- catmech_contamination_agreement(W[calibration], gold[calibration], 0)
  outcome_names <- all.vars(formula[[2L]])
  if (!all(setdiff(outcome_names, value) %in% predictors))
    flag("outcome_not_in_predictors", "The analysis outcome is not in predictors. The formula for k assumes the imputation model uses it, and omitting it attenuates the coefficient.")
  if (sum(calibration) < 100L)
    flag("small_pilot", paste0("Small pilot: k rests on ", sum(calibration), " rows, so its interval is wide. Below 100 rows, plan on the upper limit."))
  if (all(agree)) {
    flag("no_disagreements", paste0("No pilot row disagrees, so k is estimated as zero. With ", sum(calibration),
      " rows, the disagreement rate could still be as high as about ", signif(3 / sum(calibration), 2), "."))
  } else if (sum(!agree) < 10L) {
    flag("few_disagreements", paste0("Only ", sum(!agree), " pilot rows disagree, so k depends heavily on them."))
  }
  selection <- catmech_audit_selection(W, known, recorded)
  if (is.finite(selection) && selection > 3.29)
    flag("pilot_looks_selective", paste0("Recorded values on pilot rows differ from those on other rows (largest z = ", round(selection, 1), "). The planner assumes a simple random pilot."))
  gross <- catmech_audit_gross(W, gold, known, recorded)
  if (gross$pilot_ratio > 5 || gross$recorded_z > 10)
    flag("gross_errors", paste0("Possible gross errors: the largest pilot discrepancy is ", round(gross$pilot_ratio, 1),
      " robust SDs of the reference, and the most extreme recorded value is ", round(gross$recorded_z, 1),
      " robust SDs from the median. These can distort a linear fit, and linear MI failed in the synthetic gross-error cell. Do not rely on this plan until catmech_recording_audit() has been reviewed; catmech_contamination() is the alternative."))

  restore <- preserve_rng()
  on.exit(restore(), add = TRUE)
  fit_width <- function(use_recorded, fit_seed) {
    # The planner raises its own no_disagreements flag for the case these would warn about.
    fit <- suppressWarnings(catmech_reference(data, value, reference, predictors, use_recorded = use_recorded, m = m, seed = fit_seed))
    row <- suppressWarnings(catmech_pool_lm(fit, formula, conf_level = conf_level))
    row <- row[row$term == term, ]
    list(width = row$upper - row$lower, within = row$within_variance, between = row$between_variance,
         full_width = 2 * stats::qt((1 + conf_level) / 2, n - ncol(analysis$D0)) * sqrt(row$within_variance),
         note = if (use_recorded) fit$spec$recorded_fit else NA_character_)
  }
  pilot_fit <- fit_width(TRUE, as.integer(seed))
  discard_fit <- fit_width(FALSE, as.integer(seed) + 2L)

  fit_rows <- which(known)
  point <- catmech_audit_imputation(fit_rows, G, G0, gold, recorded)
  if (is.null(point)) stop("The pilot's imputation design is rank deficient.", call. = FALSE)
  point0 <- catmech_audit_imputation(fit_rows, NULL, G0, gold, recorded)
  if (is.null(point0)) stop("The pilot's imputation design without the recorded value is rank deficient.", call. = FALSE)
  k <- catmech_audit_k(analysis, point$mu, point$tau2, known, gold)
  k0 <- catmech_audit_k(analysis, point0$mu, point0$tau2, known, gold)
  classical <- catmech_audit_k_classical(analysis, fit_rows, G, G0, gold, recorded, point, known)
  classical0 <- catmech_audit_k_classical(analysis, fit_rows, NULL, G0, gold, recorded, point0, known)

  set.seed(as.integer(seed) + 3L)
  records <- pmin(n, ceiling(fractions * n - 1e-8))
  draws <- catmech_audit_bootstrap(analysis, fit_rows, G, G0, gold, recorded, known, bootstrap, m, fractions, records)
  limits <- c((1 - conf_level) / 2, (1 + conf_level) / 2)
  usable <- sum(is.finite(draws$k))
  if (usable < .8 * bootstrap)
    flag("bootstrap_lost_support", paste0("Only ", usable, " of ", bootstrap, " pilot resamples could be fitted, usually because a category is rare. The interval for k is conditional on those resamples."))
  k_interval <- catmech_audit_bca(draws$k, k, fit_rows, function(rows) {
    p <- catmech_audit_imputation(rows, G, G0, gold, recorded)
    if (is.null(p)) NA_real_ else catmech_audit_k(analysis, p$mu, p$tau2, known, gold)
  }, limits)
  k0_interval <- catmech_audit_bca(draws$k0, k0, fit_rows, function(rows) {
    p <- catmech_audit_imputation(rows, NULL, G0, gold, recorded)
    if (is.null(p)) NA_real_ else catmech_audit_k(analysis, p$mu, p$tau2, known, gold)
  }, limits)
  relative <- function(z, f) sqrt(1 + z * (1 - f) / f)
  predictive <- apply(draws$predictive, 2L, stats::quantile, probs = limits, na.rm = TRUE)
  table <- data.frame(fraction = fractions, records = records,
    additional_records = pmax(0L, records - sum(known)), small_audit = records < needed,
    relative_width = relative(k, fractions), relative_lower = relative(k_interval[1], fractions),
    relative_upper = relative(k_interval[2], fractions),
    predictive_lower = predictive[1L, ], predictive_upper = predictive[2L, ],
    width = pilot_fit$full_width * relative(k, fractions),
    width_lower = pilot_fit$full_width * relative(k_interval[1], fractions),
    width_upper = pilot_fit$full_width * relative(k_interval[2], fractions),
    discard_relative_width = relative(k0, fractions), discard_relative_lower = relative(k0_interval[1], fractions),
    discard_relative_upper = relative(k0_interval[2], fractions),
    discard_width = discard_fit$full_width * relative(k0, fractions), row.names = NULL)
  required <- function(z) z / (target_ratio^2 - 1 + z)
  reaching <- table$fraction[table$relative_width <= target_ratio]
  needed_records <- function(z) as.integer(pmin(n, ceiling(required(z) * n - 1e-8)))
  answer <- data.frame(target_ratio = target_ratio, fraction = required(k), records = needed_records(k),
    records_lower = needed_records(k_interval[1]), records_upper = needed_records(k_interval[2]),
    smallest_grid_fraction = if (length(reaching)) min(reaching) else NA_real_, row.names = NULL)
  pilot_fraction <- sum(known) / n
  check <- function(fit) (1 + 1 / m) * fit$between / fit$within * pilot_fraction / (1 - pilot_fraction)
  k_table <- data.frame(approach = c("recorded", "discard"), estimate = c(k, k0),
    lower = c(k_interval[1], k0_interval[1]), upper = c(k_interval[2], k0_interval[2]),
    classical = c(classical, classical0), rubin_check = c(check(pilot_fit), check(discard_fit)),
    row.names = NULL)
  assumptions <- c("The pilot is a simple random sample of these records.",
    "The analysis is linear MI with Rubin's rules, using the same formula and predictors.",
    "The reference is exact and measures the same quantity as the recorded value.",
    "The relative widths come from a large-sample formula; small audits run wider and vary more.")
  list(status = if (length(flags)) "review_flags" else "no_flags", flags = flags, k = k_table, table = table,
    answer = answer,
    pilot = data.frame(records = n, reference_rows = sum(known), calibration_rows = sum(calibration),
      disagreements = sum(!agree), width = pilot_fit$width, full_width = pilot_fit$full_width,
      within_variance = pilot_fit$within, between_variance = pilot_fit$between,
      discard_width = discard_fit$width, discard_full_width = discard_fit$full_width,
      largest_discrepancy_ratio = gross$pilot_ratio, largest_recorded_z = gross$recorded_z,
      selection_z = selection, bootstrap_usable = usable, recorded_fit = pilot_fit$note, row.names = NULL),
    assumptions = assumptions)
}

# The analysis as functions of the audited value. Linear in the value means the design is
# D0 + y * D1 and the offset-adjusted response O0 + y * O1, row by row,
# which gives every expectation below in closed form.
catmech_audit_analysis <- function(data, value, formula, term, fill) {
  at <- function(y) {
    d <- data
    d[[value]] <- y
    frame <- stats::model.frame(formula, d, na.action = stats::na.fail)
    response <- as.vector(stats::model.response(frame, "numeric"))
    offset <- stats::model.offset(frame)
    if (!is.null(offset)) response <- response - offset
    list(D = stats::model.matrix(formula, frame), O = response)
  }
  n <- nrow(data)
  parts <- tryCatch({
    zero <- suppressWarnings(at(rep(0, n)))
    one <- suppressWarnings(at(rep(1, n)))
    test <- suppressWarnings(at(fill))
    list(zero = zero, one = one, test = test)
  }, error = function(e) {
    stop("The analysis could not be evaluated. Its variables other than the audited value must be complete, and it must be linear in that value: ",
         conditionMessage(e), call. = FALSE)
  })
  D0 <- parts$zero$D
  O0 <- parts$zero$O
  D1 <- parts$one$D - D0
  O1 <- parts$one$O - O0
  linear <- all(is.finite(D0), is.finite(O0), is.finite(D1), is.finite(O1))
  if (linear) {
    scale <- max(1, abs(parts$test$D), abs(parts$test$O))
    linear <- max(abs(parts$test$D - (D0 + fill * D1))) <= 1e-8 * scale &&
      max(abs(parts$test$O - (O0 + fill * O1))) <= 1e-8 * scale
  }
  if (!linear) stop("The formula must be linear in the audited value. Transform the column first, then pass the transformed column as value.", call. = FALSE)
  if (all(D1 == 0) && all(O1 == 0)) stop("The formula does not use the audited value.", call. = FALSE)
  j <- match(term, colnames(D0))
  if (is.na(j))
    stop("term is not a coefficient of the analysis. Coefficients: ", paste(colnames(D0), collapse = ", "), call. = FALSE)
  if (nrow(D0) <= ncol(D0)) stop("The analysis has no residual degrees of freedom.", call. = FALSE)
  list(D0 = D0, D1 = D1, O0 = O0, O1 = O1, j = j, n = n)
}

catmech_audit_ols <- function(X, y) {
  if (nrow(X) <= ncol(X)) return(NULL)
  q <- qr(X, tol = 1e-7)
  if (q$rank != ncol(X)) return(NULL)
  beta <- qr.coef(q, y)
  residuals <- y - as.vector(X %*% beta)
  list(beta = beta, residuals = residuals, df = nrow(X) - ncol(X), variance = sum(residuals^2) / (nrow(X) - ncol(X)))
}

# Imputation means and variances for every record, from the given reference rows. Records
# without a recorded value, and all records when G is NULL, use the model without it.
catmech_audit_imputation <- function(rows, G, G0, gold, recorded, bootstrap = FALSE) {
  n <- length(gold)
  needs_plain <- is.null(G) || !all(recorded)
  spread <- function(fit, X) if (bootstrap) sum(fit$residuals^2) / max(fit$df - ncol(X), 1) else fit$variance
  mu <- numeric(n)
  tau2 <- numeric(n)
  residuals <- NULL
  if (needs_plain) {
    plain <- catmech_audit_ols(G0[rows, , drop = FALSE], gold[rows])
    if (is.null(plain)) return(NULL)
    mu <- as.vector(G0 %*% plain$beta)
    tau2 <- rep(spread(plain, G0), n)
    residuals <- plain$residuals
  }
  if (!is.null(G)) {
    use <- rows[recorded[rows]]
    W <- G[, ".catmech_recorded"]
    if (!length(use)) return(NULL)
    if (all(catmech_contamination_agreement(W[use], gold[use], 0))) {
      mu[recorded] <- W[recorded]
      tau2[recorded] <- 0
      residuals <- numeric(length(use))
    } else {
      fit <- catmech_audit_ols(G[use, , drop = FALSE], gold[use])
      if (is.null(fit)) return(NULL)
      mu[recorded] <- as.vector(G[recorded, , drop = FALSE] %*% fit$beta)
      tau2[recorded] <- spread(fit, G)
      residuals <- fit$residuals
    }
  }
  list(mu = mu, tau2 = tau2, residuals = residuals)
}

# k under the working models. Known rows enter the analysis fit at their reference value,
# other rows through the mean and variance of their imputation; the numerator averages the
# imputation variance of the influence function over all records, using auxiliaries only.
catmech_audit_k <- function(a, mu, tau2, known, gold, detail = FALSE) {
  y <- mu
  y[known] <- gold[known]
  extra <- tau2
  extra[known] <- 0
  D <- a$D0 + y * a$D1
  O <- a$O0 + y * a$O1
  root <- sqrt(extra)
  DtD <- crossprod(D) + crossprod(a$D1 * root)
  DtO <- crossprod(D, O) + crossprod(a$D1 * root, a$O1 * root)
  inverse <- tryCatch(chol2inv(chol(DtD)), error = function(e) NULL)
  if (is.null(inverse)) return(if (detail) NULL else NA_real_)
  beta <- as.vector(inverse %*% DtO)
  sse <- sum((O - as.vector(D %*% beta))^2) + sum(extra * (a$O1 - as.vector(a$D1 %*% beta))^2)
  s2 <- sse / (a$n - ncol(D))
  weight <- a$n * inverse[a$j, ]
  A0 <- as.vector(a$D0 %*% weight)
  A1 <- as.vector(a$D1 %*% weight)
  E0 <- a$O0 - as.vector(a$D0 %*% beta)
  E1 <- a$O1 - as.vector(a$D1 %*% beta)
  u <- A0 + A1 * mu
  v <- E0 + E1 * mu
  variance <- tau2 * (u * E1 + v * A1)^2 + 2 * (A1 * E1)^2 * tau2^2
  value <- mean(variance) / (s2 * weight[a$j])
  if (detail) list(k = value, A0 = A0, A1 = A1, E0 = E0, E1 = E1, s2 = s2) else value
}

# The same ratio with empirical variances on pilot rows, using leave-one-out imputation means.
catmech_audit_k_classical <- function(a, rows, G, G0, gold, recorded, point, known) {
  parts <- catmech_audit_k(a, point$mu, point$tau2, known, gold, detail = TRUE)
  if (is.null(parts)) return(NA_real_)
  loo <- rep(NA_real_, length(gold))
  leave_out <- function(X, use) {
    q <- qr(X[use, , drop = FALSE], tol = 1e-7)
    h <- rowSums(qr.Q(q)^2)
    fit <- as.vector(X[use, , drop = FALSE] %*% qr.coef(q, gold[use]))
    gold[use] - (gold[use] - fit) / pmax(1 - h, 1e-8)
  }
  plain <- rows[!recorded[rows]]
  if (is.null(G)) loo[rows] <- leave_out(G0, rows) else {
    use <- rows[recorded[rows]]
    loo[use] <- leave_out(G, use)
    if (length(plain)) loo[rows] <- ifelse(recorded[rows], loo[rows], leave_out(G0, rows))
  }
  y <- gold[rows]
  psi <- (parts$A0[rows] + parts$A1[rows] * y) * (parts$E0[rows] + parts$E1[rows] * y)
  m <- loo[rows]
  expected <- (parts$A0[rows] + parts$A1[rows] * m) * (parts$E0[rows] + parts$E1[rows] * m) +
    parts$A1[rows] * parts$E1[rows] * point$tau2[rows]
  mean((psi - expected)^2) / mean(psi^2)
}

catmech_audit_selection <- function(W, known, recorded) {
  use <- recorded
  if (sum(known & use) < 2L || sum(!known & use) < 2L) return(NA_real_)
  x <- W[use]
  pilot <- known[use]
  z <- function(v) {
    spread <- stats::var(v) / sum(pilot) * (1 - sum(pilot) / length(v))
    if (!is.finite(spread) || spread <= 0) return(0)
    abs(mean(v[pilot]) - mean(v)) / sqrt(spread)
  }
  max(z(x), z((x - mean(x))^2))
}

catmech_audit_gross <- function(W, gold, known, recorded) {
  scale <- function(v) {
    s <- stats::mad(v)
    if (!is.finite(s) || s <= 0) s <- stats::sd(v)
    if (!is.finite(s) || s <= 0) NA_real_ else s
  }
  both <- known & recorded
  spread <- scale(gold[known])
  recorded_spread <- scale(W[recorded])
  list(pilot_ratio = if (is.na(spread)) 0 else max(abs(W[both] - gold[both])) / spread,
       recorded_z = if (is.na(recorded_spread)) 0 else max(abs(W[recorded] - stats::median(W[recorded]))) / recorded_spread)
}

# Bootstrap of k, resampling reference rows. Refitting q coefficients on a resample with
# repeated rows leaves about n - 2q residual degrees of freedom, so replicates use that
# divisor; with n - q they ran about q / (n - q) low. The predictive draw then resamples an
# audit of each planned size from the replicate's residuals and adds what the fit itself
# adds at that size: the posterior variance draw and Rubin's noisy between-imputation variance.
catmech_audit_bootstrap <- function(a, rows, G, G0, gold, recorded, known, B, m, fractions, records) {
  k <- k0 <- rep(NA_real_, B)
  predictive <- matrix(NA_real_, B, length(fractions))
  r <- length(rows)
  nu <- pmax(records - ncol(G), 3)
  for (b in seq_len(B)) {
    sampled <- rows[sample.int(r, r, replace = TRUE)]
    p <- catmech_audit_imputation(sampled, G, G0, gold, recorded, bootstrap = TRUE)
    if (!is.null(p)) {
      k[b] <- catmech_audit_k(a, p$mu, p$tau2, known, gold)
      if (is.finite(k[b])) {
        res2 <- p$residuals^2
        total <- sum(res2)
        scale <- if (total > 0) vapply(records, function(size) {
          mean(sample(res2, size, replace = TRUE)) / mean(res2)
        }, numeric(1)) else rep(0, length(records))
        noise <- stats::rchisq(length(records), m - 1L) / (m - 1L)
        inflation <- (1 + 1 / m) * nu / (nu - 2)
        predictive[b, ] <- sqrt(1 + k[b] * scale * noise * inflation * (1 - fractions) / fractions)
      }
    }
    p0 <- catmech_audit_imputation(sampled, NULL, G0, gold, recorded, bootstrap = TRUE)
    if (!is.null(p0)) k0[b] <- catmech_audit_k(a, p0$mu, p0$tau2, known, gold)
  }
  list(k = k, k0 = k0, predictive = predictive)
}

# Bias-corrected and accelerated limits. The acceleration comes from a jackknife over at most
# 200 random groups of reference rows, and k stays at its limit when the draws are degenerate.
catmech_audit_bca <- function(draws, point, rows, estimate, limits) {
  draws <- draws[is.finite(draws)]
  if (!length(draws) || !is.finite(point) || length(unique(draws)) < 2L) return(c(point, point))
  groups <- split(sample(rows), rep_len(seq_len(min(length(rows), 200L)), length(rows)))
  jack <- vapply(groups, function(out) estimate(setdiff(rows, out)), numeric(1))
  jack <- jack[is.finite(jack)]
  gap <- mean(jack) - jack
  denominator <- 6 * sum(gap^2)^1.5
  acceleration <- if (denominator > 0) sum(gap^3) / denominator else 0
  bias <- stats::qnorm(min(max(mean(draws < point), 1e-4), 1 - 1e-4))
  z <- stats::qnorm(limits)
  unname(stats::quantile(draws, stats::pnorm(bias + (bias + z) / (1 - acceleration * (bias + z)))))
}

# k on data where every row has a reference: the pool-level value a pilot estimates.
catmech_audit_pool_k <- function(data, value, reference, predictors, formula, term) {
  gold <- data[[reference]]
  W <- data[[value]]
  known <- !is.na(gold)
  if (!all(known)) stop("Every row needs a reference value.", call. = FALSE)
  recorded <- !is.na(W)
  a <- catmech_audit_analysis(data, value, formula, term, gold)
  G0 <- catmech_model_matrix(data, predictors)
  G <- cbind(G0, .catmech_recorded = ifelse(recorded, W, 0))
  rows <- which(known)
  p <- catmech_audit_imputation(rows, G, G0, gold, recorded)
  p0 <- catmech_audit_imputation(rows, NULL, G0, gold, recorded)
  c(k = catmech_audit_k(a, p$mu, p$tau2, known, gold), discard_k = catmech_audit_k(a, p0$mu, p0$tau2, known, gold),
    classical = catmech_audit_k_classical(a, rows, G, G0, gold, recorded, p, known),
    discard_classical = catmech_audit_k_classical(a, rows, NULL, G0, gold, recorded, p0, known))
}
