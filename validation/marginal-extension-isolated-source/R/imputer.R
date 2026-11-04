#' Category-constrained predictive mean matching with sensitivity replacement
#'
#' A univariate imputation method for an exactly coarsened continuous variable.
#' Missing values and declared threshold conflicts are reconstructed. An observed,
#' flagged value that agrees with its category is reconstructed with the supplied
#' sensitivity probability and otherwise retained exactly.
#'
#' @param y Numeric working response, as supplied by `mice`.
#' @param ry Logical vector identifying originally observed responses available
#'   for fitting. Flagged and category-inconsistent observations are excluded
#'   regardless of this vector.
#' @param x Numeric predictor matrix with one row per response. All entries must
#'   be finite; use a zero-column matrix for an intercept-only model.
#' @param wy Logical vector identifying cells to complete or reconstruct. Defaults
#'   to `!ry`. All missing values and all threshold conflicts must be selected.
#' @param ... Additional arguments supplied by `mice`, such as predictor types.
#' @param category Observed category labels, one per row, with no missing values.
#' @param bounds Named list with one entry per category. Each entry contains
#'   scalar numeric `lower` and `upper` endpoints and scalar logical
#'   `lower_inclusive` and `upper_inclusive` flags.
#' @param flag Optional quality flags, one per row. Missing flags do not mark a
#'   value as suspect; this is a substantive trust assumption.
#' @param flag_yes Label that marks a suspect observation.
#' @param corruption_probability Sensitivity replacement probability for flagged,
#'   observed, in-bound values, as one number or one number per row in `[0, 1]`.
#'   It is not an estimated or posterior probability of measurement error.
#'   Must be explicit when flag is supplied; NULL means zero only without flags.
#' @param original_y Original response before any imputation or reconstruction.
#'   Pass an immutable copy when using this method inside `mice`.
#' @param trusted Optional logical mask restricting the fitting and donor sample.
#'   It cannot admit flagged, missing, or category-inconsistent observations.
#' @param min_donors Minimum number of eligible donors required to use predictive
#'   mean matching for a recipient.
#' @param n_donors Number of nearest eligible donors from which to sample.
#' @param fallback Policy when fewer than `min_donors` are available: `"error"`
#'   (default), or an explicitly requested truncated `"normal"` draw.
#'
#' @details
#' The fitting model uses a rank-reduced Gaussian linear regression with the
#' usual reference-prior variance and coefficient draws. Donor distances compare
#' fitted donor predictions with drawn recipient predictions (type-1 matching).
#' Rank reduction is allowed only when reconstructed recipients satisfy the same
#' predictor dependencies as the trusted fitting rows; otherwise prediction is
#' unidentified and the method errors.
#' A positive residual degree of freedom and positive residual variance are
#' required whenever reconstruction is requested. A degenerate model produces an
#' error, rather than an invented residual variance.
#'
#' Donors always come from original, trusted observations in the recipient's
#' category. The method never falls back to contaminated observations. Sparse
#' donor support and model misspecification can still invalidate inference.
#' Constrained PMM and the retain/replace rule do not establish proper multiple
#' imputation or valid confidence-interval coverage. Study those properties for
#' the intended analysis. The normal fallback is a stronger extrapolation
#' assumption and is recorded in the audit.
#'
#' @return A numeric vector containing only the `wy` cells, in row order. Its
#'   `audit` attribute is a data frame with the original row, reason,
#'   reconstruction decision, actual change, donor row, fallback use, and eligible
#'   donor count. For a missing original value, `changed` is `NA`.
#' @export
mice.impute.catmech_threshold_corruption <- function(
    y, ry, x, wy = NULL, ..., category, bounds, flag = NULL, flag_yes = "Yes",
    corruption_probability = NULL, original_y = y, trusted = NULL,
    min_donors = 1L, n_donors = 5L, fallback = c("error", "normal")) {
  fallback <- match.arg(fallback)
  n <- length(y)
  if (!is.numeric(y) || n < 1L || !is.numeric(original_y) ||
      length(original_y) != n || any(!is.na(original_y) & !is.finite(original_y))) {
    stop("y and original_y must be numeric vectors of the same positive length; observed original_y must be finite.",
         call. = FALSE)
  }
  check_mask <- function(z, name) {
    if (!is.logical(z) || length(z) != n || anyNA(z)) {
      stop(name, " must be a nonmissing logical vector of length y.", call. = FALSE)
    }
  }
  check_mask(ry, "ry")
  if (is.null(wy)) wy <- !ry
  check_mask(wy, "wy")
  if (is.null(trusted)) trusted <- ry
  check_mask(trusted, "trusted")
  x <- as.matrix(x)
  if (!is.numeric(x) || length(dim(x)) != 2L || nrow(x) != n ||
      any(!is.finite(x))) {
    stop("x must be a numeric matrix with one row per response and only finite predictors.",
         call. = FALSE)
  }
  check_count <- function(z, name) {
    if (!is.numeric(z) || length(z) != 1L || !is.finite(z) ||
        z < 1 || z != floor(z)) {
      stop(name, " must be a positive integer.", call. = FALSE)
    }
  }
  check_count(min_donors, "min_donors")
  check_count(n_donors, "n_donors")
  q <- corruption_probability
  if (is.null(q)) {
    if (!is.null(flag)) stop("Supply corruption_probability explicitly when flag is provided.", call. = FALSE)
    q <- 0
  }
  if (!is.numeric(q) || !(length(q) %in% c(1L, n)) ||
      any(!is.finite(q)) || any(q < 0 | q > 1)) {
    stop("corruption_probability must be numeric, have length one or length y, and lie in [0, 1].",
         call. = FALSE)
  }
  q <- rep_len(q, n)
  category_chr <- as.character(category)
  if (!is.list(bounds) || is.null(names(bounds)) || anyNA(names(bounds)) ||
      any(names(bounds) == "") || anyDuplicated(names(bounds))) {
    stop("bounds must be a list with unique, nonempty category names.", call. = FALSE)
  }
  if (length(category_chr) != n || anyNA(category_chr) ||
      any(!category_chr %in% names(bounds))) {
    stop("category must have one known, nonmissing bounds label per row.", call. = FALSE)
  }
  valid_bound <- function(b) {
    is.list(b) && is.numeric(b$lower) && length(b$lower) == 1L &&
      !is.na(b$lower) && is.numeric(b$upper) && length(b$upper) == 1L &&
      !is.na(b$upper) && b$lower < b$upper &&
      is.logical(b$lower_inclusive) && length(b$lower_inclusive) == 1L &&
      !is.na(b$lower_inclusive) && is.logical(b$upper_inclusive) &&
      length(b$upper_inclusive) == 1L && !is.na(b$upper_inclusive)
  }
  if (!all(vapply(bounds, valid_bound, logical(1)))) {
    stop("Each bounds entry must define increasing scalar endpoints and logical inclusion flags.",
         call. = FALSE)
  }
  if (!is.character(flag_yes) || length(flag_yes) != 1L ||
      is.na(flag_yes) || !nzchar(flag_yes)) {
    stop("flag_yes must be one nonempty label.", call. = FALSE)
  }
  if (is.null(flag)) {
    suspect <- rep(FALSE, n)
  } else {
    if (length(flag) != n) stop("flag must have length y.", call. = FALSE)
    suspect <- catmech_suspect_mask(flag, flag_yes)
  }
  in_bounds <- vapply(seq_len(n), function(i) {
    isTRUE(value_in_bounds(original_y[i], bounds[[category_chr[i]]]))
  }, logical(1))
  observed <- !is.na(original_y)
  conflict <- observed & !in_bounds
  if (any(!wy & (!observed | conflict))) {
    stop("All missing target values and threshold conflicts must be selected by wy.",
         call. = FALSE)
  }
  fit_rows <- trusted & ry & observed & !suspect & in_bounds
  selected <- which(wy)
  reason <- rep("retained", length(selected))
  reason[!observed[selected]] <- "missing"
  reason[conflict[selected]] <- "threshold_conflict"
  stochastic <- which(observed[selected] & in_bounds[selected] & suspect[selected])
  if (length(stochastic)) {
    replace_draw <- stats::rbinom(length(stochastic), 1L, q[selected[stochastic]]) == 1L
    reason[stochastic[replace_draw]] <- "flag_draw"
  }
  reconstructed <- reason != "retained"
  out <- original_y[selected]
  eligible_counts <- vapply(selected, function(i) {
    sum(fit_rows & category_chr == category_chr[i] &
          value_in_bounds(original_y, bounds[[category_chr[i]]]), na.rm = TRUE)
  }, integer(1))
  audit <- data.frame(
    row = selected, reason = reason, reconstructed = reconstructed,
    changed = rep(FALSE, length(selected)), donor_row = rep(NA_integer_, length(selected)),
    used_fallback = rep(FALSE, length(selected)), eligible_donors = eligible_counts,
    stringsAsFactors = FALSE
  )
  if (!any(reconstructed)) {
    attr(out, "audit") <- audit
    return(out)
  }
  if (fallback == "error" && any(eligible_counts[reconstructed] < min_donors)) {
    bad <- selected[reconstructed & eligible_counts < min_donors]
    stop("Insufficient trusted within-category donors for row(s): ",
         paste(bad, collapse = ", "),
         ". Reduce min_donors only with justification or explicitly choose fallback = 'normal'.",
         call. = FALSE)
  }
  prediction_rows <- fit_rows
  prediction_rows[selected[reconstructed]] <- TRUE
  draw <- catmech_bayes_lin_draw(original_y, fit_rows, x, prediction_rows)
  for (j in which(reconstructed)) {
    i <- selected[j]
    bound <- bounds[[category_chr[i]]]
    eligible <- fit_rows & category_chr == category_chr[i] &
      value_in_bounds(original_y, bound)
    eligible[is.na(eligible)] <- FALSE
    if (eligible_counts[j] >= min_donors) {
      donor <- pick_pmm_donor(draw$pred_hat, eligible, draw$pred_draw[i], n_donors)
      out[j] <- original_y[donor]
      audit$donor_row[j] <- donor
    } else {
      out[j] <- catmech_truncated_draw(draw$pred_draw[i], draw$sigma, bound)
      audit$used_fallback[j] <- TRUE
    }
    if (!is.finite(out[j]) || !isTRUE(value_in_bounds(out[j], bound))) {
      stop("Internal imputation error: a reconstructed value violates its category.",
           call. = FALSE)
    }
    audit$changed[j] <- if (observed[i]) out[j] != original_y[i] else NA
  }
  attr(out, "audit") <- audit
  out
}

# Rank reduction avoids assigning arbitrary variance to unidentified coefficients.
catmech_bayes_lin_draw <- function(y, fit_rows, x,
                                  prediction_rows = rep(TRUE, length(y))) {
  if (!is.logical(prediction_rows) || length(prediction_rows) != length(y) ||
      anyNA(prediction_rows)) {
    stop("prediction_rows must be a nonmissing logical mask.", call. = FALSE)
  }
  X_all <- cbind(`(Intercept)` = 1, x)
  X <- X_all[fit_rows, , drop = FALSE]
  yobs <- y[fit_rows]
  if (nrow(X) < 2L) {
    stop("At least two trusted observations are needed for the regression draw.",
         call. = FALSE)
  }
  decomposition <- qr(X, tol = 1e-7)
  rank <- decomposition$rank
  columns <- decomposition$pivot[seq_len(rank)]
  dropped <- setdiff(seq_len(ncol(X)), columns)
  reduced_qr <- qr(X[, columns, drop = FALSE])
  if (length(dropped)) {
    relation <- qr.coef(reduced_qr, X[, dropped, drop = FALSE])
    implied <- X_all[prediction_rows, columns, drop = FALSE] %*% relation
    actual <- X_all[prediction_rows, dropped, drop = FALSE]
    tolerance <- 1e-7 * pmax(1, abs(actual), abs(implied))
    if (any(!is.finite(implied)) || any(abs(actual - implied) > tolerance)) {
      stop("Recipient design is outside the trusted regression column span; remove unsupported predictors or supply trusted data with appropriate support.",
           call. = FALSE)
    }
  }
  X <- X[, columns, drop = FALSE]
  X_all <- X_all[, columns, drop = FALSE]
  df <- nrow(X) - rank
  if (df <= 0L) {
    stop("Insufficient trusted observations for positive residual degrees of freedom.",
         call. = FALSE)
  }
  beta_hat <- as.vector(qr.coef(reduced_qr, yobs))
  residual <- as.vector(yobs - X %*% beta_hat)
  sse <- sum(residual^2)
  # Treat roundoff from an exact fit as zero instead of fabricating uncertainty.
  roundoff_floor <- .Machine$double.eps^2 * max(1, sum(yobs^2)) * nrow(X)
  if (!is.finite(sse) || sse <= roundoff_floor) {
    stop("Trusted regression has zero residual variance; no stochastic regression draw is supported.",
         call. = FALSE)
  }
  sigma <- sqrt(sse / stats::rchisq(1L, df = df))
  R <- qr.R(reduced_qr)
  # X[, pivot] = Q R, so backsolve produces a draw in pivoted coefficient order.
  noise_pivot <- backsolve(R, stats::rnorm(rank))
  noise <- numeric(rank)
  noise[reduced_qr$pivot] <- noise_pivot
  beta_draw <- beta_hat + sigma * noise
  pred_hat <- as.vector(X_all %*% beta_hat)
  pred_draw <- as.vector(X_all %*% beta_draw)
  if (!is.finite(sigma) || any(!is.finite(pred_hat)) || any(!is.finite(pred_draw))) {
    stop("Regression draw was nonfinite; inspect scale and trusted sample support.", call. = FALSE)
  }
  list(pred_hat = pred_hat, pred_draw = pred_draw, sigma = sigma,
       rank = rank, df = df)
}

pick_pmm_donor <- function(pred, donor_rows, target_pred, n_donors = 5L) {
  donor_idx <- which(donor_rows)
  if (!length(donor_idx)) stop("No eligible donor.", call. = FALSE)
  # Randomize ties so row order does not select an arbitrary subset of tied donors.
  shuffled <- donor_idx[sample.int(length(donor_idx))]
  distance <- abs(pred[shuffled] - target_pred)
  nearest <- shuffled[order(distance)][seq_len(min(n_donors, length(shuffled)))]
  nearest[sample.int(length(nearest), 1L)]
}

catmech_truncated_draw <- function(mean, sd, bound) {
  # Open endpoints have probability zero mathematically, but rounding can land
  # on one. Reject that value rather than moving a bound by a scale-dependent eps.
  for (attempt in seq_len(100L)) {
    value <- truncnorm::rtruncnorm(1L, a = bound$lower, b = bound$upper,
                                  mean = mean, sd = sd)
    if (is.finite(value) && isTRUE(value_in_bounds(value, bound))) return(value)
  }
  stop("Could not draw a representable value inside the requested interval.",
       call. = FALSE)
}
