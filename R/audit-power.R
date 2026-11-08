#' Rough audit-size planning aid (not validated)
#'
#' Evaluates, over a grid of audit fractions, how often the expected interval width stays
#' within a target ratio of a full audit's width (\eqn{\rho \le \rho_{target}}), using the
#' two-phase relation \eqn{\rho = \sqrt{1 + k (1 - f) / f}}. The factor \eqn{k} comes from a
#' fitted [catmech_audit_plan()] or from assumed error parameters.
#'
#' \strong{Status: not validated.} No study has checked this function's recommendations, and
#' the planner it builds on, [catmech_audit_plan()], returned NO-GO under a rule fixed in
#' advance. Treat the output as a rough guide, not a sample-size guarantee.
#'
#' Without `plan`, \eqn{k = q / (1 - q)} with \eqn{q} = `error_rate` (1 - `correlation`^2)
#' (1 - `r_squared`), where `correlation` is the correlation between recorded and true values
#' over all records. Multiplying by the error rate treats the unexplained variation as if it
#' sat in erroneous records the analysis can identify and correct. That is optimistic: on
#' ClinicalTrials.gov and NHANES the contamination model could not tell which records were
#' wrong (AUC 0.54 and 0.50). For linear imputation given the recorded value, the package's
#' recommended start, the double-sampling regression estimator has
#' \eqn{k \approx (1 - \rho^2)(1 - R^2)} at the same correlation \eqn{\rho}, about
#' 1 / `error_rate` times larger when errors are rare. The result reports the audit fraction
#' that this \eqn{k} needs (`linear_mi`), and the printout shows both. That figure ignores the
#' uncertainty in \eqn{k}, so it can fall below the main recommendation, for example at very
#' high error rates.
#'
#' "Power" here is the share of draws of \eqn{k}, from an assumed log-normal spread around
#' its estimate, whose width ratio meets the target. It is not the power of a hypothesis test.
#'
#' @param plan Optional fitted `catmech_audit_plan` object. If supplied, \eqn{k} is drawn from
#'   a log-normal distribution fitted to the plan's estimate and interval for \eqn{k}.
#' @param N Total cohort size (number of patient records). Required if `plan` is `NULL`.
#' @param target_ratio Maximum acceptable confidence interval width relative to an ideal
#'   full audit (default 1.25, meaning within 25\% of the full-audit width).
#' @param power_target Target share of plausible \eqn{k} values whose width ratio meets
#'   `target_ratio` (default 0.80); see Details for why this is not the power of a test.
#' @param error_rate Anticipated error or corruption rate \eqn{p \in (0, 1)} (default 0.15).
#'   Ignored if `plan` is supplied.
#' @param correlation Anticipated correlation between recorded and true values over all
#'   records (default 0.85). Ignored if `plan` is supplied.
#' @param r_squared Anticipated auxiliary predictor \eqn{R^2} (default 0.30).
#'   Ignored if `plan` is supplied.
#' @param fractions Numeric vector of audit fractions to evaluate, strictly between 0 and 1.
#'   Defaults to `seq(0.02, 0.40, by = 0.02)`.
#' @param n_sim Number of Monte Carlo draws of \eqn{k} (default 2000).
#' @param seed Optional integer seed for reproducible draws.
#' @return An object of class `catmech_audit_power` containing:
#'   \item{table}{Data frame of evaluated audit fractions, audit sizes, expected width ratios
#'     with 95\% ranges, and `power`, the share of \eqn{k} draws meeting `target_ratio`.}
#'   \item{recommended}{Smallest evaluated audit fraction and size whose `power` reaches
#'     `power_target`, or the largest fraction if none does.}
#'   \item{target_ratio}{The target width ratio evaluated.}
#'   \item{power_target}{The target share of \eqn{k} draws.}
#'   \item{N}{Total cohort size.}
#'   \item{k_summary}{Point estimate and uncertainty bounds for the two-phase variance factor \eqn{k}.}
#'   \item{linear_mi}{Without `plan`: the linear-imputation \eqn{k} described above and the
#'     audit fraction and size it needs to meet `target_ratio`; `NULL` with `plan`.}
#'   \item{status}{Always `"not_validated"`.}
#' @examples
#' pwr <- catmech_audit_power(N = 500, target_ratio = 1.3, n_sim = 200, seed = 1)
#' print(pwr)
#' @export
catmech_audit_power <- function(plan = NULL, N = NULL, target_ratio = 1.25,
                                power_target = 0.80, error_rate = 0.15,
                                correlation = 0.85, r_squared = 0.30,
                                fractions = seq(0.02, 0.40, by = 0.02),
                                n_sim = 2000L, seed = 42L) {
  if (!is.numeric(target_ratio) || length(target_ratio) != 1L || !is.finite(target_ratio) || target_ratio <= 1)
    stop("target_ratio must be a single number strictly greater than 1.", call. = FALSE)
  if (!is.numeric(power_target) || length(power_target) != 1L || !is.finite(power_target) ||
      power_target <= 0 || power_target >= 1)
    stop("power_target must be strictly between 0 and 1.", call. = FALSE)
  if (!is.numeric(fractions) || length(fractions) == 0L || any(!is.finite(fractions)) ||
      any(fractions <= 0) || any(fractions > 1))
    stop("fractions must be numeric values strictly between 0 and 1.", call. = FALSE)
  fractions <- sort(unique(fractions))

  if (!is.null(seed)) {
    if (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed))
      stop("seed must be a single integer.", call. = FALSE)
    old_seed <- if (exists(".Random.seed", envir = .GlobalEnv)) get(".Random.seed", envir = .GlobalEnv) else NULL
    set.seed(as.integer(seed))
    on.exit({
      if (is.null(old_seed)) {
        if (exists(".Random.seed", envir = .GlobalEnv)) rm(".Random.seed", envir = .GlobalEnv)
      } else {
        assign(".Random.seed", old_seed, envir = .GlobalEnv)
      }
    })
  }

  linear_mi <- NULL
  if (!is.null(plan)) {
    if (!inherits(plan, "catmech_audit_plan"))
      stop("plan must be a catmech_audit_plan object.", call. = FALSE)
    col_key <- if ("approach" %in% names(plan$k)) "approach" else "strategy"
    k_row <- if (col_key %in% names(plan$k)) {
      plan$k[plan$k[[col_key]] == "recorded", , drop = FALSE]
    } else {
      plan$k[1L, , drop = FALSE]
    }
    if (nrow(k_row) == 0L) k_row <- plan$k[1L, , drop = FALSE]
    k_est <- k_row$estimate[1L]
    k_lo <- max(1e-6, k_row$lower[1L])
    k_hi <- max(k_lo * 1.01, k_row$upper[1L])
    N <- if (is.data.frame(plan$pilot) && "records" %in% names(plan$pilot)) {
      as.integer(plan$pilot$records[1L])
    } else if (is.list(plan$pilot) && !is.null(plan$pilot$fit$original)) {
      nrow(plan$pilot$fit$original)
    } else if (!is.null(N)) {
      as.integer(N)
    } else {
      stop("Unable to determine cohort size N from plan.", call. = FALSE)
    }
    log_se <- (log(k_hi) - log(k_lo)) / (2 * 1.96)
    k_draws <- stats::rlnorm(n_sim, meanlog = log(k_est), sdlog = max(0.05, log_se))
    k_summary <- list(estimate = k_est, lower = k_lo, upper = k_hi, source = "empirical_plan")
  } else {
    if (is.null(N) || !is.numeric(N) || length(N) != 1L || !is.finite(N) || N < 10)
      stop("N must be a single integer >= 10 when plan is not supplied.", call. = FALSE)
    N <- as.integer(N)
    if (!is.numeric(error_rate) || length(error_rate) != 1L || error_rate <= 0 || error_rate >= 1)
      stop("error_rate must be strictly between 0 and 1.", call. = FALSE)
    if (!is.numeric(correlation) || length(correlation) != 1L || correlation < 0 || correlation >= 1)
      stop("correlation must be in [0, 1).", call. = FALSE)
    if (!is.numeric(r_squared) || length(r_squared) != 1L || r_squared < 0 || r_squared >= 1)
      stop("r_squared must be in [0, 1).", call. = FALSE)

    # Theoretical two-phase variance factor k from information loss:
    # k = Var(psi - E[psi|aux]) / Var(psi)
    unexplained_noise <- error_rate * (1 - correlation^2) * (1 - r_squared)
    k_est <- max(0.01, unexplained_noise / (1 - unexplained_noise))
    # Pilot sampling variance for k with pilot size ~ 5% of cohort
    n_pilot <- max(20L, as.integer(round(0.05 * N)))
    cv_k <- sqrt(2 / (n_pilot * error_rate))
    log_se <- sqrt(log(1 + cv_k^2))
    k_draws <- stats::rlnorm(n_sim, meanlog = log(k_est), sdlog = log_se)
    k_summary <- list(estimate = k_est,
                      lower = stats::quantile(k_draws, 0.025),
                      upper = stats::quantile(k_draws, 0.975),
                      source = "prior_specification")
    # Double-sampling regression estimator for linear imputation given the recorded value,
    # at the same overall correlation.
    k_linear <- (1 - correlation^2) * (1 - r_squared)
    f_linear <- min(1, k_linear / (target_ratio^2 - 1 + k_linear))
    linear_mi <- list(k = k_linear, fraction = f_linear,
                      n_audit = as.integer(ceiling(f_linear * N)),
                      basis = "double_sampling_regression_estimator")
  }

  rows <- lapply(fractions, function(f) {
    n_audit <- as.integer(ceiling(f * N))
    # Relative width ratio: rho = sqrt(1 + k * (1 - f) / f)
    ratios <- sqrt(1 + k_draws * (1 - f) / f)
    pwr <- mean(ratios <= target_ratio)
    data.frame(
      fraction = f,
      n_audit = n_audit,
      expected_ratio = mean(ratios),
      ratio_lower = unname(stats::quantile(ratios, 0.025)),
      ratio_upper = unname(stats::quantile(ratios, 0.975)),
      power = pwr,
      stringsAsFactors = FALSE
    )
  })
  tbl <- do.call(rbind, rows)

  # Determine recommended sample size
  eligible <- which(tbl$power >= power_target)
  rec <- if (length(eligible) > 0L) {
    tbl[eligible[1L], , drop = FALSE]
  } else {
    # If target power not reached in grid, take highest fraction
    tbl[nrow(tbl), , drop = FALSE]
  }

  res <- list(
    table = tbl,
    recommended = rec,
    target_ratio = target_ratio,
    power_target = power_target,
    N = N,
    k_summary = k_summary,
    linear_mi = linear_mi,
    status = "not_validated"
  )
  class(res) <- "catmech_audit_power"
  res
}

#' Print an audit power calculation
#'
#' @param x A `catmech_audit_power` object.
#' @param digits Number of decimal digits to print (default 3).
#' @param ... Unused.
#' @return `x`, invisibly.
#' @export
print.catmech_audit_power <- function(x, digits = 3L, ...) {
  cat("-- Catmech Audit-Size Planning Aid (not validated) --\n")
  cat(sprintf("Cohort Size (N):    %d\n", x$N))
  cat(sprintf("Target Ratio:       <= %.2f (within %.0f%% of full audit width)\n",
              x$target_ratio, (x$target_ratio - 1) * 100))
  cat(sprintf("Target share of k:  %.0f%%\n", x$power_target * 100))
  cat(sprintf("Two-Phase Factor k: %.3f [%.3f, %.3f] (%s)\n\n",
              x$k_summary$estimate, x$k_summary$lower, x$k_summary$upper, x$k_summary$source))

  tbl_print <- x$table
  for (col in c("fraction", "expected_ratio", "ratio_lower", "ratio_upper", "power")) {
    tbl_print[[col]] <- round(tbl_print[[col]], digits)
  }
  print(tbl_print, row.names = FALSE)

  rec <- x$recommended
  if (rec$power >= x$power_target) {
    cat(sprintf("\nRecommended Audit Size: n = %d records (%.1f%% audit fraction)\n",
                rec$n_audit, rec$fraction * 100))
    cat(sprintf("Meets the target in %.1f%% of k draws (expected ratio: %.3f)\n",
                rec$power * 100, rec$expected_ratio))
  } else {
    cat(sprintf("\nWarning: no evaluated fraction meets the target in %.0f%% of k draws (max: %.1f%% at n = %d).\n",
                x$power_target * 100, rec$power * 100, rec$n_audit))
  }
  if (!is.null(x$linear_mi)) {
    cat(sprintf("Linear imputation (double-sampling regression, k = %.3f): about %.1f%% (n = %d) for ratio <= %.2f, ignoring uncertainty in k.\n",
                x$linear_mi$k, 100 * x$linear_mi$fraction, x$linear_mi$n_audit, x$target_ratio))
  }
  cat("Power here is the share of plausible k values meeting the target, not the power of a test.\n")
  cat("Not validated; see ?catmech_audit_power for the assumptions.\n")
  invisible(x)
}

#' Plot an audit power curve
#'
#' Displays, against audit size, the share of plausible \eqn{k} values whose width ratio meets
#' the target, with reference lines. Not validated; see [catmech_audit_power()].
#'
#' @param x A `catmech_audit_power` object.
#' @param col Line color for the power curve (default `"#1f77b4"`).
#' @param lwd Line width (default 2).
#' @param ... Additional arguments passed to [graphics::plot()].
#' @return `x`, invisibly.
#' @examples
#' pwr <- catmech_audit_power(N = 500, target_ratio = 1.3, n_sim = 200, seed = 1)
#' plot(pwr)
#' @export
plot.catmech_audit_power <- function(x, col = "#1f77b4", lwd = 2, ...) {
  old_mar <- graphics::par("mar")
  on.exit(graphics::par(mar = old_mar))
  new_mar <- old_mar
  new_mar[4L] <- max(new_mar[4L], 4.5)
  graphics::par(mar = new_mar)

  tbl <- x$table
  graphics::plot(
    tbl$n_audit, tbl$power,
    type = "b", pch = 19, col = col, lwd = lwd,
    xlab = sprintf("Audit Sample Size (n out of %d)", x$N),
    ylab = sprintf("Share of k draws with ratio <= %.2f", x$target_ratio),
    ylim = c(0, 1.05),
    main = "Audit-Size Curve (not validated)",
    panel.first = graphics::grid(col = "gray90", lty = "dotted"),
    ...
  )

  # Power target line
  graphics::abline(h = x$power_target, col = "darkred", lty = 2, lwd = 1.5)
  graphics::mtext(sprintf("Target share: %.0f%%", x$power_target * 100),
                  side = 4, at = x$power_target, col = "darkred", cex = 0.8, las = 1)

  # Recommended marker
  rec <- x$recommended
  if (rec$power >= x$power_target) {
    graphics::points(rec$n_audit, rec$power, col = "darkgreen", pch = 17, cex = 1.8)
    graphics::abline(v = rec$n_audit, col = "darkgreen", lty = 3)
    graphics::text(rec$n_audit, rec$power - 0.08,
                   sprintf("Recommended:\nn = %d (%.0f%% of k)", rec$n_audit, rec$power * 100),
                   col = "darkgreen", cex = 0.85, font = 2)
  }

  invisible(x)
}
