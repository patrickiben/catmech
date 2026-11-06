#' Prospective sample size and statistical power calculation for clinical data audits
#'
#' Computes the statistical power and sample size required to achieve a target confidence
#' interval precision ratio (\eqn{\rho \le \rho_{target}}) in a two-phase clinical audit.
#' Can evaluate power prospectively from prior clinical error parameters or empirically
#' from a fitted [catmech_audit_plan()].
#'
#' @param plan Optional fitted `catmech_audit_plan` object. If supplied, the empirical
#'   distribution of \eqn{k} is used to evaluate power.
#' @param N Total cohort size (number of patient records). Required if `plan` is `NULL`.
#' @param target_ratio Maximum acceptable confidence interval width relative to an ideal
#'   full audit (default 1.25, meaning within 25\% of the full-audit width).
#' @param power_target Desired statistical power (default 0.80).
#' @param error_rate Anticipated error or corruption rate \eqn{p \in (0, 1)} (default 0.15).
#'   Ignored if `plan` is supplied.
#' @param correlation Anticipated correlation between recorded and true values (default 0.85).
#'   Ignored if `plan` is supplied.
#' @param r_squared Anticipated auxiliary predictor \eqn{R^2} (default 0.30).
#'   Ignored if `plan` is supplied.
#' @param fractions Numeric vector of audit fractions to evaluate, strictly between 0 and 1.
#'   Defaults to `seq(0.02, 0.40, by = 0.02)`.
#' @param n_sim Number of Monte Carlo draws for evaluating the sampling distribution
#'   of the precision ratio (default 2000).
#' @param seed Optional integer seed for reproducible power evaluation.
#' @return An object of class `catmech_audit_power` containing:
#'   \item{table}{Data frame of evaluated audit fractions, audit sample sizes, expected
#'     relative widths, and statistical power.}
#'   \item{recommended}{Recommended minimum audit fraction and sample size reaching `power_target`.}
#'   \item{target_ratio}{The target precision ratio evaluated.}
#'   \item{power_target}{The target statistical power.}
#'   \item{N}{Total cohort size.}
#'   \item{k_summary}{Point estimate and uncertainty bounds for the two-phase variance factor \eqn{k}.}
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

  if (!is.null(plan)) {
    if (!inherits(plan, "catmech_audit_plan"))
      stop("plan must be a catmech_audit_plan object.", call. = FALSE)
    k_row <- plan$k[plan$k$strategy == "recorded", , drop = FALSE]
    if (nrow(k_row) == 0L) k_row <- plan$k[1L, , drop = FALSE]
    k_est <- k_row$estimate[1L]
    k_lo <- max(1e-6, k_row$lower[1L])
    k_hi <- max(k_lo * 1.01, k_row$upper[1L])
    N <- nrow(plan$pilot$fit$original)
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
    k_summary = k_summary
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
  cat("-- Catmech Audit Power & Sample Size Calculator --\n")
  cat(sprintf("Cohort Size (N):    %d\n", x$N))
  cat(sprintf("Target Ratio:       <= %.2f (within %.0f%% of full audit width)\n",
              x$target_ratio, (x$target_ratio - 1) * 100))
  cat(sprintf("Target Power:       %.0f%%\n", x$power_target * 100))
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
    cat(sprintf("Achieves %.1f%% power (expected ratio: %.3f)\n",
                rec$power * 100, rec$expected_ratio))
  } else {
    cat(sprintf("\nWarning: Evaluated fractions do not reach %.0f%% power (max: %.1f%% at n = %d).\n",
                x$power_target * 100, rec$power * 100, rec$n_audit))
  }
  invisible(x)
}

#' Plot an audit power curve
#'
#' Displays statistical power versus audit sample size with reference threshold lines.
#'
#' @param x A `catmech_audit_power` object.
#' @param col Line color for the power curve (default `"#1f77b4"`).
#' @param lwd Line width (default 2).
#' @param ... Additional arguments passed to [graphics::plot()].
#' @return `x`, invisibly.
#' @export
plot.catmech_audit_power <- function(x, col = "#1f77b4", lwd = 2, ...) {
  tbl <- x$table
  graphics::plot(
    tbl$n_audit, tbl$power,
    type = "b", pch = 19, col = col, lwd = lwd,
    xlab = sprintf("Audit Sample Size (n out of %d)", x$N),
    ylab = sprintf("Statistical Power (P[Ratio <= %.2f])", x$target_ratio),
    ylim = c(0, 1.05),
    main = "Catmech Audit Power Curve",
    panel.first = graphics::grid(col = "gray90", lty = "dotted"),
    ...
  )

  # Power target line
  graphics::abline(h = x$power_target, col = "darkred", lty = 2, lwd = 1.5)
  graphics::mtext(sprintf("Target Power: %.0f%%", x$power_target * 100),
                  side = 4, at = x$power_target, col = "darkred", cex = 0.8, las = 1)

  # Recommended marker
  rec <- x$recommended
  if (rec$power >= x$power_target) {
    graphics::points(rec$n_audit, rec$power, col = "darkgreen", pch = 17, cex = 1.8)
    graphics::abline(v = rec$n_audit, col = "darkgreen", lty = 3)
    graphics::text(rec$n_audit, rec$power - 0.08,
                   sprintf("Recommended:\nn = %d (%.0f%% pwr)", rec$n_audit, rec$power * 100),
                   col = "darkgreen", cex = 0.85, font = 2)
  }

  invisible(x)
}
