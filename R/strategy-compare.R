#' Compare audit analysis strategies
#'
#' Fits a user-specified regression model across multiple audit analysis strategies
#' on a dataset with recorded values and an audit reference subsample:
#' \enumerate{
#'   \item \strong{Validation subsample}: Audited records only (discards un-audited recorded values).
#'   \item \strong{Naive analysis}: Trusts un-audited recorded values as correct without audit adjustment.
#'   \item \strong{Linear Multiple Imputation}: Calibrates un-audited values given recorded values ([catmech_reference(use_recorded = TRUE)]).
#'   \item \strong{Contamination Model}: Imputes un-audited values under a reference-calibrated error model ([catmech_contamination()]).
#' }
#'
#' @param data A data frame containing analysis variables, recorded values, and audit references.
#' @param value Character string naming the recorded value column.
#' @param reference Character string naming the audit reference column (with `NA` for un-audited rows).
#' @param formula Model formula for the downstream analysis (e.g. `outcome ~ lab_value + age`).
#' @param predictors Character vector of predictor column names to include in the imputation models.
#'   If `NULL`, every variable in `formula` other than `value` and `reference` is used, including
#'   the outcome. When the audited value is an analysis predictor, leaving the outcome out of the
#'   imputation model biases its coefficient toward zero; a supplied vector without the outcome
#'   draws a warning from the pooling step.
#' @param family Error distribution and link function for GLM (default [stats::gaussian()] for linear models;
#'   use [stats::binomial()] for logistic regression).
#' @param m Number of completed datasets to impute (default 10L).
#' @param seed Random seed for reproducible imputation (default 123L).
#' @param conf_level Confidence level (default 0.95).
#' @param agreement_tol Tolerance for agreement audit in contamination model (default 1.0).
#' @param cluster Optional complete cluster column. Every strategy then uses a cluster-robust
#'   variance with complete-data degrees of freedom of clusters minus one: the delete-one-cluster
#'   jackknife for linear models, as in [catmech_pool_lm()], and the CR1S sandwich for other
#'   families, as in [catmech_pool_glm()]. It is not passed to `catmech_reference(cluster = )`,
#'   whose clustered imputation biased slopes toward zero in simulations when the audited value
#'   was a predictor (`validation/VARIANCE_MI_SYNTHESIS.md`); a formula term naming the same
#'   column is still used as an ordinary predictor.
#' @param ... Currently unused.
#' @return An object of class `c("catmech_strategy_comparison", "data.frame")` containing:
#'   \item{strategy}{Strategy label (`validation_subsample`, `naive_trust_recorded`, `linear_mi_recorded`, `contamination_model`).}
#'   \item{term}{Model coefficient term.}
#'   \item{estimate}{Point estimate.}
#'   \item{std.error}{Standard error.}
#'   \item{df}{Degrees of freedom.}
#'   \item{lower, upper}{Confidence interval bounds.}
#'   \item{ci_width}{Confidence interval width (upper - lower).}
#'   \item{rel_efficiency}{Relative precision (inverse variance ratio) relative to the validation subsample.}
#'   \item{n_total, n_audited}{Total observations and audited records count.}
#'   The imputation predictors used are stored in the `imputation_predictors` attribute.
#'   If the contamination model fails, its rows are `NA` and a warning gives the reason.
#' @examples
#' set.seed(123)
#' n <- 50
#' df <- data.frame(
#'   outcome = stats::rnorm(n, 10, 2),
#'   lab = stats::rnorm(n, 5, 1),
#'   age = stats::rnorm(n, 50, 10)
#' )
#' df$ref <- df$lab + stats::rnorm(n, 0, 0.2)
#' df$ref[16:n] <- NA
#' comp <- catmech_compare_strategies(df, value = "lab", reference = "ref",
#'                                    formula = outcome ~ lab + age, m = 2, seed = 1)
#' print(comp)
#' @export
catmech_compare_strategies <- function(data, value, reference, formula,
                                       predictors = NULL,
                                       family = stats::gaussian(),
                                       m = 10L, seed = 123L, conf_level = 0.95,
                                       agreement_tol = 1.0, cluster = NULL, ...) {
  if (!is.data.frame(data)) stop("data must be a data frame.", call. = FALSE)
  if (!is.character(value) || length(value) != 1L || !value %in% names(data))
    stop("value must name an existing column in data.", call. = FALSE)
  if (!is.character(reference) || length(reference) != 1L || !reference %in% names(data))
    stop("reference must name an existing column in data.", call. = FALSE)
  if (!inherits(formula, "formula")) stop("formula must be an R formula.", call. = FALSE)

  if (is.character(family)) family <- get(family, mode = "function", envir = parent.frame())()
  if (is.function(family)) family <- family()
  if (!inherits(family, "family")) stop("family must be a valid glm family object.", call. = FALSE)
  is_gaussian <- identical(family$family, "gaussian") && identical(family$link, "identity")

  audited_idx <- which(!is.na(data[[reference]]))
  n_total <- nrow(data)
  n_audited <- length(audited_idx)

  if (n_audited < 5L)
    stop("At least 5 audited records are required to compare strategies.", call. = FALSE)

  # Detect imputation predictors if not supplied. The outcome stays in: when the audited
  # value is an analysis predictor, imputing it without the outcome attenuates its slope.
  if (is.null(predictors)) {
    rhs_vars <- all.vars(stats::delete.response(stats::terms(formula, data = data)))
    lhs_vars <- if (length(formula) == 3L) all.vars(formula[[2L]]) else character()
    cand_preds <- intersect(setdiff(unique(c(lhs_vars, rhs_vars)), c(value, reference)), names(data))
    predictors <- if (length(cand_preds) > 0L) cand_preds else NULL
  }

  ids <- NULL
  if (!is.null(cluster)) {
    if (!is.character(cluster) || length(cluster) != 1L || is.na(cluster) || !cluster %in% names(data))
      stop("cluster must name a column in data.", call. = FALSE)
    ids <- data[[cluster]]
    if (anyNA(ids) || any(!nzchar(trimws(as.character(ids)))))
      stop("cluster must not contain missing or blank values.", call. = FALSE)
  }
  # Complete-data rows: a cluster-robust variance on the same footing as the pooled strategies.
  direct_rows <- function(fit, rows, label) {
    cf <- stats::coef(summary(fit))
    se <- cf[, 2L]
    df <- if (is_gaussian) fit$df.residual else Inf
    if (!is.null(ids)) {
      used <- if (is.null(fit$na.action)) rows else rows[-as.integer(fit$na.action)]
      g <- ids[used]
      if (length(unique(g)) < 2L)
        stop("A cluster-robust variance needs at least two clusters in every strategy.", call. = FALSE)
      v <- if (is_gaussian) catmech_cluster_jackknife(fit, g) else catmech_glm_cluster_vcov(fit, g)
      se <- sqrt(diag(v))[rownames(cf)]
      df <- length(unique(g)) - 1
      if (is_gaussian) {
        alone <- catmech_cluster_alone(fit, g)
        if (length(alone))
          warning("In the ", label, " strategy a single cluster determines these coefficients, so their cluster-robust intervals can be far too narrow: ",
                  paste(alone, collapse = ", "), ".", call. = FALSE)
      }
    }
    critical <- if (is.finite(df)) stats::qt((1 + conf_level) / 2, df) else stats::qnorm((1 + conf_level) / 2)
    data.frame(strategy = label, term = rownames(cf), estimate = cf[, 1L], std.error = se,
      df = df, lower = cf[, 1L] - critical * se, upper = cf[, 1L] + critical * se,
      stringsAsFactors = FALSE)
  }

  results <- list()

  # Strategy 1: Validation subsample only (audited records)
  d_sub <- data[audited_idx, , drop = FALSE]
  d_sub[[value]] <- d_sub[[reference]]
  fit_sub <- if (is_gaussian) {
    stats::lm(formula, data = d_sub, na.action = stats::na.fail)
  } else {
    stats::glm(formula, data = d_sub, family = family, na.action = stats::na.fail)
  }
  sub_df <- direct_rows(fit_sub, audited_idx, "validation_subsample")
  terms_sub <- sub_df$term
  results[["validation_subsample"]] <- sub_df

  # Strategy 2: Naive analysis (trust recorded values)
  d_naive <- data
  fit_naive <- if (is_gaussian) {
    stats::lm(formula, data = d_naive, na.action = stats::na.omit)
  } else {
    stats::glm(formula, data = d_naive, family = family, na.action = stats::na.omit)
  }
  results[["naive_trust_recorded"]] <- direct_rows(fit_naive, seq_len(n_total), "naive_trust_recorded")

  # Strategy 3: Linear MI given recorded value
  res_linear <- catmech_reference(
    data = data, value = value, reference = reference,
    predictors = predictors, use_recorded = TRUE,
    m = m, seed = seed
  )
  pooled_linear <- if (is_gaussian) {
    catmech_pool_lm(res_linear, formula, conf_level = conf_level, cluster = cluster)
  } else {
    catmech_pool_glm(res_linear, formula, family = family, conf_level = conf_level, cluster = cluster)
  }
  linear_df <- data.frame(
    strategy = "linear_mi_recorded",
    term = pooled_linear$term,
    estimate = pooled_linear$estimate,
    std.error = pooled_linear$std.error,
    df = pooled_linear$df,
    lower = pooled_linear$lower,
    upper = pooled_linear$upper,
    stringsAsFactors = FALSE
  )
  results[["linear_mi_recorded"]] <- linear_df

  # Strategy 4: Contamination model
  contam_df <- tryCatch({
    res_contam <- catmech_contamination(
      data = data, value = value, reference = reference,
      predictors = predictors, m = m, seed = seed,
      agreement_tol = agreement_tol
    )
    pooled_contam <- if (is_gaussian) {
      catmech_pool_lm(res_contam, formula, conf_level = conf_level, cluster = cluster)
    } else {
      catmech_pool_glm(res_contam, formula, family = family, conf_level = conf_level, cluster = cluster)
    }
    data.frame(
      strategy = "contamination_model",
      term = pooled_contam$term,
      estimate = pooled_contam$estimate,
      std.error = pooled_contam$std.error,
      df = pooled_contam$df,
      lower = pooled_contam$lower,
      upper = pooled_contam$upper,
      stringsAsFactors = FALSE
    )
  }, error = function(e) {
    warning("The contamination model failed, so its rows are NA: ", conditionMessage(e), call. = FALSE)
    data.frame(
      strategy = "contamination_model",
      term = terms_sub,
      estimate = rep(NA_real_, length(terms_sub)),
      std.error = rep(NA_real_, length(terms_sub)),
      df = rep(NA_real_, length(terms_sub)),
      lower = rep(NA_real_, length(terms_sub)),
      upper = rep(NA_real_, length(terms_sub)),
      stringsAsFactors = FALSE
    )
  })
  results[["contamination_model"]] <- contam_df

  # Combine and calculate efficiency
  comb <- do.call(rbind, results)
  comb$ci_width <- comb$upper - comb$lower
  comb$n_total <- n_total
  comb$n_audited <- n_audited

  # Relative efficiency: variance of validation subsample / variance of strategy
  # (Values > 1 indicate precision gains over discarding un-audited rows)
  se_sub_map <- stats::setNames(sub_df$std.error, sub_df$term)
  comb$rel_efficiency <- (se_sub_map[comb$term] / comb$std.error)^2

  rownames(comb) <- NULL
  structure(comb, class = c("catmech_strategy_comparison", "data.frame"),
            imputation_predictors = if (is.null(predictors)) character() else predictors,
            cluster = cluster)
}

#' Print strategy comparison results
#'
#' @param x A `catmech_strategy_comparison` object.
#' @param digits Number of decimal digits to print (default 4).
#' @param ... Unused.
#' @return `x`, invisibly.
#' @export
print.catmech_strategy_comparison <- function(x, digits = 4L, ...) {
  n_tot <- if (!is.null(x$n_total[1L])) x$n_total[1L] else NA_integer_
  n_aud <- if (!is.null(x$n_audited[1L])) x$n_audited[1L] else NA_integer_
  cat("-- Catmech Audit Strategy Comparison --\n")
  cat(sprintf("Cohort Size (N): %d | Audited Records (n): %d (%.1f%% audit fraction)\n",
              n_tot, n_aud, 100 * n_aud / n_tot))
  preds <- attr(x, "imputation_predictors")
  if (!is.null(preds))
    cat("Imputation predictors: ", if (length(preds)) paste(preds, collapse = ", ") else "none", "\n", sep = "")
  if (!is.null(attr(x, "cluster")))
    cat("Cluster-robust variances by: ", attr(x, "cluster"), "\n", sep = "")
  cat("\n")

  tbl <- as.data.frame(x)
  num_cols <- c("estimate", "std.error", "df", "lower", "upper", "ci_width", "rel_efficiency")
  for (col in num_cols) {
    if (col %in% names(tbl)) {
      tbl[[col]] <- ifelse(is.finite(tbl[[col]]), round(tbl[[col]], digits), tbl[[col]])
    }
  }
  disp_cols <- c("strategy", "term", "estimate", "std.error", "ci_width", "rel_efficiency")
  print(tbl[, disp_cols, drop = FALSE], row.names = FALSE)
  invisible(x)
}

#' Plot audit strategy comparison (Forest Plot)
#'
#' Draws a comparative forest plot showing point estimates and confidence intervals
#' across all evaluated audit strategies for a chosen coefficient term.
#'
#' @param x A `catmech_strategy_comparison` object.
#' @param term Character: which coefficient term to plot. If `NULL`, the first non-intercept
#'   term is plotted.
#' @param colors Named vector of colors for the strategies.
#' @param lwd Line width (default 2).
#' @param pch Point symbol (default 19).
#' @param cex Point size (default 1.3).
#' @param main Plot title.
#' @param xlab X-axis label.
#' @param ... Additional arguments passed to [graphics::plot()].
#' @return `x`, invisibly.
#' @examples
#' set.seed(123)
#' n <- 50
#' df <- data.frame(
#'   outcome = stats::rnorm(n, 10, 2),
#'   lab = stats::rnorm(n, 5, 1),
#'   age = stats::rnorm(n, 50, 10)
#' )
#' df$ref <- df$lab + stats::rnorm(n, 0, 0.2)
#' df$ref[16:n] <- NA
#' comp <- catmech_compare_strategies(
#'   df, value = "lab", reference = "ref",
#'   formula = outcome ~ lab + age, m = 2, seed = 1
#' )
#' plot(comp, term = "lab")
#' @export
plot.catmech_strategy_comparison <- function(x, term = NULL,
                                            colors = c(
                                              validation_subsample = "#898781",
                                              naive_trust_recorded = "#d03b3b",
                                              linear_mi_recorded = "#2a78d6",
                                              contamination_model = "#eb6834"
                                            ),
                                            lwd = 2, pch = 19, cex = 1.3,
                                            main = NULL, xlab = NULL, ...) {
  tbl <- as.data.frame(x)
  avail_terms <- unique(tbl$term)
  if (is.null(term)) {
    non_int <- setdiff(avail_terms, c("(Intercept)", "Intercept"))
    term <- if (length(non_int) > 0L) non_int[1L] else avail_terms[1L]
  }
  tbl_sub <- tbl[tbl$term == term, , drop = FALSE]
  if (nrow(tbl_sub) == 0L) stop(sprintf("Term '%s' not found in comparison results.", term), call. = FALSE)

  k <- nrow(tbl_sub)
  y_pos <- rev(seq_len(k))

  xlim <- range(c(tbl_sub$lower, tbl_sub$upper, 0), na.rm = TRUE)
  x_padding <- max(diff(xlim) * 0.15, 0.05)
  xlim <- c(xlim[1L] - x_padding, xlim[2L] + x_padding)

  strat_labels <- c(
    validation_subsample = "Validation Subsample (Discard un-audited)",
    naive_trust_recorded = "Naive Analysis (Trust recorded)",
    linear_mi_recorded = "Linear MI (use_recorded = TRUE)",
    contamination_model = "Contamination Model (Catmech)"
  )
  display_names <- ifelse(tbl_sub$strategy %in% names(strat_labels),
                          strat_labels[tbl_sub$strategy], tbl_sub$strategy)

  strat_cols <- ifelse(tbl_sub$strategy %in% names(colors),
                       colors[tbl_sub$strategy], "#1f77b4")

  if (is.null(xlab)) xlab <- sprintf("Estimate: %s (95%% CI)", term)
  if (is.null(main)) main <- sprintf("Audit Strategy Comparison: %s", term)

  old_mar <- graphics::par("mar")
  on.exit(graphics::par(mar = old_mar))
  graphics::par(mar = c(5, 14, 4, 7))

  graphics::plot(
    tbl_sub$estimate, y_pos,
    type = "n",
    xlim = xlim, ylim = c(0.5, k + 0.5),
    xlab = xlab, ylab = "",
    yaxt = "n", main = main,
    panel.first = graphics::grid(col = "gray90", lty = "dotted"),
    ...
  )

  graphics::abline(v = 0, lty = 2, col = "gray50", lwd = 1.2)

  for (i in seq_len(k)) {
    if (is.finite(tbl_sub$lower[i]) && is.finite(tbl_sub$upper[i])) {
      col_i <- strat_cols[i]
      graphics::segments(tbl_sub$lower[i], y_pos[i], tbl_sub$upper[i], y_pos[i], col = col_i, lwd = lwd)
      tick_h <- 0.12
      graphics::segments(tbl_sub$lower[i], y_pos[i] - tick_h, tbl_sub$lower[i], y_pos[i] + tick_h, col = col_i, lwd = lwd)
      graphics::segments(tbl_sub$upper[i], y_pos[i] - tick_h, tbl_sub$upper[i], y_pos[i] + tick_h, col = col_i, lwd = lwd)
    }
  }

  graphics::points(tbl_sub$estimate, y_pos, pch = pch, col = strat_cols, cex = cex)
  graphics::axis(2, at = y_pos, labels = display_names, las = 1, cex.axis = 0.8)

  annot <- ifelse(
    is.finite(tbl_sub$estimate),
    sprintf("%.3f [%.3f, %.3f]", tbl_sub$estimate, tbl_sub$lower, tbl_sub$upper),
    "NA [NA, NA]"
  )
  graphics::mtext(annot, side = 4, at = y_pos, las = 1, line = 0.5, cex = 0.8, col = "gray30")

  invisible(x)
}
