#' Audited predictive mean matching after masking quality flags
#'
#' Masks every flagged target and uses standard type-1 predictive mean matching
#' from `mice`. Categories may be selected as ordinary complete predictors;
#' they never force interval membership or cause trusted observations to change.
#'
#' @param data A data frame with unique column names. Original data are retained.
#' @param value Name of the numeric measurement column.
#' @param predictors Explicit complete predictor names, excluding the target and
#'   flag. Include downstream outcomes when appropriate. A category can be an
#'   ordinary factor predictor; no category likelihood is imposed.
#' @param flag Optional quality-flag column name.
#' @param flag_yes,flag_no Distinct suspect and non-suspect labels. Missing flags
#'   are treated as non-suspect, a substantive trust assumption.
#' @param donors Positive number of nearest donors, default five. When any
#'   imputation is needed, fewer trusted observations than `donors` is an error.
#' @param m Number of completions, at least two.
#' @param seed Nonnegative integer seed. The caller's RNG state is restored.
#' @return A `catmech_result`. The `mids` adapter stores a masked target as its
#'   observed data; the unchanged raw measurements remain in `original`.
#' @details
#' Calls `mice::mice.impute.pmm` with `matchtype = 1` and `ridge = 1e-5`, using
#' only originally trusted observations as donors. The trusted design must have
#' full rank and positive residual degrees of freedom and variance whenever
#' imputation is required. Completed values lie in the observed donor support;
#' missing extremes, nonignorable donor selection and misspecification can bias
#' inference. PMM is not a general robustness or coverage guarantee.
#'
#' The delegated method returns donor values, not donor row identifiers. The
#' audit records `donor_row` only when a value matches exactly one trusted row.
#' For ties, `donor_row` is `NA`, and `donor_value_rows` lists all value-compatible
#' trusted rows. This list does not identify the actual nearest-neighbor pool.
#' `eligible_donors` counts all trusted rows; `candidate_pool_size` records the
#' requested nearest-donor count used by the delegated matching algorithm.
#' @export
catmech_pmm <- function(data, value, predictors, flag = NULL,
                        flag_yes = "Yes", flag_no = "No", donors = 5L,
                        m = 20L, seed = 1L) {
  catmech_reference_inputs(data, value, predictors, m, seed, donors,
                           excluded = c(value, flag))
  scalar_name <- function(z) is.character(z) && length(z) == 1L &&
    !is.na(z) && nzchar(z)
  if ((!is.null(flag) && (!scalar_name(flag) || !flag %in% names(data))) ||
      !scalar_name(flag_yes) || !scalar_name(flag_no) || flag_yes == flag_no)
    stop("flag must name an existing column and flag labels must be distinct nonempty strings.", call. = FALSE)
  if (!is.null(flag) && flag == value)
    stop("The flag and target must be different columns.", call. = FALSE)
  f <- if (is.null(flag)) rep(flag_no, nrow(data)) else as.character(data[[flag]])
  if (any(!is.na(f) & !f %in% c(flag_no, flag_yes)))
    stop("Unknown flag labels in the data.", call. = FALSE)
  masked <- is.na(data[[value]]) | catmech_suspect_mask(f, flag_yes)
  observed <- data[[value]]
  observed[masked] <- NA_real_
  X <- if (any(masked)) catmech_model_matrix(data, predictors) else NULL
  restore <- preserve_rng()
  on.exit(restore(), add = TRUE)
  set.seed(as.integer(seed))
  draws <- catmech_reference_pmm_draws(observed, X, donors, m)
  spec <- structure(list(version = "pmm1", engine = "trusted_pmm", value = value,
    predictors = predictors,
    flag = flag, flag_yes = flag_yes, flag_no = flag_no, flag_policy = "replace",
    missing_flag_policy = "trust_observed", donors = as.integer(donors),
    matchtype = 1L, ridge = 1e-5, category_policy = "ordinary_predictor_only_if_selected",
    missingness_assumption = "ignorable_trusted_selection_given_predictors"),
    class = "catmech_pmm_spec")
  selected <- which(masked)
  audits <- lapply(seq_len(m), function(j) {
    row <- catmech_reference_audit(data[[value]], draws[, j], selected, j,
      ifelse(is.na(data[[value]][selected]), "missing", "flag_masked"),
      reconstructed = rep(TRUE, length(selected)))
    catmech_reference_donor_audit(row, observed, draws[selected, j], donors)
  })
  diagnostics <- do.call(rbind, lapply(seq_len(m), function(j) data.frame(
    imputation = j, requested = length(selected), reconstructed = length(selected),
    changed = sum(audits[[j]]$changed, na.rm = TRUE),
    trusted_observed = sum(!masked), fallback_draws = 0L)))
  result <- catmech_reference_assemble(data, value, predictors, observed, draws, audits,
    diagnostics, spec, m, seed, "conditional_predictive_matching",
    catmech_reference_independent_diagnostics("mice_pmm"))
  result$mids$call <- match.call()
  result
}

#' Complete a target using an independently measured reference subset
#'
#' Treats the reference column as the only observed truth for the target. Every
#' non-reference target is unavailable, regardless of its recorded value or flags.
#' This changes the trust assumption; it does not infer a gold standard from a
#' noisy measurement or correct missed errors without external information.
#'
#' @param data A data frame; raw target measurements remain in `original`.
#' @param value Name of the numeric target to replace in completed datasets.
#' @param reference Name of a numeric gold-standard column, observed only where
#'   independently measured. Known reference values are preserved exactly in the
#'   completed target, including when they disagree with recorded measurements.
#' @param predictors Explicit complete predictor names, excluding `value` and
#'   `reference`. Include the scientific outcome and relevant design variables.
#'   A raw measurement may be used only through a separate explicitly selected,
#'   complete auxiliary column. No missing proxy is filled automatically.
#' @param method `"normal"` (default) for conditional Gaussian posterior draws,
#'   `"pmm"` for standard predictive mean matching using reference donors, or
#'   `"pmm_abb"` for the same matching on a fresh bootstrap of the donors at each
#'   imputation. See details.
#' @param donors Nearest-donor count for matching, default five; ignored by the
#'   normal method. Fewer reference observations than donors is an error when needed.
#' @param m Number of completions, at least two.
#' @param seed Nonnegative integer seed, with caller RNG restored.
#' @param use_recorded For `method = "normal"`, also use the recorded value as a linear
#'   predictor. Records with a recorded value are imputed from a model fitted on
#'   references that also have one, and records without one from a model on the
#'   predictors fitted on all references. It needs rows with both values. If every such row
#'   agrees exactly, those records keep their recorded value with no imputation uncertainty,
#'   as `spec$recorded_fit` notes and a warning says. Intervals then treat those records as
#'   correct, which can be badly wrong when a small audit simply missed the errors.
#' @param cluster With `use_recorded = TRUE`, the name of a complete column identifying
#'   audit clusters, such as a site or facility. The imputation model then adds a random
#'   intercept for the cluster, instead of treating rows as independent. If the audit is
#'   stratified, put the strata in `predictors`, which needs some stratum with at least two
#'   audited clusters. Do not use it when the audited value will be a predictor in the
#'   analysis. See details.
#' @return A `catmech_result`. Its completed target contains exact gold references
#'   or imputations; the original reference column is unchanged. The `mids`
#'   observed target is the reference values with other rows missing. The manifest
#'   records `imputed_rows` explicitly, including an empty vector for full reference.
#' @details
#' Reference selection must be ignorable given the selected observed predictors,
#' with adequate support in the reference subset. Reference values must measure
#' the same quantity, time point and units, and are assumed exact. Unknown
#' reference error, unmeasured selection and transport failure are not repaired.
#' This is standard missing-data imputation using externally supplied validation
#' information, not automatic detection of errors in unvalidated records.
#'
#' With partial reference data, normal imputation delegates to `catmech_model`
#' with categories ignored and without a quality-flag policy. PMM delegates to
#' the same audited standard matching engine as `catmech_pmm`. Both require a
#' full-rank reference design with positive residual degrees of freedom and
#' variance. Complete reference data need no stochastic fit or donor support.
#'
#' Partial-reference PMM is experimental sensitivity output. A held-out package
#' benchmark found substantial interval undercoverage even with a randomly
#' selected, accurate reference subset. As a conservative package policy,
#' `catmech_pool_lm` requires `allow_sensitivity = TRUE` whenever reference PMM
#' imputes any target cells, and labels its calculations `sensitivity_only`.
#' This policy is not a general theorem about PMM, and the opt-in does not
#' establish inferential validity. `catmech_sensitivity` already supplies an
#' explicit sensitivity context. Complete reference data remain eligible for
#' ordinary complete-data pooling because no target cells are imputed.
#'
#' The likely cause of that undercoverage is that standard matching conditions on
#' one fixed set of donor values, which attenuates between-imputation variance
#' (Gaffert, Meinfelder and Bosch). `"pmm_abb"` resamples donors with replacement
#' for each imputation, an approximate Bayesian bootstrap (Rubin and Schenker, 1986;
#' Siddique and Belin, 2008), before the same delegated type-1 matching. It targets
#' that variance and is not gated. It still cannot leave donor support, and it does
#' not make matching more precise than the normal method when the normal model holds.
#'
#' With `use_recorded = TRUE`, reference selection may depend on the recorded value, and
#' linearity in the recorded value becomes an assumption. It was best or tied among
#' reference-using methods at 25% audits on two real datasets, and failed under synthetic
#' gross errors, which distort a linear fit;
#' [catmech_contamination()] is the alternative there.
#'
#' With `cluster`, the imputation model is the recorded-value regression plus a random
#' intercept for the cluster. Its variance components come from moments (Henderson, 1953,
#' method 3) on the audited rows and are drawn in each imputation, and coefficients are
#' drawn from the generalized least squares fit given them. An unaudited cluster gets a new
#' intercept, an audited cluster with unaudited rows gets one conditional on its residuals,
#' and every imputed row gets normal noise. If the audited rows cannot separate the two
#' variances, as with one row per cluster, the intercept variance is zero, with a warning. If
#' the predictors absorb every audited cluster, as with one audited cluster per stratum, there
#' is no estimate of between-cluster variance for unaudited clusters, and the fit is refused.
#' Fewer than five audited clusters are refused. Pool with [catmech_pool_lm()], which uses a
#' cluster-robust variance for such a fit. This is the usual random-intercept imputation model
#' (Schafer and Yucel, 2002; Grund and colleagues, 2018), with a moment-based draw of the
#' variance components in place of a full posterior. That draw is an approximation;
#' balanced clusters alone do not make it an exact posterior draw with general predictors.
#' Its justification is weaker with one or two large sites. It is not a new method. It
#' assumes one normal intercept per cluster, so cluster-specific slopes and non-constant
#' variance are not modeled. The same model is available in other software, including `mice`.
#'
#' Within a cluster, unaudited rows cannot require stochastic imputation through both
#' the recorded-present and recorded-missing paths: the current separate fits do not
#' provide a shared cluster draw across those paths, so that request is refused.
#' Different clusters may use different paths. Audited rows and recorded values kept
#' deterministically after complete calibration agreement do not trigger this refusal.
#' Resolve one pattern with reference measurements or use a joint model supporting both.
#'
#' `cluster` is experimental. Version 0.6.0 replaced an earlier version that resampled audited
#' clusters, in which resampling within strata understated variance. The replacement's planned
#' test came out PARTIAL under a rule fixed in advance: on a real facility recount at a 25%
#' audit, coverage was within its limits on eight of nine items, with 89.0% on the ninth, but
#' error on the three rare items ran 4% to 5% above the option without `cluster`, and the
#' accuracy conditions failed. Treat `cluster` as a check, not a default. Read the `df` column,
#' since very wide intervals remain possible with few audited clusters. See
#' `validation/CLUSTERED_FIX_SYNTHESIS.md`.
#'
#' In later simulations where the audited value was a predictor in the analysis and errors
#' shifted whole clusters, `cluster` biased that slope toward zero and covered far less than
#' a fit without it. The cause has not been established. [catmech_pool_lm()] warns about
#' this use. See `validation/VARIANCE_MI_SYNTHESIS.md`.
#' @export
catmech_reference <- function(data, value, reference, predictors,
                              method = c("normal", "pmm", "pmm_abb"), donors = 5L,
                              m = 20L, seed = 1L, use_recorded = FALSE,
                              cluster = NULL) {
  method <- match.arg(method)
  if (!is.logical(use_recorded) || length(use_recorded) != 1L || is.na(use_recorded))
    stop("use_recorded must be TRUE or FALSE.", call. = FALSE)
  if (use_recorded && method != "normal")
    stop("use_recorded applies only to method = 'normal'.", call. = FALSE)
  catmech_reference_inputs(data, value, predictors, m, seed, donors,
                           excluded = c(value, reference))
  if (!is.null(cluster)) {
    if (!use_recorded) stop("cluster applies only with use_recorded = TRUE.", call. = FALSE)
    catmech_reference_check_cluster(data, cluster, c(value, reference, predictors))
  }
  if (!is.character(reference) || length(reference) != 1L || is.na(reference) ||
      !nzchar(reference) || !reference %in% names(data) || reference == value)
    stop("reference must name an existing column different from the target.", call. = FALSE)
  gold <- data[[reference]]
  if (!is.numeric(gold) || any(!is.na(gold) & !is.finite(gold)))
    stop("reference must contain independently measured finite numeric values or NA.", call. = FALSE)
  known <- !is.na(gold)
  if (!any(known)) stop("At least one independently observed reference value is required.", call. = FALSE)
  X <- if (all(known) || use_recorded) NULL else catmech_model_matrix(data, predictors)
  recorded_fit <- NULL
  cluster_uncertainty <- "random_intercept_drawn_each_imputation"
  restore <- preserve_rng()
  on.exit(restore(), add = TRUE)
  set.seed(as.integer(seed))
  if (all(known)) {
    draws <- matrix(rep(gold, m), nrow = nrow(data), ncol = m)
    model_diagnostics <- catmech_reference_independent_diagnostics("reference_complete")
  } else if (use_recorded) {
    sampled <- if (is.null(cluster)) catmech_reference_recorded_draws(data, value, gold, predictors, m, seed) else
      catmech_reference_clustered_draws(data, value, gold, predictors, m, cluster)
    draws <- sampled$draws
    recorded_fit <- sampled$note
    if (!is.null(sampled$cluster_uncertainty)) cluster_uncertainty <- sampled$cluster_uncertainty
    if (identical(recorded_fit, "recorded_kept_all_calibration_rows_agree")) {
      checked <- sum(known & !is.na(data[[value]]))
      warning(sprintf(paste0("All %d audited rows agree with their recorded values, so %d unaudited recorded values were kept ",
        "with no imputation uncertainty. An audit of this size cannot rule out an error rate of about %.1f%%, ",
        "and intervals here treat those records as correct."),
        checked, sum(!known & !is.na(data[[value]])), 100 * min(1, 3 / checked)), call. = FALSE)
    }
    model_diagnostics <- catmech_reference_independent_diagnostics("reference_normal_recorded")
  } else if (method == "normal") {
    work <- data
    work[[value]] <- gold
    fitted <- catmech_model(work, value = value, predictors = predictors,
                            category_mode = "ignore", m = m, seed = seed)
    draws <- vapply(fitted$completed, function(z) z[[value]], numeric(nrow(data)))
    model_diagnostics <- fitted$model_diagnostics
  } else {
    draws <- catmech_reference_pmm_draws(gold, X, donors, m, bootstrap = method == "pmm_abb")
    model_diagnostics <- catmech_reference_independent_diagnostics(paste0("reference_mice_", method))
  }
  engine <- if (use_recorded) "reference_normal_recorded" else paste0("reference_", method)
  spec <- structure(list(version = "reference1", engine = engine, value = value,
    predictors = predictors,
    reference = reference, reference_policy = "observed_gold_only_all_others_unavailable",
    raw_target_policy = "never_trusted_automatically",
    reference_selection_assumption = "ignorable_given_complete_selected_predictors",
    reference_measurement_assumption = "exact_same_quantity_time_and_units",
    method = method, donors = if (method != "normal") as.integer(donors) else NULL,
    matchtype = if (method != "normal") 1L else NULL,
    ridge = if (method != "normal") 1e-5 else NULL,
    prior = if (method == "normal") "p(beta,sigma_squared) proportional to 1/sigma_squared" else NULL),
    class = "catmech_reference_spec")
  if (method == "pmm_abb") spec$donor_resampling <- "approximate_bayesian_bootstrap_each_imputation"
  if (use_recorded) {
    spec$use_recorded <- TRUE
    spec$raw_target_policy <- "linear_predictor_never_trusted_automatically"
    spec$reference_selection_assumption <- "ignorable_given_complete_selected_predictors_and_recorded_value"
    spec$recorded_fit <- recorded_fit
  }
  if (!is.null(cluster)) {
    spec$cluster <- cluster
    spec$cluster_uncertainty <- cluster_uncertainty
    spec$audited_clusters <- length(unique(data[[cluster]][known]))
  }
  audits <- lapply(seq_len(m), function(j) {
    row <- catmech_reference_audit(data[[value]], draws[, j], seq_len(nrow(data)), j,
      ifelse(known, "gold_reference", "reference_missing"), reconstructed = !known)
    row$reference_observed <- known
    if (method != "normal") {
      info <- catmech_reference_donor_audit(row[!known, , drop = FALSE], gold, draws[!known, j], donors)
      row$donor_row <- NA_integer_
      row$eligible_donors <- ifelse(known, NA_integer_, sum(known))
      row$donor_value_match_count <- NA_integer_
      row$donor_value_rows <- NA_character_
      row$donor_row_identified <- FALSE
      row$candidate_pool_size <- NA_integer_
      fields <- c("donor_row", "eligible_donors", "donor_value_match_count",
                  "donor_value_rows", "donor_row_identified", "candidate_pool_size")
      row[!known, fields] <- info[fields]
    }
    row
  })
  diagnostics <- do.call(rbind, lapply(seq_len(m), function(j) data.frame(
    imputation = j, requested = sum(!known), reconstructed = sum(!known),
    changed = sum(audits[[j]]$changed, na.rm = TRUE),
    reference_observed = sum(known),
    reference_numeric_replacements = sum(known & !is.na(data[[value]]) & gold != data[[value]], na.rm = TRUE),
    reference_fills = sum(known & is.na(data[[value]])), fallback_draws = 0L)))
  result <- catmech_reference_assemble(data, value, predictors, gold, draws, audits,
    diagnostics, spec, m, seed,
    if (method == "normal") "conditional_reference_model" else if (method == "pmm" && any(!known))
      "sensitivity_only" else "conditional_reference_matching",
    model_diagnostics)
  result$manifest$reference_rows <- which(known)
  result$manifest$reference_column <- reference
  result$mids$call <- match.call()
  result
}

catmech_reference_inputs <- function(data, value, predictors, m, seed, donors, excluded) {
  if (!is.data.frame(data) || nrow(data) < 2L || anyNA(names(data)) ||
      anyDuplicated(names(data)) || any(!nzchar(names(data))))
    stop("data must have at least two rows and unique nonempty column names.", call. = FALSE)
  if (any(c(".imp", ".id") %in% names(data)))
    stop(".imp and .id are reserved column names.", call. = FALSE)
  if (any(!vapply(data, function(z) is.atomic(z) && is.null(dim(z)), logical(1))))
    stop("List and matrix columns are not supported.", call. = FALSE)
  if (!is.character(value) || length(value) != 1L || is.na(value) || !nzchar(value) ||
      !value %in% names(data)) stop("value must name an existing numeric target.", call. = FALSE)
  if (!is.numeric(data[[value]]) || any(!is.na(data[[value]]) & !is.finite(data[[value]])))
    stop("The target must contain finite numeric values or NA.", call. = FALSE)
  if (!is.character(predictors) || anyNA(predictors) || anyDuplicated(predictors) ||
      any(!nzchar(predictors)) || any(!predictors %in% names(data)))
    stop("predictors must be explicit, unique existing column names.", call. = FALSE)
  if (any(predictors %in% excluded))
    stop("predictors must exclude the target, reference or flag role columns.", call. = FALSE)
  check_count(m, "m", 2L)
  check_count(seed, "seed", 0L)
  check_count(donors, "donors", 1L)
  px <- data[predictors]
  if (anyNA(px) || any(!vapply(px, function(z) {
    is.numeric(z) || is.factor(z) || is.character(z) || is.logical(z)
  }, logical(1))))
    stop("Predictors must be complete numeric, logical, character or factor columns.", call. = FALSE)
  if (any(vapply(px, function(z) is.numeric(z) && any(!is.finite(z)), logical(1))))
    stop("Predictors must be finite.", call. = FALSE)
  invisible(TRUE)
}

catmech_reference_pmm_draws <- function(observed, X, donors, m, bootstrap = FALSE) {
  trusted <- !is.na(observed)
  masked <- !trusted
  n <- length(observed)
  if (!any(masked)) return(matrix(rep(observed, m), nrow = n, ncol = m))
  if (sum(trusted) < donors)
    stop("Fewer trusted/reference observations than requested donors; no silent donor-count reduction is allowed.", call. = FALSE)
  design <- catmech_model_design(X[trusted, , drop = FALSE])
  catmech_model_fit(design, observed[trusted])
  draws <- matrix(NA_real_, n, m)
  x <- X[, -1L, drop = FALSE]
  for (j in seq_len(m)) {
    out <- observed
    out[masked] <- if (bootstrap) catmech_reference_abb_draw(observed, trusted, x, masked, donors) else
      mice::mice.impute.pmm(observed, trusted, x, wy = masked, donors = donors, matchtype = 1L, ridge = 1e-5)
    if (any(!is.finite(out)) || any(!out[masked] %in% observed[trusted]))
      stop("PMM produced a nonfinite value or a value outside trusted donor support.", call. = FALSE)
    draws[, j] <- out
  }
  draws
}

catmech_reference_recorded_draws <- function(data, value, gold, predictors, m, seed) {
  W <- data[[value]]
  known <- !is.na(gold)
  recorded <- !is.na(W)
  calibration <- known & recorded
  if (!any(calibration))
    stop("use_recorded requires rows with both recorded and reference values.", call. = FALSE)
  draws <- matrix(rep(gold, m), nrow = nrow(data), ncol = m)
  note <- "linear_given_recorded"
  complete_with <- function(rows, extra, fit_seed) {
    work <- data[rows, predictors, drop = FALSE]
    work[[value]] <- gold[rows]
    names_all <- predictors
    if (!is.null(extra)) {
      auxiliary <- ".catmech_recorded"
      while (auxiliary %in% c(names(data), value)) auxiliary <- paste0(auxiliary, "_")
      work[[auxiliary]] <- extra[rows]
      names_all <- c(predictors, auxiliary)
    }
    fit <- catmech_model(work, value, names_all, category_mode = "ignore", m = m, seed = fit_seed)
    vapply(fit$completed, function(z) z[[value]], numeric(sum(rows)))
  }
  if (any(recorded & !known)) {
    if (all(catmech_contamination_agreement(W[calibration], gold[calibration], 0))) {
      draws[recorded & !known, ] <- W[recorded & !known]
      note <- "recorded_kept_all_calibration_rows_agree"
    } else {
      draws[recorded, ] <- complete_with(recorded, W, as.integer(seed))
    }
  }
  if (any(!recorded & !known)) {
    rows <- known | !recorded
    completed <- complete_with(rows, NULL, as.integer(seed) + 1L)
    draws[!recorded & !known, ] <- completed[!recorded[rows] & !known[rows], , drop = FALSE]
  }
  list(draws = draws, note = note)
}

catmech_reference_check_cluster <- function(data, cluster, excluded) {
  if (!is.character(cluster) || length(cluster) != 1L || is.na(cluster) || !cluster %in% names(data))
    stop("cluster must name an existing column.", call. = FALSE)
  if (cluster %in% excluded)
    stop("cluster must differ from the target, reference and predictor columns.", call. = FALSE)
  if (anyNA(data[[cluster]])) stop("cluster must not contain missing values.", call. = FALSE)
  if (any(!nzchar(trimws(as.character(data[[cluster]])))))
    stop("cluster must not contain blank values.", call. = FALSE)
  invisible(TRUE)
}

# Moment (Henderson III) variance components for a random intercept, from the pooled fit
# and the fit with cluster indicators. They separate only with residual degrees of freedom
# both within and between clusters.
catmech_reference_cluster_moments <- function(X, Y, cl) {
  n <- nrow(X); p <- ncol(X)
  pooled <- stats::lm.fit(X, Y)
  sse_pooled <- sum(pooled$residuals^2)
  if (!is.finite(sse_pooled) || sse_pooled <= .Machine$double.eps^2 * max(1, sum(Y^2)) * length(Y))
    stop("The audited rows leave no residual variance.", call. = FALSE)
  D <- stats::model.matrix(~ cl - 1)
  within <- stats::lm.fit(cbind(X, D), Y)
  df_within <- n - within$rank; df_between <- within$rank - p
  sse_within <- sum(within$residuals^2)
  separable <- df_within > 0L && df_between > 0L && sse_within > 0
  list(sse_pooled = sse_pooled, df_pooled = n - p, sse_within = sse_within, df_within = df_within,
       ss_between = max(0, sse_pooled - sse_within), df_between = df_between, separable = separable,
       trace = if (separable) n - sum(crossprod(qr.Q(pooled$qr), D)^2) else NA_real_)
}

# Generalized least squares given the variance components, by quasi-demeaning within clusters.
catmech_reference_cluster_gls <- function(X, Y, cl, sigma2, tau2) {
  sizes <- tabulate(cl)
  lambda <- 1 - sqrt(sigma2 / (sigma2 + sizes * tau2))
  xbar <- rowsum(X, cl, reorder = TRUE) / sizes
  ybar <- as.vector(rowsum(Y, cl, reorder = TRUE)) / sizes
  stats::lm.fit(X - (lambda * xbar)[as.integer(cl), , drop = FALSE], Y - (lambda * ybar)[as.integer(cl)])
}

# One draw per imputation: variance components, then coefficients given them, then a
# cluster intercept and row noise for each imputed row.
catmech_reference_cluster_fill <- function(x, y, fit_rows, fill_rows, ids, m) {
  if (qr(x[fit_rows, , drop = FALSE], tol = 1e-7)$rank != ncol(x) || length(fit_rows) <= ncol(x))
    stop("The trusted regression design must have full column rank; remove aliased or unsupported predictors.", call. = FALSE)
  if (length(unique(ids[fit_rows])) < 5L)
    stop("Fewer than five audited clusters: cluster uncertainty is refused.", call. = FALSE)
  X <- x[fit_rows, , drop = FALSE]; Y <- y[fit_rows]; p <- ncol(X)
  cl <- factor(as.character(ids[fit_rows])); sizes <- tabulate(cl)
  moments <- catmech_reference_cluster_moments(X, Y, cl)
  fill_clusters <- as.character(ids[fill_rows]); fill_units <- unique(fill_clusters)
  if (moments$df_between == 0L && moments$df_within > 0L && any(!fill_units %in% levels(cl)))
    stop("The predictors absorb every audited cluster, as with one audited cluster per stratum, so between-cluster variance cannot be estimated for unaudited clusters. Audit more clusters per stratum, or remove the strata from predictors.", call. = FALSE)
  if (!moments$separable)
    warning("The audited rows cannot separate within-cluster and between-cluster variance, as with one row per cluster, so the intercept variance is zero and this fit has no cluster effect.", call. = FALSE)
  out <- matrix(NA_real_, length(fill_rows), m)
  for (j in seq_len(m)) {
    draw <- NULL
    for (attempt in seq_len(100L)) {
      if (moments$separable) {
        sigma2 <- moments$sse_within / stats::rchisq(1L, moments$df_within)
        tau2 <- max(0, moments$df_between * (moments$ss_between / stats::rchisq(1L, moments$df_between) - sigma2) / moments$trace)
      } else {
        sigma2 <- moments$sse_pooled / stats::rchisq(1L, moments$df_pooled)
        tau2 <- 0
      }
      gls <- catmech_reference_cluster_gls(X, Y, cl, sigma2, tau2)
      if (gls$rank != p) next
      beta <- gls$coefficients + as.vector(backsolve(qr.R(gls$qr), stats::rnorm(p))) * sqrt(sigma2)
      effect <- stats::rnorm(length(fill_units), 0, sqrt(tau2)); names(effect) <- fill_units
      seen <- fill_units[fill_units %in% levels(cl)]
      if (length(seen) && tau2 > 0) {
        k <- sizes[match(seen, levels(cl))]
        residual_mean <- as.vector(rowsum(Y - as.vector(X %*% beta), cl, reorder = TRUE)) / sizes
        effect[seen] <- k * tau2 / (sigma2 + k * tau2) * residual_mean[match(seen, levels(cl))] +
          stats::rnorm(length(seen), 0, sqrt(tau2 * sigma2 / (sigma2 + k * tau2)))
      }
      draw <- as.vector(x[fill_rows, , drop = FALSE] %*% beta) + unname(effect[fill_clusters]) + stats::rnorm(length(fill_rows), 0, sqrt(sigma2))
      break
    }
    if (is.null(draw))
      stop("Audited clusters lost design support in imputation; audit more clusters.", call. = FALSE)
    out[, j] <- draw
  }
  attr(out, "separable") <- moments$separable
  out
}

catmech_reference_clustered_draws <- function(data, value, gold, predictors, m, cluster) {
  W <- data[[value]]
  known <- !is.na(gold)
  recorded <- !is.na(W)
  calibration <- known & recorded
  if (!any(calibration))
    stop("use_recorded requires rows with both recorded and reference values.", call. = FALSE)
  draws <- matrix(rep(gold, m), nrow = nrow(data), ncol = m)
  note <- "linear_given_recorded"
  uncertainty <- "random_intercept_drawn_each_imputation"
  ids <- data[[cluster]]
  keep_recorded <- any(recorded & !known) &&
    all(catmech_contamination_agreement(W[calibration], gold[calibration], 0))
  if (!keep_recorded && any(as.character(ids[recorded & !known]) %in%
                           as.character(ids[!recorded & !known])))
    stop(paste0("A cluster cannot require stochastic imputation for both recorded-present and recorded-missing rows: ",
      "a shared cluster draw across these paths is not implemented. Supply reference values for one pattern ",
      "in each affected cluster, or use a joint model that supports both."), call. = FALSE)
  X <- catmech_model_matrix(data, predictors)
  if (any(recorded & !known)) {
    if (length(unique(ids[calibration])) < 5L)
      stop("Fewer than five audited clusters: cluster uncertainty is refused.", call. = FALSE)
    if (keep_recorded) {
      draws[recorded & !known, ] <- W[recorded & !known]
      note <- "recorded_kept_all_calibration_rows_agree"
      uncertainty <- "none_recorded_kept_all_calibration_rows_agree"
    } else {
      with_recorded <- cbind(X, .catmech_recorded = ifelse(recorded, W, 0))
      filled <- catmech_reference_cluster_fill(with_recorded, gold, which(calibration),
        which(recorded & !known), ids, m)
      if (!isTRUE(attr(filled, "separable"))) uncertainty <- "intercept_variance_zero_not_separable"
      draws[recorded & !known, ] <- filled
    }
  }
  if (any(!recorded & !known)) {
    filled <- catmech_reference_cluster_fill(X, gold, which(known), which(!recorded & !known), ids, m)
    if (!isTRUE(attr(filled, "separable"))) uncertainty <- "intercept_variance_zero_not_separable"
    else if (keep_recorded) uncertainty <- "random_intercept_drawn_each_imputation"
    draws[!recorded & !known, ] <- filled
  }
  list(draws = draws, note = note, cluster_uncertainty = uncertainty)
}

# Type-1 matching on a bootstrap of the donors. A resample that loses design support
# is redrawn, which conditions on estimability rather than changing the donor rule.
catmech_reference_abb_draw <- function(observed, trusted, x, masked, donors) {
  pool <- which(trusted)
  recipients <- sum(masked)
  for (attempt in seq_len(100L)) {
    chosen <- pool[sample.int(length(pool), length(pool), replace = TRUE)]
    x_chosen <- x[chosen, , drop = FALSE]
    if (length(unique(observed[chosen])) < 2L || qr(cbind(1, x_chosen), tol = 1e-7)$rank != ncol(x) + 1L) next
    ry <- c(rep(TRUE, length(chosen)), rep(FALSE, recipients))
    return(mice::mice.impute.pmm(c(observed[chosen], rep(NA_real_, recipients)), ry,
      rbind(x_chosen, x[masked, , drop = FALSE]), wy = !ry, donors = donors, matchtype = 1L, ridge = 1e-5))
  }
  stop("Bootstrap donor samples repeatedly lost design support; obtain more references or reduce predictors.", call. = FALSE)
}

catmech_reference_audit <- function(original, completed, rows, imputation, reason, reconstructed) {
  data.frame(row = rows, reason = reason, reconstructed = reconstructed,
    changed = ifelse(is.na(original[rows]), NA, completed[rows] != original[rows]),
    donor_row = rep(NA_integer_, length(rows)), used_fallback = rep(FALSE, length(rows)),
    eligible_donors = rep(NA_integer_, length(rows)), imputation = rep(imputation, length(rows)),
    category_derived = rep(FALSE, length(rows)), stringsAsFactors = FALSE)
}

catmech_reference_donor_audit <- function(audit, observed, values, donors) {
  trusted <- which(!is.na(observed))
  possible <- lapply(values, function(v) trusted[observed[trusted] == v])
  count <- lengths(possible)
  if (any(count == 0L)) stop("Returned PMM value has no trusted source row.", call. = FALSE)
  audit$donor_row <- vapply(possible, function(rows) if (length(rows) == 1L) rows else NA_integer_, integer(1))
  audit$donor_value_match_count <- count
  audit$donor_value_rows <- vapply(possible, paste, collapse = ";", FUN.VALUE = character(1))
  audit$donor_row_identified <- count == 1L
  audit$eligible_donors <- rep(length(trusted), nrow(audit))
  audit$candidate_pool_size <- rep(as.integer(donors), nrow(audit))
  audit
}

catmech_reference_independent_diagnostics <- function(engine) {
  list(engine = engine, status = "not_applicable", convergence_pass = TRUE,
    summary = data.frame(parameter = character(), rhat = numeric(),
                         ess_bulk = numeric(), ess_tail = numeric()), traces = NULL,
    reason = "No MCMC chain. This status does not establish calibration or model adequacy.")
}

catmech_reference_assemble <- function(data, value, predictors, observed, draws,
                                      audits, diagnostics, spec, m, seed,
                                      claim_level, model_diagnostics) {
  n <- nrow(data)
  if (!is.matrix(draws) || !identical(dim(draws), c(n, as.integer(m))) || any(!is.finite(draws)))
    stop("Internal error: invalid completion matrix.", call. = FALSE)
  completed <- lapply(seq_len(m), function(j) {
    out <- data
    out[[value]] <- draws[, j]
    out
  })
  audit <- do.call(rbind, audits)
  rownames(audit) <- NULL
  masked <- is.na(observed)
  base <- data[unique(c(predictors, value))]
  base[[value]] <- observed
  if (ncol(base) == 1L) {
    auxiliary <- ".catmech_intercept"
    while (auxiliary %in% names(base)) auxiliary <- paste0(auxiliary, "_")
    base[[auxiliary]] <- 1
  }
  where <- matrix(FALSE, n, ncol(base), dimnames = list(NULL, names(base)))
  where[, value] <- masked
  mids <- mice::mice(base, m = m, where = where,
    method = stats::setNames(rep("", ncol(base)), names(base)),
    predictorMatrix = matrix(0L, ncol(base), ncol(base), dimnames = list(names(base), names(base))),
    maxit = 0L, remove.constant = FALSE, remove.collinear = FALSE, printFlag = FALSE)
  target_imputations <- as.data.frame(draws[masked, , drop = FALSE])
  names(target_imputations) <- as.character(seq_len(m))
  rownames(target_imputations) <- rownames(base)[masked]
  mids$imp[[value]] <- target_imputations
  mids$method[value] <- paste0("catmech_", spec$engine, "_external_draws")
  manifest <- list(schema_version = spec$version, engine = spec$engine,
    assumption_hash = catmech_hash(spec),
    data_hash = digest::digest(serialize(data, NULL, version = 2), algo = "sha256", serialize = FALSE),
    data_hash_format = "R serialization version 2; environment-dependent data fingerprint",
    predictors = predictors, m = as.integer(m), seed = as.integer(seed), rng_kind = RNGkind(),
    imputed_rows = which(masked),
    adapter_input = "trusted/reference target values with all unavailable targets masked",
    package_version = as.character(utils::packageVersion("catmech")),
    r_version = as.character(getRversion()), mice_version = as.character(utils::packageVersion("mice")),
    claim_level = claim_level, identification = "conditional_on_trust_reference_and_selection_assumptions",
    assumptions = unclass(spec))
  structure(list(original = data, completed = completed, mids = mids, audit = audit,
    diagnostics = diagnostics, spec = spec, manifest = manifest,
    model_diagnostics = model_diagnostics), class = "catmech_result")
}
