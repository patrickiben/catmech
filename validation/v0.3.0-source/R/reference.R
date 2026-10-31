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
#'   or `"pmm"` for standard predictive mean matching using reference donors.
#' @param donors Nearest-donor count for PMM, default five; ignored by the normal
#'   method. Fewer reference observations than donors is an error when needed.
#' @param m Number of completions, at least two.
#' @param seed Nonnegative integer seed, with caller RNG restored.
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
#' @export
catmech_reference <- function(data, value, reference, predictors,
                              method = c("normal", "pmm"), donors = 5L,
                              m = 20L, seed = 1L) {
  method <- match.arg(method)
  catmech_reference_inputs(data, value, predictors, m, seed, donors,
                           excluded = c(value, reference))
  if (!is.character(reference) || length(reference) != 1L || is.na(reference) ||
      !nzchar(reference) || !reference %in% names(data) || reference == value)
    stop("reference must name an existing column different from the target.", call. = FALSE)
  gold <- data[[reference]]
  if (!is.numeric(gold) || any(!is.na(gold) & !is.finite(gold)))
    stop("reference must contain independently measured finite numeric values or NA.", call. = FALSE)
  known <- !is.na(gold)
  if (!any(known)) stop("At least one independently observed reference value is required.", call. = FALSE)
  X <- if (all(known)) NULL else catmech_model_matrix(data, predictors)
  restore <- preserve_rng()
  on.exit(restore(), add = TRUE)
  set.seed(as.integer(seed))
  if (all(known)) {
    draws <- matrix(rep(gold, m), nrow = nrow(data), ncol = m)
    model_diagnostics <- catmech_reference_independent_diagnostics("reference_complete")
  } else if (method == "normal") {
    work <- data
    work[[value]] <- gold
    fitted <- catmech_model(work, value = value, predictors = predictors,
                            category_mode = "ignore", m = m, seed = seed)
    draws <- vapply(fitted$completed, function(z) z[[value]], numeric(nrow(data)))
    model_diagnostics <- fitted$model_diagnostics
  } else {
    draws <- catmech_reference_pmm_draws(gold, X, donors, m)
    model_diagnostics <- catmech_reference_independent_diagnostics("reference_mice_pmm")
  }
  engine <- paste0("reference_", method)
  spec <- structure(list(version = "reference1", engine = engine, value = value,
    predictors = predictors,
    reference = reference, reference_policy = "observed_gold_only_all_others_unavailable",
    raw_target_policy = "never_trusted_automatically",
    reference_selection_assumption = "ignorable_given_complete_selected_predictors",
    reference_measurement_assumption = "exact_same_quantity_time_and_units",
    method = method, donors = if (method == "pmm") as.integer(donors) else NULL,
    matchtype = if (method == "pmm") 1L else NULL,
    ridge = if (method == "pmm") 1e-5 else NULL,
    prior = if (method == "normal") "p(beta,sigma_squared) proportional to 1/sigma_squared" else NULL),
    class = "catmech_reference_spec")
  audits <- lapply(seq_len(m), function(j) {
    row <- catmech_reference_audit(data[[value]], draws[, j], seq_len(nrow(data)), j,
      ifelse(known, "gold_reference", "reference_missing"), reconstructed = !known)
    row$reference_observed <- known
    if (method == "pmm") {
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
    if (method == "normal") "conditional_reference_model" else "conditional_reference_matching",
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

catmech_reference_pmm_draws <- function(observed, X, donors, m) {
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
    out[masked] <- mice::mice.impute.pmm(observed, trusted, x, wy = masked,
      donors = donors, matchtype = 1L, ridge = 1e-5)
    if (any(!is.finite(out)) || any(!out[masked] %in% observed[trusted]))
      stop("PMM produced a nonfinite value or a value outside trusted donor support.", call. = FALSE)
    draws[, j] <- out
  }
  draws
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
