#' Generate auditable sensitivity completions
#'
#' Supports one continuous target with complete predictors. Each completion
#' uses an independent parameter and donor draw; no iterative chain is needed
#' for this restricted setting. Categories must be observed, except that a
#' missing category can be derived from an unflagged observed value.
#' Predictors must be selected explicitly to avoid accidental use of identifiers
#' or simulation truth. Include downstream outcomes when appropriate for the
#' intended analysis. Completed data and Rubin pooling are experimental and do
#' not establish identification, proper MI, or interval coverage.
#' @param data A data.frame. The original is stored without modification.
#' @param spec A [catmech_spec()] or equivalent validated list.
#' @param predictors Character vector of complete predictor column names.
#'   Use character(0) for an intercept-only model. The category, target and flag
#'   cannot be predictors. Category information restricts the donor support.
#' @param m Number of completed datasets, at least 2.
#' @param seed Nonnegative integer random seed. The caller's RNG is restored.
#' @return A catmech_result containing original, completed, mids, audit,
#'   diagnostics, spec and manifest. mids is an analysis adapter created from
#'   the supplied completions, not evidence of mice chain convergence.
#' @export
#' @examples
#' d <- simulate_catmech_lab(150, seed = 12)
#' d <- d[setdiff(names(d), c("lab_true", "corrupted"))]
#' s <- catmech_spec("lab_value", "lab_category", c("Low", "Medium", "High"),
#'   c(-Inf, 5.5, 7, Inf), flag = "outlier_flag", corruption_probability = 1,
#'   conflict_policy = "reconstruct")
#' fit <- catmech_impute(d, s, c("age", "sex", "inflammation", "outcome"), m = 3)
#' fit$diagnostics
catmech_impute <- function(data, spec, predictors, m = 20L, seed = 1L) {
  spec <- validate_spec(spec)
  check_count(m, "m", 2L)
  check_count(seed, "seed", 0L)
  if (!is.data.frame(data) || nrow(data) < 2L || anyNA(names(data)) || anyDuplicated(names(data)) ||
      any(!nzchar(names(data)))) stop("data must have at least two rows and unique nonempty column names.", call. = FALSE)
  if (any(c(".imp", ".id") %in% names(data))) stop(".imp and .id are reserved column names.", call. = FALSE)
  if (!is.character(predictors) || anyNA(predictors) || anyDuplicated(predictors))
    stop("predictors must be an explicit character vector of unique column names.", call. = FALSE)
  special <- c(spec$value, spec$category, spec$flag)
  needed <- c(special, predictors)
  if (length(setdiff(needed, names(data)))) stop("Missing data columns: ",
    paste(setdiff(needed, names(data)), collapse = ", "), call. = FALSE)
  if (any(predictors %in% special))
    stop("Predictors must exclude the value, category and flag columns.", call. = FALSE)
  if (any(!vapply(data, function(z) is.atomic(z) && is.null(dim(z)), logical(1))))
    stop("List and matrix columns are not supported.", call. = FALSE)
  y <- data[[spec$value]]
  if (!is.numeric(y) || any(!is.finite(y) & !is.na(y)))
    stop("The value column must be numeric, with finite values or NA.", call. = FALSE)
  category <- as.character(data[[spec$category]])
  if (any(!is.na(category) & !category %in% spec$levels))
    stop("Unknown category labels in the data.", call. = FALSE)
  flag <- if (is.null(spec$flag)) rep(spec$flag_no, nrow(data)) else as.character(data[[spec$flag]])
  if (any(!is.na(flag) & !flag %in% c(spec$flag_no, spec$flag_yes)))
    stop("Unknown flag labels in the data.", call. = FALSE)
  suspect <- catmech_suspect_mask(flag, spec$flag_yes)
  derived <- is.na(category)
  if (any(derived & (is.na(y) | suspect)))
    stop("A missing category cannot be derived from a missing or flagged measurement.", call. = FALSE)
  category[derived] <- as.character(catmech_derive_category(y[derived], spec$thresholds, spec$levels))
  if (anyNA(category)) stop("Values used to derive categories fall outside the specified domain.", call. = FALSE)
  bounds <- make_interval_bounds(spec$thresholds, spec$levels)
  in_bounds <- vapply(seq_along(y), function(i) value_in_bounds(y[i], bounds[[category[i]]]), logical(1))
  conflict <- !is.na(y) & !in_bounds
  if (any(conflict) && spec$conflict_policy == "error")
    stop(sum(conflict), " category/value conflict(s). Resolve them or explicitly set conflict_policy='reconstruct'.", call. = FALSE)
  where_y <- is.na(y) | suspect | conflict
  trusted <- !is.na(y) & !suspect & !conflict
  px <- data[predictors]
  if (length(predictors)) {
    if (anyNA(px) || any(!vapply(px, function(z) is.numeric(z) || is.factor(z) || is.character(z) || is.logical(z), logical(1))))
      stop("Predictors must be complete numeric, logical, character or factor columns.", call. = FALSE)
    if (any(vapply(px, function(z) is.numeric(z) && any(!is.finite(z)), logical(1))))
      stop("Predictors must be finite.", call. = FALSE)
    for (nm in names(px)) if (is.character(px[[nm]]) || is.factor(px[[nm]])) px[[nm]] <- droplevels(factor(px[[nm]]))
    constant <- vapply(px, function(z) length(unique(z)) < 2L, logical(1))
    if (any(constant)) stop("Remove constant predictors: ", paste(names(px)[constant], collapse=", "), call.=FALSE)
    x <- stats::model.matrix(~ ., data = px)[, -1L, drop = FALSE]
  } else x <- matrix(numeric(0), nrow = nrow(data), ncol = 0L)
  restore <- preserve_rng()
  on.exit(restore(), add = TRUE)
  set.seed(as.integer(seed))
  completed <- vector("list", m)
  audits <- vector("list", m)
  for (j in seq_len(m)) {
    draw <- mice.impute.catmech_threshold_corruption(y, trusted, x, wy = where_y,
      category = category, bounds = bounds, flag = flag, flag_yes = spec$flag_yes,
      corruption_probability = spec$corruption_probability, original_y = y,
      trusted = trusted, min_donors = spec$min_donors,
      n_donors = spec$n_donors, fallback = spec$fallback)
    d <- data
    d[[spec$value]][where_y] <- as.numeric(draw)
    d[[spec$category]] <- factor(category, levels = spec$levels)
    a <- attr(draw, "audit")
    if (is.null(a)) stop("Internal error: missing audit record.", call. = FALSE)
    a$imputation <- rep(j, nrow(a))
    a$category_derived <- derived[a$row]
    audits[[j]] <- a
    ok <- vapply(seq_along(y), function(i) value_in_bounds(d[[spec$value]][i], bounds[[category[i]]]), logical(1))
    if (anyNA(d[[spec$value]]) || !all(ok)) stop("Internal error: completion violates declared intervals.", call. = FALSE)
    completed[[j]] <- d
  }
  audit <- do.call(rbind, audits)
  rownames(audit) <- NULL
  diagnostics <- do.call(rbind, lapply(seq_len(m), function(j) {
    a <- audits[[j]]
    data.frame(imputation = j, requested = sum(where_y), reconstructed = sum(a$reconstructed),
      changed = sum(a$changed, na.rm = TRUE), threshold_conflicts = sum(conflict),
      fallback_draws = sum(a$used_fallback), threshold_violations = 0L)
  }))
  # Restrict the adapter to the scientific variables explicitly selected by the
  # caller. Unused metadata columns remain intact in original/completed.
  adapter_names <- unique(c(predictors, spec$value, spec$category))
  base <- data[adapter_names]
  base[[spec$category]] <- factor(as.character(base[[spec$category]]), levels = spec$levels)
  long <- do.call(rbind, lapply(0:m, function(j) {
    z <- if (j == 0L) base else completed[[j]][adapter_names]
    z$.imp <- j
    z$.id <- as.character(seq_len(nrow(z)))
    z
  }))
  where <- matrix(FALSE, nrow(data), length(adapter_names), dimnames = list(NULL, adapter_names))
  where[, spec$value] <- where_y
  where[, spec$category] <- derived
  mids <- mice::as.mids(long, where = where)
  # as.mids initializes a default method vector: label this adapter truthfully.
  mids$method[] <- ""
  mids$method[spec$value] <- "catmech_external_draws"
  mids$call <- match.call()
  manifest <- list(schema_version = "1", assumption_hash = catmech_hash(spec),
    data_hash = digest::digest(serialize(data, NULL, version = 2), algo = "sha256", serialize = FALSE),
    data_hash_format = "R serialization version 2; environment-dependent data fingerprint",
    predictors = predictors, m = as.integer(m), seed = as.integer(seed),
    rng_kind = RNGkind(), package_version = as.character(utils::packageVersion("catmech")),
    r_version = as.character(getRversion()), mice_version = as.character(utils::packageVersion("mice")),
    claim_level = "assumption_based_sensitivity", identification = "not_identified",
    category_derived_rows = which(derived), assumptions = unclass(spec))
  structure(list(original = data, completed = completed, mids = mids, audit = audit,
    diagnostics = diagnostics, spec = spec, manifest = manifest), class = "catmech_result")
}

preserve_rng <- function() {
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  old <- if (had_seed) get(".Random.seed", envir = .GlobalEnv) else NULL
  kind <- RNGkind()
  function() {
    do.call(RNGkind, as.list(kind))
    if (had_seed) assign(".Random.seed", old, envir = .GlobalEnv)
    else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
      rm(".Random.seed", envir = .GlobalEnv)
  }
}

#' Write a self-contained result bundle
#' @param result A catmech_result.
#' @param directory New or empty output directory. Existing files are never overwritten.
#' @return The normalized output path, invisibly. Files include a lossless RDS,
#'   long CSV, per-cell audit, diagnostics, JSON manifest and session information.
#' @export
write_catmech <- function(result, directory) {
  if (!inherits(result, "catmech_result")) stop("Expected catmech_result.", call. = FALSE)
  if (dir.exists(directory) && length(list.files(directory, all.files = TRUE, no.. = TRUE)))
    stop("Output directory must be new or empty.", call. = FALSE)
  if (!dir.exists(directory) && !dir.create(directory, recursive = TRUE)) stop("Could not create output directory.", call. = FALSE)
  long <- do.call(rbind, lapply(0:length(result$completed), function(j) {
    z <- if (j == 0L) result$original else result$completed[[j]]
    z$.imp <- j
    z$.id <- seq_len(nrow(z))
    z
  }))
  saveRDS(result, file.path(directory, "result.rds"), version = 2)
  utils::write.csv(long, file.path(directory, "completed_long.csv"), row.names = FALSE)
  utils::write.csv(result$audit, file.path(directory, "audit.csv"), row.names = FALSE)
  utils::write.csv(result$diagnostics, file.path(directory, "diagnostics.csv"), row.names = FALSE)
  manifest <- result$manifest
  json_safe <- function(x) {
    if (is.list(x)) return(lapply(x, json_safe))
    if (is.numeric(x) && any(!is.finite(x))) return(lapply(x, function(z) {
      if (is.na(z)) return(NULL)
      if (is.infinite(z)) return(if (z > 0) "Inf" else "-Inf")
      z
    }))
    x
  }
  manifest <- json_safe(manifest)
  jsonlite::write_json(manifest, file.path(directory, "manifest.json"), auto_unbox=TRUE,
    pretty=TRUE, null="null", digits=NA)
  writeLines(canonical_json(result$spec), file.path(directory, "assumptions.canonical.json"))
  writeLines(utils::capture.output(utils::sessionInfo()), file.path(directory, "sessionInfo.txt"))
  invisible(normalizePath(directory))
}
