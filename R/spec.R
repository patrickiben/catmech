#' Declare assumptions for interval-constrained sensitivity completion
#'
#' The interval convention is [a,b] for the first category and (a,b] for
#' subsequent categories. Finite endpoints may define a restricted domain.
#' The probability is a sensitivity weight for replacing flagged, in-bound
#' observed values, not an estimated posterior corruption probability.
#' @param value,category Names of the continuous measurement and interval category.
#' @param levels Ordered, unique category labels.
#' @param thresholds Strictly increasing numeric endpoints, length levels + 1.
#' @param flag Optional quality-flag column name; NULL means no flags.
#' @param flag_yes,flag_no Labels for suspect and non-suspect flags. Missing flags
#'   do not mark an observed value as suspect; this is an explicit assumption.
#' @param corruption_probability Fixed scalar replacement probability in [0,1].
#'   Must be supplied explicitly when a flag column is configured. With no flag,
#'   NULL is normalized to zero.
#' @param conflict_policy Either "error" or "reconstruct". The latter declares
#'   the category authoritative for all observed category/value conflicts.
#' @param min_donors Minimum eligible donors before using the fallback policy.
#' @param n_donors Number of nearest eligible donors from which to sample.
#' @param fallback "error" (default) or "normal" for an explicit, experimental
#'   truncated-normal fallback when too few category donors exist.
#' @return A validated catmech_spec list.
#' @export
#' @examples
#' spec <- catmech_spec("value", "band", c("low", "high"), c(-Inf, 0, Inf))
#' spec$conflict_policy
catmech_spec <- function(value, category, levels, thresholds, flag = NULL,
                         flag_yes = "Yes", flag_no = "No",
                         corruption_probability = NULL,
                         conflict_policy = c("error", "reconstruct"),
                         min_donors = 1L, n_donors = 5L,
                         fallback = c("error", "normal")) {
  conflict_policy <- match.arg(conflict_policy)
  fallback <- match.arg(fallback)
  scalar_name <- function(x) is.character(x) && length(x) == 1L && !is.na(x) && nzchar(x)
  if (!scalar_name(value) || !scalar_name(category) ||
      (!is.null(flag) && !scalar_name(flag))) stop("Column names must be nonempty strings.", call. = FALSE)
  if (anyDuplicated(c(value, category, flag))) stop("value, category and flag must name different columns.", call. = FALSE)
  if (!scalar_name(flag_yes) || !scalar_name(flag_no) || flag_yes == flag_no)
    stop("flag_yes and flag_no must be distinct strings.", call. = FALSE)
  if (!is.character(levels) || !length(levels) || anyNA(levels) ||
      any(!nzchar(levels)) || anyDuplicated(levels))
    stop("levels must be unique nonempty character labels.", call. = FALSE)
  levels <- unname(levels)
  thresholds <- normalize_thresholds(thresholds)
  make_interval_bounds(thresholds, levels)
  if (is.null(corruption_probability)) {
    if (!is.null(flag)) stop("Supply corruption_probability explicitly when a flag is configured; it is a sensitivity assumption.", call. = FALSE)
    corruption_probability <- 0
  }
  check_probability(corruption_probability, "corruption_probability")
  check_count(min_donors, "min_donors")
  check_count(n_donors, "n_donors")
  structure(list(version = "1", value = value, category = category,
    levels = levels, thresholds = thresholds, flag = flag,
    flag_yes = flag_yes, flag_no = flag_no,
    corruption_probability = as.numeric(corruption_probability),
    conflict_policy = conflict_policy, min_donors = as.integer(min_donors),
    n_donors = as.integer(n_donors), fallback = fallback), class = "catmech_spec")
}

check_probability <- function(x, name) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x < 0 || x > 1)
    stop(name, " must be a finite scalar in [0,1].", call. = FALSE)
}

check_count <- function(x, name, minimum = 1L) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) ||
      x != floor(x) || x < minimum || x > .Machine$integer.max)
    stop(name, " must be an integer >= ", minimum, ".", call. = FALSE)
}

normalize_thresholds <- function(x) {
  if (!is.atomic(x) && !is.list(x)) stop("Invalid thresholds.", call. = FALSE)
  vapply(as.list(x), function(z) {
    if (length(z) != 1L || is.na(z)) stop("thresholds cannot contain missing endpoints.", call. = FALSE)
    if (is.character(z)) {
      if (z %in% c("Inf", "+Inf", "inf", ".inf", "+.inf")) return(Inf)
      if (z %in% c("-Inf", "-inf", "-.inf")) return(-Inf)
      stop("Use numeric thresholds or explicit infinity strings.", call. = FALSE)
    }
    if (!is.numeric(z)) stop("thresholds must be numeric.", call. = FALSE)
    as.numeric(z)
  }, numeric(1), USE.NAMES = FALSE)
}

make_interval_bounds <- function(thresholds, levels) {
  thresholds <- normalize_thresholds(thresholds)
  if (length(thresholds) != length(levels) + 1L || anyNA(thresholds) ||
      any(diff(thresholds) <= 0) || any(!is.finite(thresholds[-c(1L,length(thresholds))])))
    stop("thresholds must increase strictly and have length levels + 1; only outer endpoints may be infinite.", call. = FALSE)
  if (thresholds[1L] == Inf || utils::tail(thresholds, 1L) == -Inf)
    stop("Invalid infinite endpoints.", call. = FALSE)
  stats::setNames(lapply(seq_along(levels), function(k) list(
    lower = thresholds[k], upper = thresholds[k+1L],
    lower_inclusive = k == 1L, upper_inclusive = TRUE)), levels)
}

value_in_bounds <- function(v, bound) {
  lower <- if (isTRUE(bound$lower_inclusive)) v >= bound$lower else v > bound$lower
  upper <- if (isTRUE(bound$upper_inclusive)) v <= bound$upper else v < bound$upper
  is.finite(v) & lower & upper
}

#' Derive exact interval categories
#' @param x Numeric values. Missing values remain missing.
#' @param thresholds,levels Endpoints and category labels as in [catmech_spec()].
#' @return A factor with the specified levels. Out-of-domain values give NA.
#' @examples
#' vals <- c(3.2, 5.5, 7.8, 12.0)
#' catmech_derive_category(
#'   vals,
#'   thresholds = c(-Inf, 5, 10, Inf),
#'   levels = c("low", "med", "high")
#' )
#' @export
catmech_derive_category <- function(x, thresholds, levels) {
  thresholds <- normalize_thresholds(thresholds)
  make_interval_bounds(thresholds, levels)
  if (!is.numeric(x) || any(!is.finite(x) & !is.na(x)))
    stop("x must contain finite numeric values or NA.", call. = FALSE)
  factor(cut(x, breaks = thresholds, labels = levels, right = TRUE,
             include.lowest = TRUE), levels = levels)
}

catmech_suspect_mask <- function(flag, flag_yes = "Yes") {
  !is.na(flag) & as.character(flag) == flag_yes
}

validate_spec <- function(spec) {
  if (!is.list(spec) || is.null(names(spec)) || anyDuplicated(names(spec)))
    stop("spec must be a named list with unique fields.", call. = FALSE)
  allowed <- c("version", names(formals(catmech_spec)))
  extra <- setdiff(names(spec), allowed)
  if (length(extra)) stop("Unsupported specification field(s): ", paste(extra, collapse=", "), call.=FALSE)
  if (!is.null(spec$version) && !identical(as.character(spec$version), "1"))
    stop("Only specification version 1 is supported.", call. = FALSE)
  spec$version <- NULL
  do.call(catmech_spec, unclass(spec))
}

#' Read a strict YAML assumption specification
#' @param path Path to a YAML file. YAML expression evaluation is disabled.
#' @return A validated catmech_spec object. Legacy prototype YAML is deliberately
#'   rejected; see inst/specs/lab.yaml for the supported schema.
#' @examples
#' spec_file <- system.file("specs", "lab.yaml", package = "catmech")
#' if (file.exists(spec_file)) {
#'   sp <- read_catmech_spec(spec_file)
#'   print(sp)
#' }
#' @export
read_catmech_spec <- function(path) {
  validate_spec(yaml::read_yaml(path, eval.expr = FALSE))
}

#' Compute a canonical assumption hash
#'
#' Named maps are sorted, sequence order is retained, and numeric values use
#' full round-trip precision. Missing values, signed infinities and strings
#' have distinct type tags. No numerical rounding or null-dropping is applied.
#' @param spec A list, usually a validated catmech_spec.
#' @return A SHA-256 hex string. This hashes assumptions, not the input dataset.
#' @examples
#' sp <- catmech_spec(
#'   value = "lab", category = "stage",
#'   levels = c("low", "med", "high"),
#'   thresholds = c(-Inf, 5, 10, Inf)
#' )
#' catmech_hash(sp)
#' @export
catmech_hash <- function(spec) {
  digest::digest(canonical_json(spec), algo = "sha256", serialize = FALSE)
}

canonical_json <- function(x) {
  encode <- function(v) {
    if (is.null(v)) return(list(type = "null"))
    if (is.list(v)) {
      nm <- names(v)
      if (!is.null(nm)) {
        if (anyNA(nm) || any(!nzchar(nm)) || anyDuplicated(nm))
          stop("Canonical maps need unique nonempty names.", call. = FALSE)
        v <- v[order(nm, method = "radix")]
        return(list(type = "map", value = lapply(v, encode)))
      }
      return(list(type = "sequence", value = lapply(v, encode)))
    }
    if (!is.atomic(v) || is.factor(v)) stop("Unsupported canonical value type.", call. = FALSE)
    if (!is.null(names(v))) return(list(type = "named_vector",
      names = encode(names(v)), value = encode(unname(v))))
    if (is.numeric(v)) {
      s <- vapply(v, function(z) {
        if (is.nan(z)) return("NaN")
        if (is.na(z)) return("NA")
        if (is.infinite(z)) return(if (z > 0) "+Inf" else "-Inf")
        if (z == 0) return("0")
        formatC(z, digits = 17L, format = "g", decimal.mark = ".")
      }, character(1), USE.NAMES = FALSE)
      return(list(type = "number", value = as.list(s)))
    }
    if (is.character(v) || is.logical(v)) return(list(type = typeof(v), value =
      lapply(v, function(z) if (is.na(z)) list(missing = TRUE) else list(value = z))))
    stop("Unsupported canonical value type.", call. = FALSE)
  }
  as.character(jsonlite::toJSON(encode(unclass(x)), auto_unbox = TRUE, null = "null", digits = NA))
}
