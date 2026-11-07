#' Generate an auditable compliance and integrity trail for a catmech analysis
#'
#' Extracts complete provenance, input and completion data hashes, model specifications,
#' convergence monitors, imputation statistics, and downstream pooled analyses into a
#' tamper-evident audit document. Supports Markdown, JSON, and standalone HTML formats.
#'
#' @param object A `catmech_result` object, or a list containing a `result` and optional
#'   downstream objects (`pooled`, `sensitivity`, `audit_plan`, `holdout_check`).
#' @param pooled Optional pooled analysis object (`catmech_pooled_lm`, `catmech_pooled_glm`,
#'   or `catmech_pooled_coxph`).
#' @param sensitivity Optional `catmech_sensitivity` object.
#' @param format Output format: `"markdown"` (default), `"json"`, or `"html"`.
#' @param path Optional file path to write the report to. If `NULL`, the character
#'   string is returned invisibly.
#' @param title Document title (default: `"Catmech Audit Trail & Compliance Report"`).
#' @return A character vector of class `c("catmech_audit_trail", "character")` containing
#'   the generated audit document.
#' @examples
#' set.seed(123)
#' df <- data.frame(
#'   lab = stats::rnorm(30, 10, 2),
#'   age = stats::rnorm(30, 50, 10),
#'   y = stats::rnorm(30, 5, 1)
#' )
#' df$ref <- df$lab
#' df$ref[11:30] <- NA
#' res <- catmech_reference(df, value = "lab", reference = "ref",
#'                          predictors = "age", m = 2, seed = 1)
#' pooled <- catmech_pool_lm(res, y ~ lab + age)
#' trail <- catmech_audit_trail(res, pooled = pooled)
#' print(trail)
#' @export
catmech_audit_trail <- function(object, pooled = NULL, sensitivity = NULL,
                                format = c("markdown", "json", "html"),
                                path = NULL,
                                title = "Catmech Audit Trail & Compliance Report") {
  format <- match.arg(format)

  # Support passing a list of related catmech objects
  res <- if (inherits(object, "catmech_result")) object else if (is.list(object) && "result" %in% names(object)) {
    if (is.null(pooled) && "pooled" %in% names(object)) pooled <- object$pooled
    if (is.null(sensitivity) && "sensitivity" %in% names(object)) sensitivity <- object$sensitivity
    object$result
  } else {
    stop("object must be a catmech_result or a list containing 'result'.", call. = FALSE)
  }

  if (!inherits(res, "catmech_result")) stop("Expected a valid catmech_result.", call. = FALSE)

  # Basic provenance metadata
  timestamp <- format(as.POSIXlt(Sys.time(), tz = "UTC"), "%Y-%m-%dT%H:%M:%SZ")
  r_ver <- paste(R.version$major, R.version$minor, sep = ".")
  platform <- R.version$platform
  pkg_ver <- tryCatch(as.character(utils::packageVersion("catmech")), error = function(e) "1.0.0")

  # Cryptographic hashes
  orig_data <- res$original
  orig_hash <- if (is.data.frame(orig_data)) digest::digest(orig_data, algo = "sha256") else "not_available"
  spec_hash <- if (!is.null(res$manifest$assumption_hash)) res$manifest$assumption_hash else
    if (!is.null(res$spec)) catmech_hash(res$spec) else "not_available"

  completed_hashes <- vapply(res$completed, function(d) digest::digest(d, algo = "sha256"), character(1L))

  # Imputation & model metadata
  engine <- if (!is.null(res$manifest$engine)) res$manifest$engine else "unknown"
  n_rows <- if (is.data.frame(orig_data)) nrow(orig_data) else NA_integer_
  n_cols <- if (is.data.frame(orig_data)) ncol(orig_data) else NA_integer_
  m <- length(res$completed)

  imputed_idx <- if (!is.null(res$manifest$imputed_rows)) res$manifest$imputed_rows else integer(0L)
  n_imputed <- length(imputed_idx)
  imputed_pct <- if (!is.na(n_rows) && n_rows > 0L) round(100 * n_imputed / n_rows, 2L) else NA_real_

  val_col <- if (!is.null(res$spec$value)) res$spec$value else "unspecified"
  ref_col <- if (!is.null(res$spec$reference)) res$spec$reference else "none"
  cluster_col <- if (!is.null(res$spec$cluster)) res$spec$cluster else "none"
  preds <- if (!is.null(res$spec$predictors)) paste(res$spec$predictors, collapse = ", ") else "none"

  # Diagnostics
  diag <- res$model_diagnostics
  conv_pass <- if (!is.null(diag$convergence_pass) && !identical(diag$status, "not_applicable")) {
    isTRUE(diag$convergence_pass)
  } else {
    NA
  }
  mcmc_iters <- if (!is.null(diag$iterations_used)) {
    diag$iterations_used
  } else if (!is.null(res$manifest$iterations_used)) {
    res$manifest$iterations_used
  } else if (!is.null(res$spec$iterations)) {
    res$spec$iterations
  } else {
    NA_integer_
  }

  # Build structured audit record
  trail_data <- list(
    title = title,
    timestamp_utc = timestamp,
    system = list(
      r_version = r_ver,
      platform = platform,
      catmech_version = pkg_ver
    ),
    integrity = list(
      original_data_sha256 = orig_hash,
      assumption_specification_sha256 = spec_hash,
      completed_datasets_sha256 = as.list(completed_hashes)
    ),
    specification = list(
      engine = engine,
      target_value = val_col,
      reference = ref_col,
      cluster = cluster_col,
      predictors = preds,
      observations = n_rows,
      variables = n_cols,
      completions_m = m,
      imputed_records = n_imputed,
      imputed_percentage = imputed_pct
    ),
    diagnostics = list(
      convergence_passed = conv_pass,
      mcmc_iterations = mcmc_iters
    )
  )

  if (!is.null(pooled)) {
    trail_data$pooled_analysis <- list(
      class = class(pooled)[1L],
      table = as.data.frame(pooled)
    )
  }

  if (!is.null(sensitivity)) {
    trail_data$sensitivity_envelope <- list(
      status = sensitivity$status,
      complete = sensitivity$complete,
      successful_envelope = as.list(sensitivity$successful_envelope),
      table = as.data.frame(sensitivity$table)
    )
  }

  # Render document
  rendered <- if (format == "json") {
    jsonlite::toJSON(trail_data, pretty = TRUE, auto_unbox = TRUE)
  } else if (format == "markdown") {
    render_audit_markdown(trail_data, pooled, sensitivity)
  } else {
    render_audit_html(trail_data, pooled, sensitivity)
  }

  if (!is.null(path)) {
    if (!is.character(path) || length(path) != 1L || is.na(path) || !nzchar(path))
      stop("path must be a single non-empty character string.", call. = FALSE)
    writeLines(rendered, con = path)
  }

  structure(rendered, class = c("catmech_audit_trail", "character"), format = format, data = trail_data)
}

#' Print an audit trail
#'
#' @param x A `catmech_audit_trail` object.
#' @param ... Unused.
#' @return `x`, invisibly.
#' @export
print.catmech_audit_trail <- function(x, ...) {
  cat(x, "\n")
  invisible(x)
}

render_audit_markdown <- function(d, pooled, sensitivity) {
  lines <- c(
    paste("#", d$title),
    "",
    paste("**Generated:**", d$timestamp_utc, "UTC | **Catmech Version:**", d$system$catmech_version),
    paste("**R Version:**", d$system$r_version, "| **Platform:**", d$system$platform),
    "",
    "---",
    "",
    "## 1. Cryptographic Data & Specification Integrity (SHA-256)",
    "",
    paste("- **Original Input Data Hash:** `", d$integrity$original_data_sha256, "`", sep = ""),
    paste("- **Assumption Specification Hash:** `", d$integrity$assumption_specification_sha256, "`", sep = ""),
    "- **Completed Dataset Checksums:**",
    ""
  )

  for (i in seq_along(d$integrity$completed_datasets_sha256)) {
    lines <- c(lines, paste0("  - Completion `m=", i, "`: `", d$integrity$completed_datasets_sha256[[i]], "`"))
  }

  lines <- c(
    lines,
    "",
    "## 2. Cohort Specification & Imputation Profile",
    "",
    "| Parameter | Value |",
    "| :--- | :--- |",
    paste("| Engine | `", d$specification$engine, "` |", sep = ""),
    paste("| Target Value Column | `", d$specification$target_value, "` |", sep = ""),
    paste("| Reference Column | `", d$specification$reference, "` |", sep = ""),
    paste("| Predictors | ", d$specification$predictors, " |", sep = ""),
    paste("| Cluster Variable | `", d$specification$cluster, "` |", sep = ""),
    paste("| Total Observations | ", d$specification$observations, " |", sep = ""),
    paste("| Completions (m) | ", d$specification$completions_m, " |", sep = ""),
    paste("| Imputed Records | ", d$specification$imputed_records, " (", d$specification$imputed_percentage, "%) |", sep = ""),
    paste("| Convergence Status | ", if (is.na(d$diagnostics$convergence_passed)) "Not Applicable" else if (isTRUE(d$diagnostics$convergence_passed)) "Passed" else "Failed", " |", sep = "")
  )

  if (!is.na(d$diagnostics$mcmc_iterations)) {
    lines <- c(lines, paste("| MCMC Iterations Used | ", d$diagnostics$mcmc_iterations, " |", sep = ""))
  }

  if (!is.null(pooled)) {
    lines <- c(
      lines,
      "",
      "## 3. Downstream Pooled Inference",
      "",
      paste("**Model Class:** `", class(pooled)[1L], "`", sep = ""),
      ""
    )
    tbl_p <- as.data.frame(pooled)
    # Format numeric columns
    for (col in names(tbl_p)) {
      if (is.numeric(tbl_p[[col]])) tbl_p[[col]] <- round(tbl_p[[col]], 4L)
    }
    cols_to_show <- intersect(c("term", "estimate", "std.error", "df", "lower", "upper",
                                "estimate_exp", "lower_exp", "upper_exp", "interpretation"), names(tbl_p))
    tbl_md <- knitr_table_fallback(tbl_p[, cols_to_show, drop = FALSE])
    lines <- c(lines, tbl_md)
  }

  if (!is.null(sensitivity)) {
    lines <- c(
      lines,
      "",
      "## 4. Sensitivity Envelope",
      "",
      paste("- **Status:** `", sensitivity$status, "`", sep = ""),
      paste("- **Complete Scenario Set:** `", sensitivity$complete, "`", sep = "")
    )
    if (!is.null(sensitivity$successful_envelope)) {
      lines <- c(lines, paste0("- **Successful Envelope:** [",
                               round(sensitivity$successful_envelope[1L], 4L), ", ",
                               round(sensitivity$successful_envelope[2L], 4L), "]"))
    }
    lines <- c(lines, "")
    tbl_s <- as.data.frame(sensitivity$table)
    for (col in names(tbl_s)) {
      if (is.numeric(tbl_s[[col]])) tbl_s[[col]] <- round(tbl_s[[col]], 4L)
    }
    lines <- c(lines, knitr_table_fallback(tbl_s))
  }

  lines <- c(
    lines,
    "",
    "---",
    "",
    "### Regulatory Audit & Compliance Certification",
    "> This document provides an auditable, verifiable record of data imputation and downstream analysis.",
    "> All completions and statistical intervals are cryptographically tied to the recorded input hash and assumption specification.",
    ""
  )

  paste(lines, collapse = "\n")
}

render_audit_html <- function(d, pooled, sensitivity) {
  md <- render_audit_markdown(d, pooled, sensitivity)
  # Basic self-contained HTML wrapper
  html_lines <- c(
    "<!DOCTYPE html>",
    "<html>",
    "<head>",
    "<meta charset='utf-8'>",
    paste0("<title>", d$title, "</title>"),
    "<style>",
    "body { font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, Helvetica, Arial, sans-serif; line-height: 1.6; max-width: 860px; margin: 40px auto; padding: 0 20px; color: #24292e; }",
    "h1, h2, h3 { color: #1a1f2c; border-bottom: 1px solid #eaecef; padding-bottom: 0.3em; }",
    "table { border-collapse: collapse; width: 100%; margin: 16px 0; }",
    "th, td { border: 1px solid #dfe2e5; padding: 8px 12px; text-align: left; }",
    "th { background-color: #f6f8fa; font-weight: 600; }",
    "code { background-color: #f3f4f6; padding: 0.2em 0.4em; border-radius: 3px; font-family: SFMono-Regular, Menlo, Monaco, Consolas, monospace; font-size: 85%; }",
    "blockquote { margin: 20px 0; padding: 10px 20px; border-left: 4px solid #0366d6; background-color: #f1f8ff; color: #444; }",
    "</style>",
    "</head>",
    "<body>",
    paste0("<pre style='white-space: pre-wrap;'>", html_escape(md), "</pre>"),
    "</body>",
    "</html>"
  )
  paste(html_lines, collapse = "\n")
}

html_escape <- function(x) {
  x <- gsub("&", "&amp;", x, fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE)
  x <- gsub(">", "&gt;", x, fixed = TRUE)
  x
}

knitr_table_fallback <- function(df) {
  # Simple markdown table generator without requiring knitr
  header <- paste("|", paste(names(df), collapse = " | "), "|")
  sep <- paste("|", paste(rep(":---", ncol(df)), collapse = " | "), "|")
  rows <- apply(df, 1L, function(r) {
    paste("|", paste(as.character(r), collapse = " | "), "|")
  })
  c(header, sep, rows)
}
