#!/usr/bin/env Rscript
# Static scientific figure from completed, source-verified coverage summaries.
# Run from any directory. Arguments: synthetic results, clinical results, output.
# The legacy curve is its rerun on the SAME held-out cohorts, not old pilot data.
args <- commandArgs(trailingOnly = TRUE)
file_arg <- grep("^--file=", commandArgs(), value = TRUE)
root <- normalizePath(file.path(dirname(sub("^--file=", "", file_arg[1L])), ".."))
synthetic_dir <- if (length(args) >= 1L) args[[1L]] else file.path(root, "validation", "coverage-confirmatory")
clinical_dir <- if (length(args) >= 2L) args[[2L]] else file.path(root, "validation", "coverage-clinical")
output_dir <- if (length(args) >= 3L) args[[3L]] else file.path(root, "validation", "figures")
if (dir.exists(output_dir) && length(list.files(output_dir, all.files = TRUE, no.. = TRUE)))
  stop("Figure output directory must be new or empty.")

read_finished <- function(directory, expected_mode) {
  required <- c("summary.csv", "RESULTS.md", "source-hashes.csv", "settings.json", "session-info.txt")
  if (!all(file.exists(file.path(directory, required))) ||
      file.exists(file.path(directory, "checkpoint-replicates.csv")))
    stop("Only completed result bundles are accepted: ", directory)
  settings <- jsonlite::read_json(file.path(directory, "settings.json"))
  if (!identical(settings$mode, expected_mode)) stop("Unexpected benchmark mode in ", directory)
  summary <- utils::read.csv(file.path(directory, "summary.csv"), stringsAsFactors = FALSE)
  fields <- c("family", "scenario", "method", "requested", "returned_intervals",
    "coverage_conditional", "coverage_mcse", "coverage_wilson_lower", "coverage_wilson_upper")
  if (!all(fields %in% names(summary))) stop("Coverage summary schema is incomplete.")
  summary[fields]
}

synthetic <- read_finished(synthetic_dir, "confirmatory")
clinical <- read_finished(clinical_dir, "clinical")
panels <- list(
  list(family = "synthetic", scenario = "false_positive_flags",
    title = "Uncertain quality flags", subtitle = "Synthetic benchmark",
    methods = c("legacy", "model_ignore", "trusted_normal")),
  list(family = "synthetic", scenario = "category_error",
    title = "Incorrect categories", subtitle = "Synthetic benchmark",
    methods = c("legacy", "model_ignore", "trusted_normal")),
  list(family = "clinical", scenario = "missed_flags",
    title = "Missed recording errors", subtitle = "Clinical-format synthetic stress",
    methods = c("model_ignore", "trusted_normal")),
  list(family = "clinical", scenario = "site_shift",
    title = "Site shift + missed errors", subtitle = "Clinical-format compound stress",
    methods = c("model_ignore", "trusted_normal"))
)
plotted <- do.call(rbind, lapply(seq_along(panels), function(i) {
  panel <- panels[[i]]
  data <- if (panel$family == "synthetic") synthetic else clinical
  data <- data[data$family == panel$family & data$scenario == panel$scenario &
                 data$method %in% panel$methods, , drop = FALSE]
  if (nrow(data) != length(panel$methods) || anyDuplicated(data$method))
    stop("Missing or duplicate result for panel: ", panel$scenario)
  data <- data[match(panel$methods, data$method), , drop = FALSE]
  data$panel <- i
  data
}))
rownames(plotted) <- NULL
ok <- plotted$returned_intervals > 0L
if (any(plotted$returned_intervals > plotted$requested) ||
    any(!is.finite(plotted$coverage_conditional[ok])) ||
    any(plotted$coverage_conditional[ok] < 0 | plotted$coverage_conditional[ok] > 1) ||
    any(plotted$coverage_wilson_lower[ok] > plotted$coverage_conditional[ok]) ||
    any(plotted$coverage_wilson_upper[ok] < plotted$coverage_conditional[ok]))
  stop("Invalid coverage/count fields.")
expected_mcse <- sqrt(plotted$coverage_conditional[ok] *
  (1 - plotted$coverage_conditional[ok]) / plotted$returned_intervals[ok])
if (any(abs(plotted$coverage_mcse[ok] - expected_mcse) > 1e-10))
  stop("Stored MCSE does not match the plotted coverage denominator.")

colors <- c(legacy = "#B35B36", model_ignore = "#176E9D", trusted_normal = "#53636B")
labels <- c(legacy = "Legacy PMM\nsensitivity", model_ignore = "Gaussian model\n(categories ignored)",
            trusted_normal = "Trusted normal MI\n(mice baseline)")
draw_figure <- function() {
  old <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old), add = TRUE)
  graphics::par(mfrow = c(2L, 2L), mar = c(4.9, 4.1, 3.2, 1.1),
    oma = c(5.0, 0.5, 4.0, 0.5), family = "sans", las = 1,
    col.axis = "#34434B", col.lab = "#34434B", fg = "#34434B", xaxs = "i", yaxs = "i")
  for (i in seq_along(panels)) {
    panel <- panels[[i]]
    d <- plotted[plotted$panel == i, , drop = FALSE]
    x <- seq_len(nrow(d))
    graphics::plot(NA, xlim = c(0.45, length(x) + 0.55), ylim = c(0, 1.075),
      axes = FALSE, xlab = "", ylab = "Coverage of 95% intervals (%)", cex.lab = 0.85)
    graphics::abline(h = seq(0, 1, .2), col = "#E8ECEE", lwd = 0.8)
    graphics::abline(h = .95, col = "#7E4747", lty = 2L, lwd = 1.25)
    graphics::axis(2, at = seq(0, 1, .2), labels = seq(0, 100, 20),
                    tick = FALSE, line = -0.25, cex.axis = .82)
    graphics::axis(1, at = x, labels = FALSE, tck = 0, col = "#B5BEC3")
    for (j in x) {
      method <- d$method[j]
      graphics::mtext(labels[[method]], side = 1, at = j, line = 1.0, cex = .79,
                       col = colors[[method]])
      graphics::mtext(sprintf("%d / %d returned", d$returned_intervals[j], d$requested[j]),
                       side = 1, at = j, line = 3.45, cex = .64, col = "#566168")
      if (d$returned_intervals[j] == 0L) {
        graphics::text(j, .50, "No returned\nintervals", cex = .80, col = colors[[method]])
        next
      }
      lo <- d$coverage_wilson_lower[j]
      hi <- d$coverage_wilson_upper[j]
      graphics::arrows(j, lo, j, hi, angle = 90, code = 3L, length = .045,
                       col = colors[[method]], lwd = 1.7)
      graphics::points(j, d$coverage_conditional[j], pch = 21L, bg = colors[[method]],
                        col = "white", cex = 1.35, lwd = 0.75)
      graphics::text(j, min(1.045, hi + .041),
        sprintf("%.1f%%", 100 * d$coverage_conditional[j]),
        cex = .87, font = 2L, col = colors[[method]])
    }
    graphics::mtext(paste0(LETTERS[i], "  ", panel$title), side = 3, line = 1.65,
                     adj = 0, font = 2L, cex = .98, col = "#1E2D36")
    graphics::mtext(panel$subtitle, side = 3, line = .45, adj = 0,
                     cex = .78, col = "#637078")
  }
  graphics::mtext("Coverage improvements have a transfer boundary", outer = TRUE,
                   side = 3, line = 2.0, adj = 0.025, cex = 1.30, font = 2L, col = "#182A36")
  graphics::mtext("Held-out synthetic evaluations: masking flagged values still leaves error when flags miss contamination.",
                   outer = TRUE, side = 3, line = .55, adj = 0.025, cex = .86, col = "#4A5C67")
  graphics::mtext("Dashed line: nominal 95%. Bars: 95% Wilson intervals for Monte Carlo coverage uncertainty.",
                   outer = TRUE, side = 1, line = .65, adj = 0.025, cex = .75, col = "#4A5C67")
  graphics::mtext("Coverage is conditional on returned intervals; counts show returned / requested. MCSEs are retained in plotted-data.csv.",
                   outer = TRUE, side = 1, line = 1.85, adj = 0.025, cex = .75, col = "#4A5C67")
  graphics::mtext("The site panel combines a shift with missed errors; it does not isolate a site effect. All data are synthetic.",
                   outer = TRUE, side = 1, line = 3.05, adj = 0.025, cex = .75, col = "#4A5C67")
}

if (!dir.exists(output_dir) && !dir.create(output_dir, recursive = TRUE)) stop("Could not create output directory.")
grDevices::pdf(file.path(output_dir, "coverage-comparison.pdf"), width = 11, height = 9,
                family = "Helvetica", useDingbats = FALSE)
draw_figure()
grDevices::dev.off()
grDevices::png(file.path(output_dir, "coverage-comparison.png"), width = 2420, height = 1980,
                res = 220, type = "cairo")
draw_figure()
grDevices::dev.off()
utils::write.csv(plotted, file.path(output_dir, "plotted-data.csv"), row.names = FALSE)
input_files <- c(file.path(synthetic_dir, c("summary.csv", "source-hashes.csv", "settings.json")),
                 file.path(clinical_dir, c("summary.csv", "source-hashes.csv", "settings.json")),
                 file.path(root, "tools", "plot-coverage.R"))
provenance <- list(
  description = "Descriptive static scientific coverage comparison from completed held-out runs",
  nominal = 0.95, point = "coverage conditional on returned intervals",
  interval = "95 percent Wilson interval for binomial Monte Carlo coverage uncertainty",
  mcse = "sqrt(coverage * (1 - coverage) / returned intervals); stored alongside points",
  legacy_comparison = "Same held-out cohort rerun, not historical pilot values",
  clinical_legacy = "Not evaluated because category-boundary conventions differ",
  clinical_site_stress = "Compound: fixed site003 shift plus random errors and 20 percent flag sensitivity; SITEID is in the analysis, so no isolated site-effect attribution",
  input_files = lapply(input_files, function(p) list(
    path = sub(paste0("^", root, "/"), "", normalizePath(p)),
    sha256 = digest::digest(file = p, algo = "sha256"))),
  r_version = R.version.string)
jsonlite::write_json(provenance, file.path(output_dir, "figure-provenance.json"),
                     auto_unbox = TRUE, pretty = TRUE)
writeLines(c("# Coverage comparison figure", "",
  "The figure compares legacy sensitivity PMM, the Gaussian model ignoring categories, and independently implemented trusted normal MI on the same held-out synthetic cohorts. The clinical-format panels expose remaining failures when contaminated measurements are missed by flags.", "",
  "The site panel is a compound stress: site 003 has a fixed 0.20 mg/dL shift, 8% of observed records have random additive errors, and flag sensitivity is only 20%. The missed-flags panel has 60% sensitivity. SITEID is included in the analysis, where a constant site shift alone can be absorbed; this comparison does not isolate the effect of site bias.", "",
  "Points show coverage among returned intervals. Error bars are 95% Wilson Monte Carlo intervals, which remain informative at zero or unit coverage; the accompanying CSV also gives each binomial MCSE and the returned/requested denominator. The dashed line is the nominal 95% target.", "",
  "The legacy constraint method was not run on clinical-format fixtures because their endpoint convention differs. The figure is descriptive evidence for the specific simulated settings, not a clinical validation or universal coverage claim.", "",
  "Artifacts: coverage-comparison.png, coverage-comparison.pdf, plotted-data.csv, figure-provenance.json.", "",
  "Reproduce with `Rscript tools/plot-coverage.R validation/coverage-confirmatory validation/coverage-clinical validation/new-figures` from the repository root. Existing figure directories are never overwritten."
), file.path(output_dir, "README.md"))
cat("Wrote figure and plotted data to", normalizePath(output_dir), "\n")
