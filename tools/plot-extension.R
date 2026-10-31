#!/usr/bin/env Rscript
# Four-panel scientific figure using only completed frozen v0.3 study outputs.
# Usage: Rscript tools/plot-extension.R [NEW_OUTPUT_DIRECTORY]
args <- commandArgs(trailingOnly = TRUE)
file_arg <- grep("^--file=", commandArgs(), value = TRUE)
root <- normalizePath(file.path(dirname(sub("^--file=", "", file_arg[1L])), ".."))
outdir <- if (length(args)) args[1L] else file.path(root, "validation", "extension-figures")
if (dir.exists(outdir) && length(list.files(outdir, all.files = TRUE, no.. = TRUE)))
  stop("Use a new or empty figure directory; existing evidence is not overwritten.")
archive <- file.path(root, "validation", "v0.3.0-source")
input_files <- character()
source_hashes <- NULL
read_finished <- function(mode) {
  directory <- file.path(root, "validation", paste0("extension-", mode))
  required <- c("summary.csv", "RESULTS.md", "source-hashes.csv", "settings.json", "session-info.txt")
  if (!all(file.exists(file.path(directory, required))) || file.exists(file.path(directory, "checkpoint.rds")))
    stop("Only completed bundles may be plotted: ", directory)
  if (!identical(jsonlite::read_json(file.path(directory, "settings.json"))$mode, mode))
    stop("Unexpected benchmark mode.")
  hashes <- utils::read.csv(file.path(directory, "source-hashes.csv"), stringsAsFactors = FALSE)
  for (i in seq_len(nrow(hashes))) {
    path <- file.path(archive, hashes$file[i])
    if (!file.exists(path) || digest::digest(file = path, algo = "sha256") != hashes$sha256[i])
      stop("Frozen numerical source integrity failed: ", hashes$file[i])
  }
  if (!is.null(source_hashes) && !identical(source_hashes, hashes)) stop("Studies used different numerical sources.")
  source_hashes <<- hashes
  input_files <<- c(input_files, file.path(directory, required))
  tab <- utils::read.csv(file.path(directory, "summary.csv"), stringsAsFactors = FALSE)
  tab$source_bundle <- paste0("extension-", mode)
  tab
}
distribution <- read_finished("distribution")
reference <- read_finished("reference")
stress <- read_finished("stress")
colors <- c(normal = "#176E9D", pmm = "#B85C32", gold_only = "#687B83", mask = "#7B5B56")
shapes <- c(normal = 21L, pmm = 24L, gold_only = 22L, mask = 21L)
plotted_rows <- list()
add_row <- function(data, scenario, method, panel, label, series, fraction = NA_real_, statistic = "coverage") {
  d <- data[data$scenario == scenario & data$method == method, , drop = FALSE]
  if (nrow(d) != 1L) stop("Missing or duplicate result: ", scenario, "/", method)
  d$panel <- panel; d$plot_label <- label; d$plot_series <- series
  d$reference_fraction <- fraction; d$statistic <- statistic
  d$point <- if (statistic == "coverage") d$coverage_conditional else d$mean_width
  d$mcse <- if (statistic == "coverage") d$coverage_mcse else d$width_mcse
  d$lower <- if (statistic == "coverage") d$coverage_wilson_lower else d$mean_width - stats::qnorm(.975) * d$width_mcse
  d$upper <- if (statistic == "coverage") d$coverage_wilson_upper else d$mean_width + stats::qnorm(.975) * d$width_mcse
  d$interval_kind <- if (statistic == "coverage") "95% Wilson Monte Carlo interval" else "mean +/- 1.96 MCSE"
  plotted_rows[[length(plotted_rows) + 1L]] <<- d
}
# The shared reference rows appear ONCE; do not imply five independent studies.
add_row(reference, "missed_flags", "v02_mask_normal", "A", "Flag masking: missed errors", "mask")
add_row(reference, "site_shift", "v02_mask_normal", "A", "Flag masking: compound stress", "mask")
add_row(reference, "missed_flags", "reference_normal", "A", "25% gold: normal MI", "normal", .25)
add_row(reference, "missed_flags", "gold_only", "A", "25% gold: gold-only analysis", "gold_only", .25)
add_row(reference, "missed_flags", "reference_pmm", "A", "25% gold: PMM (experimental)", "pmm", .25)
for (s in c("heavy_tails", "strong_mar")) for (method in c("normal_no_cat", "pmm_no_cat")) {
  series <- if (method == "normal_no_cat") "normal" else "pmm"
  label <- paste(if (s == "heavy_tails") "Heavy tails" else "Strong MAR", "/", if (series == "normal") "Normal MI" else "PMM")
  add_row(distribution, s, method, "B", label, series)
}
for (fraction in c(.10, .25, .50)) for (method in c("reference_normal", "gold_only", "reference_pmm")) {
  data <- if (fraction == .25) reference else stress
  scenario <- if (fraction == .25) "missed_flags" else paste0("missed_flags_gold", fraction * 100)
  series <- switch(method, reference_normal = "normal", gold_only = "gold_only", reference_pmm = "pmm")
  add_row(data, scenario, method, "C", paste0(fraction * 100, "% gold"), series, fraction, "width")
}
add_row(reference, "clean", "reference_normal", "D", "Random, accurate gold\n1,000 / 1,000", "normal", .25)
add_row(stress, "nonignorable_selection", "reference_normal", "D", "Selective gold (MNAR)\n1,000 / 1,000", "normal", .25)
add_row(stress, "noisy_reference", "reference_normal", "D", "Noisy reference assay\n500 / 500", "normal", .25)
plotted <- do.call(rbind, plotted_rows); rownames(plotted) <- NULL
if (any(plotted$returned_intervals != plotted$requested) || any(plotted$other_errors != 0L) ||
    any(plotted$refused != 0L) || any(!is.finite(plotted$point)) ||
    any(plotted$lower > plotted$point) || any(plotted$upper < plotted$point))
  stop("Unexpected failures or invalid plot data. Figure must be reviewed before changing denominator conventions.")
coverage_rows <- plotted$statistic == "coverage"
if (any(abs(plotted$coverage_mcse[coverage_rows] - sqrt(plotted$point[coverage_rows] *
    (1 - plotted$point[coverage_rows]) / plotted$returned_intervals[coverage_rows])) > 1e-12))
  stop("Stored coverage MCSE disagrees with denominator.")
input_files <- c(input_files, file.path(root, "tools", "plot-extension.R"))
hash_inputs <- function() vapply(input_files, function(p) digest::digest(file = p, algo = "sha256"), character(1))
input_hashes <- hash_inputs()

panel_heading <- function(title, subtitle) {
  left <- graphics::grconvertX(.025, from = "nfc", to = "user")
  graphics::mtext(title, side = 3, line = 2.05, at = left, adj = 0, cex = 1.03, font = 2L, col = "#1E303B")
  graphics::mtext(subtitle, side = 3, line = .78, at = left, adj = 0, cex = .77, col = "#596C77")
}
draw_coverage <- function(panel, title, subtitle, limits = c(0, 1.075), ypos = NULL) {
  d <- plotted[plotted$panel == panel, ]
  y <- if (is.null(ypos)) rev(seq_len(nrow(d))) else ypos
  graphics::par(mar = c(3.9, 15.5, 3.65, 1.05))
  graphics::plot(NA, xlim = limits, ylim = c(.45, max(y) + .63), axes = FALSE,
    xlab = "Coverage of 95% intervals (%)", ylab = "", cex.lab = .87)
  ticks <- if (limits[1] == 0) seq(0, 1, .2) else seq(.85, 1, .05)
  graphics::abline(v = ticks, col = "#E7ECEF", lwd = .8)
  graphics::abline(v = .95, col = "#8D5555", lty = 2L, lwd = 1.2)
  graphics::axis(1, at = ticks, labels = round(ticks * 100), cex.axis = .84, col = "#B6C0C5", col.axis = "#42535E", tck = -.018)
  graphics::axis(2, at = y, labels = d$plot_label, tick = FALSE, las = 1L, cex.axis = .80, line = .15, col.axis = "#34454F")
  for (i in seq_len(nrow(d))) {
    col <- colors[[d$plot_series[i]]]
    graphics::segments(d$lower[i], y[i], d$upper[i], y[i], col = col, lwd = 1.8)
    graphics::segments(c(d$lower[i], d$upper[i]), y[i] - .075, c(d$lower[i], d$upper[i]), y[i] + .075, col = col, lwd = 1.5)
    graphics::points(d$point[i], y[i], pch = shapes[[d$plot_series[i]]], cex = 1.15, bg = col, col = "white", lwd = .7)
    graphics::text(d$point[i], y[i] + .25, sprintf("%.1f%%", 100 * d$point[i]), cex = .80, font = 2L, col = col)
  }
  panel_heading(title, subtitle)
}
draw_width <- function() {
  graphics::par(mar = c(3.9, 6.1, 3.65, 1.05))
  graphics::plot(NA, xlim = c(5, 56), ylim = c(0, .65), axes = FALSE,
    xlab = "Randomly measured reference fraction (%)", ylab = "Mean coefficient-interval width", cex.lab = .87)
  graphics::abline(h = seq(0, .6, .1), col = "#E7ECEF", lwd = .8)
  graphics::axis(1, at = c(10, 25, 50), labels = c("10", "25", "50"), cex.axis = .84, col = "#B6C0C5", tck = -.018)
  graphics::axis(2, at = seq(0, .6, .1), labels = sprintf("%.1f", seq(0, .6, .1)), cex.axis = .84, col = "#B6C0C5", tck = -.018)
  for (series in c("gold_only", "normal", "pmm")) {
    d <- plotted[plotted$panel == "C" & plotted$plot_series == series, ]; d <- d[order(d$reference_fraction), ]
    x <- d$reference_fraction * 100; col <- colors[[series]]
    graphics::lines(x, d$point, col = col, lwd = 1.6, lty = switch(series, normal = 1L, gold_only = 2L, pmm = 3L))
    graphics::arrows(x, d$lower, x, d$upper, angle = 90, code = 3L, length = .045, col = col, lwd = 1.4)
    graphics::points(x, d$point, pch = shapes[[series]], bg = col, col = "white", cex = 1.10, lwd = .65)
  }
  graphics::legend("topright", inset = c(.02, .015), legend = c("Normal MI (default)", "Gold-only analysis", "PMM (experimental)"),
    col = colors[c("normal", "gold_only", "pmm")], pch = shapes[c("normal", "gold_only", "pmm")],
    pt.bg = colors[c("normal", "gold_only", "pmm")], lty = c(1L, 2L, 3L), lwd = 1.4,
    cex = .77, bty = "n", text.col = "#34454F", seg.len = 1.6)
  panel_heading("C  More references buy precision", "10% and 50%: 250 / 250; 25%: 1,000 / 1,000")
}
draw_figure <- function() {
  old <- graphics::par(no.readonly = TRUE); on.exit(graphics::par(old), add = TRUE)
  graphics::par(mfrow = c(2L, 2L), oma = c(5.2, .4, 4.1, .4), family = "sans", las = 1,
    cex = 1, fg = "#34454F", col.lab = "#34454F", col.axis = "#42535E", xaxs = "i", yaxs = "i")
  draw_coverage("A", "A  New measurements change what is identifiable", "Shared reference results shown once; 1,000 / 1,000 per row")
  draw_coverage("B", "B  Matching is not uniformly safer", "Same predictors, no category; 1,000 / 1,000 per row", limits = c(.85, 1.015), ypos = c(5, 4, 2, 1))
  draw_width()
  draw_coverage("D", "D  References must be credible", "Normal MI; 25% expected reference sample")
  graphics::mtext("Measured references help; assumptions still matter", outer = TRUE, side = 3,
    line = 2.0, adj = .018, cex = 1.43, font = 2L, col = "#182D3B")
  graphics::mtext("Held-out synthetic evidence | 20 imputations | Accurate reference assays provide additional information",
    outer = TRUE, side = 3, line = .45, adj = .018, cex = .88, col = "#536977")
  graphics::mtext("Coverage bars: 95% Wilson Monte Carlo intervals; dashed line: nominal 95%. Width bars: mean +/- 1.96 MCSE.",
    outer = TRUE, side = 1, line = .75, adj = .018, cex = .77, col = "#536977")
  graphics::mtext("Counts are returned / requested. Panel B expands the coverage axis. Panel A reuses the same 1,000 reference cohorts.",
    outer = TRUE, side = 1, line = 1.95, adj = .018, cex = .77, col = "#536977")
  graphics::mtext("Gold must be accurate and representative. Partial-reference PMM is experimental / sensitivity-only. All data are synthetic.",
    outer = TRUE, side = 1, line = 3.15, adj = .018, cex = .77, col = "#536977")
}
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
grDevices::pdf(file.path(outdir, "extension-evidence.pdf"), width = 14, height = 10.5, pointsize = 12, family = "Helvetica", useDingbats = FALSE)
draw_figure(); grDevices::dev.off()
grDevices::png(file.path(outdir, "extension-evidence.png"), width = 3080, height = 2310, res = 220, type = "cairo", pointsize = 12)
draw_figure(); grDevices::dev.off()
if (!identical(input_hashes, hash_inputs())) stop("An input changed while plotting; do not use this output.")
utils::write.csv(plotted, file.path(outdir, "plotted-data.csv"), row.names = FALSE)
provenance <- list(description = "Four-panel static figure from completed, source-verified frozen v0.3 numerical studies",
  numerical_source = "validation/v0.3.0-source", frozen_source_files_verified = nrow(source_hashes), nominal_coverage = .95,
  coverage_interval = "95% Wilson interval for Monte Carlo binomial uncertainty", width_interval = "mean +/- qnorm(.975) * MCSE of mean width",
  counts = "returned / requested; all plotted performance requests returned intervals", category_comparison = "Panel B excludes category predictors in both methods",
  reference_pairing = "Panel A reference rows shown once; same 1000 cohorts, gold subsets and imputation seeds across five recording settings",
  reference_information = "References are additional exact measured values on selected rows; W and categories are not used in reference imputation",
  site_stress = "Compound fixed site003 shift plus random errors and 20% flag sensitivity; SITEID included; no isolated site-effect attribution",
  panel_c_pairing = "10% and 50% use 250 cohorts, 25% uses 1000; reference subsets are nested within common cohorts",
  panel_d_stress = "MNAR selection uses hidden residuals and cohort-specific intercept; noisy gold adds independent assay error SD .08",
  reference_pmm = "Unchanged numerical evidence of undercoverage; v0.3.1 partial-reference pooling requires explicit sensitivity-only override",
  input_files = lapply(seq_along(input_files), function(i) list(path = substring(input_files[i], nchar(root) + 2L), sha256 = unname(input_hashes[i]))),
  output_files = lapply(c("extension-evidence.pdf", "extension-evidence.png", "plotted-data.csv"), function(p) list(path = p, sha256 = digest::digest(file = file.path(outdir, p), algo = "sha256"))),
  r_version = R.version.string)
jsonlite::write_json(provenance, file.path(outdir, "figure-provenance.json"), pretty = TRUE, auto_unbox = TRUE)
writeLines(c("# Measured references and remaining assumptions", "",
  "A: a new random 25% accurate reference sample supports normal MI under its model when flag masking misses corruption. Shared reference results are shown once, alongside two recording-error scenarios from the same 1000 cohorts. The compound stress combines a site shift, random errors and 20% flag sensitivity; SITEID is included, so it does not isolate a site effect.", "",
  "B: holding predictor information fixed and excluding categories, PMM improves this heavy-tail example but loses coverage under strong MAR. Its horizontal coverage axis is expanded; the other coverage panels start at zero.", "",
  "C: more random references improve precision. Width is not a validity measure: experimental reference PMM has narrower intervals but confirmed undercoverage (74.8%, 85.1%, 90.8% at the tested 10%, 25%, 50% fractions). The 10% and 50% rows use 250 cohorts; the 25% row uses 1000. Subsets are nested within common cohorts.", "",
  "D: normal reference MI fails when validation selection is nonignorable or supplied references are noisy. The three cases use a 25% expected reference fraction. Counts retain the 500-replicate noisy-assay design.", "",
  "Coverage error bars are 95% Wilson intervals for Monte Carlo uncertainty, not individual scientific confidence intervals. Mean-width bars are mean +/- 1.96 MCSE. Exact MCSEs, counts, labels and every plotted source row are preserved in plotted-data.csv. Dashed vertical lines show nominal 95% coverage.", "",
  "All results are synthetic, conditional on the stated designs, and use additional measured information for reference methods. Partial-reference PMM is experimental and its v0.3.1 pooled analysis requires explicit sensitivity-only opt-in. The figure does not certify general coverage or identify error without assumptions.", "",
  "The script accepts completed bundles only and verifies all 12 frozen numerical source/plan files against validation/v0.3.0-source. figure-provenance.json records SHA-256 hashes of inputs, script and generated figure/data files.", "",
  "Reproduce with `Rscript tools/plot-extension.R validation/new-extension-figures` from the repository root. Existing output directories are never overwritten. See ../EXTENSION_SYNTHESIS.md and ../EXTENSION_PLAN.md for the full evidence and prospective design.", "",
  "Artifacts: extension-evidence.png, extension-evidence.pdf, plotted-data.csv, figure-provenance.json."
), file.path(outdir, "README.md"))
cat("Wrote verified-input figure bundle:", normalizePath(outdir), "\n")
