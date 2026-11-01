#!/usr/bin/env Rscript
# Six-panel figure for the 25% real-data audit cells, from the completed study only.
# Usage: Rscript tools/plot-real-data.R [NEW_OUTPUT_DIRECTORY]
args <- commandArgs(trailingOnly = TRUE)
file_arg <- grep("^--file=", commandArgs(), value = TRUE)
root <- normalizePath(file.path(dirname(sub("^--file=", "", file_arg[1L])), ".."))
outdir <- if (length(args)) args[1L] else file.path(root, "validation", "real-data-figures")
study <- file.path(root, "validation", "real-data-study")
if (dir.exists(outdir) && length(list.files(outdir, all.files = TRUE, no.. = TRUE)))
  stop("Use a new or empty figure directory; existing evidence is not overwritten.")
if (!file.exists(file.path(study, "session-info.txt")) || file.exists(file.path(study, "checkpoint.rds")))
  stop("Only a completed study may be plotted.")
tab <- utils::read.csv(file.path(study, "summary.csv"), stringsAsFactors = FALSE)
tab <- tab[tab$scenario == "audit_25pct", ]
methods <- c(oracle = "Oracle (all references)", contamination = "Contamination model", cole_mi = "Linear MI with recorded value",
             raking = "Generalized raking", reference_normal = "Reference normal (v0.3)", naive_recorded = "Recorded as truth")
colors <- c(oracle = "#1E303B", contamination = "#C2410C", cole_mi = "#7C5CC4", raking = "#3F8F4F", reference_normal = "#176E9D", naive_recorded = "#687B83")
shapes <- c(oracle = 18L, contamination = 21L, cole_mi = 23L, raking = 25L, reference_normal = 22L, naive_recorded = 24L)
datasets <- c(ctgov = "ClinicalTrials.gov: enrollment against participant-flow total", nhanes = "NHANES: self-reported against measured weight")
panel <- function(dataset, statistic, show_labels, title) {
  d <- tab[tab$family == dataset, ]
  d <- d[match(names(methods), d$method), ]
  y <- rev(seq_along(methods))
  if (statistic == "coverage") {
    point <- d$coverage_conditional; lower <- d$coverage_wilson_lower; upper <- d$coverage_wilson_upper
    limits <- c(.86, 1); ticks <- seq(.86, 1, .02); label <- "Coverage of 95% intervals (%)"; nominal <- .95
  } else if (statistic == "width") {
    point <- d$mean_width; half <- stats::qnorm(.975) * d$width_mcse; lower <- point - half; upper <- point + half
    limits <- c(0, 1.12 * max(upper)); ticks <- pretty(c(0, max(upper)), 4); label <- "Mean interval width"; nominal <- NA
  } else {
    # Standard errors for RMSE from the delta method on squared errors are in the study tables; points only here.
    point <- d$rmse; lower <- upper <- point
    limits <- c(0, 1.12 * max(point)); ticks <- pretty(c(0, max(point)), 4); label <- "RMSE"; nominal <- NA
  }
  ticks <- ticks[ticks >= limits[1] & ticks <= limits[2]]
  graphics::par(mar = c(3.9, if (show_labels) 14.5 else 1.5, 3.2, 1.2))
  graphics::plot(NA, xlim = limits, ylim = c(.4, length(methods) + .6), axes = FALSE, xlab = label, ylab = "", cex.lab = .87)
  graphics::abline(v = ticks, col = "#E7ECEF", lwd = .8)
  if (!is.na(nominal)) graphics::abline(v = nominal, col = "#8D5555", lty = 2L, lwd = 1.2)
  graphics::axis(1, at = ticks, labels = if (statistic == "coverage") round(ticks * 100) else ticks, cex.axis = .82, col = "#B6C0C5", col.axis = "#42535E", tck = -.018)
  if (show_labels) graphics::axis(2, at = y, labels = methods, tick = FALSE, las = 1L, cex.axis = .82, line = .1, col.axis = "#34454F")
  for (i in seq_along(methods)) {
    col <- colors[[names(methods)[i]]]
    graphics::segments(lower[i], y[i], upper[i], y[i], col = col, lwd = 1.6)
    graphics::points(point[i], y[i], pch = shapes[[names(methods)[i]]], cex = 1.15, bg = col, col = if (names(methods)[i] == "oracle") col else "white", lwd = .6)
    graphics::text(upper[i] + .018 * diff(limits), y[i], if (statistic == "coverage") sprintf("%.1f", 100 * point[i]) else if (statistic == "width") sprintf("%.3f", point[i]) else sprintf("%.4f", point[i]),
                   adj = 0, cex = .70, col = "#1E303B")
  }
  graphics::mtext(title, side = 3, line = .9, at = graphics::grconvertX(.02, from = "nfc", to = "user"), adj = 0, cex = .86, font = 2L, col = "#1E303B")
}
draw_figure <- function() {
  old <- graphics::par(no.readonly = TRUE); on.exit(graphics::par(old), add = TRUE)
  graphics::par(oma = c(3.6, .4, 4.3, .4), family = "sans", las = 1, cex = 1, fg = "#34454F", col.lab = "#34454F", xaxs = "i", yaxs = "i")
  graphics::layout(matrix(1:6, 2L, 3L, byrow = TRUE), widths = c(1.65, 1, 1))
  for (name in names(datasets)) {
    panel(name, "coverage", TRUE, datasets[[name]])
    panel(name, "width", FALSE, "")
    panel(name, "rmse", FALSE, "")
  }
  graphics::mtext("Real recording differences: a useful gain, but not the best tool", outer = TRUE, side = 3,
    line = 2.2, adj = .012, cex = 1.3, font = 2L, col = "#182D3B")
  graphics::mtext("Random 25% audit | 1,000 cohorts of 600 drawn from each evaluation pool | Target is the whole pool's coefficient using references",
    outer = TRUE, side = 3, line = .75, adj = .012, cex = .85, col = "#536977")
  graphics::mtext("Coverage bars: 95% Wilson Monte Carlo intervals; dashed line: nominal 95%. Width bars: mean +/- 1.96 MCSE. Oracle coverage on ClinicalTrials.gov is below 95%, consistent with non-constant residual variance.",
    outer = TRUE, side = 1, line = 1.1, adj = .012, cex = .72, col = "#536977")
  graphics::mtext("Conditional on returned intervals. Neither dataset is a clinical source-data audit. Both values on ClinicalTrials.gov are sponsor-entered.",
    outer = TRUE, side = 1, line = 2.2, adj = .012, cex = .72, col = "#536977")
}
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
grDevices::pdf(file.path(outdir, "real-data-evidence.pdf"), width = 14, height = 7.6, pointsize = 12, family = "Helvetica", useDingbats = FALSE)
draw_figure(); grDevices::dev.off()
grDevices::png(file.path(outdir, "real-data-evidence.png"), width = 3080, height = 1672, res = 220, type = "cairo", pointsize = 12)
draw_figure(); grDevices::dev.off()
utils::write.csv(tab[tab$method %in% names(methods), ], file.path(outdir, "plotted-data.csv"), row.names = FALSE)
cat("Wrote figure bundle:", normalizePath(outdir), "\n")
