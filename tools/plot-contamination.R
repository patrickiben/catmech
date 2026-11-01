#!/usr/bin/env Rscript
# Four-panel scientific figure using only the completed frozen v0.4 study.
# Usage: Rscript tools/plot-contamination.R [NEW_OUTPUT_DIRECTORY] [STUDY_DIRECTORY]
args <- commandArgs(trailingOnly = TRUE)
file_arg <- grep("^--file=", commandArgs(), value = TRUE)
root <- normalizePath(file.path(dirname(sub("^--file=", "", file_arg[1L])), ".."))
outdir <- if (length(args)) args[1L] else file.path(root, "validation", "contamination-figures")
study <- if (length(args) >= 2L) args[2L] else file.path(root, "validation", "contamination-study")
if (dir.exists(outdir) && length(list.files(outdir, all.files = TRUE, no.. = TRUE)))
  stop("Use a new or empty figure directory; existing evidence is not overwritten.")
required <- c("summary.csv", "RESULTS.md", "source-hashes.csv", "settings.json", "session-info.txt")
if (!all(file.exists(file.path(study, required))) || file.exists(file.path(study, "checkpoint.rds")))
  stop("Only completed bundles may be plotted: ", study)
tab <- utils::read.csv(file.path(study, "summary.csv"), stringsAsFactors = FALSE)
# Palette checked with the dataviz validator: lightness, chroma, CVD and contrast pass.
colors <- c(contamination = "#C2410C", normal = "#176E9D", cole = "#7C5CC4", raking = "#3F8F4F", gold_only = "#687B83")
shapes <- c(contamination = 21L, normal = 22L, cole = 23L, gold_only = 22L, raking = 25L)
names_method <- c(contamination = "contamination", reference_normal = "normal", cole_mi = "cole",
                  gold_only = "gold_only", raking = "raking")
labels <- c(contamination = "Contamination model", normal = "Reference normal (v0.3)", cole = "Linear MI with recorded value",
            raking = "Generalized raking", gold_only = "Gold-only regression")
plotted_rows <- list()
row_for <- function(scenario, method, panel, statistic = "coverage", fraction = NA_real_) {
  d <- tab[tab$scenario == scenario & tab$method == method, , drop = FALSE]
  if (nrow(d) != 1L) stop("Missing or duplicate result: ", scenario, "/", method)
  d$panel <- panel; d$statistic <- statistic; d$reference_fraction <- fraction
  d$point <- if (statistic == "coverage") d$coverage_conditional else d$mean_width
  d$lower <- if (statistic == "coverage") d$coverage_wilson_lower else d$mean_width - stats::qnorm(.975) * d$width_mcse
  d$upper <- if (statistic == "coverage") d$coverage_wilson_upper else d$mean_width + stats::qnorm(.975) * d$width_mcse
  plotted_rows[[length(plotted_rows) + 1L]] <<- d
  d
}
panel_heading <- function(title, subtitle) {
  left <- graphics::grconvertX(.025, from = "nfc", to = "user")
  graphics::mtext(title, side = 3, line = 2.05, at = left, adj = 0, cex = 1.03, font = 2L, col = "#1E303B")
  graphics::mtext(subtitle, side = 3, line = .78, at = left, adj = 0, cex = .77, col = "#596C77")
}
dodged_coverage <- function(panel, rows, methods, limits, title, subtitle, ticks, label_method = "contamination") {
  graphics::par(mar = c(3.9, 13.5, 3.65, 1.05))
  y <- rev(seq_along(rows))
  offsets <- seq(.26, -.26, length.out = length(methods))
  graphics::plot(NA, xlim = limits, ylim = c(.45, length(rows) + .62), axes = FALSE,
                 xlab = "Coverage of 95% intervals (%)", ylab = "", cex.lab = .87)
  graphics::abline(v = ticks, col = "#E7ECEF", lwd = .8)
  graphics::abline(v = .95, col = "#8D5555", lty = 2L, lwd = 1.2)
  graphics::axis(1, at = ticks, labels = round(ticks * 100), cex.axis = .84, col = "#B6C0C5", col.axis = "#42535E", tck = -.018)
  graphics::axis(2, at = y, labels = names(rows), tick = FALSE, las = 1L, cex.axis = .80, line = .15, col.axis = "#34454F")
  for (i in seq_along(rows)) for (j in seq_along(methods)) {
    method <- methods[j]
    if (!nrow(tab[tab$scenario == rows[[i]] & tab$method == method, , drop = FALSE])) next
    d <- row_for(rows[[i]], method, panel)
    series <- names_method[[method]]
    col <- colors[[series]]
    yy <- y[i] + offsets[j]
    if (d$upper >= limits[1]) graphics::segments(max(d$lower, limits[1]), yy, d$upper, yy, col = col, lwd = 1.5)
    if (d$point >= limits[1]) graphics::points(d$point, yy, pch = shapes[[series]], cex = 1.05, bg = col, col = "white", lwd = .6)
    if (method == label_method) graphics::text(min(d$upper + .006, limits[2]), yy, sprintf("%.1f", 100 * d$point), adj = 0, cex = .72, col = "#1E303B", font = 2L)
  }
  panel_heading(title, subtitle)
}
draw_width <- function() {
  graphics::par(mar = c(3.9, 6.1, 3.65, 1.05))
  graphics::plot(NA, xlim = c(2, 56), ylim = c(0, .95), axes = FALSE,
    xlab = "Randomly measured reference fraction (%)", ylab = "Mean coefficient-interval width", cex.lab = .87)
  graphics::abline(h = seq(0, .9, .1), col = "#E7ECEF", lwd = .8)
  graphics::axis(1, at = c(5, 10, 25, 50), labels = c("5", "10", "25", "50"), cex.axis = .84, col = "#B6C0C5", tck = -.018)
  graphics::axis(2, at = seq(0, .9, .1), labels = sprintf("%.1f", seq(0, .9, .1)), cex.axis = .84, col = "#B6C0C5", tck = -.018)
  oracle <- tab$mean_width[tab$scenario == "missed_flags" & tab$method == "oracle"]
  graphics::abline(h = oracle, col = "#8D5555", lty = 2L, lwd = 1.1)
  graphics::text(55.5, oracle - .028, "Full-truth oracle", adj = 1, cex = .72, col = "#596C77")
  fractions <- c(reference_05pct = 5, reference_10pct = 10, missed_flags = 25, reference_50pct = 50)
  for (method in c("gold_only", "reference_normal", "contamination")) {
    d <- do.call(rbind, lapply(names(fractions), function(s) row_for(s, method, "A", "width", fractions[[s]] / 100)))
    series <- names_method[[method]]
    col <- colors[[series]]
    x <- d$reference_fraction * 100
    graphics::lines(x, d$point, col = col, lwd = 1.6, lty = if (series == "gold_only") 2L else 1L)
    graphics::arrows(x, d$lower, x, d$upper, angle = 90, code = 3L, length = .035, col = col, lwd = 1.2)
    graphics::points(x, d$point, pch = shapes[[series]], bg = col, col = "white", cex = 1.05, lwd = .6)
  }
  series <- c("gold_only", "normal", "contamination")
  graphics::legend("topright", inset = c(.02, .01), legend = labels[series], col = colors[series], pch = shapes[series],
    pt.bg = colors[series], lty = c(2L, 1L, 1L), lwd = 1.4, cex = .77, bty = "n", text.col = "#34454F", seg.len = 1.6)
  panel_heading("A  Using recorded values buys precision", "Missed flags. 25%: 1,000 cohorts; other fractions: 150 cohorts")
}
primary <- c("Clean records" = "clean", "Flagged errors" = "flagged_errors", "Missed flags" = "missed_flags",
             "Site shift (declared)" = "site_shift_declared")
four <- c("contamination", "reference_normal", "cole_mi", "raking")
credibility <- c("Targeted audit" = "targeted_audit", "Selective, clean records" = "nonignorable_clean",
                 "Selective, missed flags" = "nonignorable_missed", "Noisy, SD declared" = "noisy_declared",
                 "Noisy, wrong SD" = "noisy_wrong_sd", "Noisy, not declared" = "noisy_undeclared")
stress <- c("Heavy-tailed target" = "heavy_tail_target", "Nonlinear target" = "nonlinear_target",
            "Unequal target variance" = "heteroskedastic_target", "Two-sided errors" = "error_bimodal",
            "Gross and small errors" = "error_gross_and_small", "Multiplicative errors" = "error_multiplicative",
            "Errors depend on outcome" = "differential_outcome", "Errors depend on value" = "value_dependent_error",
            "Flag computed from value" = "flag_from_recorded", "Assay noise, tolerance" = "assay_noise_tolerance",
            "Assay noise, no tolerance" = "assay_noise_no_tolerance", "Recorded missing on value" = "recorded_missing_mnar",
            "Site shift, undeclared" = "site_shift_undeclared")
draw_figure <- function() {
  plotted_rows <<- list()
  old <- graphics::par(no.readonly = TRUE); on.exit(graphics::par(old), add = TRUE)
  graphics::par(mfrow = c(2L, 2L), oma = c(6.9, .4, 4.1, .4), family = "sans", las = 1,
    cex = 1, fg = "#34454F", col.lab = "#34454F", col.axis = "#42535E", xaxs = "i", yaxs = "i")
  draw_width()
  dodged_coverage("B", primary, four, c(.85, 1.015), "B  Coverage where the model holds",
                  "Same data for every method; 1,000 cohorts per row. Errors match the model here, which favors it", seq(.85, 1, .05))
  dodged_coverage("C", credibility, c("contamination", "reference_normal"), c(0, 1.09), "C  Selective or noisy references",
                  "25% expected references; 150 cohorts per row. Blank rows ran the contamination model alone", seq(0, 1, .2))
  dodged_coverage("D", stress, c("contamination", "cole_mi"), c(.4, 1.06), "D  Assumptions deliberately broken",
                  "Missed-flag fixture unless named; 150 cohorts per row. Coverage below 40% leaves the axis", seq(.4, 1, .1))
  graphics::mtext("Using recorded values without trusting them", outer = TRUE, side = 3,
    line = 2.0, adj = .018, cex = 1.43, font = 2L, col = "#182D3B")
  graphics::mtext("Held-out synthetic evidence | 20 imputations | A reference subset calibrates a model of recording errors",
    outer = TRUE, side = 3, line = .45, adj = .018, cex = .88, col = "#536977")
  series <- c("contamination", "normal", "cole", "raking")
  graphics::legend(x = graphics::grconvertX(.018, from = "ndc", to = "user"), y = graphics::grconvertY(.083, from = "ndc", to = "user"),
    legend = labels[series], horiz = TRUE, xpd = NA, bty = "n", pch = shapes[series], pt.bg = colors[series],
    col = colors[series], cex = .80, text.col = "#34454F")
  graphics::mtext("Coverage bars: 95% Wilson Monte Carlo intervals; dashed line: nominal 95%. Width bars: mean +/- 1.96 MCSE. Coverage and width are conditional on returned intervals.",
    outer = TRUE, side = 1, line = 3.85, adj = .018, cex = .74, col = "#536977")
  graphics::mtext("Refusals and failures are in the study tables. Cells sharing a fixture reuse cohorts and are not independent replications. All data are synthetic.",
    outer = TRUE, side = 1, line = 5.05, adj = .018, cex = .74, col = "#536977")
}
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
grDevices::pdf(file.path(outdir, "contamination-evidence.pdf"), width = 14, height = 11, pointsize = 12, family = "Helvetica", useDingbats = FALSE)
draw_figure(); grDevices::dev.off()
grDevices::png(file.path(outdir, "contamination-evidence.png"), width = 3080, height = 2420, res = 220, type = "cairo", pointsize = 12)
draw_figure(); grDevices::dev.off()
plotted <- do.call(rbind, plotted_rows)
if (any(!is.finite(plotted$point)) || any(plotted$lower > plotted$point) || any(plotted$upper < plotted$point))
  stop("Invalid plot data.")
utils::write.csv(plotted, file.path(outdir, "plotted-data.csv"), row.names = FALSE)
hashes <- utils::read.csv(file.path(study, "source-hashes.csv"), stringsAsFactors = FALSE)
jsonlite::write_json(list(description = "Four-panel static figure from the completed frozen v0.4 study",
  study_summary_sha256 = digest::digest(file = file.path(study, "summary.csv"), algo = "sha256"),
  plot_script_sha256 = digest::digest(file = file.path(root, "tools/plot-contamination.R"), algo = "sha256"),
  source_files_recorded = nrow(hashes), plotted_rows = nrow(plotted),
  coverage_interval = "95% Wilson interval for Monte Carlo binomial uncertainty, conditional on returned intervals",
  width_interval = "mean +/- qnorm(.975) * MCSE of mean width", r_version = R.version.string),
  file.path(outdir, "figure-provenance.json"), pretty = TRUE, auto_unbox = TRUE)
cat("Wrote figure bundle:", normalizePath(outdir), "\n")
