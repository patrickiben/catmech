#!/usr/bin/env Rscript
# Plot completed confirmation evidence only; refuse partial checkpoints.
args <- commandArgs(trailingOnly = TRUE)
script <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1L])
root <- normalizePath(file.path(dirname(script), ".."), mustWork = TRUE)
input <- if (length(args)) args[1L] else file.path(root, "validation/joint-calibration-confirmation")
output <- if (length(args) > 1L) args[2L] else file.path(root, "validation/joint-calibration-figures")
if (!file.exists(file.path(input, "COMPLETE.json")))
  stop("Confirmation is incomplete: no result files were read and no figures were written.")
complete <- jsonlite::read_json(file.path(input, "COMPLETE.json"), simplifyVector = TRUE)
if (!identical(complete$status, "complete") || !identical(complete$mode, "confirmation"))
  stop("Figures require a completed confirmation stage.")
files <- c("COMPLETE.json", "settings.json", "records.csv", "fit-records.csv", "summary.csv",
  "prior-triples.csv", "prior-summary.csv", "availability.csv", "source-hashes-before.csv", "source-hashes-after.csv")
if (!all(file.exists(file.path(input, files)))) stop("Completed confirmation artifacts are missing.")
read_csv <- function(name) utils::read.csv(file.path(input, name), stringsAsFactors = FALSE)
config <- jsonlite::read_json(file.path(input, "settings.json"), simplifyVector = TRUE)
records <- read_csv("records.csv"); fits <- read_csv("fit-records.csv")
summary <- read_csv("summary.csv"); triples <- read_csv("prior-triples.csv")
prior_summary <- read_csv("prior-summary.csv"); availability <- read_csv("availability.csv")
stopifnot(identical(config$mode, "confirmation"), isTRUE(config$accuracy_exported),
  identical(read_csv("source-hashes-before.csv"), read_csv("source-hashes-after.csv")),
  nrow(records) == complete$analysis_requests, nrow(fits) == complete$fitted_workflow_requests,
  !anyDuplicated(records[c("cell", "replicate", "method")]),
  !anyDuplicated(fits[c("cell", "replicate", "fit_id")]),
  all(records$source_fingerprint == config$source_fingerprint),
  all(fits$source_fingerprint == config$source_fingerprint),
  identical(config$source_fingerprint, complete$source_fingerprint))
cells <- c("gaussian_balanced", "gaussian_unbalanced", "near_zero_icc", "observed_mar", "null_slope", "student_t5")
prior_cells <- c("gaussian_balanced", "near_zero_icc", "observed_mar")
methods <- c("old_cv3", "new_cv3", "direct")
stopifnot(identical(config$cells, c(cells, "contextual_x")),
  identical(config$prior_cells, prior_cells))
check <- function(actual, expected, label) {
  actual <- as.numeric(actual); expected <- as.numeric(expected)
  good <- length(actual) == length(expected) && identical(is.na(actual), is.na(expected))
  if (good) {
    observed <- !is.na(expected)
    good <- all(is.finite(actual[observed]) & is.finite(expected[observed]) &
      abs(actual[observed] - expected[observed]) <= 1e-12 + 1e-11 * abs(expected[observed]))
  }
  if (!good) stop("Plotted-data self-check failed: ", label)
}
wilson <- function(k, n) {
  if (!n) return(c(NA_real_, NA_real_))
  z <- stats::qnorm(.975); p <- k / n; denominator <- 1 + z^2 / n
  center <- (p + z^2 / (2 * n)) / denominator
  half <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / denominator
  c(max(0, center - half), min(1, center + half))
}
se_mean <- function(x) if (length(x) > 1L) stats::sd(x) / sqrt(length(x)) else NA_real_
mean_or_na <- function(x) if (length(x)) mean(x) else NA_real_
plot_data <- list()
add <- function(panel, cell, method, estimate, interval, n, events = NA_integer_, uncertainty) {
  plot_data[[length(plot_data) + 1L]] <<- data.frame(panel = panel, cell = cell,
    method = method, estimate = estimate, lower = interval[1L], upper = interval[2L],
    n = n, events = events, uncertainty = uncertainty, stringsAsFactors = FALSE)
}
# Recompute every coverage point and Wilson interval from individual requests,
# then compare to the saved summary rather than trusting it as plotting input.
for (cell in c(cells, "contextual_x")) for (method in methods) {
  rows <- records[records$cell == cell & records$method == method, ]
  expected_n <- if (cell == "contextual_x") config$contextual_reps else config$primary_reps
  stopifnot(nrow(rows) == expected_n, setequal(rows$replicate, seq_len(expected_n)),
    all(rows$status %in% c("ok", "error")), all(rows$covered %in% c(0, 1)),
    all(rows$covered[rows$status != "ok"] == 0))
  n <- sum(rows$status == "ok"); k <- sum(rows$covered)
  p <- if (n) k / n else NA_real_; interval <- wilson(k, n)
  saved <- summary[summary$cell == cell & summary$method == method, ]
  stopifnot(nrow(saved) == 1L)
  check(unlist(saved[c("requests", "intervals", "failures", "coverage_conditional",
    "coverage_conditional_lower", "coverage_conditional_upper")]),
    c(nrow(rows), n, nrow(rows) - n, p, interval), paste(cell, method, "coverage"))
  add(if (cell == "contextual_x") "B" else "A", cell, method,
    p, interval, n, k, "95% Wilson Monte Carlo interval; denominator = returned intervals")
}
# Reconstruct sensitivity from all three raw widths. Both analyses must have
# all three intervals, giving the same cohort denominator for their comparison.
for (cell in prior_cells) {
  saved_triples <- triples[triples$cell == cell, ]
  stopifnot(nrow(saved_triples) == config$prior_reps,
    !anyDuplicated(saved_triples$replicate), setequal(saved_triples$replicate, seq_len(config$prior_reps)))
  direct <- cv3 <- rep(NA_real_, config$prior_reps)
  for (replicate in seq_len(config$prior_reps)) {
    r <- records[records$cell == cell & records$replicate == replicate, ]
    saved <- saved_triples[saved_triples$replicate == replicate, ]
    for (method in c("new_cv3", "direct")) {
      selected <- r[match(paste0(method, c("", "_prior_quarter", "_prior_four")), r$method), ]
      stopifnot(nrow(selected) == 3L, !anyNA(selected$method))
      good <- all(selected$status == "ok")
      prefix <- if (method == "direct") "direct" else "cv3"
      stopifnot(identical(good, saved[[paste0(prefix, "_complete")]]))
      sensitivity <- if (good) mean(abs(log(selected$width[-1L] / selected$width[1L]))) else NA_real_
      check(saved[[paste0(prefix, "_sensitivity")]], sensitivity, paste(cell, replicate, method, "S"))
      if (method == "direct") direct[replicate] <- sensitivity else cv3[replicate] <- sensitivity
    }
    check(saved$difference, direct[replicate] - cv3[replicate], paste(cell, replicate, "paired S"))
  }
  common <- is.finite(direct) & is.finite(cv3)
  d <- direct[common]; c <- cv3[common]
  saved <- prior_summary[prior_summary$cell == cell, ]; stopifnot(nrow(saved) == 1L)
  check(unlist(saved[c("requested_triples", "direct_triples", "cv3_triples", "complete_pairs",
    "mean_direct_sensitivity", "mean_cv3_sensitivity", "difference", "mcse")]),
    c(config$prior_reps, sum(is.finite(direct)), sum(is.finite(cv3)), sum(common),
      mean_or_na(d), mean_or_na(c), mean_or_na(d - c), se_mean(d - c)), paste(cell, "prior summary"))
  for (method in c("new_cv3", "direct")) {
    values <- if (method == "direct") d else c
    average <- mean_or_na(values)
    add("C", cell, method, average, average + c(-1, 1) * stats::qnorm(.975) * se_mean(values),
      length(values), uncertainty = "95% normal Monte Carlo interval for mean S; common complete prior triples")
  }
}
# Fit failures are distinct from analysis-stage failures; use fit-records.csv.
for (cell in cells) for (id in c("old", "new")) {
  rows <- fits[fits$cell == cell & fits$fit_id == id, ]
  stopifnot(nrow(rows) == config$primary_reps, setequal(rows$replicate, seq_len(config$primary_reps)),
    all(rows$fit_status %in% c("ok", "error")))
  k <- sum(rows$fit_status != "ok"); n <- nrow(rows)
  saved <- availability[availability$cell == cell, ]; stopifnot(nrow(saved) == 1L)
  check(saved$requests, n, paste(cell, "fit requests"))
  check(saved[[paste0(id, "_failures")]], k, paste(cell, id, "fit failures"))
  check(saved[[paste0(id, "_failure_rate")]], k / n, paste(cell, id, "fit refusal rate"))
  add("D", cell, paste0(id, "_cv3"), k / n, wilson(k, n), n, k,
    "95% Wilson Monte Carlo interval; denominator = all fit requests")
}
plot_data <- do.call(rbind, plot_data)
labels <- c(gaussian_balanced = "Gaussian balanced", gaussian_unbalanced = "Gaussian unbalanced",
  near_zero_icc = "Near-zero ICC", observed_mar = "Observed-data MAR", null_slope = "Null slope", student_t5 = "Student-t5")
colors <- c(old_cv3 = "#697386", new_cv3 = "#2468A0", direct = "#B95716")
draw_points <- function(rows, y, multiplier = 1) {
  color <- unname(colors[rows$method]); good <- is.finite(rows$estimate)
  ci <- good & is.finite(rows$lower) & is.finite(rows$upper)
  graphics::segments(multiplier * rows$lower[ci], y[ci], multiplier * rows$upper[ci], y[ci],
    col = color[ci], lwd = 1.4)
  graphics::segments(multiplier * rows$lower[ci], y[ci] - .045, multiplier * rows$lower[ci], y[ci] + .045, col = color[ci])
  graphics::segments(multiplier * rows$upper[ci], y[ci] - .045, multiplier * rows$upper[ci], y[ci] + .045, col = color[ci])
  graphics::points(multiplier * rows$estimate[good], y[good], pch = 16, cex = 1.0, col = color[good])
}
panel_axes <- function(xlim, ylim, yticks, ylabels, title, subtitle, xlab, reference = NULL, xticks = NULL) {
  graphics::plot(NA_real_, NA_real_, type = "n", xlim = xlim, ylim = ylim,
    xlab = xlab, ylab = "", axes = FALSE, xaxs = "i", yaxs = "i")
  ticks <- if (is.null(xticks)) pretty(xlim) else xticks
  ticks <- ticks[ticks >= xlim[1L] & ticks <= xlim[2L]]
  graphics::abline(v = ticks, col = "#E8EBEF", lwd = .8)
  if (!is.null(reference)) graphics::abline(v = reference, col = "#59616E", lty = 2L, lwd = 1.2)
  graphics::axis(1L, at = ticks, col = "#9BA4B0", col.axis = "#273142", cex.axis = .85)
  graphics::axis(2L, at = yticks, labels = ylabels, las = 1L, tick = FALSE,
    col.axis = "#273142", cex.axis = .88)
  graphics::title(main = title, adj = 0, line = 2.1, cex.main = 1.08, col.main = "#172B45")
  graphics::mtext(subtitle, side = 3L, adj = 0, line = .5, cex = .76, col = "#526073")
}
render <- function() {
  graphics::par(mfrow = c(2L, 2L), mar = c(4.4, 8.9, 3.8, 1.6), oma = c(4.0, 0, 5.2, 0),
    family = "sans", cex = 1, col.lab = "#273142", mgp = c(2.7, .65, 0))
  z <- plot_data[plot_data$panel == "A", ]
  lower <- max(0, min(c(.9, z$lower[is.finite(z$lower)])) * 100 - 1)
  panel_axes(c(lower, 100), c(.45, 6.55), 6:1, labels[cells],
    "A  Coverage in the six original scenarios",
    sprintf("%d cohorts per scenario; coverage among returned intervals", config$primary_reps),
    "Coverage (%)", 95)
  draw_points(z, 7 - match(z$cell, cells) + c(old_cv3 = .18, new_cv3 = 0, direct = -.18)[z$method], 100)
  z <- plot_data[plot_data$panel == "B", ]
  panel_axes(c(0, 100), c(.5, 3.5), 3:1, c("Old CV3 CI", "New CV3 CI", "Posterior credible"),
    "B  Contextual-X assumption stress",
    sprintf("%d cohorts; marginal target 0.9, structural slope 0.7", config$contextual_reps),
    "Coverage among returned intervals (%)", 95)
  draw_points(z, 4 - match(z$method, methods), 100)
  z <- plot_data[plot_data$panel == "C", ]
  stopifnot(all(z$lower[is.finite(z$lower)] >= 0), all(z$upper[is.finite(z$upper)] <= .15))
  counts <- z$n[match(prior_cells, z$cell)]
  panel_axes(c(0, .15), c(.5, 3.5), 3:1,
    paste0(c("Balanced", "Near-zero ICC", "Observed MAR"), " (n=", counts, ")"),
    "C  Relative sensitivity of interval width",
    "Common complete triples; covariance-prior scales 0.25, 1 and 4",
    "S = mean absolute log width ratio (smaller = less change)", xticks = c(0, .05, .10, .15))
  draw_points(z, 4 - match(z$cell, prior_cells) + c(new_cv3 = .1, direct = -.1)[z$method])
  z <- plot_data[plot_data$panel == "D", ]
  upper <- max(c(2, z$upper[is.finite(z$upper)] * 100)) * 1.06
  panel_axes(c(-.025 * upper, upper), c(.45, 6.55), 6:1, labels[cells],
    "D  Numerical fit refusals",
    sprintf("Old / longer-run sampler; all %d fit requests per scenario", config$primary_reps),
    "Fit refusals (%)")
  draw_points(z, 7 - match(z$cell, cells) + c(old_cv3 = .12, new_cv3 = -.12)[z$method], 100)
  graphics::mtext("Joint reference model: confirmation evidence", side = 3L, outer = TRUE,
    line = 3.3, font = 2L, cex = 1.35, col = "#172B45")
  graphics::mtext("Prespecified synthetic scenarios; model-based credible intervals and CV3/Rubin confidence intervals", side = 3L,
    outer = TRUE, line = 1.9, cex = .9, col = "#526073")
  graphics::mtext("Error bars show 95% Monte Carlo uncertainty: Wilson intervals for proportions; normal intervals for mean S.",
    side = 1L, outer = TRUE, line = 1.3, cex = .8, col = "#526073")
  graphics::mtext("Coverage conditions on a returned interval. Contextual X deliberately violates the model. These results do not establish universal validity.",
    side = 1L, outer = TRUE, line = 2.5, cex = .8, col = "#526073")
  graphics::par(fig = c(0, 1, .904, .941), mar = rep(0, 4), oma = rep(0, 4), new = TRUE)
  graphics::plot.new()
  graphics::legend("center", legend = c("Old sampler + CV3 CI", "Longer-run sampler + CV3 CI", "Direct posterior credible interval"),
    col = colors, pch = 16, horiz = TRUE, bty = "n", cex = .9, x.intersp = .7)
}
dir.create(output, recursive = TRUE, showWarnings = FALSE)
utils::write.csv(plot_data, file.path(output, "plotted-data.csv"), row.names = FALSE, na = "")
grDevices::pdf(file.path(output, "joint-calibration-evidence.pdf"), width = 14, height = 10,
  family = "sans", useDingbats = FALSE)
render(); invisible(grDevices::dev.off())
grDevices::png(file.path(output, "joint-calibration-evidence.png"), width = 2800, height = 2000, res = 200,
  type = if (capabilities("cairo")) "cairo" else "quartz")
render(); invisible(grDevices::dev.off())
jsonlite::write_json(list(stage = "confirmation", source_fingerprint = complete$source_fingerprint,
  script_sha256 = digest::digest(file = script, algo = "sha256"),
  input_files = stats::setNames(lapply(file.path(input, files), function(f)
    digest::digest(file = f, algo = "sha256")), files),
  points = nrow(plot_data), self_checks = "Every plotted value recomputed from raw records and checked against saved summaries; prior triples share a common denominator.",
  uncertainty = "95% Monte Carlo intervals, not uncertainty intervals for an individual fitted coefficient."),
  file.path(output, "figure-provenance.json"), auto_unbox = TRUE, pretty = TRUE)
writeLines(c("# Joint calibration confirmation figure", "",
  "Generated only after the confirmation COMPLETE marker. The PDF is vector output; the PNG is 2800 by 2000 pixels.", "",
  "Panels A and B show coverage among returned intervals, with Wilson Monte Carlo intervals and a 95% reference line. B is a deliberately misspecified contextual-X stress, separate from the original scenarios.", "",
  "Panel C compares the mean absolute log width ratio across the prespecified prior scales on cohorts where both analyses returned all three intervals. Error bars are normal Monte Carlo intervals for each mean; the plotted-data file records the common denominator. These individual error bars do not replace the paired sensitivity analysis.", "",
  "Panel D shows model-fit refusals over every request, with Wilson intervals; it is not the count of analysis-stage refusals. Full analysis refusal information remains in the study records.", "",
  "Each plotted scalar is independently recomputed from saved raw records and checked against its saved summary before output. plotted-data.csv and figure-provenance.json retain numerical and source provenance.", "",
  "The direct intervals are conditional model-based Bayesian credible intervals. CV3 intervals use cluster jackknife variance and Rubin pooling. Neither the figure nor diagnostic success establishes universal validity."),
  file.path(output, "README.md"))
cat(sprintf("Wrote PDF, PNG and %d self-checked plotted points from completed confirmation.\n", nrow(plot_data)))
