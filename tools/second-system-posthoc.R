#!/usr/bin/env Rscript
# Readings of the second-system test made after its results. Post hoc: they enter no rule.
# Usage: Rscript tools/second-system-posthoc.R [OUTPUT_FILE]
args <- commandArgs(trailingOnly = TRUE)
file_arg <- grep("^--file=", commandArgs(), value = TRUE)
root <- normalizePath(file.path(dirname(sub("^--file=", "", file_arg[1])), ".."))
output <- if (length(args)) args[1] else file.path(root, "validation", "second-system-posthoc.txt")
second <- parse(file.path(root, "tools/second-system-benchmark.R"))
named <- vapply(second, function(x) if (is.call(x) && identical(x[[1]], as.name("<-")) && is.name(x[[2]])) as.character(x[[2]]) else "", character(1))
SEEDS <- jsonlite::read_json(file.path(root, "validation/second-system-study/settings.json"), simplifyVector = TRUE)$seeds
for (n in c("ITEMS", "draw_cohort")) eval(second[[which(named == n)[1]]], envir = globalenv())
pool_all <- utils::read.csv(file.path(root, "validation/second-system/zambia-pool.csv"), stringsAsFactors = FALSE)
saved <- utils::read.csv(gzfile(file.path(root, "validation/second-system-study/replicates.csv.gz")), stringsAsFactors = FALSE)
saved <- saved[saved$run == "primary", ]
slope <- function(d) unname(stats::coef(stats::lm(verified ~ covariate, data = d))[2])
lines <- c("Post hoc: made after the results, enters no rule", "",
           "THE MOST INFLUENTIAL DISTRICT, verified values: its clinics of the pool; the pool slope with and without it; share of cohorts that include it; coverage of the slope at clinic_25 when it is in the cohort and when it is not")
for (i in seq_along(ITEMS)) {
  pool <- pool_all[pool_all$item == ITEMS[i] & pool_all$both_zero == 0L, ]
  full <- slope(pool)
  deleted <- vapply(sort(unique(pool$district)), function(g) slope(pool[pool$district != g, ]), numeric(1))
  top <- names(which.max(abs(deleted - full)))
  by_district <- split(seq_len(nrow(pool)), pool$district)
  present <- vapply(1:1000, function(r) top %in% unique(pool$district[draw_cohort(pool, i, r)$rows]), logical(1))
  parts <- vapply(c("oracle", "linear_mi", "random_intercept"), function(m) {
    s <- saved[saved$item == ITEMS[i] & saved$audit == "clinic_25" & saved$method == m & saved$status == "ok", ]
    p <- present[s$replicate]
    sprintf("%s %.1f, %.1f", m, 100 * mean(s$cover[p]), 100 * mean(s$cover[!p]))
  }, character(1))
  lines <- c(lines, sprintf("  %s: %d of %d clinics; slope %.1f, without it %.1f; in %.1f%% of cohorts; %s", ITEMS[i], length(by_district[[top]]), nrow(pool), full, deleted[[top]], 100 * mean(present), paste(parts, collapse = "; ")))
}
s <- saved[saved$method == "linear_mi" & saved$status == "ok", ]
wide <- s[order(-s$width), ][1:3, ]
lines <- c(lines, "", "WIDEST LINEAR MI INTERVALS, against the median width for that item and audit")
for (j in seq_len(nrow(wide))) {
  same <- s[s$item == wide$item[j] & s$audit == wide$audit[j], ]
  lines <- c(lines, sprintf("  %s %s replicate %d: width %.3g, median %.3g, pooled df %.3g", wide$item[j], wide$audit[j], wide$replicate[j], wide$width[j], stats::median(same$width), wide$df[j]))
}
tiny <- tapply(s$df < 1, s$audit, mean)
lines <- c(lines, paste0("  Share of linear MI requests with pooled df under 1: ", paste(sprintf("%s %.2f%%", names(tiny), 100 * tiny), collapse = ", ")))
kept <- saved[saved$method == "linear_mi" & saved$kept > 0, ]
lines <- c(lines, "", sprintf("REQUESTS WHERE EVERY AUDITED ROW AGREED, linear MI: %d", nrow(kept)))
writeLines(lines, output)
cat(lines, sep = "\n"); cat("\n")
