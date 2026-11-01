#!/usr/bin/env Rscript
# Prospective test of catmech_audit_plan; rule and settings locked in validation/AUDIT_PLAN_PLAN.md.
# Usage: Rscript tools/audit-plan-benchmark.R pilot|study [NEW_DIRECTORY] [PILOT_REPS]
args <- commandArgs(trailingOnly = TRUE)
mode <- if (length(args)) match.arg(args[1], c("pilot", "study")) else "pilot"
file_arg <- grep("^--file=", commandArgs(), value = TRUE)
root <- normalizePath(file.path(dirname(sub("^--file=", "", file_arg[1])), ".."))
outdir <- if (length(args) >= 2L) args[2] else file.path(root, "validation", if (mode == "pilot") "audit-plan-pilot" else "audit-plan-study")
if (dir.exists(outdir) && length(list.files(outdir, all.files = TRUE, no.. = TRUE)))
  stop("Use a new or empty output directory; existing evidence is never overwritten.")
hash_sources <- function() {
  paths <- c(sort(list.files(file.path(root, "R"), "\\.R$", full.names = TRUE)),
    file.path(root, c("DESCRIPTION", "NAMESPACE", "tools/audit-plan-benchmark.R", "tools/contamination-fixtures.R", "tools/lab-fixture.R",
                      "validation/real-data/ctgov-pool.csv", "validation/real-data/nhanes-pool.csv")),
    if (mode == "study") file.path(root, "validation/AUDIT_PLAN_PLAN.md"))
  data.frame(file = substring(paths, nchar(root) + 2L),
    sha256 = vapply(paths, function(p) digest::digest(file = p, algo = "sha256"), character(1)), row.names = NULL)
}
hashes_before <- hash_sources()
pkgload::load_all(root, quiet = TRUE)
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
workers <- suppressWarnings(as.integer(Sys.getenv("CATMECH_BENCHMARK_WORKERS", "14")))
if (length(workers) != 1L || is.na(workers) || workers < 1L || workers > 16L) stop("Invalid worker count (1 to 16).")
if (.Platform$OS.type == "windows") workers <- 1L
pilot_reps <- if (length(args) >= 3L) as.integer(args[3]) else 40L
is_study <- mode == "study"
config <- list(mode = mode, workers = workers, m = 20L, n = 600L, bootstrap = 500L,
  real_reps = if (is_study) 1500L else pilot_reps, synthetic_reps = if (is_study) 300L else min(pilot_reps, 40L),
  retest_reps = if (is_study) 200L else min(pilot_reps, 40L),
  split = if (is_study) "evaluation" else "development",
  target_records = c(30L, 60L, 150L, 300L), pilot_records = c(60L, 150L), target_ratio = 1.25,
  cohort_base = if (is_study) 600100000 else 700100000, order_base = if (is_study) 610100000 else 710100000,
  imputation_base = if (is_study) 620000000 else 720000000, plan_10_base = if (is_study) 630000000 else 730000000,
  plan_25_base = if (is_study) 640000000 else 740000000, retest_base = if (is_study) 650000000 else 750000000,
  offsets = list(ctgov = 0, nhanes = 2000000, missed_flags = 4000000, error_gross_and_small = 6000000),
  fixture_master_seed = if (is_study) 430271L else 430272L, fixture_phase = if (is_study) "heldout" else "development",
  fixture_selection_base = if (is_study) 670100000 else 770100000, fixture_noise_base = if (is_study) 676100000 else 776100000,
  fixture_mechanism_base = if (is_study) 678100000 else 778100000, fixture_imputation_base = if (is_study) 690000000 else 790000000)
jsonlite::write_json(config, file.path(outdir, "settings.json"), pretty = TRUE, auto_unbox = TRUE)
fractions <- config$target_records / config$n

load_pool <- function(dataset) {
  if (dataset == "ctgov") {
    d <- utils::read.csv(file.path(root, "validation/real-data/ctgov-pool.csv"), stringsAsFactors = FALSE)
    d <- d[d$split == config$split, ]
    data.frame(outcome = log(d$duration_months), truth = log(d$reference_n), lab_value = log(d$recorded_n),
      phase = factor(d$phase, levels = c("not_applicable", "early", "phase2", "phase3_4")),
      sponsor = factor(d$sponsor, levels = c("other", "industry")), log_sites = log(d$sites), start_year = d$start_year)
  } else {
    d <- utils::read.csv(file.path(root, "validation/real-data/nhanes-pool.csv"), stringsAsFactors = FALSE)
    d <- d[d$split == config$split, ]
    data.frame(outcome = d$sbp, truth = d$reference_kg, lab_value = d$recorded_kg, age = d$age,
      sex = factor(d$sex, levels = c("female", "male")),
      race = factor(d$race, levels = c("white", "black", "mexican_american", "other_hispanic", "other")), height = d$height)
  }
}
pools <- list(ctgov = list(data = load_pool("ctgov"), covariates = c("phase", "sponsor", "log_sites", "start_year")),
              nhanes = list(data = load_pool("nhanes"), covariates = c("age", "sex", "race", "height")))
analysis_formula <- function(covariates) stats::reformulate(c("lab_value", covariates), response = "outcome")

# Synthetic boundary cells reuse the contamination fixtures with new seed streams. Only the
# stream bases change, so the generators and their error mechanisms stay as they were.
fixture_env <- new.env(parent = globalenv())
sys.source(file.path(root, "tools/lab-fixture.R"), envir = fixture_env)
sys.source(file.path(root, "tools/contamination-fixtures.R"), envir = fixture_env)
assign("contamination_stream_bases", function(phase) list(selection = config$fixture_selection_base, noise = config$fixture_noise_base,
  mechanism = config$fixture_mechanism_base, imputation = config$fixture_imputation_base), envir = fixture_env)
synthetic_cohort <- function(cell, r) {
  s <- list(fixture = "missed_flags", mechanism = if (cell == "missed_flags") "none" else cell, fraction = .25,
            selection = "random", noise = 0, agreement_tol = 0)
  p <- fixture_env$make_contamination_case(s, r, config$fixture_phase, config$fixture_master_seed, fixture_env$make_lab_fixture)
  d <- p$input[c("age", "sex", "inflammation", "SITEID", "outcome", "lab_value")]
  d$truth <- p$truth
  d
}

# The whole-pool values a pilot estimates, with every row treated as audited.
pool_table <- do.call(rbind, lapply(names(pools), function(name) {
  z <- pools[[name]]
  fit <- stats::lm(analysis_formula(z$covariates), data = transform(z$data, lab_value = truth))
  values <- catmech_audit_pool_k(z$data, "lab_value", "truth", c(z$covariates, "outcome"), analysis_formula(z$covariates), "lab_value")
  data.frame(dataset = name, split = config$split, pool_rows = nrow(z$data), target_coefficient = unname(stats::coef(fit)["lab_value"]),
    k = values[["k"]], discard_k = values[["discard_k"]], classical = values[["classical"]], discard_classical = values[["discard_classical"]])
}))
utils::write.csv(pool_table, file.path(outdir, "pool-k.csv"), row.names = FALSE)

classify <- function(e) {
  if (grepl("full column rank|more observations than coefficients|lost design support|rank deficient|constant predictor|Remove constant|needs at least|every coefficient|estimable|zero residual variance",
            conditionMessage(e), ignore.case = TRUE)) "refused" else "error"
}
capture <- function(expr) {
  warnings <- character()
  value <- tryCatch(withCallingHandlers(expr, warning = function(w) {
    warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning")
  }), error = identity)
  list(value = value, warnings = paste(unique(warnings), collapse = " | "))
}
fit_target <- function(d, truth, rows, predictors, formula, seed) {
  input <- d
  input$gold <- NA_real_
  input$gold[rows] <- truth[rows]
  run <- capture({
    fit <- catmech_reference(input, "lab_value", "gold", predictors, use_recorded = TRUE, m = config$m, seed = seed)
    row <- catmech_pool_lm(fit, formula)
    row <- row[row$term == "lab_value", ]
    data.frame(status = "ok", message = "", estimate = row$estimate, lower = row$lower, upper = row$upper, width = row$upper - row$lower,
               within_variance = row$within_variance, between_variance = row$between_variance, df = row$df, recorded_fit = fit$spec$recorded_fit)
  })
  if (inherits(run$value, "error"))
    data.frame(status = classify(run$value), message = conditionMessage(run$value), estimate = NA_real_, lower = NA_real_, upper = NA_real_, width = NA_real_,
               within_variance = NA_real_, between_variance = NA_real_, df = NA_real_, recorded_fit = NA_character_, warnings = run$warnings) else
    cbind(run$value, warnings = run$warnings)
}
run_cohort <- function(cell, r) {
  offset <- config$offsets[[cell]]
  if (cell %in% names(pools)) {
    z <- pools[[cell]]
    set.seed(as.integer(config$cohort_base + offset + r))
    d <- z$data[sample.int(nrow(z$data), config$n, replace = TRUE), ]
    covariates <- z$covariates
  } else {
    d <- synthetic_cohort(cell, r)
    covariates <- c("age", "sex", "inflammation", "SITEID")
  }
  rownames(d) <- NULL
  set.seed(as.integer(config$order_base + offset + r))
  order_rows <- sample.int(config$n)
  truth <- d$truth
  d$truth <- NULL
  predictors <- c(covariates, "outcome")
  formula <- analysis_formula(covariates)
  oracle <- stats::lm(formula, data = transform(d, lab_value = truth))
  full_width <- unname(diff(stats::confint(oracle)["lab_value", ]))
  targets <- do.call(rbind, lapply(seq_along(config$target_records), function(i) {
    rows <- order_rows[seq_len(config$target_records[i])]
    out <- fit_target(d, truth, rows, predictors, formula, as.integer(config$imputation_base + offset + 2L * r))
    retest <- if (r <= config$retest_reps) fit_target(d, truth, rows, predictors, formula, as.integer(config$retest_base + offset + 2L * r))$width else NA_real_
    cbind(data.frame(cell = cell, replicate = r, target_fraction = fractions[i], target_records = config$target_records[i], full_width = full_width), out, retest_width = retest)
  }))
  plans <- do.call(rbind, lapply(config$pilot_records, function(size) {
    rows <- order_rows[seq_len(size)]
    input <- d
    input$gold <- NA_real_
    input$gold[rows] <- truth[rows]
    seed <- as.integer((if (size == 60L) config$plan_10_base else config$plan_25_base) + offset + 4L * r)
    start <- Sys.time()
    run <- capture(catmech_audit_plan(input, "lab_value", "gold", predictors, formula, "lab_value", fractions = fractions,
      target_ratio = config$target_ratio, bootstrap = config$bootstrap, m = config$m, seed = seed))
    elapsed <- as.numeric(difftime(Sys.time(), start, units = "secs"))
    base <- data.frame(cell = cell, replicate = r, pilot_records = size, pilot_fraction = size / config$n, target_fraction = fractions, plan_seed = seed, elapsed = elapsed)
    if (inherits(run$value, "error")) {
      return(cbind(base, status = classify(run$value), message = c(conditionMessage(run$value), rep("", nrow(base) - 1L)), warnings = c(run$warnings, rep("", nrow(base) - 1L)), flags = "",
        k = NA_real_, k_lower = NA_real_, k_upper = NA_real_, classical = NA_real_, rubin_check = NA_real_, discard_k = NA_real_, discard_lower = NA_real_, discard_upper = NA_real_,
        relative_width = NA_real_, relative_lower = NA_real_, relative_upper = NA_real_, predictive_lower = NA_real_, predictive_upper = NA_real_,
        predicted_width = NA_real_, width_lower = NA_real_, width_upper = NA_real_, discard_relative_width = NA_real_, discard_width = NA_real_,
        pilot_width = NA_real_, pilot_full_width = NA_real_, disagreements = NA_real_, bootstrap_usable = NA_real_, answer_records = NA_real_, smallest_grid_fraction = NA_real_))
    }
    plan <- run$value
    tb <- plan$table
    cbind(base, status = "ok", message = "", warnings = c(run$warnings, rep("", nrow(base) - 1L)), flags = paste(plan$flags, collapse = "|"),
      k = plan$k$estimate[1], k_lower = plan$k$lower[1], k_upper = plan$k$upper[1], classical = plan$k$classical[1], rubin_check = plan$k$rubin_check[1],
      discard_k = plan$k$estimate[2], discard_lower = plan$k$lower[2], discard_upper = plan$k$upper[2],
      relative_width = tb$relative_width, relative_lower = tb$relative_lower, relative_upper = tb$relative_upper,
      predictive_lower = tb$predictive_lower, predictive_upper = tb$predictive_upper, predicted_width = tb$width, width_lower = tb$width_lower, width_upper = tb$width_upper,
      discard_relative_width = tb$discard_relative_width, discard_width = tb$discard_width, pilot_width = plan$pilot$width, pilot_full_width = plan$pilot$full_width,
      disagreements = plan$pilot$disagreements, bootstrap_usable = plan$pilot$bootstrap_usable, answer_records = plan$answer$records, smallest_grid_fraction = plan$answer$smallest_grid_fraction)
  }))
  list(targets = targets, plans = plans)
}

cells <- list(list(cell = "ctgov", reps = config$real_reps), list(cell = "nhanes", reps = config$real_reps),
              list(cell = "missed_flags", reps = config$synthetic_reps), list(cell = "error_gross_and_small", reps = config$synthetic_reps))
target_rows <- list(); plan_rows <- list(); start <- Sys.time()
for (s in cells) {
  for (chunk in split(seq_len(s$reps), ceiling(seq_len(s$reps) / 140L))) {
    batch <- parallel::mclapply(chunk, function(r) run_cohort(s$cell, r), mc.cores = workers, mc.preschedule = TRUE, mc.set.seed = FALSE)
    if (any(vapply(batch, inherits, logical(1), what = "try-error")))
      stop("A worker failed; checkpoint is not verified evidence: ", as.character(batch[[which(vapply(batch, inherits, logical(1), what = "try-error"))[1]]]))
    target_rows <- c(target_rows, lapply(batch, `[[`, "targets")); plan_rows <- c(plan_rows, lapply(batch, `[[`, "plans"))
    saveRDS(list(targets = target_rows, plans = plan_rows), file.path(outdir, "checkpoint.rds"))
    cat(s$cell, max(chunk), "/", s$reps, "; elapsed", round(as.numeric(difftime(Sys.time(), start, units = "secs"))), "seconds\n"); flush.console()
  }
}
if (!identical(hashes_before, hash_sources())) stop("Source, plan or pool changed during run. Checkpoint retained; evidence is not verified.")
targets <- do.call(rbind, target_rows)
plans <- do.call(rbind, plan_rows)
targets$ratio <- targets$width / targets$full_width
utils::write.csv(targets, file.path(outdir, "targets.csv"), row.names = FALSE)
utils::write.csv(plans, file.path(outdir, "plans.csv"), row.names = FALSE)

# Pairs each cohort's prediction with its actual width. A failed target has no actual width and
# leaves the cell, with a count; a failed plan where the target succeeded counts as infinite error.
target_key <- targets[c("cell", "replicate", "target_fraction", "status", "width", "full_width", "ratio", "retest_width")]
names(target_key)[names(target_key) == "status"] <- "target_status"
joined <- merge(plans, target_key, by = c("cell", "replicate", "target_fraction"))
joined <- joined[abs(joined$target_fraction - joined$pilot_fraction) > 1e-9, ]
joined <- joined[order(match(joined$cell, names(config$offsets)), joined$replicate, joined$pilot_records, joined$target_fraction), ]
order_statistic_interval <- function(x) {
  x <- sort(x); n <- length(x)
  low <- max(1L, floor(n / 2 - 1.96 * sqrt(n) / 2)); high <- min(n, ceiling(n / 2 + 1 + 1.96 * sqrt(n) / 2))
  c(x[low], x[high])
}
summary_table <- do.call(rbind, lapply(split(joined, interaction(joined$cell, joined$pilot_records, joined$target_fraction, drop = TRUE)), function(g) {
  usable <- g[g$target_status == "ok", ]
  error <- ifelse(usable$status == "ok", abs(usable$predicted_width - usable$width) / usable$width, Inf)
  signed <- (usable$predicted_width - usable$width) / usable$width
  planned <- usable[usable$status == "ok", ]
  ci <- order_statistic_interval(error)
  center <- stats::median(usable$width)
  retest <- usable[is.finite(usable$retest_width), ]
  data.frame(cell = g$cell[1], pilot_records = g$pilot_records[1], target_fraction = g$target_fraction[1], requested = nrow(g),
    target_failures = sum(g$target_status != "ok"), plan_failures = sum(usable$status != "ok"), compared = nrow(usable),
    median_error = stats::median(error), median_error_lower = ci[1], median_error_upper = ci[2],
    p90_error = unname(stats::quantile(error, .9, type = 7)), median_signed_error = stats::median(signed[usable$status == "ok"]),
    interval_coverage = mean(planned$ratio >= planned$relative_lower & planned$ratio <= planned$relative_upper),
    predictive_coverage = mean(planned$ratio >= planned$predictive_lower & planned$ratio <= planned$predictive_upper),
    floor_median_error = stats::median(abs(center - usable$width) / usable$width),
    retest_median_error = if (nrow(retest)) stats::median(abs(retest$retest_width - retest$width) / retest$width) else NA_real_,
    mean_actual_ratio = mean(usable$ratio), median_actual_ratio = stats::median(usable$ratio), mean_predicted_ratio = mean(planned$relative_width))
}))
rownames(summary_table) <- NULL
utils::write.csv(summary_table, file.path(outdir, "summary.csv"), row.names = FALSE)

# How often the planner picks the truth's smallest fraction for within 25% of the full-audit width.
truth_pick <- do.call(rbind, lapply(split(targets, interaction(targets$cell, targets$replicate, drop = TRUE)), function(g) {
  reached <- g$target_fraction[g$status == "ok" & g$ratio <= config$target_ratio]
  data.frame(cell = g$cell[1], replicate = g$replicate[1], true_pick = if (length(reached)) min(reached) else NA_real_)
}))
plan_pick <- plans[!duplicated(plans[c("cell", "replicate", "pilot_records")]), c("cell", "replicate", "pilot_records", "status", "smallest_grid_fraction")]
picks <- merge(plan_pick, truth_pick, by = c("cell", "replicate"))
picks$same <- picks$status == "ok" & ((is.na(picks$true_pick) & is.na(picks$smallest_grid_fraction)) |
  (!is.na(picks$true_pick) & !is.na(picks$smallest_grid_fraction) & abs(picks$true_pick - picks$smallest_grid_fraction) < 1e-9))
label <- function(x) ifelse(is.na(x), "none", format(x, trim = TRUE))
picks_table <- do.call(rbind, lapply(split(picks, interaction(picks$cell, picks$pilot_records, drop = TRUE)), function(g) {
  data.frame(cell = g$cell[1], pilot_records = g$pilot_records[1], cohorts = nrow(g), same_pick = mean(g$same),
    plan_refused = sum(g$status != "ok"), truth_05 = mean(g$true_pick %in% .05), truth_10 = mean(g$true_pick %in% .10), truth_25 = mean(g$true_pick %in% .25),
    truth_50 = mean(g$true_pick %in% .50), truth_none = mean(is.na(g$true_pick)),
    plan_05 = mean(g$smallest_grid_fraction %in% .05), plan_10 = mean(g$smallest_grid_fraction %in% .10), plan_25 = mean(g$smallest_grid_fraction %in% .25),
    plan_50 = mean(g$smallest_grid_fraction %in% .50), plan_none = mean(g$status == "ok" & is.na(g$smallest_grid_fraction)))
}))
rownames(picks_table) <- NULL
utils::write.csv(picks_table, file.path(outdir, "fraction-picks.csv"), row.names = FALSE)

# Pilot estimates of k against the whole-pool value, and how often each flag was raised.
first <- plans[!duplicated(plans[c("cell", "replicate", "pilot_records")]), ]
pool_lookup <- function(cell, field) if (cell %in% pool_table$dataset) pool_table[[field]][pool_table$dataset == cell] else NA_real_
k_table <- do.call(rbind, lapply(split(first, interaction(first$cell, first$pilot_records, drop = TRUE)), function(g) {
  ok <- g[g$status == "ok", ]
  pool <- pool_lookup(g$cell[1], "k")
  q <- function(v, p) unname(stats::quantile(v, p, type = 7))
  data.frame(cell = g$cell[1], pilot_records = g$pilot_records[1], cohorts = nrow(g), plans_returned = nrow(ok), pool_k = pool, pool_classical = pool_lookup(g$cell[1], "classical"),
    pool_discard_k = pool_lookup(g$cell[1], "discard_k"),
    k_mean = mean(ok$k), k_median = stats::median(ok$k), k_q05 = q(ok$k, .05), k_q25 = q(ok$k, .25), k_q75 = q(ok$k, .75), k_q95 = q(ok$k, .95),
    interval_covers_pool = if (is.na(pool)) NA_real_ else mean(ok$k_lower <= pool & pool <= ok$k_upper),
    median_interval_lower = stats::median(ok$k_lower), median_interval_upper = stats::median(ok$k_upper),
    classical_median = stats::median(ok$classical), rubin_check_median = stats::median(ok$rubin_check),
    discard_k_median = stats::median(ok$discard_k), discard_k_q05 = q(ok$discard_k, .05), discard_k_q95 = q(ok$discard_k, .95),
    median_disagreements = stats::median(ok$disagreements), bootstrap_usable_min = min(ok$bootstrap_usable),
    flag_small_pilot = mean(grepl("small_pilot", ok$flags)), flag_few_disagreements = mean(grepl("few_disagreements", ok$flags)),
    flag_no_disagreements = mean(grepl("no_disagreements", ok$flags)), flag_gross_errors = mean(grepl("gross_errors", ok$flags)),
    flag_selective = mean(grepl("pilot_looks_selective", ok$flags)), flag_bootstrap = mean(grepl("bootstrap_lost_support", ok$flags)))
}))
rownames(k_table) <- NULL
utils::write.csv(k_table, file.path(outdir, "k-summary.csv"), row.names = FALSE)

# The verdict rule from the plan, on the real datasets only.
holds <- function(size) all(vapply(c("ctgov", "nhanes"), function(name) {
  rows <- summary_table[summary_table$cell == name & summary_table$pilot_records == size, ]
  nrow(rows) == 3L && all(rows$median_error <= .10)
}, logical(1)))
verdict <- if (holds(60L)) "GO" else if (holds(150L)) "PARTIAL" else "NO-GO"
utils::write.csv(data.frame(verdict = verdict, holds_10pct_pilot = holds(60L), holds_25pct_pilot = holds(150L)), file.path(outdir, "verdict.csv"), row.names = FALSE)
status_table <- rbind(
  cbind(kind = "target", as.data.frame(table(cell = targets$cell, status = targets$status), stringsAsFactors = FALSE)),
  cbind(kind = "plan", as.data.frame(table(cell = first$cell, status = first$status), stringsAsFactors = FALSE)))
utils::write.csv(status_table, file.path(outdir, "request-status.csv"), row.names = FALSE)
failures <- rbind(
  cbind(kind = "target", targets[targets$status != "ok", c("cell", "replicate", "status", "message")]),
  cbind(kind = "plan", first[first$status != "ok", c("cell", "replicate", "status", "message")]))
utils::write.csv(failures, file.path(outdir, "failures.csv"), row.names = FALSE)
utils::write.csv(hashes_before, file.path(outdir, "source-hashes.csv"), row.names = FALSE)
writeLines(c(paste("Completed UTC:", format(Sys.time(), tz = "UTC", usetz = TRUE)),
  paste("Elapsed seconds:", round(as.numeric(difftime(Sys.time(), start, units = "secs")), 1)), capture.output(sessionInfo())), file.path(outdir, "session-info.txt"))
unlink(file.path(outdir, "checkpoint.rds"))
cat("Completed", mode, ": verdict", verdict, "\n")
print(summary_table[c("cell", "pilot_records", "target_fraction", "compared", "target_failures", "plan_failures", "median_error", "p90_error", "floor_median_error", "interval_coverage", "predictive_coverage")], digits = 3, row.names = FALSE)
cat("Output:", normalizePath(outdir), "\n")
