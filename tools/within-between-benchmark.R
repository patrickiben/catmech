#!/usr/bin/env Rscript
# Prospective test of separate within-cluster and between-cluster slopes in the option; settings locked in validation/WITHIN_BETWEEN_PLAN.md.
# Guard cells reuse the earlier studies' cohorts, so the candidate pairs with their saved results. SNAP and the new generator's cells run every method.
# Usage: Rscript tools/within-between-benchmark.R control|pilot|study [NEW_DIRECTORY] [PILOT_REPS]
args <- commandArgs(trailingOnly = TRUE)
mode <- if (length(args)) match.arg(args[1], c("control", "pilot", "study")) else "pilot"
file_arg <- grep("^--file=", commandArgs(), value = TRUE)
root <- normalizePath(file.path(dirname(sub("^--file=", "", file_arg[1])), ".."))
outdir <- if (length(args) >= 2L) args[2] else file.path(root, "validation", paste0("within-between-", mode))
if (dir.exists(outdir) && length(list.files(outdir, all.files = TRUE, no.. = TRUE)))
  stop("Use a new or empty output directory; existing evidence is never overwritten.")
hash_sources <- function() {
  paths <- c(sort(list.files(file.path(root, "R"), "\\.R$", full.names = TRUE)),
    file.path(root, c("DESCRIPTION", "NAMESPACE", "tools/within-between-benchmark.R", "tools/within-between-score.py", "tools/within-between-prepare.py", "tools/independent-setups/ind_mice.R", "tools/variance-mi-benchmark.R", "tools/clustered-fix-benchmark.R", "tools/second-system-benchmark.R",
                      "validation/WITHIN_BETWEEN_PLAN.md", "validation/within-between/snap-pool.csv", "validation/within-between/generator/generator.R", "validation/within-between/generator/grid.csv",
                      "validation/real-audit/rwanda-pool.csv", "validation/second-system/zambia-pool.csv", "validation/variance-mi-study/replicates.csv.gz", "validation/second-system-study/replicates.csv.gz", "validation/independent-setups-study/replicates.csv.gz")))
  data.frame(file = substring(paths, nchar(root) + 2L),
    sha256 = vapply(paths, function(p) digest::digest(file = p, algo = "sha256"), character(1)), row.names = NULL)
}
hashes_before <- hash_sources()
pkgload::load_all(root, quiet = TRUE)
suppressPackageStartupMessages(library(mice))
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
pilot <- mode == "pilot"
pilot_reps <- if (length(args) >= 3L) as.integer(args[3]) else 3L
M <- 20L
config <- jsonlite::read_json(file.path(root, "validation/clustered-fix-study/settings.json"), simplifyVector = TRUE)
earlier <- jsonlite::read_json(file.path(root, "validation/variance-mi-study/settings.json"), simplifyVector = TRUE)
study_config <- config
if (pilot) {
  config$held_out <- list(cohort = 910100000, selection = 920100000, imputation = 930000000)
  config$development <- list(cohort = 940100000, selection = 950100000, imputation = 960000000)
  config$simulation <- list(cohort = 970100000, imputation = 980000000)
}
SEEDS <- if (pilot) list(covariate_cohort = 990100000, cohort = 1110000000, audit = 1120000000, snap_cohort = 1500000000, snap_audit = 1510000000, fresh = 1520000000) else
  list(covariate_cohort = earlier$seeds$covariate_cohort, cohort = 2010000000, audit = 2020000000, snap_cohort = 2100000000, snap_audit = 2110000000, fresh = 2120000000)
BASE <- if (pilot) c(guard = 300000000, linear_mi = 310000000, random_intercept = 320000000, within_between = 330000000, ind_mice = 340000000) else
  c(guard = 1400000000, linear_mi = 1300000000, random_intercept = 1310000000, within_between = 1320000000, ind_mice = 1330000000)
settings <- list(mode = mode, m = M, reps_real_25 = if (pilot) pilot_reps else 2000L, reps_other = if (pilot) pilot_reps else 500L, reps_zambia = if (pilot) pilot_reps else 1000L,
                 reps_snap = if (pilot) pilot_reps else 1000L, reps_fresh = if (pilot) pilot_reps else 500L, base = as.list(BASE), seeds = SEEDS)
jsonlite::write_json(settings, file.path(outdir, "settings.json"), pretty = TRUE, auto_unbox = TRUE)

named_from <- function(path) {
  e <- parse(path)
  list(e = e, names = vapply(e, function(x) if (is.call(x) && identical(x[[1]], as.name("<-")) && is.name(x[[2]])) as.character(x[[2]]) else "", character(1)))
}
fix <- named_from(file.path(root, "tools/clustered-fix-benchmark.R"))
for (e in fix$e[which(fix$names == "rwanda")[1]:which(fix$names == "real_case")[1]]) eval(e, envir = globalenv())
variance <- named_from(file.path(root, "tools/variance-mi-benchmark.R"))
eval(variance$e[[which(variance$names == "simulate_covariate_cell")[1]]], envir = globalenv())
second <- named_from(file.path(root, "tools/second-system-benchmark.R"))
for (n in c("draw_cohort", "draw_audit", "AUDITS", "ITEMS")) eval(second$e[[which(second$names == n)[1]]], envir = globalenv())
ind_env <- new.env(parent = globalenv()); sys.source(file.path(root, "tools/independent-setups/ind_mice.R"), envir = ind_env)
fresh_env <- new.env(parent = globalenv()); sys.source(file.path(root, "validation/within-between/generator/generator.R"), envir = fresh_env)
fresh_grid <- utils::read.csv(file.path(root, "validation/within-between/generator/grid.csv"), stringsAsFactors = FALSE)

# The candidate: the option, with cluster means added to its predictors.
candidate_means <- function(p) {
  d <- p$input; ids <- d$cluster_id
  columns <- c("lab_value", setdiff(p$predictors, "stratum"))
  means <- list()
  for (v in columns) {
    centered <- stats::ave(d[[v]], ids)
    if (max(abs(d[[v]] - centered)) > 1e-10 * (1 + max(abs(d[[v]])))) means[[paste0(".mean_", v)]] <- centered
  }
  first <- !duplicated(ids) & !is.na(ids)
  seen <- ids %in% unique(ids[p$audited])
  rows <- which(first & seen)
  base <- matrix(1, length(rows), 1L)
  if ("stratum" %in% p$predictors) base <- cbind(base, stats::model.matrix(~ stratum - 1, data = droplevels(d[rows, , drop = FALSE])))
  for (v in setdiff(columns, "lab_value")) if (is.null(means[[paste0(".mean_", v)]])) base <- cbind(base, d[[v]][rows])
  if (!length(means)) return(list())
  Z <- cbind(base, do.call(cbind, lapply(means, function(z) z[rows])))
  q <- qr(Z, tol = 1e-7)
  kept <- q$pivot[seq_len(q$rank)]
  kept_means <- names(means)[kept[kept > ncol(base)] - ncol(base)]
  if (length(rows) - q$rank < 3L) return(list())
  means[kept_means]
}
as_row <- function(estimate, se, df, target, fallback = NA) {
  half <- stats::qt(.975, df) * se
  list(estimate = estimate, se = se, df = df, width = 2 * half, cover = as.numeric(abs(estimate - target) <= half), fallback = fallback)
}
one_value <- function(p, values) {
  d <- p$input; d$lab_value <- values
  term <- if (is.null(p$term)) "covariate" else p$term
  fit <- stats::lm(p$formula, data = d, na.action = stats::na.fail)
  vc <- catmech_cluster_jackknife(fit, d$cluster_id)
  j <- match(term, names(stats::coef(fit)))
  as_row(unname(stats::coef(fit)[j]), sqrt(vc[j, j]), length(unique(d$cluster_id)) - 1L, p$target)
}
imputed <- function(p, predictors, seed, cluster, input = p$input, fallback = NA) {
  term <- if (is.null(p$term)) "covariate" else p$term
  options <- list(input, "lab_value", "gold_value", predictors, method = "normal", m = M, seed = seed, use_recorded = TRUE)
  if (cluster) options$cluster <- "cluster_id"
  fit <- do.call(catmech_reference, options)
  row <- catmech_pool_lm(fit, p$formula, cluster = "cluster_id")
  row <- row[row$term == term, ]
  as_row(row$estimate, row$std.error, row$df, p$target, fallback)
}
run_method <- function(method, p, seed, fraction) {
  switch(method,
    oracle = one_value(p, p$truth),
    recorded = one_value(p, p$input$lab_value),
    linear_mi = imputed(p, p$predictors, seed, FALSE),
    random_intercept = imputed(p, p$predictors, seed, TRUE),
    within_between = {
      means <- candidate_means(p)
      input <- p$input
      for (n in names(means)) input[[n]] <- means[[n]]
      imputed(p, c(p$predictors, names(means)), seed, TRUE, input, fallback = as.numeric(!length(means)))
    },
    ind_mice = {
      got <- ind_env$fit_method(list(input = p$input, formula = p$formula, term = if (is.null(p$term)) "covariate" else p$term, audited = p$audited, fraction = fraction, predictors = p$predictors), seed)
      as_row(got$estimate, got$se, got$df, p$target)
    })
}
limited <- function(expr) {
  job <- parallel::mcparallel(tryCatch(expr, error = identity), silent = TRUE)
  got <- parallel::mccollect(job, wait = FALSE, timeout = 60)
  if (is.null(got)) {
    tools::pskill(job$pid, tools::SIGKILL); parallel::mccollect(job, wait = TRUE)
    return(simpleError("exceeded the 60 second limit"))
  }
  parallel::mccollect(job, wait = TRUE)
  out <- got[[1L]]
  if (is.null(out)) simpleError("the child returned nothing") else out
}
call_method <- function(method, p, seed, fraction, labels, r) {
  started <- Sys.time()
  got <- limited(suppressWarnings(suppressMessages(run_method(method, p, seed, fraction))))
  seconds <- as.numeric(difftime(Sys.time(), started, units = "secs"))
  bad <- inherits(got, "error") || !all(is.finite(c(got$estimate, got$se, got$df))) || got$se <= 0 || got$df <= 0
  data.frame(labels, replicate = r, method = method, status = if (bad) "failed" else "ok", message = if (inherits(got, "error")) substr(conditionMessage(got), 1L, 100L) else "",
             estimate = if (bad) NA_real_ else got$estimate, se = if (bad) NA_real_ else got$se, df = if (bad) NA_real_ else got$df, width = if (bad) NA_real_ else got$width,
             cover = if (bad) NA_real_ else got$cover, target = p$target, fallback = if (bad || is.null(got$fallback)) NA_real_ else got$fallback, seconds = round(seconds, 4), stringsAsFactors = FALSE)
}
run_all <- function(n, f) {
  batch <- parallel::mclapply(seq_len(n), f, mc.cores = 14L, mc.preschedule = TRUE, mc.set.seed = FALSE)
  bad <- vapply(batch, inherits, logical(1), what = "try-error")
  if (any(bad)) stop("A worker failed: ", conditionMessage(attr(batch[[which(bad)[1]]], "condition")))
  do.call(rbind, batch)
}
write_gz <- function(x, name) { out <- gzfile(file.path(outdir, name), "wt"); utils::write.csv(x, out, row.names = FALSE); close(out) }

# Cells in the earlier studies' order, so cell numbers match them.
grid <- utils::read.csv(file.path(root, "validation/clustered-fix-study/simulation-grid.csv"), stringsAsFactors = FALSE)
covariate_grid <- expand.grid(design = "rwanda", errors = c("shift", "under", "row_only"), fraction = config$fractions, stringsAsFactors = FALSE)
covariate_grid$rows <- 3L
guard_cells <- list()
for (z in real) for (j in seq_along(config$fractions)) guard_cells[[length(guard_cells) + 1L]] <- list(kind = "real", z = z, fraction_index = j, level = z$level, analysis = z$id, fraction = config$fractions[j],
  reps = if (j == 1L) settings$reps_real_25 else settings$reps_other, number = 100L + length(guard_cells) + 1L)
for (i in seq_len(nrow(grid))) guard_cells[[length(guard_cells) + 1L]] <- list(kind = "simulation", spec = as.list(grid[i, ]), index = i, level = "simulation", analysis = grid$label[i], fraction = grid$fraction[i], reps = settings$reps_other, number = i)
for (i in seq_len(nrow(covariate_grid))) guard_cells[[length(guard_cells) + 1L]] <- list(kind = "covariate", spec = as.list(covariate_grid[i, ]), index = i, level = "imputed_covariate", analysis = paste0("covariate_", covariate_grid$errors[i]), fraction = covariate_grid$fraction[i], reps = settings$reps_other, number = 200L + i)
guard_cohort <- function(s, r) {
  if (s$kind == "real") real_case(s$z, s$fraction_index, r) else if (s$kind == "simulation")
    simulate_cell(s$spec, as.integer(config$simulation$cohort + s$index * 1e5 + r), s$fraction) else
    simulate_covariate_cell(s$spec, as.integer(SEEDS$covariate_cohort + s$index * 1e5 + r), s$fraction)
}

# Control: rerunning the saved ind_mice requests must reproduce them, which shows that the guard cohorts are the same.
if (mode != "pilot") {
  con <- gzfile(file.path(root, "validation/independent-setups-study/replicates.csv.gz"), "rt")
  saved <- utils::read.csv(con, stringsAsFactors = FALSE); saved <- saved[saved$method == "ind_mice" & saved$fraction == .25 & saved$replicate <= 10L & saved$level %in% c("held_out", "development", "unused"), ]
  real_cells <- Filter(function(s) s$kind == "real" && s$fraction_index == 1L, guard_cells)
  replay <- do.call(rbind, lapply(real_cells, function(s) run_all(10L, function(r) {
    p <- guard_cohort(s, r)
    call_method("ind_mice", p, as.integer(1260000000 + s$number * 1e5 + r), s$fraction, data.frame(level = s$level, analysis = s$analysis, fraction = s$fraction), r)
  })))
  both <- merge(replay, saved, by = c("analysis", "replicate"), suffixes = c("", ".saved"))
  stopifnot(nrow(both) == 110L, all(abs(both$estimate - both$estimate.saved) < 1e-8), all(abs(both$width - both$width.saved) < 1e-8))
  writeLines(sprintf("Reran %d saved ind_mice requests on real items at 25%%. Largest absolute difference in estimate: %.2e; in width: %.2e.", nrow(both), max(abs(both$estimate - both$estimate.saved)), max(abs(both$width - both$width.saved))), file.path(outdir, "control.txt"))
  cat("Control passed\n")
  if (mode == "control") quit(save = "no")
}

start <- Sys.time()
guard <- do.call(rbind, lapply(guard_cells, function(s) {
  got <- run_all(s$reps, function(r) {
    p <- guard_cohort(s, r)
    call_method("within_between", p, as.integer(BASE[["guard"]] + s$number * 1e5 + r), s$fraction, data.frame(level = s$level, analysis = s$analysis, fraction = s$fraction), r)
  })
  cat("guard", s$analysis, s$fraction, "; elapsed", round(as.numeric(difftime(Sys.time(), start, units = "secs"))), "seconds\n"); flush.console()
  got
}))
write_gz(guard, "guard.csv.gz")

pool_all <- utils::read.csv(file.path(root, "validation/second-system/zambia-pool.csv"), stringsAsFactors = FALSE)
zambia <- do.call(rbind, lapply(seq_along(ITEMS), function(i) {
  pool <- pool_all[pool_all$item == ITEMS[i] & pool_all$both_zero == 0L, ]
  target <- unname(stats::coef(stats::lm(verified ~ covariate, data = pool))[2])
  out <- do.call(rbind, lapply(seq_along(AUDITS), function(a) {
    run_all(settings$reps_zambia, function(r) {
      co <- draw_cohort(pool, i, r)
      audited <- draw_audit(co, AUDITS[[a]], a, i, r)
      truth <- pool$verified[co$rows]; gold <- truth; gold[!audited] <- NA_real_
      p <- list(input = data.frame(covariate = pool$covariate[co$rows], cluster_id = co$ids, lab_value = pool$reported[co$rows], gold_value = gold), truth = truth, target = target, formula = lab_value ~ covariate, term = "covariate", audited = audited, predictors = "covariate")
      call_method("within_between", p, as.integer(BASE[["guard"]] + (400L + i * 10L + a) * 1e5 + r), .25, data.frame(item = ITEMS[i], audit = names(AUDITS)[a]), r)
    })
  }))
  cat("zambia", ITEMS[i], "; elapsed", round(as.numeric(difftime(Sys.time(), start, units = "secs"))), "seconds\n"); flush.console()
  out
}))
write_gz(zambia, "zambia.csv.gz")

# Deciding data.
METHODS <- c("oracle", "recorded", "linear_mi", "random_intercept", "within_between", "ind_mice")
snap <- utils::read.csv(file.path(root, "validation/within-between/snap-pool.csv"), stringsAsFactors = FALSE)
rows_by_state <- split(seq_len(nrow(snap)), snap$state)
snap$balance <- 1 / lengths(rows_by_state)[as.character(snap$state)]
SNAP_ANALYSES <- list(
  benefit_on_size = list(role = "outcome", target = unname(stats::coef(stats::lm(reference ~ size, data = snap, weights = balance))["size"]), covariate = "size"),
  benefit_on_gross = list(role = "outcome", target = unname(stats::coef(stats::lm(reference ~ gross, data = snap, weights = balance))["gross"]), covariate = "gross"),
  gross_on_benefit = list(role = "predictor", target = unname(stats::coef(stats::lm(gross ~ reference + size, data = snap, weights = balance))["reference"]), covariate = "size"))
utils::write.csv(data.frame(analysis = names(SNAP_ANALYSES), target = vapply(SNAP_ANALYSES, `[[`, numeric(1), "target")), file.path(outdir, "snap-targets.csv"), row.names = FALSE)
SNAP_K <- c(10L, 40L); SNAP_AUDITED <- c(13L, 26L)
snap_case <- function(a, k_index, audit_index, r) {
  spec <- SNAP_ANALYSES[[a]]
  set.seed(as.integer(SEEDS$snap_cohort + ((a - 1L) * 2L + k_index) * 1e6 + r))
  picked <- sample.int(length(rows_by_state), length(rows_by_state), replace = TRUE)
  rows <- unlist(lapply(picked, function(s) { u <- rows_by_state[[s]]; u[sample.int(length(u), SNAP_K[k_index])] }), use.names = FALSE)
  ids <- rep(seq_along(picked), each = SNAP_K[k_index])
  set.seed(as.integer(SEEDS$snap_audit + ((a - 1L) * 2L + k_index) * 1e6 + audit_index * 1e5 + r))
  audited <- ids %in% sample.int(length(picked), SNAP_AUDITED[audit_index])
  d <- snap[rows, ]
  gold <- d$reference; gold[!audited] <- NA_real_
  input <- data.frame(covariate = d[[spec$covariate]], cluster_id = ids, lab_value = d$recorded, gold_value = gold)
  if (spec$role == "predictor") {
    input$outcome <- d$gross
    list(input = input, truth = d$reference, target = spec$target, audited = audited, formula = outcome ~ lab_value + covariate, predictors = c("covariate", "outcome"), term = "lab_value")
  } else {
    list(input = input, truth = d$reference, target = spec$target, audited = audited, formula = lab_value ~ covariate, predictors = "covariate", term = "covariate")
  }
}
results <- list()
for (a in seq_along(SNAP_ANALYSES)) for (k in seq_along(SNAP_K)) for (u in seq_along(SNAP_AUDITED)) {
  number <- ((a - 1L) * 2L + k - 1L) * 2L + u
  got <- run_all(settings$reps_snap, function(r) {
    p <- snap_case(a, k, u, r)
    do.call(rbind, lapply(METHODS, function(m) call_method(m, p, if (m %in% names(BASE)) as.integer(BASE[[m]] + number * 1e5 + r) else 1L, SNAP_AUDITED[u] / 53,
      data.frame(source = "snap", analysis = paste0(names(SNAP_ANALYSES)[a], "_k", SNAP_K[k]), role = SNAP_ANALYSES[[a]]$role, fraction = c(.25, .5)[u]), r)))
  })
  results[[length(results) + 1L]] <- got
  cat("snap", names(SNAP_ANALYSES)[a], SNAP_K[k], SNAP_AUDITED[u], "; elapsed", round(as.numeric(difftime(Sys.time(), start, units = "secs"))), "seconds\n"); flush.console()
}
for (i in seq_len(nrow(fresh_grid))) for (u in 1:2) {
  fraction <- c(.25, .5)[u]; number <- 100L + (i - 1L) * 2L + u
  got <- run_all(settings$reps_fresh, function(r) {
    p <- fresh_env$simulate_cell(as.list(fresh_grid[i, ]), as.integer(SEEDS$fresh + number * 1e5 + r), fraction)
    do.call(rbind, lapply(METHODS, function(m) call_method(m, p, if (m %in% names(BASE)) as.integer(BASE[[m]] + number * 1e5 + r) else 1L, fraction,
      data.frame(source = "fresh", analysis = fresh_grid$label[i], role = fresh_grid$role[i], fraction = fraction), r)))
  })
  results[[length(results) + 1L]] <- got
  cat("fresh", fresh_grid$label[i], fraction, "; elapsed", round(as.numeric(difftime(Sys.time(), start, units = "secs"))), "seconds\n"); flush.console()
}
results <- do.call(rbind, results)
write_gz(results, "deciding.csv.gz")
if (!identical(hashes_before, hash_sources())) stop("Source, plan, pool or reused results changed during the run; evidence is not verified.")
utils::write.csv(hashes_before, file.path(outdir, "source-hashes.csv"), row.names = FALSE)
writeLines(c(paste("Completed UTC:", format(Sys.time(), tz = "UTC", usetz = TRUE)), paste("Elapsed seconds:", round(as.numeric(difftime(Sys.time(), start, units = "secs")), 1)), capture.output(sessionInfo())), file.path(outdir, "session-info.txt"))
cat("Completed", mode, ":", nrow(guard), "guard,", nrow(zambia), "Zambia,", nrow(results), "deciding; statuses", paste(names(table(results$status)), table(results$status), collapse = ", "), "\n")
