#!/usr/bin/env Rscript
# Prospective test of comparators set up by fresh agents; settings locked in validation/INDEPENDENT_SETUPS_PLAN.md.
# The setups in tools/independent-setups/ are the agents' code, unchanged. This harness only hands them cohorts.
# Usage: Rscript tools/independent-setups-benchmark.R pilot|study [NEW_DIRECTORY] [PILOT_REPS]
args <- commandArgs(trailingOnly = TRUE)
mode <- if (length(args)) match.arg(args[1], c("pilot", "study")) else "pilot"
file_arg <- grep("^--file=", commandArgs(), value = TRUE)
root <- normalizePath(file.path(dirname(sub("^--file=", "", file_arg[1])), ".."))
outdir <- if (length(args) >= 2L) args[2] else file.path(root, "validation", paste0("independent-setups-", mode))
if (dir.exists(outdir) && length(list.files(outdir, all.files = TRUE, no.. = TRUE)))
  stop("Use a new or empty output directory; existing evidence is never overwritten.")
METHODS <- c("ind_mice", "ind_survey", "ind_design")
hash_sources <- function() {
  paths <- c(sort(list.files(file.path(root, "R"), "\\.R$", full.names = TRUE)),
    file.path(root, c("DESCRIPTION", "NAMESPACE", "tools/independent-setups-benchmark.R", "tools/independent-setups-score.py", "tools/variance-mi-benchmark.R", "tools/clustered-fix-benchmark.R", "tools/second-system-benchmark.R",
                      file.path("tools/independent-setups", paste0(METHODS, ".R")), "validation/INDEPENDENT_SETUPS_PLAN.md", "validation/real-audit/rwanda-pool.csv", "validation/second-system/zambia-pool.csv", "validation/variance-mi-study/replicates.csv.gz")))
  data.frame(file = substring(paths, nchar(root) + 2L),
    sha256 = vapply(paths, function(p) digest::digest(file = p, algo = "sha256"), character(1)), row.names = NULL)
}
hashes_before <- hash_sources()
pkgload::load_all(root, quiet = TRUE)
suppressPackageStartupMessages({library(survey); library(mice)})
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
pilot <- mode == "pilot"
pilot_reps <- if (length(args) >= 3L) as.integer(args[3]) else 3L
config <- jsonlite::read_json(file.path(root, "validation/clustered-fix-study/settings.json"), simplifyVector = TRUE)
earlier <- jsonlite::read_json(file.path(root, "validation/variance-mi-study/settings.json"), simplifyVector = TRUE)
if (pilot) {
  config$held_out <- list(cohort = 910100000, selection = 920100000, imputation = 930000000)
  config$development <- list(cohort = 940100000, selection = 950100000, imputation = 960000000)
  config$simulation <- list(cohort = 970100000, imputation = 980000000)
}
BASE <- if (pilot) c(ind_mice = 760000000, ind_survey = 770000000, ind_design = 780000000) else c(ind_mice = 1260000000, ind_survey = 1270000000, ind_design = 1280000000)
SEEDS <- if (pilot) list(covariate_cohort = 990100000, cohort = 1110000000, audit = 1120000000) else
  list(covariate_cohort = earlier$seeds$covariate_cohort, cohort = 2010000000, audit = 2020000000)
settings <- list(mode = mode, reps_real_25 = if (pilot) pilot_reps else 2000L, reps_other = if (pilot) pilot_reps else 500L, reps_zambia = if (pilot) pilot_reps else 1000L, base = as.list(BASE), seeds = SEEDS)
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
setups <- lapply(METHODS, function(m) { env <- new.env(parent = globalenv()); sys.source(file.path(root, "tools/independent-setups", paste0(m, ".R")), envir = env); env$fit_method })
names(setups) <- METHODS
grid <- utils::read.csv(file.path(root, "validation/clustered-fix-study/simulation-grid.csv"), stringsAsFactors = FALSE)
covariate_grid <- expand.grid(design = "rwanda", errors = c("shift", "under", "row_only"), fraction = config$fractions, stringsAsFactors = FALSE)
covariate_grid$rows <- 3L
cells <- list()
for (z in real) for (j in seq_along(config$fractions)) cells[[length(cells) + 1L]] <- list(kind = "real", z = z, fraction_index = j, level = z$level, analysis = z$id, fraction = config$fractions[j],
  reps = if (j == 1L) settings$reps_real_25 else settings$reps_other, number = 100L + length(cells) + 1L)
for (i in seq_len(nrow(grid))) cells[[length(cells) + 1L]] <- list(kind = "simulation", spec = as.list(grid[i, ]), index = i, level = "simulation", analysis = grid$label[i], fraction = grid$fraction[i], reps = settings$reps_other, number = i)
for (i in seq_len(nrow(covariate_grid))) cells[[length(cells) + 1L]] <- list(kind = "covariate", spec = as.list(covariate_grid[i, ]), index = i, level = "imputed_covariate", analysis = paste0("covariate_", covariate_grid$errors[i]), fraction = covariate_grid$fraction[i], reps = settings$reps_other, number = 200L + i)
cohort <- function(s, r) {
  if (s$kind == "real") real_case(s$z, s$fraction_index, r) else if (s$kind == "simulation")
    simulate_cell(s$spec, as.integer(config$simulation$cohort + s$index * 1e5 + r), s$fraction) else
    simulate_covariate_cell(s$spec, as.integer(SEEDS$covariate_cohort + s$index * 1e5 + r), s$fraction)
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
# The setups get only the fields the brief lists, never the truth or target.
call_setups <- function(p, fraction, target, number, r, labels) {
  shown <- list(input = p$input, formula = p$formula, term = if (is.null(p$term)) "covariate" else p$term, audited = p$audited, fraction = fraction, predictors = p$predictors)
  do.call(rbind, lapply(METHODS, function(m) {
    started <- Sys.time()
    seed <- as.integer(BASE[[m]] + number * 1e5 + r)
    got <- limited(suppressWarnings(suppressMessages(setups[[m]](shown, seed))))
    seconds <- as.numeric(difftime(Sys.time(), started, units = "secs"))
    bad <- inherits(got, "error") || !all(c("estimate", "se", "df") %in% names(got)) || !all(is.finite(unlist(got[c("estimate", "se", "df")]))) || got$se <= 0 || got$df <= 0
    half <- if (bad) NA_real_ else stats::qt(.975, got$df) * got$se
    data.frame(labels, replicate = r, method = m, status = if (bad) "failed" else "ok", message = if (inherits(got, "error")) substr(conditionMessage(got), 1L, 100L) else "",
               estimate = if (bad) NA_real_ else got$estimate, se = if (bad) NA_real_ else got$se, df = if (bad) NA_real_ else got$df, width = if (bad) NA_real_ else 2 * half,
               cover = if (bad) NA_real_ else as.numeric(abs(got$estimate - target) <= half), target = target, seconds = round(seconds, 4), stringsAsFactors = FALSE)
  }))
}
run_all <- function(n, f) {
  batch <- parallel::mclapply(seq_len(n), f, mc.cores = 14L, mc.preschedule = TRUE, mc.set.seed = FALSE)
  bad <- vapply(batch, inherits, logical(1), what = "try-error")
  if (any(bad)) stop("A worker failed: ", conditionMessage(attr(batch[[which(bad)[1]]], "condition")))
  do.call(rbind, batch)
}
write_gz <- function(x, name) { out <- gzfile(file.path(outdir, name), "wt"); utils::write.csv(x, out, row.names = FALSE); close(out) }
start <- Sys.time()
results <- do.call(rbind, lapply(cells, function(s) {
  got <- run_all(s$reps, function(r) {
    p <- cohort(s, r)
    call_setups(p, s$fraction, p$target, s$number, r, data.frame(level = s$level, analysis = s$analysis, fraction = s$fraction))
  })
  cat(s$level, s$analysis, s$fraction, "; elapsed", round(as.numeric(difftime(Sys.time(), start, units = "secs"))), "seconds\n"); flush.console()
  got
}))
write_gz(results, "replicates.csv.gz")
# Description only: whole-cluster audits on the second health system, with that study's cohorts and audits.
pool_all <- utils::read.csv(file.path(root, "validation/second-system/zambia-pool.csv"), stringsAsFactors = FALSE)
zambia <- do.call(rbind, lapply(seq_along(ITEMS), function(i) {
  pool <- pool_all[pool_all$item == ITEMS[i] & pool_all$both_zero == 0L, ]
  target <- unname(stats::coef(stats::lm(verified ~ covariate, data = pool))[2])
  do.call(rbind, lapply(which(names(AUDITS) %in% c("district_33", "district_50")), function(a) {
    got <- run_all(settings$reps_zambia, function(r) {
      co <- draw_cohort(pool, i, r)
      audited <- draw_audit(co, AUDITS[[a]], a, i, r)
      truth <- pool$verified[co$rows]; truth[!audited] <- NA_real_
      p <- list(input = data.frame(covariate = pool$covariate[co$rows], cluster_id = co$ids, lab_value = pool$reported[co$rows], gold_value = truth), formula = lab_value ~ covariate, term = "covariate", audited = audited, predictors = "covariate")
      call_setups(p, AUDITS[[a]]$clusters / 18, target, 300L + i * 2L + a, r, data.frame(item = ITEMS[i], audit = names(AUDITS)[a]))
    })
    got
  }))
}))
write_gz(zambia, "zambia.csv.gz")
if (!identical(hashes_before, hash_sources())) stop("Source, plan, pool or reused results changed during the run; evidence is not verified.")
utils::write.csv(hashes_before, file.path(outdir, "source-hashes.csv"), row.names = FALSE)
writeLines(c(paste("Completed UTC:", format(Sys.time(), tz = "UTC", usetz = TRUE)), paste("Elapsed seconds:", round(as.numeric(difftime(Sys.time(), start, units = "secs")), 1)), capture.output(sessionInfo())), file.path(outdir, "session-info.txt"))
cat("Completed", mode, ":", nrow(results), "+", nrow(zambia), "requests;", paste(names(table(results$status)), table(results$status), collapse = ", "), "\n")
