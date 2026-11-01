#!/usr/bin/env Rscript
# use_recorded was added after the frozen studies, so this checks it against what those
# studies ran as linear MI, on the same real-data cohorts and audit subsets. Different
# code draws different numbers, so agreement is judged on operating characteristics.
# Usage: Rscript tools/check-recorded-option.R [NEW_OUTPUT_DIRECTORY]
args <- commandArgs(trailingOnly = TRUE)
file_arg <- grep("^--file=", commandArgs(), value = TRUE)
root <- normalizePath(file.path(dirname(sub("^--file=", "", file_arg[1])), ".."))
outdir <- if (length(args)) args[1] else file.path(root, "validation", "recorded-option-check")
if (dir.exists(outdir) && length(list.files(outdir, all.files = TRUE, no.. = TRUE)))
  stop("Use a new or empty output directory; existing evidence is never overwritten.")
pkgload::load_all(root, quiet = TRUE)
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
workers <- suppressWarnings(as.integer(Sys.getenv("CATMECH_BENCHMARK_WORKERS", "14")))
settings <- jsonlite::read_json(file.path(root, "validation/real-data-study/settings.json"))
main <- utils::read.csv(file.path(root, "validation/real-data-study/replicates.csv"), stringsAsFactors = FALSE)
targets <- utils::read.csv(file.path(root, "validation/real-data-study/pool-targets.csv"), stringsAsFactors = FALSE)
ctgov <- utils::read.csv(file.path(root, "validation/real-data/ctgov-pool.csv"), stringsAsFactors = FALSE)
ctgov <- ctgov[ctgov$split == "evaluation", ]
nhanes <- utils::read.csv(file.path(root, "validation/real-data/nhanes-pool.csv"), stringsAsFactors = FALSE)
nhanes <- nhanes[nhanes$split == "evaluation", ]
pools <- list(
  ctgov = list(offset = 0, covariates = c("phase", "sponsor", "log_sites", "start_year"),
    data = data.frame(outcome = log(ctgov$duration_months), truth = log(ctgov$reference_n), lab_value = log(ctgov$recorded_n),
      phase = factor(ctgov$phase, levels = c("not_applicable", "early", "phase2", "phase3_4")),
      sponsor = factor(ctgov$sponsor, levels = c("other", "industry")), log_sites = log(ctgov$sites), start_year = ctgov$start_year)),
  nhanes = list(offset = settings$nhanes_offset, covariates = c("age", "sex", "race", "height"),
    data = data.frame(outcome = nhanes$sbp, truth = nhanes$reference_kg, lab_value = nhanes$recorded_kg, age = nhanes$age,
      sex = factor(nhanes$sex, levels = c("female", "male")),
      race = factor(nhanes$race, levels = c("white", "black", "mexican_american", "other_hispanic", "other")), height = nhanes$height)))
run_one <- function(dataset, fraction, r) {
  z <- pools[[dataset]]
  set.seed(as.integer(settings$cohort_base + z$offset + r))
  d <- z$data[sample.int(nrow(z$data), settings$n, replace = TRUE), ]
  rownames(d) <- NULL
  set.seed(as.integer(settings$selection_base + z$offset + r))
  u <- stats::runif(nrow(d))
  d$gold <- ifelse(u < fraction, d$truth, NA_real_)
  input <- d[setdiff(names(d), "truth")]
  fit <- tryCatch(catmech_reference(input, "lab_value", "gold", c(z$covariates, "outcome"), use_recorded = TRUE,
                                    m = settings$m, seed = as.integer(settings$imputation_base + z$offset + r)), error = identity)
  if (inherits(fit, "error")) return(data.frame(dataset = dataset, replicate = r, status = "refused", estimate = NA_real_, lower = NA_real_, upper = NA_real_))
  pooled <- catmech_pool_lm(fit, stats::reformulate(c("lab_value", z$covariates), response = "outcome"))
  row <- pooled[pooled$term == "lab_value", ]
  data.frame(dataset = dataset, replicate = r, status = "ok", estimate = row$estimate, lower = row$lower, upper = row$upper)
}
cells <- list(list(name = "audit_25pct", fraction = .25, reps = 1000L), list(name = "audit_10pct", fraction = .10, reps = 200L), list(name = "audit_50pct", fraction = .50, reps = 200L))
out <- list()
for (cell in cells) for (dataset in names(pools)) {
  batch <- parallel::mclapply(seq_len(cell$reps), function(r) run_one(dataset, cell$fraction, r), mc.cores = workers, mc.preschedule = TRUE, mc.set.seed = FALSE)
  failed <- vapply(batch, inherits, logical(1), what = "try-error")
  if (any(failed)) stop(as.character(batch[[which(failed)[1]]]))
  res <- do.call(rbind, batch)
  target <- targets$target_coefficient[targets$dataset == dataset]
  good <- res[res$status == "ok", ]
  old <- main[main$family == dataset & main$scenario == cell$name & main$method == "cole_mi" & main$status == "ok", ]
  summarize <- function(est, low, up) c(bias = mean(est - target), rmse = sqrt(mean((est - target)^2)), coverage = mean(low <= target & target <= up), width = mean(up - low))
  new <- summarize(good$estimate, good$lower, good$upper)
  was <- summarize(old$estimate, old$lower, old$upper)
  out[[length(out) + 1L]] <- data.frame(dataset = dataset, cell = cell$name, returned = nrow(good), requested = cell$reps,
    use_recorded_bias = new[["bias"]], study_linear_mi_bias = was[["bias"]], use_recorded_rmse = new[["rmse"]], study_linear_mi_rmse = was[["rmse"]],
    use_recorded_coverage = new[["coverage"]], study_linear_mi_coverage = was[["coverage"]], use_recorded_width = new[["width"]], study_linear_mi_width = was[["width"]])
  utils::write.csv(res, file.path(outdir, paste0("replicates-", dataset, "-", cell$name, ".csv")), row.names = FALSE)
}
table <- do.call(rbind, out)
utils::write.csv(table, file.path(outdir, "comparison.csv"), row.names = FALSE)
print(table, digits = 3, row.names = FALSE)
