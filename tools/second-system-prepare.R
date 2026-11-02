#!/usr/bin/env Rscript
# Builds the item pools for the second-system test from the unchanged raw file, as defined in validation/SECOND_SYSTEM_PLAN.md.
# Usage: Rscript tools/second-system-prepare.R [OUTPUT_DIRECTORY]
args <- commandArgs(trailingOnly = TRUE)
file_arg <- grep("^--file=", commandArgs(), value = TRUE)
root <- normalizePath(file.path(dirname(sub("^--file=", "", file_arg[1])), ".."))
output <- if (length(args)) args[1] else file.path(root, "validation", "second-system")
target <- file.path(output, "zambia-pool.csv")
if (file.exists(target)) stop("An existing pool is never overwritten.")
raw <- as.data.frame(readxl::read_excel(file.path(root, "validation/second-system/raw/Zambia_data.xlsx")))
items <- c("opd", "deliveries", "anc", "postnatal", "immunization", "ipt", "fpmethods", "counselling", "exposedbabies")
pool <- do.call(rbind, lapply(items, function(item) {
  data.frame(item = item, clinic = raw$id, district = raw$district, province = raw$province, catchment = raw$q9,
             reported = raw[[paste0(item, "_qtr1_hia")]], verified = raw[[paste0(item, "_qtr1_ts")]], stringsAsFactors = FALSE)
}))
pool <- pool[stats::complete.cases(pool$reported, pool$verified, pool$catchment), ]
pool$both_zero <- as.integer(pool$reported == 0 & pool$verified == 0)
pool$covariate <- log10(pool$catchment)
rownames(pool) <- NULL
dir.create(output, recursive = TRUE, showWarnings = FALSE)
utils::write.csv(pool, target, row.names = FALSE)
print(table(pool$item, pool$both_zero))
