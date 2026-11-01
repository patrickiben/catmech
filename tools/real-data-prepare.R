#!/usr/bin/env Rscript
# Builds nhanes-pool.csv from the downloaded CDC transport files, as defined in
# validation/REAL_DATA_PLAN.md, and adds exclusion counts to the fetch manifest.
# Usage: Rscript tools/real-data-prepare.R RAW_DIRECTORY OUTPUT_DIRECTORY
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("Usage: real-data-prepare.R RAW_DIRECTORY OUTPUT_DIRECTORY")
raw <- args[1]; output <- args[2]
target <- file.path(output, "nhanes-pool.csv")
if (file.exists(target)) stop("An existing pool is never overwritten.")
read_file <- function(name) foreign::read.xport(file.path(raw, paste0(name, ".xpt")))
demo <- read_file("DEMO_J")[c("SEQN", "RIDAGEYR", "RIAGENDR", "RIDRETH1", "RIDEXPRG")]
body <- read_file("BMX_J")[c("SEQN", "BMXWT", "BMXHT")]
weight <- read_file("WHQ_J")[c("SEQN", "WHD020")]
pressure <- read_file("BPX_J")[c("SEQN", "BPXSY1", "BPXSY2", "BPXSY3", "BPXSY4")]
d <- Reduce(function(a, b) merge(a, b, by = "SEQN", all.x = TRUE), list(demo, body, weight, pressure))
counts <- c(participants = nrow(d))
d <- d[d$RIDAGEYR >= 20, ]; counts["adults_20_plus"] <- nrow(d)
d <- d[is.na(d$RIDEXPRG) | d$RIDEXPRG != 1, ]; counts["not_pregnant"] <- nrow(d)
readings <- as.matrix(d[c("BPXSY1", "BPXSY2", "BPXSY3", "BPXSY4")])
d$sbp <- ifelse(rowSums(!is.na(readings)) > 0, rowMeans(readings, na.rm = TRUE), NA_real_)
d <- d[!is.na(d$BMXWT) & !is.na(d$BMXHT) & !is.na(d$sbp), ]; counts["complete_analysis_values"] <- nrow(d)
recorded <- d$WHD020
recorded[recorded %in% c(7777, 9999)] <- NA_real_
counts["self_report_missing"] <- sum(is.na(recorded))
split_of <- function(id) {
  # The same rule as the Python fetcher: SHA-256 of the identifier, taken modulo 5.
  # The whole hash is reduced one hex digit at a time, since it exceeds double precision.
  hash <- digest::digest(as.character(id), algo = "sha256", serialize = FALSE)
  value <- 0L
  for (digit in strsplit(hash, "")[[1]]) value <- (value * 16L + strtoi(digit, 16L)) %% 5L
  if (value == 0L) "development" else "evaluation"
}
pool <- data.frame(seqn = d$SEQN, split = vapply(d$SEQN, split_of, character(1)),
  recorded_kg = round(recorded * 0.45359237, 3), reference_kg = d$BMXWT, sbp = round(d$sbp, 3), age = d$RIDAGEYR,
  sex = ifelse(d$RIAGENDR == 1, "male", "female"),
  race = c("mexican_american", "other_hispanic", "white", "black", "other")[d$RIDRETH1], height = d$BMXHT)
if (anyNA(pool[setdiff(names(pool), "recorded_kg")])) stop("Unexpected missing analysis value.")
utils::write.csv(pool, target, row.names = FALSE)
manifest_path <- file.path(output, "fetch-manifest.json")
manifest <- jsonlite::read_json(manifest_path)
manifest$nhanes_pool <- c(as.list(counts), list(rows_retained = nrow(pool),
  pool_sha256 = digest::digest(file = target, algo = "sha256")))
jsonlite::write_json(manifest, manifest_path, pretty = TRUE, auto_unbox = TRUE)
print(counts)
print(table(pool$split))
