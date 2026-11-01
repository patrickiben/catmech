#!/usr/bin/env Rscript
# Builds the two pools defined in validation/REAL_AUDIT_PLAN.md from the downloaded
# source files and writes a manifest with sources, licenses and hashes. Needs the haven package, which is not a
# package dependency, to read the Stata file.
# Usage: Rscript tools/real-audit-prepare.R RAW_DIRECTORY OUTPUT_DIRECTORY
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("Usage: real-audit-prepare.R RAW_DIRECTORY OUTPUT_DIRECTORY")
raw <- args[1]; output <- args[2]
dir.create(output, recursive = TRUE, showWarnings = FALSE)
rwanda_target <- file.path(output, "rwanda-pool.csv")
effect_target <- file.path(output, "effect-size-pool.csv")
if (file.exists(rwanda_target) || file.exists(effect_target)) stop("An existing pool is never overwritten.")
rwanda_file <- file.path(raw, "pone.0235823.s001.dta")
effect_file <- file.path(raw, "codebook-primary-studies-complete.csv")

rwanda <- as.data.frame(lapply(haven::read_dta(rwanda_file), function(x) as.vector(haven::zap_labels(x))))
if (nrow(rwanda) != 228L || length(unique(rwanda$facility_id)) != 76L) stop("Unexpected Rwanda file shape.")
rwanda$facility_type <- c("hospital", "health_center")[rwanda$facility_cat]
if (anyNA(rwanda$facility_type)) stop("Unexpected facility category.")
utils::write.csv(rwanda, rwanda_target, row.names = FALSE, na = "")

effect <- utils::read.table(effect_file, header = TRUE, stringsAsFactors = FALSE, fill = TRUE, comment.char = "", encoding = "latin1")
if (nrow(effect) != 500L) stop("Unexpected effect-size file shape.")
num <- function(x) suppressWarnings(as.numeric(x))
pool <- data.frame(row = seq_len(nrow(effect)), meta = effect$meta, efftype = effect$efftype,
  family = ifelse(effect$efftype %in% c("g", "d"), "smd", "r"), stratum_sampled = effect$typestudy,
  n = num(effect$n), recorded = num(effect$effest), reference = num(effect$effestnew), info = num(effect$info), stringsAsFactors = FALSE)
counts <- table(factor(pool$info, levels = 0:3))
if (!identical(as.vector(counts), c(276L, 74L, 54L, 96L))) stop("Effect-size categories do not match the published counts.")
if (anyNA(pool[c("n", "recorded", "reference")])) stop("Unexpected missing effect-size value.")
utils::write.csv(pool, effect_target, row.names = FALSE)

source_entry <- function(path, url, license, citation)
  list(file = basename(path), url = url, license = license, citation = citation, bytes = file.size(path),
       sha256 = digest::digest(file = path, algo = "sha256"))
manifest <- list(retrieved_date = "2026-11-01",
  rwanda = c(source_entry(rwanda_file, "https://journals.plos.org/plosone/article/file?type=supplementary&id=10.1371/journal.pone.0235823.s001",
      "CC BY (article license)", "Health management information system (HMIS) data verification: A case study in four districts in Rwanda. PLOS ONE 2020. doi:10.1371/journal.pone.0235823"),
    list(rows = nrow(rwanda), facilities = length(unique(rwanda$facility_id)), pool_sha256 = digest::digest(file = rwanda_target, algo = "sha256"))),
  effect_sizes = c(source_entry(effect_file, "https://osf.io/7nsmd/ (codebooks/codebook-primary-studies-complete.csv in the linked storage)",
      "CC0 1.0 Universal", "Maassen et al. Reproducibility of individual effect sizes in meta-analyses in psychology. PLOS ONE 2020. doi:10.1371/journal.pone.0233107"),
    list(rows = nrow(pool), meta_analyses = length(unique(pool$meta)), pool_sha256 = digest::digest(file = effect_target, algo = "sha256"))))
jsonlite::write_json(manifest, file.path(output, "fetch-manifest.json"), pretty = TRUE, auto_unbox = TRUE)
cat("Rwanda rows", nrow(rwanda), "facilities", length(unique(rwanda$facility_id)), "health centers",
    length(unique(rwanda$facility_id[rwanda$facility_type == "health_center"])), "\n")
cat("Effect-size rows", nrow(pool), "meta-analyses", length(unique(pool$meta)), "with independent reference", sum(pool$info %in% 0:1), "\n")
