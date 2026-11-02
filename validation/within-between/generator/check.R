# check.R -- verifies generator.R and grid.csv, and writes check.txt.
#
# For every cell, at BOTH fractions: runs cohorts on seeds separate from the
# target run, fits OLS on full true values (the formula with lab_value replaced
# by truth), and compares the mean coefficient with the stored target.
# Fits nothing else. The only other outputs are descriptive summaries.
#
# Usage: Rscript check.R [cohorts_per_fraction]

args <- commandArgs(FALSE)
script <- sub("^--file=", "", args[grep("^--file=", args)])
if (length(script)) setwd(dirname(normalizePath(script)))
source("generator.R")
library(parallel)

targs <- commandArgs(TRUE)
N_PER <- if (length(targs) >= 1) as.integer(targs[1]) else 2000L
CORES <- max(1L, detectCores() - 2L)
SEED_BASE <- c(`0.25` = 1200000000L, `0.5` = 1600000000L)
SEED_STEP <- 1000000L

grid <- read_grid("grid.csv")
t_start <- Sys.time()
out <- character()
say <- function(...) out <<- c(out, sprintf(...))

full_truth_coef <- function(sim) {
  d <- sim$input
  d$lab_value <- sim$truth
  unname(coef(lm(sim$formula, data = d))[sim$term])
}

# One-way ANOVA ICC, a plain descriptive.
icc <- function(v, cl, sizes) {
  K <- length(sizes)
  N <- sum(sizes)
  vbar <- as.vector(rowsum(v, cl, reorder = FALSE)) / sizes
  ssb <- sum(sizes * (vbar - mean(v))^2)
  ssw <- sum((v - vbar[cl])^2)
  if (ssw + ssb < 1e-12) return(NA_real_)
  msb <- ssb / (K - 1)
  msw <- ssw / max(N - K, 1)
  n0 <- (N - sum(sizes^2) / N) / (K - 1)
  (msb - msw) / (msb + (n0 - 1) * msw)
}

summarise_cohort <- function(cell, seed, fraction) {
  s <- simulate_cell(cell, seed, fraction)
  d <- s$input
  cl <- d$cluster_id
  sizes <- tabulate(cl)
  others <- setdiff(names(d), "gold_value")
  aud_k <- as.vector(rowsum(as.numeric(s$audited), cl, reorder = FALSE))
  err <- d$lab_value - s$truth
  a <- s$audited
  c(coef = full_truth_coef(s),
    na_outside = sum(vapply(d[others], function(x) sum(is.na(x)), numeric(1))),
    gold_wrong = sum(abs(d$gold_value[a] - s$truth[a]) > 1e-12) + sum(is.na(d$gold_value[a])),
    gold_leak = sum(!is.na(d$gold_value[!a])),
    partial_audit = sum(aud_k != 0 & aud_k != sizes),
    clusters = length(sizes),
    aud_clusters = sum(aud_k > 0),
    aud_share_clusters = mean(aud_k > 0),
    aud_share_rows = mean(a),
    rows = nrow(d),
    size_min = min(sizes),
    size_max = max(sizes),
    size_cv = sd(sizes) / mean(sizes),
    differ = mean(abs(err) > 1e-9),
    cor_lab_truth = if (sd(d$lab_value) > 0) cor(d$lab_value, s$truth) else NA_real_,
    err_sd = sd(err),
    truth_sd = sd(s$truth),
    icc_truth = icc(s$truth, cl, sizes),
    icc_err = if (sd(err) > 0) icc(err, cl, sizes) else NA_real_,
    icc_x = icc(d$covariate, cl, sizes),
    min_strata_aud = if ("stratum" %in% names(d))
      min(tabulate(unique(data.frame(cl, h = as.integer(d$stratum))[a, ])$h, nlevels(d$stratum))) else NA_real_,
    aud_x_degenerate = as.numeric(length(unique(d$covariate[a])) < 2),
    aud_gold_degenerate = as.numeric(length(unique(d$gold_value[a])) < 2))
}

# ---------- single-cohort structure tests, and generation time (serial) ----------
struct <- lapply(seq_len(nrow(grid)), function(i) {
  cell <- as.list(grid[i, ])
  seed <- 700000000L + i
  a <- simulate_cell(cell, seed, 0.25)
  b <- simulate_cell(cell, seed, 0.5)
  a2 <- simulate_cell(cell, seed, 0.25)
  d <- a$input
  must <- c("covariate", "cluster_id", "lab_value", "gold_value")
  ok_types <- all(must %in% names(d)) && is.numeric(d$covariate) && is.integer(d$cluster_id) &&
    is.numeric(d$lab_value) && is.numeric(d$gold_value)
  if ("stratum" %in% names(d)) ok_types <- ok_types && is.factor(d$stratum)
  if ("outcome" %in% names(d)) ok_types <- ok_types && is.numeric(d$outcome)
  ok_form <- all(all.vars(a$formula) %in% names(d)) && all(a$predictors %in% names(d))
  role_ok <- if (cell$role == "outcome") {
    a$imputed_outcome && a$term == "covariate" && all.vars(a$formula)[1] == "lab_value"
  } else {
    !a$imputed_outcome && a$term == "lab_value" && all.vars(a$formula)[1] == "outcome"
  }
  strata_ok <- (cell$n_strata > 1) == ("stratum" %in% names(d)) &&
    (cell$n_strata > 1) == ("stratum" %in% all.vars(a$formula)) &&
    (cell$n_strata > 1) == ("stratum" %in% a$predictors)
  outcome_ok <- (cell$role == "predictor") == ("outcome" %in% names(d)) &&
    (cell$role == "predictor") == ("outcome" %in% a$predictors)
  same_data <- identical(a$truth, b$truth) && identical(d$lab_value, b$input$lab_value) &&
    identical(d$covariate, b$input$covariate)
  repro <- identical(a$input, a2$input) && identical(a$audited, a2$audited)
  t1 <- system.time(for (r in 1:50) simulate_cell(cell, 800000000L + r, 0.5))[["elapsed"]] / 50
  c(types = ok_types, formula = ok_form, role = role_ok, strata = strata_ok, outcome = outcome_ok,
    same_across_fractions = same_data, reproducible = repro, target_matches = isTRUE(all.equal(a$target, cell$target)),
    sec_per_cohort = t1)
})
struct <- do.call(rbind, struct)

# ---------- cohort runs at both fractions ----------
run_cell <- function(i, fraction) {
  cell <- as.list(grid[i, ])
  base <- SEED_BASE[[as.character(fraction)]] + i * SEED_STEP
  splits <- splitIndices(N_PER, CORES)
  res <- mclapply(splits, function(idx) {
    t(vapply(idx, function(r) summarise_cohort(cell, base + r, fraction),
             numeric(23)))
  }, mc.cores = CORES)
  do.call(rbind, res)
}
fractions <- c(0.25, 0.5)
results <- vector("list", nrow(grid))
for (i in seq_len(nrow(grid))) {
  results[[i]] <- lapply(fractions, function(f) run_cell(i, f))
}

# ---------- report ----------
say("check.txt -- written by check.R; reflects the code and grid hashed below")
say("run date: %s", format(Sys.time(), "%Y-%m-%d %H:%M %Z"))
say("%s", R.version.string)
gen_md5 <- unname(tools::md5sum("generator.R"))
grid_md5 <- unname(tools::md5sum("grid.csv"))
tr <- if (file.exists("target_run.txt")) readLines("target_run.txt") else character()
tr_md5 <- sub("generator_md5 ", "", tr[grepl("^generator_md5", tr)])
say("generator.R md5: %s", gen_md5)
say("grid.csv md5:    %s", grid_md5)
say("target run used this generator.R: %s", if (length(tr_md5) && tr_md5 == gen_md5) "YES" else "NO")
say("cohorts per cell: %d at each of 2 fractions = %d, on seeds separate from the target run",
    N_PER, 2L * N_PER)
say("")
say("GRID: %d cells; outcome %d, predictor %d (%.0f%% predictor)", nrow(grid), sum(grid$role == "outcome"),
    sum(grid$role == "predictor"), 100 * mean(grid$role == "predictor"))
say("Clusters range %d to %d. With strata: %d cells. Balanced: %d cells. Unbalanced: %d cells.",
    min(grid$K), max(grid$K), sum(grid$n_strata > 1), sum(grid$size_type == "balanced"), sum(grid$size_type != "balanced"))
say("Deliberately error-free cells (label contains 'errorfree'): %s",
    paste(grid$label[grepl("errorfree", grid$label)], collapse = ", "))
say("")

say("TARGET CHECK: mean full-truth OLS coefficient vs stored target; flag if |diff| > 3 x check MC SE")
say("%-24s %-9s %-6s %11s %11s %9s %6s %8s %4s %5s", "label", "role", "method", "target", "mean_ols", "mc_se", "z",
    "tgt_se", "NA", "flag")
n_flag <- 0L
tab <- list()
for (i in seq_len(nrow(grid))) {
  m <- do.call(rbind, results[[i]])
  b <- m[, "coef"]
  mb <- mean(b, na.rm = TRUE)
  se <- sd(b, na.rm = TRUE) / sqrt(sum(!is.na(b)))
  z <- (mb - grid$target[i]) / se
  flag <- abs(z) > 3 || any(is.na(b))
  n_flag <- n_flag + flag
  tab[[i]] <- c(z = z)
  say("%-24s %-9s %-6s %11.5f %11.5f %9.5f %6.2f %8.5f %4d %5s", grid$label[i], grid$role[i],
      if (grid$target_method[i] == "closed_form") "closed" else "MC", grid$target[i], mb, se, z,
      if (grid$target_method[i] == "closed_form") 0 else grid$target_mc_se[i], sum(is.na(b)), if (flag) "FLAG" else "")
}
zs <- unlist(tab)
say("")
say("FLAGGED CELLS: %d of %d", n_flag, nrow(grid))
say("z scores: sd %.2f, largest |z| %.2f; about %.1f cells beyond 2 expected by chance, %d seen",
    sd(zs), max(abs(zs)), 0.0455 * length(zs), sum(abs(zs) > 2))
closed <- grid$target_method == "closed_form"
z2 <- (grid$target_mc[closed] - grid$target[closed]) / grid$target_mc_se[closed]
say("Closed forms against the separate target run (%d cells): largest |z| %.2f, %d beyond 3",
    sum(closed), max(abs(z2)), sum(abs(z2) > 3))
say("")

say("STRUCTURE: counts summed over all cohorts in each cell; all should be 0")
say("%-24s %8s %8s %8s %8s %8s %8s", "label", "NAother", "goldbad", "goldleak", "partial", "strat<1", "coefNA")
bad_struct <- 0
for (i in seq_len(nrow(grid))) {
  m <- do.call(rbind, results[[i]])
  sa <- if (grid$n_strata[i] > 1) sum(m[, "min_strata_aud"] < 1) else 0
  v <- c(sum(m[, "na_outside"]), sum(m[, "gold_wrong"]), sum(m[, "gold_leak"]), sum(m[, "partial_audit"]), sa,
         sum(is.na(m[, "coef"])))
  bad_struct <- bad_struct + sum(v)
  say("%-24s %8d %8d %8d %8d %8d %8d", grid$label[i], v[1], v[2], v[3], v[4], v[5], v[6])
}
say("Total structural failures: %d", bad_struct)
say("")

say("SINGLE-COHORT TESTS: column types, formula, role fields, strata fields, outcome fields; same truth,")
say("recorded values and covariate at both fractions for one seed; same seed reproduces; target stored.")
fails <- which(rowSums(struct[, 1:8, drop = FALSE] != 1) > 0)
say("Cells failing any: %d%s", length(fails), if (length(fails)) paste0(" (", paste(grid$label[fails], collapse = ", "), ")") else "")
say("Generation time per cohort: mean %.2f ms, slowest cell %.2f ms (%s), serial, excluding OLS",
    1000 * mean(struct[, "sec_per_cohort"]), 1000 * max(struct[, "sec_per_cohort"]),
    grid$label[which.max(struct[, "sec_per_cohort"])])
say("")

say("AUDIT AND ROWS: audited share of clusters and audited clusters (min-max), rows per cohort (min/mean/max), cluster size range and mean CV")
say("%-24s %4s | %6s %7s | %6s %7s | %5s %6s %5s | %6s %6s", "label", "K", "shr25", "n25", "shr50", "n50", "rmin",
    "rmean", "rmax", "size", "cv")
over_cap <- 0
for (i in seq_len(nrow(grid))) {
  m25 <- results[[i]][[1]]
  m50 <- results[[i]][[2]]
  m <- rbind(m25, m50)
  over_cap <- over_cap + sum(m[, "rows"] > 1500)
  say("%-24s %4d | %6.3f %3d-%3d | %6.3f %3d-%3d | %5d %6.0f %5d | %2d-%2d %6.2f", grid$label[i], grid$K[i],
      mean(m25[, "aud_share_clusters"]), min(m25[, "aud_clusters"]), max(m25[, "aud_clusters"]),
      mean(m50[, "aud_share_clusters"]), min(m50[, "aud_clusters"]), max(m50[, "aud_clusters"]),
      min(m[, "rows"]), mean(m[, "rows"]), max(m[, "rows"]),
      min(m[, "size_min"]), max(m[, "size_max"]), mean(m[, "size_cv"]))
}
say("Cohorts over 1,500 rows: %d", over_cap)
say("")

say("DESCRIPTIVES: share of rows where recorded differs from true (cell mean, lowest cohort); correlation of recorded with true;")
say("error SD / true SD (median); ICC of true value, of error, of covariate (median); audited rows share; share of cohorts")
say("where audited rows have a constant covariate (at fraction 0.25, the harder one)")
say("%-24s %6s %6s %6s %8s %6s %6s %6s %6s %6s", "label", "differ", "minDif", "cor", "errSD/SD", "iccT", "iccE", "iccX", "rows25", "xDeg")
no_err <- character()
for (i in seq_len(nrow(grid))) {
  m25 <- results[[i]][[1]]
  m <- rbind(m25, results[[i]][[2]])
  ef <- grepl("errorfree", grid$label[i])
  if (!ef && max(m[, "differ"]) == 0) no_err <- c(no_err, grid$label[i])
  if (!ef && min(m[, "differ"]) == 0) no_err <- c(no_err, paste0(grid$label[i], " (some cohorts)"))
  if (ef && max(m[, "differ"]) != 0) no_err <- c(no_err, paste0(grid$label[i], " (should be error-free)"))
  say("%-24s %6.3f %6.3f %6.3f %8.3f %6.2f %6.2f %6.2f %6.3f %6.3f", grid$label[i], mean(m[, "differ"]), min(m[, "differ"]),
      median(m[, "cor_lab_truth"], na.rm = TRUE), median(m[, "err_sd"] / m[, "truth_sd"]),
      median(m[, "icc_truth"]), if (ef) NA_real_ else median(m[, "icc_err"], na.rm = TRUE), median(m[, "icc_x"]),
      mean(m25[, "aud_share_rows"]), mean(m25[, "aud_x_degenerate"]))
}
say("Recorded-vs-true problems (non-error-free cell with an all-equal cohort, or error-free cell that differs): %d%s",
    length(no_err), if (length(no_err)) paste0(" -- ", paste(no_err, collapse = "; ")) else "")
say("")
elapsed <- as.numeric(difftime(Sys.time(), t_start, units = "secs"))
say("OVERALL: %d flagged, %d structural failures, %d single-cohort failures, %d over the row cap, %d recorded-vs-true problems",
    n_flag, bad_struct, length(fails), over_cap, length(no_err))
say("check.R elapsed: %.0f seconds on %d cores", elapsed, CORES)
writeLines(out, "check.txt")
cat(out[grepl("^(FLAGGED|OVERALL|Total|Recorded|Generation|Closed|z scores)", out)], sep = "\n")
cat(out[grepl("FLAG$", out)], sep = "\n")
