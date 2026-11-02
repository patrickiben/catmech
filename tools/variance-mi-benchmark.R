#!/usr/bin/env Rscript
# Prospective test of corrected cluster-robust variances through imputation; settings locked in validation/VARIANCE_MI_PLAN.md.
# Imputations come from the package's own catmech_reference(). The variances are computed here, because the corrected ones are not yet in the package.
# Usage: Rscript tools/variance-mi-benchmark.R pilot|study [NEW_DIRECTORY] [PILOT_REPS]
args <- commandArgs(trailingOnly = TRUE)
mode <- if (length(args)) match.arg(args[1], c("pilot", "study")) else "pilot"
file_arg <- grep("^--file=", commandArgs(), value = TRUE)
root <- normalizePath(file.path(dirname(sub("^--file=", "", file_arg[1])), ".."))
outdir <- if (length(args) >= 2L) args[2] else file.path(root, "validation", paste0("variance-mi-", mode))
if (dir.exists(outdir) && length(list.files(outdir, all.files = TRUE, no.. = TRUE)))
  stop("Use a new or empty output directory; existing evidence is never overwritten.")
hash_sources <- function() {
  paths <- c(sort(list.files(file.path(root, "R"), "\\.R$", full.names = TRUE)),
    file.path(root, c("DESCRIPTION", "NAMESPACE", "tools/variance-mi-benchmark.R", "tools/variance-mi-score.py", "tools/clustered-fix-benchmark.R", "validation/VARIANCE_MI_PLAN.md", "validation/real-audit/rwanda-pool.csv")))
  data.frame(file = substring(paths, nchar(root) + 2L),
    sha256 = vapply(paths, function(p) digest::digest(file = p, algo = "sha256"), character(1)), row.names = NULL)
}
hashes_before <- hash_sources()
pkgload::load_all(root, quiet = TRUE)
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
pilot <- mode == "pilot"
pilot_reps <- if (length(args) >= 3L) as.integer(args[3]) else 3L
config <- jsonlite::read_json(file.path(root, "validation/clustered-fix-study/settings.json"), simplifyVector = TRUE)
if (pilot) {
  config$held_out <- list(cohort = 1810100000, selection = 1820100000, imputation = 1830000000)
  config$development <- list(cohort = 1840100000, selection = 1850100000, imputation = 1860000000)
  config$simulation <- list(cohort = 1800100000, imputation = 1870000000)
}
expressions <- parse(file.path(root, "tools/clustered-fix-benchmark.R"))
assigned <- vapply(expressions, function(e) if (is.call(e) && identical(e[[1]], as.name("<-")) && is.name(e[[2]])) as.character(e[[2]]) else "", character(1))
first <- which(assigned == "rwanda")[1]; last <- which(assigned == "real_case")[1]
for (e in expressions[first:last]) eval(e, envir = globalenv())
grid <- utils::read.csv(file.path(root, "validation/clustered-fix-study/simulation-grid.csv"), stringsAsFactors = FALSE)
SEEDS <- list(covariate_cohort = if (pilot) 1880100000 else 1700100000, covariate_imputation = if (pilot) 1890000000 else 1710000000,
              covariate_wild = 1720000000, wild = if (pilot) 1900000000 else 1730000000)
settings <- list(mode = mode, m = 20L, reps_real_25 = if (pilot) pilot_reps else 2000L, reps_other = if (pilot) pilot_reps else 500L, wild_replications = 999L, seeds = SEEDS)
jsonlite::write_json(settings, file.path(outdir, "settings.json"), pretty = TRUE, auto_unbox = TRUE)
BOOT <- settings$wild_replications

sym_power <- function(A, power) {
  e <- eigen((A + t(A)) / 2, symmetric = TRUE)
  keep <- e$values > 1e-10
  e$vectors[, keep, drop = FALSE] %*% (t(e$vectors[, keep, drop = FALSE]) * e$values[keep]^power)
}
# For a fixed design, each variance's cluster scores are linear in the response: score_g = row g of a matrix times y.
analysis_parts <- function(X, ids, k) {
  n <- nrow(X); p <- ncol(X)
  units <- split(seq_len(n), ids); G <- length(units)
  XtXi <- chol2inv(chol(crossprod(X)))
  b <- XtXi[, k]; Xb <- as.vector(X %*% b)
  S <- W <- C <- matrix(0, G, n)
  for (g in seq_len(G)) {
    u <- units[[g]]; Xg <- X[u, , drop = FALSE]
    I_H <- diag(length(u)) - Xg %*% XtXi %*% t(Xg)
    M_col <- function(v) {
      out <- -as.vector(X %*% (XtXi %*% crossprod(Xg, v)))
      out[u] <- out[u] + v
      out
    }
    S[g, ] <- M_col(Xb[u])
    W[g, ] <- M_col(as.vector(sym_power(I_H, -.5) %*% Xb[u]))
    C[g, ] <- M_col(as.vector(sym_power(I_H, -1) %*% Xb[u]))
  }
  WW <- tcrossprod(W)
  list(Xb = Xb, S = S, W = W, C = C, f = (G / (G - 1)) * ((n - 1) / (n - p)), j = (G - 1) / G, df = sum(diag(WW))^2 / sum(WW^2), G = G)
}
one_analysis <- function(parts, Y) {
  list(q = as.vector(crossprod(parts$Xb, Y)), cr1s = parts$f * colSums((parts$S %*% Y)^2), cr2 = colSums((parts$W %*% Y)^2), cv3 = parts$j * colSums((parts$C %*% Y)^2))
}
# Rubin's rules with Barnard-Rubin degrees of freedom, with the small floor on lambda that mice uses.
pool_one <- function(q, u, dfcom) {
  m <- length(q)
  if (m == 1L) return(list(estimate = q, se = sqrt(u), df = dfcom))
  ubar <- mean(u); b <- stats::var(q); total <- ubar + (1 + 1 / m) * b
  lambda <- max((1 + 1 / m) * b / total, 1e-4)
  dfold <- (m - 1) / lambda^2; dfobs <- (dfcom + 1) / (dfcom + 3) * dfcom * (1 - lambda)
  list(estimate = mean(q), se = sqrt(total), df = dfold * dfobs / (dfold + dfobs))
}
# Exploratory: wild reference distribution for the pooled statistic, from the same weights in every imputation.
wild_pooled_p <- function(X, Y, ids, k, target, seed, observed, f) {
  Z <- X[, -k, drop = FALSE]; x <- X[, k]
  qz <- qr(Z)
  x_tilde <- qr.resid(qz, x)
  U <- qr.resid(qz, Y - x * target)
  groups <- as.integer(factor(ids)); G <- max(groups)
  XtXi <- chol2inv(chol(crossprod(X)))
  L <- rowsum(x_tilde * X, groups, reorder = TRUE)
  set.seed(seed)
  xi <- matrix(sample(c(-1, 1), BOOT * G, replace = TRUE), BOOT, G)
  numerator <- numeric(BOOT); variance <- numeric(BOOT)
  m <- ncol(Y)
  for (j in seq_len(m)) {
    u <- U[, j]
    s <- as.vector(rowsum(x_tilde * u, groups, reorder = TRUE))
    A <- L %*% XtXi %*% t(rowsum(u * X, groups, reorder = TRUE))
    numerator <- numerator + as.vector(xi %*% s)
    cluster <- sweep(xi, 2L, s, `*`) - xi %*% t(A)
    variance <- variance + rowSums(cluster^2)
  }
  stat <- (numerator / m) / sqrt(f * variance / m)
  (1 + sum(abs(stat) >= abs(observed))) / (BOOT + 1)
}

simulate_covariate_cell <- function(cell, seed, fraction) {
  p <- simulate_cell(cell, seed, fraction)
  d <- p$input
  id <- d$cluster_id; h <- as.integer(d$stratum)
  d$outcome <- 1 + .5 * p$truth + .3 * d$covariate + stats::rnorm(nlevels(d$stratum))[h] + rep(stats::rnorm(max(id), 0, sqrt(5)), each = cell$rows) + stats::rnorm(nrow(d), 0, sqrt(20))
  list(input = d, truth = p$truth, target = .5, audited = p$audited, formula = outcome ~ lab_value + covariate + stratum, predictors = c("covariate", "stratum", "outcome"), term = "lab_value", imputed_outcome = FALSE)
}

IMPUTERS <- c("oracle", "linear_mi", "random_intercept")
run_imputer <- function(imputer, p, seed, wild_seed) {
  input <- p$input
  completed_y <- if (imputer == "oracle") matrix(p$truth, ncol = 1L) else {
    fit <- if (imputer == "linear_mi") {
      catmech_reference(input, "lab_value", "gold_value", p$predictors, method = "normal", m = settings$m, seed = seed, use_recorded = TRUE)
    } else {
      catmech_reference(input, "lab_value", "gold_value", p$predictors, method = "normal", m = settings$m, seed = seed, use_recorded = TRUE, cluster = "cluster_id")
    }
    vapply(fit$completed, function(z) z$lab_value, numeric(nrow(input)))
  }
  term <- if (is.null(p$term)) "covariate" else p$term
  imputed_outcome <- is.null(p$imputed_outcome) || p$imputed_outcome
  if (imputed_outcome) {
    d <- input; d$lab_value <- completed_y[, 1L]
    X <- stats::model.matrix(p$formula, data = d)
    k <- match(term, colnames(X))
    parts <- analysis_parts(X, d$cluster_id, k)
    a <- one_analysis(parts, completed_y)
    df_cr2 <- parts$df
    f <- parts$f; G <- parts$G
  } else {
    pieces <- lapply(seq_len(ncol(completed_y)), function(j) {
      d <- input; d$lab_value <- completed_y[, j]
      X <- stats::model.matrix(p$formula, data = d)
      parts <- analysis_parts(X, d$cluster_id, match(term, colnames(X)))
      c(one_analysis(parts, matrix(d$outcome, ncol = 1L)), list(df = parts$df, G = parts$G))
    })
    a <- list(q = vapply(pieces, `[[`, numeric(1), "q"), cr1s = vapply(pieces, `[[`, numeric(1), "cr1s"), cr2 = vapply(pieces, `[[`, numeric(1), "cr2"), cv3 = vapply(pieces, `[[`, numeric(1), "cv3"))
    df_cr2 <- mean(vapply(pieces, `[[`, numeric(1), "df")); G <- pieces[[1L]]$G
  }
  P <- list(cr1s = pool_one(a$q, a$cr1s, G - 1), cr2 = pool_one(a$q, a$cr2, df_cr2), cv3 = pool_one(a$q, a$cv3, G - 1))
  row <- list(estimate = P$cr1s$estimate, target = p$target, clusters = G, df_complete_cr2 = df_cr2, between = if (ncol(completed_y) > 1L) stats::var(a$q) else 0)
  for (v in names(P)) {
    half <- stats::qt(.975, P[[v]]$df) * P[[v]]$se
    row[[paste0("se_", v)]] <- P[[v]]$se; row[[paste0("df_", v)]] <- P[[v]]$df; row[[paste0("width_", v)]] <- 2 * half
    row[[paste0("cover_", v)]] <- as.numeric(abs(P[[v]]$estimate - p$target) <= half)
  }
  row$wcr_p <- if (imputed_outcome) wild_pooled_p(X, completed_y, input$cluster_id, k, p$target, wild_seed, (P$cr1s$estimate - p$target) / P$cr1s$se, f) else NA_real_
  row
}
NA_ROW <- list(estimate = NA_real_, target = NA_real_, clusters = NA_real_, df_complete_cr2 = NA_real_, between = NA_real_)
for (v in c("cr1s", "cr2", "cv3")) for (f in c("se_", "df_", "width_", "cover_")) NA_ROW[[paste0(f, v)]] <- NA_real_
NA_ROW$wcr_p <- NA_real_

cells <- list()
for (z in real) for (j in seq_along(config$fractions)) cells[[length(cells) + 1L]] <- list(kind = "real", z = z, fraction_index = j, level = z$level, analysis = z$id, fraction = config$fractions[j],
  reps = if (j == 1L) settings$reps_real_25 else settings$reps_other, wild_cell = 100L + length(cells) + 1L)
for (i in seq_len(nrow(grid))) cells[[length(cells) + 1L]] <- list(kind = "simulation", spec = as.list(grid[i, ]), index = i, level = "simulation", analysis = grid$label[i], fraction = grid$fraction[i], reps = settings$reps_other, wild_cell = i)
covariate_grid <- expand.grid(design = "rwanda", errors = c("shift", "under", "row_only"), fraction = config$fractions, stringsAsFactors = FALSE)
covariate_grid$rows <- 3L
for (i in seq_len(nrow(covariate_grid))) cells[[length(cells) + 1L]] <- list(kind = "covariate", spec = as.list(covariate_grid[i, ]), index = i, level = "imputed_covariate", analysis = paste0("covariate_", covariate_grid$errors[i]), fraction = covariate_grid$fraction[i], reps = settings$reps_other, wild_cell = i)

one_replicate <- function(s, r) {
  if (s$kind == "real") {
    p <- real_case(s$z, s$fraction_index, r)
    seed <- as.integer(s$z$seeds$imputation + s$z$index * 1e6 + s$fraction_index * 1e5 + r)
    wild_seed <- as.integer(SEEDS$wild + s$wild_cell * 1e5 + r)
  } else if (s$kind == "simulation") {
    p <- simulate_cell(s$spec, as.integer(config$simulation$cohort + s$index * 1e5 + r), s$fraction)
    seed <- as.integer(config$simulation$imputation + s$index * 1e5 + r)
    wild_seed <- as.integer(SEEDS$wild + s$wild_cell * 1e5 + r)
  } else {
    p <- simulate_covariate_cell(s$spec, as.integer(SEEDS$covariate_cohort + s$index * 1e5 + r), s$fraction)
    seed <- as.integer(SEEDS$covariate_imputation + s$index * 1e5 + r)
    wild_seed <- as.integer(SEEDS$covariate_wild + s$index * 1e5 + r)
  }
  do.call(rbind, lapply(IMPUTERS, function(imputer) {
    warnings <- character()
    got <- tryCatch(withCallingHandlers(run_imputer(imputer, p, seed, wild_seed), warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning") }), error = identity)
    failed <- inherits(got, "error")
    status <- if (!failed) "ok" else if (grepl("Fewer than five|lost design support|full column rank|more observations than coefficients|constant predictor|no residual variance", conditionMessage(got))) "refused" else "error"
    data.frame(level = s$level, analysis = s$analysis, fraction = s$fraction, replicate = r, imputer = imputer, status = status, message = if (failed) conditionMessage(got) else "",
               warnings = paste(unique(warnings), collapse = " | "), as.data.frame(if (failed) NA_ROW else got), stringsAsFactors = FALSE)
  }))
}
start <- Sys.time(); records <- list()
for (s in cells) {
  batch <- parallel::mclapply(seq_len(s$reps), function(r) one_replicate(s, r), mc.cores = 14L, mc.preschedule = TRUE, mc.set.seed = FALSE)
  bad <- vapply(batch, inherits, logical(1), what = "try-error")
  if (any(bad)) stop("A worker failed on ", s$analysis, ": ", conditionMessage(attr(batch[[which(bad)[1]]], "condition")))
  records <- c(records, batch)
  cat(s$level, s$analysis, s$fraction, "; elapsed", round(as.numeric(difftime(Sys.time(), start, units = "secs"))), "seconds\n"); flush.console()
}
if (!identical(hashes_before, hash_sources())) stop("Source, plan or pool changed during run; evidence is not verified.")
results <- do.call(rbind, records)
out <- gzfile(file.path(outdir, "replicates.csv.gz"), "wt")
utils::write.csv(results, out, row.names = FALSE)
close(out)
utils::write.csv(hashes_before, file.path(outdir, "source-hashes.csv"), row.names = FALSE)
writeLines(c(paste("Completed UTC:", format(Sys.time(), tz = "UTC", usetz = TRUE)), paste("Elapsed seconds:", round(as.numeric(difftime(Sys.time(), start, units = "secs")), 1)), capture.output(sessionInfo())), file.path(outdir, "session-info.txt"))
cat("Completed", mode, ":", nrow(results), "requests;", paste(names(table(results$status)), table(results$status), collapse = ", "), "\n")
