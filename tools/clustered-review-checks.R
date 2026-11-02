#!/usr/bin/env Rscript
# Unplanned checks raised by the independent review of the clustered-audit test. Post hoc: they enter no rule.
# Same cohorts, audits and imputation seeds as the study. Only the resampling inside the imputation changes.
# Usage: Rscript tools/clustered-review-checks.R [NEW_DIRECTORY]
args <- commandArgs(trailingOnly = TRUE)
file_arg <- grep("^--file=", commandArgs(), value = TRUE)
root <- normalizePath(file.path(dirname(sub("^--file=", "", file_arg[1])), ".."))
outdir <- if (length(args)) args[1] else file.path(root, "validation", "clustered-review")
if (dir.exists(outdir) && length(list.files(outdir, all.files = TRUE, no.. = TRUE)))
  stop("Use a new or empty output directory; existing evidence is never overwritten.")
pkgload::load_all(root, quiet = TRUE)
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
workers <- suppressWarnings(as.integer(Sys.getenv("CATMECH_BENCHMARK_WORKERS", "14")))
config <- list(m = 20L, fractions = c(.25, .50), reps = 500L,
  held_out = list(cohort = 860100000, selection = 870100000, imputation = 880000000),
  development = list(cohort = 400100000, selection = 410100000, imputation = 420000000))

rwanda <- utils::read.csv(file.path(root, "validation/real-audit/rwanda-pool.csv"), stringsAsFactors = FALSE)
rwanda <- rwanda[rwanda$facility_type == "health_center", ]
add <- function(...) { parts <- list(...); out <- parts[[1]]; for (p in parts[-1]) out <- out + p; out }
HELD_OUT <- list(
  lbw = list(w = with(rwanda, add(lbw_notpremature_h, lbw_premature_h)), r = rwanda$lbw_r, x = rwanda$hf_deliveries_h),
  resus = list(w = rwanda$newborns_nocryres_succ_h, r = rwanda$nocry_res_success_r, x = rwanda$hf_deliveries_h),
  stillbirths = list(w = with(rwanda, add(sbmacerated_h, sbfresh_h)), r = with(rwanda, add(sbmacerated_r, sbfresh_r)), x = rwanda$hf_deliveries_h),
  deliveries = list(w = rwanda$hf_deliveries_h, r = rwanda$hf_deliveries_r, x = rwanda$ancnewreg_h),
  anc1_early = list(w = rwanda$anc1vis1trim_h, r = rwanda$anc1vis_16weeks_r, x = rwanda$ancnewreg_h))
DEVELOPMENT <- list(
  pnc = list(w = rwanda$pnc_newborn_24h_h, r = rwanda$pnc_newborn_24h_r, x = rwanda$hf_deliveries_h),
  syphilis = list(w = rwanda$anc_syphilis_test_h, r = rwanda$anc_syphilis_test_r, x = rwanda$ancnewreg_h),
  iron = list(w = rwanda$womenreciron_h, r = rwanda$womenreciron_r, x = rwanda$ancnewreg_h),
  anc4 = list(w = rwanda$anc_4visits_h, r = rwanda$anc_4visits_r, x = rwanda$ancnewreg_h))
real_analysis <- function(id, values, index, level, seeds) {
  keep <- stats::complete.cases(values$w, values$r, values$x)
  pool <- data.frame(unit = paste0("f", rwanda$facility_id[keep]), stratum = factor(paste0("c", rwanda$hosp_id[keep])),
                     lab_value = values$w[keep], lab_true = values$r[keep], covariate = values$x[keep], stringsAsFactors = FALSE)
  list(id = id, index = index, level = level, seeds = seeds, pool = pool, rows_by_unit = split(seq_len(nrow(pool)), pool$unit),
       units_by_stratum = lapply(split(pool$unit, pool$stratum), unique))
}
real <- c(Map(real_analysis, names(HELD_OUT), HELD_OUT, seq_along(HELD_OUT), "held_out", list(config$held_out)),
          Map(real_analysis, names(DEVELOPMENT), DEVELOPMENT, seq_along(DEVELOPMENT), "development", list(config$development)))
real_case <- function(z, fraction_index, r) {
  set.seed(as.integer(z$seeds$cohort + z$index * 1e6 + r))
  drawn <- lapply(z$units_by_stratum, function(u) u[sample.int(length(u), length(u), replace = TRUE)])
  set.seed(as.integer(z$seeds$selection + z$index * 1e6 + fraction_index * 1e5 + r))
  fraction <- config$fractions[fraction_index]
  picked <- lapply(drawn, function(u) seq_along(u) %in% sample.int(length(u), max(1L, as.integer(round(fraction * length(u))))))
  parts <- list(); position <- 0L
  for (s in names(drawn)) for (j in seq_along(drawn[[s]])) {
    position <- position + 1L
    rows <- z$rows_by_unit[[drawn[[s]][j]]]
    parts[[position]] <- data.frame(row = rows, cluster_id = position, audited = picked[[s]][j])
  }
  index <- do.call(rbind, parts)
  d <- z$pool[index$row, ]; rownames(d) <- NULL
  gold <- d$lab_true; gold[!index$audited] <- NA_real_
  list(input = data.frame(covariate = d$covariate, stratum = d$stratum, cluster_id = index$cluster_id, row_id = seq_len(nrow(d)), lab_value = d$lab_value, gold_value = gold),
       formula = lab_value ~ covariate + stratum, predictors = c("covariate", "stratum"))
}

# Rescaled bootstrap, not a package feature: n_h - 1 clusters drawn in each stratum, weighted n_h / (n_h - 1).
# A stratum with one audited cluster keeps it. In memory only, in place of the package's within-stratum resample.
namespace <- asNamespace("catmech")
package_fill <- get("catmech_reference_cluster_fill", envir = namespace)
rescaled_fill <- function(x, y, fit_rows, fill_rows, ids, strata, m) {
  by_cluster <- split(fit_rows, as.character(ids[fit_rows]))
  groups <- split(names(by_cluster), as.character(strata[vapply(by_cluster, `[`, integer(1), 1L)]))
  out <- matrix(NA_real_, length(fill_rows), m)
  for (j in seq_len(m)) {
    draw <- NULL
    for (attempt in seq_len(200L)) {
      rows <- integer(); w <- numeric()
      for (u in groups) {
        nh <- length(u)
        pick <- if (nh == 1L) u else u[sample.int(nh, nh - 1L, replace = TRUE)]
        got <- by_cluster[pick]
        rows <- c(rows, unlist(got, use.names = FALSE)); w <- c(w, rep(if (nh == 1L) 1 else nh / (nh - 1), sum(lengths(got))))
      }
      fit <- stats::lm.wfit(x[rows, , drop = FALSE], y[rows], w)
      if (fit$rank != ncol(x) || sum(w) <= ncol(x)) next
      variance <- sum(w * fit$residuals^2) / (sum(w) - ncol(x))
      if (!is.finite(variance) || variance <= 0) next
      draw <- as.vector(x[fill_rows, , drop = FALSE] %*% fit$coefficients) + stats::rnorm(length(fill_rows), 0, sqrt(variance))
      break
    }
    if (is.null(draw)) stop("lost design support")
    out[, j] <- draw
  }
  out
}
set_fill <- function(f) {
  unlockBinding("catmech_reference_cluster_fill", namespace)
  assign("catmech_reference_cluster_fill", f, envir = namespace)
  lockBinding("catmech_reference_cluster_fill", namespace)
}
VARIANTS <- list(control = list(cluster = "cluster_id", strata = "stratum", fill = package_fill),
                 rows = list(cluster = "row_id", strata = "stratum", fill = package_fill),
                 no_strata = list(cluster = "cluster_id", strata = NULL, fill = package_fill),
                 rescaled = list(cluster = "cluster_id", strata = "stratum", fill = rescaled_fill))
run_variant <- function(v, p, seed) {
  set_fill(v$fill)
  on.exit(set_fill(package_fill), add = TRUE)
  fit <- catmech_reference(p$input, "lab_value", "gold_value", p$predictors, method = "normal", m = config$m, seed = seed,
                           use_recorded = TRUE, cluster = v$cluster, cluster_strata = v$strata)
  pooled <- catmech_pool_lm(fit, p$formula, cluster = "cluster_id")
  pooled[pooled$term == "covariate", ]
}
icc_anova <- function(e, id) {
  id <- factor(id); k <- tapply(e, id, length); g <- nlevels(id); n <- length(e)
  means <- tapply(e, id, mean)
  between <- sum(k * (means - mean(e))^2) / (g - 1); within <- sum((e - means[id])^2) / (n - g)
  k0 <- (n - sum(k^2) / n) / (g - 1)
  value <- (between - within) / (between + (k0 - 1) * within)
  if (is.finite(value)) value else NA_real_
}
audit_icc <- function(p) {
  audited <- !is.na(p$input$gold_value)
  d <- p$input[audited, ]
  fit <- tryCatch(stats::lm(gold_value ~ lab_value + covariate + stratum, data = d), error = function(e) NULL)
  if (is.null(fit) || sum(stats::residuals(fit)^2) < 1e-18) return(NA_real_)
  icc_anova(stats::residuals(fit), d$cluster_id)
}
one_request <- function(z, fraction_index, r) {
  p <- real_case(z, fraction_index, r)
  seed <- as.integer(z$seeds$imputation + z$index * 1e6 + fraction_index * 1e5 + r)
  rows <- lapply(names(VARIANTS), function(name) {
    result <- tryCatch(run_variant(VARIANTS[[name]], p, seed), error = identity)
    failed <- inherits(result, "error")
    data.frame(analysis = z$id, fraction = config$fractions[fraction_index], replicate = r, variant = name, status = if (failed) "refused" else "ok",
               estimate = if (failed) NA_real_ else result$estimate, se = if (failed) NA_real_ else result$std.error, df = if (failed) NA_real_ else result$df,
               lower = if (failed) NA_real_ else result$lower, upper = if (failed) NA_real_ else result$upper,
               between_variance = if (failed) NA_real_ else result$between_variance, audited_icc = NA_real_, stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows)
  out$audited_icc <- audit_icc(p)
  out
}
started <- Sys.time()
records <- list()
for (z in real) for (j in seq_along(config$fractions)) {
  got <- parallel::mclapply(seq_len(config$reps), function(r) one_request(z, j, r), mc.cores = workers, mc.preschedule = TRUE, mc.set.seed = FALSE)
  bad <- vapply(got, inherits, logical(1), what = "try-error")
  if (any(bad)) stop("Worker failure on ", z$id, ": ", conditionMessage(attr(got[[which(bad)[1]]], "condition")))
  records <- c(records, got)
  cat(z$id, config$fractions[j], round(as.numeric(difftime(Sys.time(), started, units = "secs"))), "seconds\n")
}
utils::write.csv(do.call(rbind, records), file.path(outdir, "variants.csv"), row.names = FALSE)

# Within-facility correlation of the imputation model's residuals, on each whole pool.
pool_icc <- do.call(rbind, lapply(real, function(z) {
  d <- z$pool
  e <- stats::residuals(stats::lm(lab_true ~ lab_value + covariate + stratum, data = d))
  by_facility <- tapply(e^2, d$unit, sum)
  keep <- !d$unit %in% names(sort(by_facility, decreasing = TRUE))[1:3]
  trimmed <- stats::residuals(stats::lm(lab_true ~ lab_value + covariate + stratum, data = d[keep, ]))
  data.frame(analysis = z$id, icc = icc_anova(e, d$unit), icc_without_three_facilities = icc_anova(trimmed, d$unit[keep]))
}))
utils::write.csv(pool_icc, file.path(outdir, "pool-icc.csv"), row.names = FALSE)

# The same quantity on very large draws from the study's simulation settings.
a <- list(b0 = 2, b1 = .6, size_log_mean = 3, size_log_sd = .5, x_noise_sd = 3, residual_variance = 100, wrong_rate = .4, rate_beta = c(1, 1.5), error_sd = 6, tracking_slope = .25)
set.seed(20261101)
g <- 200000L; each <- 3L
simulation_icc <- list()
for (errors in c("iid", "persistent", "tracking")) for (icc in c(0, .2, .5)) {
  id <- rep(seq_len(g), each = each)
  size <- rep(stats::rlnorm(g, a$size_log_mean, a$size_log_sd), each = each)
  x <- size + stats::rnorm(g * each, 0, a$x_noise_sd)
  u <- rep(stats::rnorm(g, 0, sqrt(icc * a$residual_variance)), each = each)
  truth <- a$b0 + a$b1 * x + u + stats::rnorm(g * each, 0, sqrt((1 - icc) * a$residual_variance))
  rate <- if (errors == "iid") rep(a$wrong_rate, g * each) else rep(stats::rbeta(g, a$rate_beta[1], a$rate_beta[2]), each = each)
  wrong <- stats::runif(g * each) < rate
  shift <- if (errors == "tracking") a$tracking_slope * (x - mean(x)) else 0
  recorded <- truth + wrong * (shift + stats::rnorm(g * each, 0, a$error_sd))
  simulation_icc[[length(simulation_icc) + 1L]] <- data.frame(errors = errors, target_icc = icc, residual_icc = icc_anova(stats::residuals(stats::lm(truth ~ x + recorded)), id))
}
utils::write.csv(do.call(rbind, simulation_icc), file.path(outdir, "simulation-icc.csv"), row.names = FALSE)
writeLines(utils::capture.output(utils::sessionInfo()), file.path(outdir, "session-info.txt"))
cat("Completed in", round(as.numeric(difftime(Sys.time(), started, units = "secs"))), "seconds\n")
