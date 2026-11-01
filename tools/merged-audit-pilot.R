#!/usr/bin/env Rscript
# Go/no-go simulation for audit allocation under interference; plan in
# validation/merged-audit-pilot/PLAN.md. Ports the kernel, exposure, data-generating
# model and estimators of the CRT manuscript's deposited code as shown in the public
# demo (crt-design-demo, commit 04ee03b). Trial outcomes are simulated.
# Usage: Rscript tools/merged-audit-pilot.R reproduce|pilot|study [NEW_DIRECTORY] [REPS]
args <- commandArgs(trailingOnly = TRUE)
mode <- if (length(args)) match.arg(args[1], c("reproduce", "pilot", "study")) else "pilot"
file_arg <- grep("^--file=", commandArgs(), value = TRUE)
root <- normalizePath(file.path(dirname(sub("^--file=", "", file_arg[1])), ".."))
outdir <- if (length(args) >= 2L) args[2] else file.path(root, "validation", "merged-audit-pilot", mode)
if (dir.exists(outdir) && length(list.files(outdir, all.files = TRUE, no.. = TRUE)))
  stop("Use a new or empty output directory; existing evidence is never overwritten.")
hash_sources <- function() {
  paths <- c(sort(list.files(file.path(root, "R"), "\\.R$", full.names = TRUE)),
    file.path(root, c("DESCRIPTION", "NAMESPACE", "tools/merged-audit-pilot.R", "validation/merged-audit-pilot/PLAN.md",
                      "validation/merged-audit-pilot/sites.csv")))
  paths <- paths[file.exists(paths)]
  data.frame(file = substring(paths, nchar(root) + 2L),
    sha256 = vapply(paths, function(p) digest::digest(file = p, algo = "sha256"), character(1)), row.names = NULL)
}
hashes_before <- hash_sources()
pkgload::load_all(root, quiet = TRUE)
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
workers <- suppressWarnings(as.integer(Sys.getenv("CATMECH_BENCHMARK_WORKERS", "14")))
if (.Platform$OS.type == "windows") workers <- 1L
reps <- if (length(args) >= 3L) as.integer(args[3]) else if (mode == "reproduce") 1000L else if (mode == "pilot") 3L else 400L
sites <- utils::read.csv(file.path(root, "validation/merged-audit-pilot/sites.csv"))
n <- nrow(sites)
effect <- log(.5)
contamination_strength <- .35
wave_one <- 8L
wave_two <- 8L
imputations <- 20L
regimes <- list(moderate = list(rate = .20, draw = function(k) stats::rnorm(k, 0, 1)),
                gross = list(rate = .10, draw = function(k) log(10) * sample(c(-1, 1), k, replace = TRUE)))
directions <- c(0, .5, 1)
policies <- c("random", "reactive", "undirected", "transmitter", "influence", "oracle")
base <- if (mode == "study") 400100000 else 500100000
jsonlite::write_json(list(mode = mode, reps = reps, workers = workers, effect = effect, c = contamination_strength,
  wave_one = wave_one, wave_two = wave_two, m = imputations, directions = directions, base_seed = base),
  file.path(outdir, "settings.json"), pretty = TRUE, auto_unbox = TRUE)

haversine <- function(la1, lo1, la2, lo2) {
  rad <- pi / 180
  dp <- (la2 - la1) * rad; dl <- (lo2 - lo1) * rad
  a <- sin(dp / 2)^2 + cos(la1 * rad) * cos(la2 * rad) * sin(dl / 2)^2
  2 * 6371 * asin(pmin(1, sqrt(a)))
}
distance <- outer(seq_len(n), seq_len(n), function(i, j) haversine(sites$lat[i], sites$lon[i], sites$lat[j], sites$lon[j]))
inverse_square <- 1 / pmax(distance, .05)^2
diag(inverse_square) <- 0
p0 <- sites$prevalence
logit_p <- stats::qlogis(p0)
trials <- sites$examined
strata <- as.numeric(sites$lat >= stats::median(sites$lat))
# K[i, j] is flow from i into j. The symmetric form is the manuscript's gravity kernel.
symmetric_kernel <- function(mass) outer(mass, mass) * inverse_square
directed_kernel <- function(mass, delta) {
  K <- outer(p0 * mass, mass * (1 - p0)) * inverse_square
  (1 - delta) * (K + t(K)) / 2 + delta * K
}
allocate <- function() {
  a <- numeric(n)
  for (s in 0:1) {
    idx <- which(strata == s)
    idx <- idx[sample.int(length(idx))]
    a[idx] <- as.numeric(seq_along(idx) <= ceiling(length(idx) / 2))
  }
  a
}
exposure <- function(K, a) {
  cross <- outer(a, a, "!=")
  kappa <- .5 * colSums(K)
  inflow <- colSums(K * cross)
  s <- 1 - exp(-inflow / kappa)
  lbar <- colSums(K * cross * logit_p) / pmax(inflow, .Machine$double.xmin)
  lbar[inflow <= 0] <- logit_p[inflow <= 0]
  list(s = s, dspill = s * (lbar - logit_p), inflow = inflow, kappa = kappa, cross = cross)
}
simulate_outcome <- function(K, a) {
  x <- exposure(K, a)
  eta0 <- logit_p + effect * a
  etac <- colSums(K * x$cross * eta0) / pmax(x$inflow, .Machine$double.xmin)
  etac[x$inflow <= 0] <- eta0[x$inflow <= 0]
  stats::rbinom(n, trials, stats::plogis((1 - contamination_strength * x$s) * eta0 + contamination_strength * x$s * etac))
}
fit_arm <- function(y, a, x = NULL) {
  design <- cbind(1, a, strata, logit_p)
  if (!is.null(x)) design <- cbind(design, x$s, a * x$s, x$dspill)
  fit <- tryCatch(suppressWarnings(stats::glm.fit(design, cbind(y, trials - y), family = stats::binomial(),
    control = list(epsilon = 1e-8, maxit = 50))), error = function(e) NULL)
  if (is.null(fit) || !fit$converged || fit$rank != ncol(design)) return(NULL)
  list(estimate = unname(fit$coefficients[2]), variance = chol2inv(fit$R)[2, 2])
}
# Sensitivity of each site's cross-arm exposure share to the recorded mass of one site.
# A site's own mass cancels from its own exposure, so only downstream sites contribute.
leverage <- function(K, a) {
  x <- exposure(K, a)
  ratio <- x$inflow / x$kappa
  share <- sweep(K, 2L, x$kappa, "/")
  gap <- abs(sweep(x$cross * 1, 2L, .5 * ratio, "-"))
  score <- sweep(share * gap, 2L, 1 - x$s, "*")
  diag(score) <- 0
  rowSums(score)
}
pool_interval <- function(estimates, variances) {
  m <- length(estimates)
  qbar <- mean(estimates); ubar <- mean(variances); b <- stats::var(estimates)
  total <- ubar + (1 + 1 / m) * b
  df <- if (b > 0) (m - 1) * (1 + ubar / ((1 + 1 / m) * b))^2 else Inf
  half <- stats::qt(.975, df) * sqrt(total)
  c(estimate = qbar, se = sqrt(total), lower = qbar - half, upper = qbar + half)
}
wald_interval <- function(fit) {
  if (is.null(fit)) return(NULL)
  half <- 1.96 * sqrt(fit$variance)
  c(estimate = fit$estimate, se = sqrt(fit$variance), lower = fit$estimate - half, upper = fit$estimate + half)
}
reproduce_replicate <- function(r) {
  set.seed(as.integer(base + 5000000 + r))
  a <- allocate()
  K <- symmetric_kernel(trials)
  y <- simulate_outcome(K, a)
  x <- exposure(K, a)
  rbind(c(method = 1, wald_interval(fit_arm(y, a))), c(method = 2, wald_interval(fit_arm(y, a, x))))
}
if (mode == "reproduce") {
  out <- do.call(rbind, parallel::mclapply(seq_len(reps), reproduce_replicate, mc.cores = workers, mc.set.seed = FALSE))
  deposited <- data.frame(estimator = c("naive", "adjusted"), bias = c(.289, .005), rmse = c(.372, .149), cover = c(.168, .963))
  tab <- do.call(rbind, lapply(1:2, function(k) {
    g <- out[out[, "method"] == k, ]
    error <- g[, "estimate"] - effect
    cover <- g[, "lower"] <= effect & effect <= g[, "upper"]
    data.frame(estimator = deposited$estimator[k], replicates = nrow(g), port_bias = mean(error), bias_mcse = stats::sd(error) / sqrt(nrow(g)),
      deposited_bias = deposited$bias[k], port_rmse = sqrt(mean(error^2)), deposited_rmse = deposited$rmse[k],
      port_coverage = mean(cover), coverage_mcse = sqrt(mean(cover) * (1 - mean(cover)) / nrow(g)), deposited_coverage = deposited$cover[k])
  }))
  utils::write.csv(tab, file.path(outdir, "port-check.csv"), row.names = FALSE)
  print(tab, digits = 3, row.names = FALSE)
  quit(save = "no")
}
site_table <- function(w, r) data.frame(w = w, r = r, elevation = as.numeric(scale(sites$elevation_m)),
                                        latitude = as.numeric(scale(sites$lat)), logit_p = logit_p)
fit_recording <- function(w, r, seed) {
  tryCatch(suppressWarnings(catmech_contamination(site_table(w, r), "w", "r", c("elevation", "latitude", "logit_p"),
    m = imputations, seed = seed, iterations = 1200L, max_iterations = 4800L)), error = identity)
}
analyze <- function(w, r, y, a, delta, seed) {
  # Plug-in analysis: audited masses where known, recorded masses elsewhere.
  plug <- ifelse(is.na(r), w, r)
  plugin <- wald_interval(fit_arm(y, a, exposure(directed_kernel(exp(plug), delta), a)))
  recording <- fit_recording(w, r, seed)
  pooled <- NULL
  status <- "ok"
  if (inherits(recording, "error")) {
    status <- if (inherits(recording, "catmech_convergence_error")) "convergence_failure" else "error"
  } else {
    fits <- lapply(recording$completed, function(d) fit_arm(y, a, exposure(directed_kernel(exp(d$w), delta), a)))
    if (any(vapply(fits, is.null, logical(1)))) status <- "analysis_failure" else
      pooled <- pool_interval(vapply(fits, `[[`, 0, "estimate"), vapply(fits, `[[`, 0, "variance"))
  }
  list(plugin = plugin, pooled = pooled, status = status, recording = recording)
}
run_replicate <- function(regime_name, delta, r) {
  regime <- regimes[[regime_name]]
  # Allocation, errors and first-wave audit are shared across directions and policies.
  set.seed(as.integer(base + r))
  a <- allocate()
  wrong <- stats::runif(n) < regime$rate
  shift <- ifelse(wrong, regime$draw(n), 0)
  truth <- log(trials)
  w <- truth + shift
  first <- sort(sample.int(n, wave_one))
  set.seed(as.integer(base + 1000000 + r))
  y <- simulate_outcome(directed_kernel(trials, delta), a)
  seed <- as.integer(base + 2000000 + r)
  r1 <- rep(NA_real_, n)
  r1[first] <- truth[first]
  rows <- list()
  add <- function(policy, result) {
    for (analysis in c("pooled", "plugin")) {
      z <- result[[analysis]]
      rows[[length(rows) + 1L]] <<- data.frame(regime = regime_name, direction = delta, replicate = r, policy = policy, analysis = analysis,
        status = if (!is.null(z)) "ok" else if (analysis == "plugin") "analysis_failure" else result$status,
        estimate = if (is.null(z)) NA_real_ else z[["estimate"]], se = if (is.null(z)) NA_real_ else z[["se"]],
        lower = if (is.null(z)) NA_real_ else z[["lower"]], upper = if (is.null(z)) NA_real_ else z[["upper"]],
        errors_found = if (is.null(result$found)) NA_real_ else result$found, stringsAsFactors = FALSE)
    }
  }
  recorded_only <- wald_interval(fit_arm(y, a, exposure(directed_kernel(exp(w), delta), a)))
  full <- wald_interval(fit_arm(y, a, exposure(directed_kernel(trials, delta), a)))
  add("recorded_only", list(plugin = recorded_only, pooled = recorded_only, status = "ok", found = 0))
  add("full_audit", list(plugin = full, pooled = full, status = "ok", found = sum(wrong)))
  first_result <- analyze(w, r1, y, a, delta, seed)
  first_result$found <- sum(wrong[first])
  add("first_wave_only", first_result)
  if (inherits(first_result$recording, "error")) {
    for (policy in policies) add(policy, list(plugin = NULL, pooled = NULL, status = "first_wave_failure"))
    return(do.call(rbind, rows))
  }
  table <- first_result$recording$contamination$corruption
  expected_error <- abs(w - table$posterior_target)
  plug <- ifelse(is.na(r1), w, r1)
  K_plug <- directed_kernel(exp(plug), delta)
  directed_leverage <- leverage(K_plug, a)
  undirected_leverage <- leverage((K_plug + t(K_plug)) / 2, a)
  base_fit <- fit_arm(y, a, exposure(K_plug, a))
  influence <- vapply(seq_len(n), function(i) {
    if (is.null(base_fit) || !is.na(r1[i])) return(0)
    moved <- plug
    moved[i] <- moved[i] + 1
    fit <- fit_arm(y, a, exposure(directed_kernel(exp(moved), delta), a))
    if (is.null(fit)) 0 else abs(fit$estimate - base_fit$estimate)
  }, numeric(1))
  set.seed(as.integer(base + 3000000 + r))
  open <- setdiff(seq_len(n), first)
  tie <- stats::runif(n) * 1e-9
  scores <- list(random = stats::runif(n), reactive = expected_error + tie,
    undirected = expected_error * undirected_leverage + tie, transmitter = expected_error * directed_leverage + tie,
    influence = expected_error * influence + tie, oracle = abs(shift) * directed_leverage + tie)
  for (policy in policies) {
    chosen <- open[order(scores[[policy]][open], decreasing = TRUE)[seq_len(wave_two)]]
    r2 <- r1
    r2[chosen] <- truth[chosen]
    result <- analyze(w, r2, y, a, delta, seed + 1L)
    result$found <- sum(wrong[c(first, chosen)])
    add(policy, result)
  }
  do.call(rbind, rows)
}
start <- Sys.time()
results <- list()
# Three directions under the moderate regime, and full direction only under the gross control.
cells <- list(c("moderate", 0), c("moderate", .5), c("moderate", 1), c("gross", 1))
for (cell in cells) {
  regime_name <- cell[1]
  delta <- as.numeric(cell[2])
  batch <- parallel::mclapply(seq_len(reps), function(r) run_replicate(regime_name, delta, r),
                              mc.cores = workers, mc.preschedule = TRUE, mc.set.seed = FALSE)
  failed <- vapply(batch, inherits, logical(1), what = "try-error")
  if (any(failed)) stop(as.character(batch[[which(failed)[1]]]))
  results <- c(results, batch)
  cat(regime_name, delta, "done; elapsed", round(as.numeric(difftime(Sys.time(), start, units = "secs"))), "seconds\n"); flush.console()
}
if (!identical(hashes_before, hash_sources())) stop("Source or plan changed during run; evidence is not verified.")
results <- do.call(rbind, results)
results$error <- results$estimate - effect
results$covered <- results$lower <= effect & effect <= results$upper
utils::write.csv(results, file.path(outdir, "replicates.csv"), row.names = FALSE)
mcse <- function(x) if (length(x) > 1L) stats::sd(x) / sqrt(length(x)) else NA_real_
summary_table <- do.call(rbind, lapply(split(results, interaction(results$regime, results$direction, results$policy, results$analysis, drop = TRUE)), function(g) {
  good <- g[g$status == "ok", ]
  data.frame(regime = g$regime[1], direction = g$direction[1], policy = g$policy[1], analysis = g$analysis[1], requested = nrow(g), returned = nrow(good),
    bias = mean(good$error), bias_mcse = mcse(good$error), mse = mean(good$error^2), mse_mcse = mcse(good$error^2), rmse = sqrt(mean(good$error^2)),
    empirical_sd = stats::sd(good$estimate), mean_se = mean(good$se), coverage = mean(good$covered),
    coverage_mcse = sqrt(mean(good$covered) * (1 - mean(good$covered)) / nrow(good)), mean_width = mean(good$upper - good$lower),
    mean_errors_found = mean(g$errors_found, na.rm = TRUE))
}))
rownames(summary_table) <- NULL
utils::write.csv(summary_table, file.path(outdir, "summary.csv"), row.names = FALSE)
contrast <- function(regime_name, delta, analysis, new, old) {
  pick <- function(policy) {
    g <- results[results$regime == regime_name & results$direction == delta & results$analysis == analysis & results$policy == policy & results$status == "ok", ]
    g[c("replicate", "error", "covered")]
  }
  both <- merge(pick(new), pick(old), by = "replicate", suffixes = c("_new", "_old"))
  d <- both$error_new^2 - both$error_old^2
  old_mse <- mean(both$error_old^2)
  data.frame(regime = regime_name, direction = delta, analysis = analysis, new = new, old = old, pairs = nrow(both),
    mse_new = mean(both$error_new^2), mse_old = old_mse, mse_ratio = mean(both$error_new^2) / old_mse, mse_difference = mean(d), difference_mcse = mcse(d),
    z = mean(d) / mcse(d), coverage_difference = mean(both$covered_new - both$covered_old))
}
pairs_wanted <- list(c("full_audit", "recorded_only"), c("transmitter", "reactive"), c("transmitter", "random"), c("undirected", "reactive"),
                     c("transmitter", "undirected"), c("influence", "transmitter"), c("influence", "reactive"), c("oracle", "transmitter"), c("oracle", "reactive"))
contrasts <- do.call(rbind, lapply(cells, function(cell) do.call(rbind, lapply(c("pooled", "plugin"), function(analysis)
  do.call(rbind, lapply(pairs_wanted, function(z) contrast(cell[1], as.numeric(cell[2]), analysis, z[1], z[2])))))))
utils::write.csv(contrasts, file.path(outdir, "paired-contrasts.csv"), row.names = FALSE)
utils::write.csv(hashes_before, file.path(outdir, "source-hashes.csv"), row.names = FALSE)
# The verdict rule fixed in PLAN.md, applied mechanically to the pooled analysis.
get_row <- function(regime_name, delta, new, old) contrasts[contrasts$regime == regime_name & contrasts$direction == delta & contrasts$analysis == "pooled" & contrasts$new == new & contrasts$old == old, ]
verdict_for <- function(regime_name = "moderate") {
  audits_matter <- with(get_row(regime_name, 1, "full_audit", "recorded_only"), mse_ratio <= .9 && z <= -2.5)
  reactive <- get_row(regime_name, 1, "transmitter", "reactive")
  random <- get_row(regime_name, 1, "transmitter", "random")
  network <- reactive$mse_ratio <= .9 && reactive$z <= -2.5 && random$mse_ratio < 1 && random$z <= -2.5
  direction_full <- get_row(regime_name, 1, "transmitter", "undirected")
  direction_half <- get_row(regime_name, .5, "transmitter", "undirected")
  direction <- direction_full$z <= -2 && direction_full$mse_difference < direction_half$mse_difference
  data.frame(regime = regime_name, audits_matter = audits_matter, network_beats_reactive_and_random = network, direction_adds_value = direction,
    verdict = if (!audits_matter || !network) "NO-GO" else if (direction) "GO" else "PARTIAL")
}
verdicts <- verdict_for()
utils::write.csv(verdicts, file.path(outdir, "verdict.csv"), row.names = FALSE)
writeLines(c(paste("Completed UTC:", format(Sys.time(), tz = "UTC", usetz = TRUE)),
  paste("Elapsed seconds:", round(as.numeric(difftime(Sys.time(), start, units = "secs")), 1)), capture.output(sessionInfo())), file.path(outdir, "session-info.txt"))
print(verdicts, row.names = FALSE)
cat("Output:", normalizePath(outdir), "\n")
