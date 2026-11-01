#!/usr/bin/env Rscript
# Independent posterior check of catmech_contamination(). A random-walk Metropolis
# sampler targets the posterior with status and target integrated out analytically,
# written here from the model statement. It shares no code with the package sampler,
# which uses a different algorithm, and its answers must agree. Known-wrong
# likelihoods must disagree. This verifies computation, not coverage or adequacy.
# Usage: Rscript tools/check-contamination-oracle.R [NEW_OUTPUT_DIRECTORY]
args <- commandArgs(trailingOnly = TRUE)
file_arg <- grep("^--file=", commandArgs(), value = TRUE)
root <- normalizePath(file.path(dirname(sub("^--file=", "", file_arg[1])), ".."))
outdir <- if (length(args)) args[1] else file.path(root, "validation", "contamination-oracle")
if (dir.exists(outdir) && length(list.files(outdir, all.files = TRUE, no.. = TRUE)))
  stop("Use a new or empty output directory; existing evidence is never overwritten.")
hash_sources <- function() {
  paths <- c(sort(list.files(file.path(root, "R"), "\\.R$", full.names = TRUE)),
    file.path(root, c("DESCRIPTION", "NAMESPACE", "tools/check-contamination-oracle.R")))
  data.frame(file = substring(paths, nchar(root) + 2L),
    sha256 = vapply(paths, function(p) digest::digest(file = p, algo = "sha256"), character(1)), row.names = NULL)
}
hashes_before <- hash_sources()
pkgload::load_all(root, quiet = TRUE)
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
workers <- if (.Platform$OS.type == "windows") 1L else 6L
gibbs_iterations <- 20000L
metropolis_iterations <- 400000L
metropolis_burnin <- 100000L

make_case <- function(name) {
  set.seed(switch(name, exact = 6101L, strata_delta = 6102L, noisy = 6103L))
  n <- 70L
  x <- rnorm(n)
  truth <- 1 + .5 * x + rnorm(n, 0, .5)
  stratum <- if (name == "strata_delta") ifelse(runif(n) < .3, "yes", "no") else rep("all", n)
  rate <- if (name == "strata_delta") ifelse(stratum == "yes", .5, .08) else rep(.2, n)
  wrong <- runif(n) < rate
  w <- truth + ifelse(wrong, rnorm(n, 1, .4), 0)
  w[runif(n) < .15] <- NA_real_
  r <- rep(NA_real_, n)
  r[seq(1L, n, by = 3L)] <- truth[seq(1L, n, by = 3L)]
  sd_r <- if (name == "noisy") .25 else 0
  if (sd_r > 0) r[!is.na(r)] <- r[!is.na(r)] + rnorm(sum(!is.na(r)), 0, sd_r)
  list(name = name, data = data.frame(w = w, r = r, x = x, stratum = stratum),
       reference_sd = sd_r, delta = if (name == "strata_delta") .3 else 0)
}

# The marginal log posterior on (beta, log sigma, logit pi, mu, log tau). `wrong_model`
# selects deliberately incorrect likelihoods for the negative controls.
make_log_posterior <- function(case, prior, wrong_model = "none") {
  d <- case$data
  X <- cbind(1, d$x)
  levels <- sort(unique(d$stratum[!is.na(d$w)]))
  s <- match(d$stratum, levels)
  k <- length(levels)
  rv <- case$reference_sd^2
  rec <- !is.na(d$w); ref <- !is.na(d$r)
  agree <- rec & ref & abs(d$w - d$r) <= 1e-8 * pmax(1, abs(d$r))
  function(theta) {
    beta <- theta[1:2]
    sigma2 <- exp(2 * theta[3])
    logit <- theta[3 + seq_len(k)]
    mu <- theta[4 + k]
    tau2 <- exp(2 * theta[5 + k])
    pi <- plogis(logit)[s]
    m <- as.vector(X %*% beta)
    out <- sum(prior$pi_a * plogis(logit, log.p = TRUE) + prior$pi_b * plogis(-logit, log.p = TRUE))
    out <- out - (prior$nu0 + 1) / 2 * log(tau2) -
      (prior$nu0 * prior$tau0^2 + prior$kappa0 * (mu - prior$mu0)^2) / (2 * tau2)
    if (rv > 0) out <- out - .5 * log(sigma2) - prior$tau0^2 / (2 * sigma2)
    density <- function(z, mean, variance) -.5 * log(2 * pi_const * variance) - (z - mean)^2 / (2 * variance)
    mix <- function(a, b) { top <- pmax(a, b); top + log(exp(a - top) + exp(b - top)) }
    open <- rec & !ref
    centre <- m[open] + case$delta
    spread <- if (wrong_model == "omit_error_variance") sigma2 else sigma2 + tau2
    shift <- if (wrong_model == "omit_error_mean") 0 else mu
    out <- out + sum(mix(log1p(-pi[open]) + density(d$w[open], centre, sigma2),
                         log(pi[open]) + density(d$w[open], centre + shift, spread)))
    if (rv == 0) {
      out <- out + sum(density(d$r[ref], m[ref], sigma2))
      both <- rec & ref
      out <- out + sum(log1p(-pi[both & agree]))
      if (wrong_model != "drop_known_errors") {
        bad <- both & !agree
        out <- out + sum(log(pi[bad]) + density(d$w[bad] - d$r[bad], mu, tau2))
      }
    } else {
      only <- ref & !rec
      out <- out + sum(density(d$r[only], m[only], sigma2 + rv))
      both <- which(ref & rec)
      a <- log1p(-pi[both]) + density(d$w[both], m[both], sigma2) + density(d$r[both], d$w[both], rv)
      # Explicit bivariate normal, written without the package's factorization.
      v11 <- sigma2 + tau2; v22 <- sigma2 + rv; v12 <- sigma2
      det <- v11 * v22 - v12^2
      e1 <- d$w[both] - m[both] - mu; e2 <- d$r[both] - m[both]
      b <- log(pi[both]) - log(2 * pi_const) - .5 * log(det) - (v22 * e1^2 - 2 * v12 * e1 * e2 + v11 * e2^2) / (2 * det)
      out <- out + sum(mix(a, b))
    }
    out
  }
}
pi_const <- base::pi

metropolis <- function(log_post, start, iterations, burnin, seed) {
  set.seed(seed)
  d <- length(start)
  chain <- matrix(NA_real_, iterations, d)
  current <- start
  lp <- log_post(current)
  # Adaptation stops at burn-in, so the retained chain is an ordinary Metropolis chain.
  proposal <- diag(.05^2, d)
  root_cov <- chol(proposal)
  accepted <- 0
  for (i in seq_len(iterations)) {
    candidate <- current + as.vector(stats::rnorm(d) %*% root_cov)
    lp_candidate <- log_post(candidate)
    if (is.finite(lp_candidate) && log(stats::runif(1L)) < lp_candidate - lp) {
      current <- candidate
      lp <- lp_candidate
      if (i > burnin) accepted <- accepted + 1
    }
    chain[i, ] <- current
    if (i <= burnin && i %% 5000L == 0L) {
      proposal <- stats::cov(chain[(i %/% 2L):i, , drop = FALSE]) * 2.38^2 / d + diag(1e-8, d)
      root_cov <- chol(proposal)
    }
  }
  list(draws = chain[(burnin + 1L):iterations, , drop = FALSE], acceptance = accepted / (iterations - burnin))
}

run_oracle <- function(case, prior, start, iterations, burnin, wrong_model = "none", base_seed = 1L) {
  log_post <- make_log_posterior(case, prior, wrong_model)
  chains <- lapply(1:2, function(j) metropolis(log_post, start + stats::rnorm(length(start), 0, .05), iterations, burnin, base_seed + j))
  list(a = chains[[1]]$draws, b = chains[[2]]$draws, acceptance = mean(c(chains[[1]]$acceptance, chains[[2]]$acceptance)))
}
natural <- function(draws, k) {
  out <- draws
  out[, 3] <- exp(draws[, 3])
  out[, 3 + seq_len(k)] <- plogis(draws[, 3 + seq_len(k)])
  out[, 5 + k] <- exp(draws[, 5 + k])
  out
}
moment_table <- function(gibbs, oracle, k, labels, case_name) {
  a <- natural(oracle$a, k); b <- natural(oracle$b, k)
  do.call(rbind, lapply(seq_along(labels), function(j) {
    g <- gibbs[, , labels[j]]
    pair <- cbind(a[, j], b[, j])
    ess_g <- posterior::ess_bulk(g); ess_o <- posterior::ess_bulk(pair)
    sd_g <- stats::sd(as.vector(g)); sd_o <- stats::sd(as.vector(pair))
    # posterior::mcse_sd accounts for kurtosis, which matters for error-size parameters.
    mcse_mean <- sqrt(posterior::mcse_mean(g)^2 + posterior::mcse_mean(pair)^2)
    mcse_sd <- sqrt(posterior::mcse_sd(g)^2 + posterior::mcse_sd(pair)^2)
    data.frame(case = case_name, parameter = labels[j], statistic = c("mean", "sd"),
      gibbs = c(mean(g), sd_g), metropolis = c(mean(pair), sd_o), mcse = c(mcse_mean, mcse_sd),
      gibbs_ess = ess_g, metropolis_ess = ess_o)
  }))
}
fit_case <- function(case, iterations = gibbs_iterations) {
  strata <- if (case$name == "strata_delta") "stratum" else NULL
  catmech_contamination(case$data, "w", "r", "x", corruption_strata = strata,
    reference_sd = case$reference_sd, selection_delta = case$delta, m = 20L,
    seed = 7700L, burnin = 1000L, iterations = iterations, max_iterations = iterations)
}
parameter_labels <- function(fit) {
  all <- dimnames(fit$model_diagnostics$traces)[[3]]
  c(all[startsWith(all, "beta:")], "sigma", all[startsWith(all, "pi:")], all[startsWith(all, "mu:")], all[startsWith(all, "tau:")])
}
# Status and target posteriors for individual rows, from the Metropolis parameter
# draws and closed-form conditionals written here.
row_checks <- function(case, fit, oracle, prior) {
  d <- case$data
  rows <- utils::head(which(!is.na(d$w) & is.na(d$r)), 6L)
  levels <- sort(unique(d$stratum[!is.na(d$w)]))
  k <- length(levels)
  draws <- rbind(oracle$a, oracle$b)
  draws <- draws[seq(1L, nrow(draws), by = 20L), , drop = FALSE]
  do.call(rbind, lapply(rows, function(i) {
    s <- match(d$stratum[i], levels)
    m <- draws[, 1] + draws[, 2] * d$x[i] + case$delta
    sigma2 <- exp(2 * draws[, 3]); tau2 <- exp(2 * draws[, 5 + k]); mu <- draws[, 4 + k]
    pi <- plogis(draws[, 3 + s])
    a <- (1 - pi) * dnorm(d$w[i], m, sqrt(sigma2))
    b <- pi * dnorm(d$w[i], m + mu, sqrt(sigma2 + tau2))
    p <- b / (a + b)
    wrong_mean <- (m / sigma2 + (d$w[i] - mu) / tau2) / (1 / sigma2 + 1 / tau2)
    table <- fit$contamination$corruption[i, ]
    data.frame(case = case$name, row = i, quantity = c("posterior_corruption", "posterior_target"),
      gibbs = c(table$posterior_corruption, table$posterior_target),
      metropolis = c(mean(p), mean((1 - p) * d$w[i] + p * wrong_mean)))
  }))
}
one_case <- function(name) {
  case <- make_case(name)
  fit <- fit_case(case)
  prior <- fit$spec$prior
  labels <- parameter_labels(fit)
  k <- sum(startsWith(labels, "pi:"))
  start <- c(colMeans(apply(fit$model_diagnostics$traces[, , labels[1:2]], 3, as.vector)),
             log(mean(fit$model_diagnostics$traces[, , "sigma"])))
  start <- c(start, stats::qlogis(vapply(labels[startsWith(labels, "pi:")], function(z) mean(fit$model_diagnostics$traces[, , z]), numeric(1))),
             mean(fit$model_diagnostics$traces[, , labels[startsWith(labels, "mu:")]]),
             log(mean(fit$model_diagnostics$traces[, , labels[startsWith(labels, "tau:")]])))
  oracle <- run_oracle(case, prior, unname(start), metropolis_iterations, metropolis_burnin, base_seed = 8800L)
  moments <- moment_table(fit$model_diagnostics$traces, oracle, k, labels, name)
  rows <- row_checks(case, fit, oracle, prior)
  controls <- do.call(rbind, lapply(if (name == "exact") c("omit_error_variance", "omit_error_mean", "drop_known_errors") else character(), function(wrong) {
    bad <- run_oracle(case, prior, unname(start), metropolis_iterations %/% 2L, metropolis_burnin %/% 2L, wrong, base_seed = 9900L)
    z <- moment_table(fit$model_diagnostics$traces, bad, k, labels, name)
    data.frame(case = name, wrong_likelihood = wrong, largest_abs_z = max(abs(z$gibbs - z$metropolis) / z$mcse),
               detected = max(abs(z$gibbs - z$metropolis) / z$mcse) > 5)
  }))
  known <- !is.na(case$data$r)
  reasons <- fit$audit$reason[fit$audit$imputation == 1L]
  contract <- data.frame(case = name,
    check = c("references_preserved_or_modeled", "retained_rows_equal_recorded", "replaced_rows_differ", "original_unchanged", "rhat_at_most_1.01"),
    pass = c(if (case$reference_sd == 0) all(fit$completed[[1]]$w[known] == case$data$r[known]) else !any(fit$completed[[1]]$w[known] == case$data$r[known]),
      all(fit$completed[[1]]$w[reasons == "recorded_retained"] == case$data$w[reasons == "recorded_retained"]),
      all(fit$completed[[1]]$w[reasons == "recorded_replaced"] != case$data$w[reasons == "recorded_replaced"]),
      identical(fit$original, case$data), all(fit$model_diagnostics$summary$rhat <= 1.01)))
  list(moments = moments, rows = rows, controls = controls, contract = contract, acceptance = oracle$acceptance)
}
start_time <- Sys.time()
results <- parallel::mclapply(c("exact", "strata_delta", "noisy"), one_case, mc.cores = min(workers, 3L), mc.set.seed = FALSE)
failed <- vapply(results, inherits, logical(1), what = "try-error")
if (any(failed)) stop("Oracle case failed: ", paste(vapply(results[failed], as.character, ""), collapse = "\n"))
moments <- do.call(rbind, lapply(results, `[[`, "moments"))
moments$abs_z <- abs(moments$gibbs - moments$metropolis) / moments$mcse
moments$pass <- moments$abs_z <= 5
rows <- do.call(rbind, lapply(results, `[[`, "rows"))
rows$difference <- rows$gibbs - rows$metropolis
rows$pass <- abs(rows$difference) <= .015
controls <- do.call(rbind, lapply(results, `[[`, "controls"))
contract <- do.call(rbind, lapply(results, `[[`, "contract"))
if (!identical(hashes_before, hash_sources())) stop("Source changed during the oracle run.")
utils::write.csv(moments, file.path(outdir, "posterior-moments.csv"), row.names = FALSE)
utils::write.csv(rows, file.path(outdir, "row-posteriors.csv"), row.names = FALSE)
utils::write.csv(controls, file.path(outdir, "negative-controls.csv"), row.names = FALSE)
utils::write.csv(contract, file.path(outdir, "data-contract-controls.csv"), row.names = FALSE)
utils::write.csv(hashes_before, file.path(outdir, "source-hashes.csv"), row.names = FALSE)
writeLines(utils::capture.output(utils::sessionInfo()), file.path(outdir, "sessionInfo.txt"))
ok <- all(moments$pass) && all(rows$pass) && all(controls$detected) && all(contract$pass)
elapsed <- round(as.numeric(difftime(Sys.time(), start_time, units = "secs")))
writeLines(c("# Independent contamination-model oracle", "",
  paste0("A random-walk Metropolis sampler on the posterior with status and target integrated out, written in `tools/check-contamination-oracle.R` from the model statement, is compared with the package's slice-within-Gibbs sampler. Three synthetic cases (n = 70) cover one stratum, two strata with a selection shift, and a declared noisy reference. Run time ", elapsed, " seconds."), "",
  paste0("- Posterior means and SDs: ", sum(moments$pass), "/", nrow(moments), " within five combined MCSEs; largest |z| ", formatC(max(moments$abs_z), digits = 2, format = "f"), "."),
  paste0("- Row-level corruption probabilities and target means: ", sum(rows$pass), "/", nrow(rows), " within 0.015; largest difference ", formatC(max(abs(rows$difference)), digits = 4, format = "f"), "."),
  paste0("- Known-wrong likelihoods detected: ", sum(controls$detected), "/", nrow(controls), "; smallest largest-|z| ", formatC(min(controls$largest_abs_z), digits = 1, format = "f"), "."),
  paste0("- Data-contract controls: ", sum(contract$pass), "/", nrow(contract), "."), "",
  "This is implementation evidence for these examples. It says nothing about repeated-sampling coverage, model adequacy or real data."),
  file.path(outdir, "README.md"))
cat(readLines(file.path(outdir, "README.md")), sep = "\n")
if (!ok) stop("Contamination oracle failed; see the CSV files.")
