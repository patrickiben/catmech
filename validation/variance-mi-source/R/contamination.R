#' Complete a target using a reference-calibrated recording-error model
#'
#' Uses the recorded measurement without trusting it. Each recorded value is
#' either an exact recording of the true target or a recording with additive
#' normal error. Rows with both a recorded value and a reference calibrate that
#' mixture, and every other recorded value receives a posterior probability of
#' being wrong. This keeps the information in correct records, which
#' [catmech_reference()] discards by default.
#'
#' @param data A data frame; raw target measurements remain in `original`.
#' @param value Name of the numeric recorded target to replace in completions.
#' @param reference Name of a numeric reference column, observed only where
#'   independently measured.
#' @param predictors Explicit complete predictor names, excluding `value` and
#'   `reference`. Include the scientific outcome and relevant design variables.
#' @param corruption_strata Optional categorical columns with separate
#'   corruption rates, usually a quality flag and any site or source effect.
#'   These must be complete on recorded rows and have calibration rows.
#' @param error_strata Optional categorical columns with a separate error mean
#'   and standard deviation, such as a site with a systematic recording shift.
#' @param agreement_tol Recorded values within this distance of the reference
#'   count as correct, for example half a rounding unit. The default requires
#'   numerical equality apart from a relative fuzz of `1e-8`.
#' @param reference_sd Declared standard deviation of independent reference
#'   measurement error. The default treats references as exact. A positive value
#'   is an assumption to vary in a sensitivity analysis; it is not estimated.
#'   It cannot be combined with `agreement_tol`.
#' @param selection_delta Sensitivity shift added to the target mean of
#'   unreferenced rows. Zero assumes ignorable reference selection. Compare
#'   several values with [catmech_sensitivity()]; the data cannot identify it.
#' @param prior Optional named list overriding `pi_a` and `pi_b` (Beta prior for
#'   each corruption rate, default 0.5 and 0.5), `mu0`, `kappa0`, `nu0` and `tau0`
#'   (normal-scaled-inverse-chi-squared prior for each error mean and variance,
#'   defaults 0, 0.1, 2 and the residual SD of the reference regression).
#' @param m Number of completions, at least two.
#' @param seed Nonnegative integer seed, with caller RNG restored.
#' @param burnin Discarded Gibbs iterations per chain.
#' @param iterations Retained iterations per chain, at least 100.
#' @param max_iterations If the convergence screen fails, the sampler restarts with
#'   twice as many iterations, up to this limit, and then refuses. The manifest
#'   records `iterations_used`. Set equal to `iterations` to disable restarts.
#' @param chains Number of Gibbs chains, at least four.
#' @return A `catmech_result`. Its `contamination` element holds `corruption`
#'   (one row per record, with `posterior_corruption`), `parameters` (posterior
#'   summaries of the recording and target models) and `recording_audit` (from
#'   [catmech_recording_audit()]).
#' @details
#' The model is `T | Z ~ Normal(Z beta, sigma^2)`, `A ~ Bernoulli(pi[stratum])`,
#' `W = T` when `A = 0`, and `W | T ~ Normal(T + mu[stratum], tau[stratum]^2)` when
#' `A = 1`. Each scan updates the corruption rates and error parameters by slice
#' sampling with `A` and `T` integrated out, draws `A` and then `T | A`, and draws the
#' regression parameters from their conjugate conditional. With exact references,
#' `A` is known on calibration rows. The reference design must have full column rank,
#' which keeps the regression posterior proper under `p(beta, sigma^2)` proportional
#' to `1/sigma^2`.
#'
#' Reference selection and missing recorded values must be ignorable given the
#' predictors, recorded value and strata. Selection may therefore depend on the
#' recorded value, as in a targeted audit of unusual records; [catmech_reference()]
#' does not support that. Errors must be nondifferential given the strata: they may
#' not depend on the outcome or true value once the strata are fixed.
#'
#' Recorded values that look unusual relative to `Z beta` are more likely to be
#' classified as errors, so a non-normal or nonlinear target model can misclassify
#' genuine extremes. Flags computed from the recorded value itself double count
#' that information. Run [catmech_holdout_check()] and compare with
#' [catmech_reference()] rather than assuming the model holds. This is a
#' contamination model used with validation data, not a new identification result.
#' Without calibration rows the call is refused.
#' @export
#' @examples
#' d <- simulate_catmech_lab(200, seed = 8, corruption_probability = .5)
#' d$gold <- NA_real_
#' d$gold[seq(1, 200, 4)] <- d$lab_true[seq(1, 200, 4)]
#' d[c("lab_true", "corrupted")] <- NULL
#' fit <- catmech_contamination(d, "lab_value", "gold",
#'   c("age", "inflammation", "outcome"), corruption_strata = "outlier_flag",
#'   m = 5, seed = 12, iterations = 200)
#' fit$contamination$parameters[1:3, ]
catmech_contamination <- function(data, value, reference, predictors,
                                  corruption_strata = NULL, error_strata = NULL,
                                  agreement_tol = 0, reference_sd = 0,
                                  selection_delta = 0, prior = list(), m = 20L,
                                  seed = 1L, burnin = 200L, iterations = 300L,
                                  max_iterations = 1200L, chains = 4L) {
  catmech_reference_inputs(data, value, predictors, m, seed, 1L,
                           excluded = c(value, reference))
  catmech_contamination_reference_column(data, value, reference)
  check_count(burnin, "burnin", 0L)
  check_count(iterations, "iterations", 100L)
  check_count(max_iterations, "max_iterations", iterations)
  check_count(chains, "chains", 4L)
  if (m > iterations * chains) stop("m cannot exceed the number of retained Gibbs draws.", call. = FALSE)
  nonnegative <- function(z) is.numeric(z) && length(z) == 1L && is.finite(z) && z >= 0
  if (!nonnegative(agreement_tol) || !nonnegative(reference_sd))
    stop("agreement_tol and reference_sd must be single nonnegative finite numbers.", call. = FALSE)
  if (agreement_tol > 0 && reference_sd > 0)
    stop("agreement_tol applies only to exact references; choose it or reference_sd.", call. = FALSE)
  if (!is.numeric(selection_delta) || length(selection_delta) != 1L || !is.finite(selection_delta))
    stop("selection_delta must be one finite number.", call. = FALSE)
  W <- data[[value]]
  R <- data[[reference]]
  known <- !is.na(R)
  recorded <- !is.na(W)
  calibration <- known & recorded
  if (!any(known)) stop("At least one independently observed reference value is required.", call. = FALSE)
  if (!any(calibration))
    stop("Calibrating the recording model requires rows with both recorded and reference values.", call. = FALSE)
  if (length(intersect(c(corruption_strata, error_strata), c(value, reference))))
    stop("Strata must not use the target or reference columns.", call. = FALSE)
  corruption <- catmech_contamination_strata(data, corruption_strata, recorded, "corruption_strata")
  errors <- catmech_contamination_strata(data, error_strata, recorded, "error_strata")
  unreferenced <- recorded & !known
  for (stratum in list(corruption, errors)) {
    needed <- unique(stratum$index[unreferenced])
    absent <- setdiff(needed, unique(stratum$index[calibration]))
    if (length(absent))
      stop("No calibration rows support stratum: ", paste(stratum$levels[absent], collapse = ", "), call. = FALSE)
  }
  agree <- catmech_contamination_agreement(W, R, agreement_tol)
  audit_by <- unique(c(corruption_strata, error_strata))
  recording_audit <- catmech_recording_audit(data, value, reference, by = audit_by,
                                             agreement_tol = agreement_tol)
  if (identical(recording_audit$status, "mostly_disagreeing") && reference_sd == 0)
    warning("Nearly all calibration rows disagree. The exact-recording component is then barely used; ",
            "check units, agreement_tol and reference_sd.", call. = FALSE)
  complete_reference <- all(known) && reference_sd == 0
  X <- if (complete_reference) NULL else catmech_model_matrix(data, predictors)
  resolved <- NULL
  iterations_used <- 0L
  restore <- preserve_rng()
  on.exit(restore(), add = TRUE)
  set.seed(as.integer(seed))
  if (complete_reference) {
    draws <- matrix(rep(R, m), nrow = nrow(data), ncol = m)
    model_diagnostics <- catmech_reference_independent_diagnostics("reference_complete")
    probability <- ifelse(recorded, as.numeric(!agree), NA_real_)
    posterior_target <- R
    parameters <- data.frame(parameter = character(), mean = numeric(), sd = numeric(),
                             lower = numeric(), upper = numeric())
  } else {
    reference_design <- catmech_model_design(X[known, , drop = FALSE])
    reference_fit <- catmech_model_fit(reference_design, R[known])
    resolved <- catmech_contamination_prior(prior, sqrt(reference_fit$sse / reference_design$df))
    problem <- list(X = X, W = W, R = R, known = known, recorded = recorded, agree = agree,
      reference_sd = reference_sd, delta = selection_delta, corruption = corruption, errors = errors,
      reference_design = reference_design, reference_fit = reference_fit,
      full_design = catmech_model_design(X))
    iterations_used <- iterations
    repeat {
      sampled <- catmech_contamination_sampler(problem, resolved, m, burnin, iterations_used, chains)
      model_diagnostics <- catmech_model_chain_diagnostics(sampled$traces)
      if (model_diagnostics$convergence_pass || iterations_used >= max_iterations) break
      iterations_used <- min(2L * iterations_used, as.integer(max_iterations))
    }
    model_diagnostics$engine <- "gibbs_contamination"
    if (!model_diagnostics$convergence_pass) {
      condition <- structure(list(
        message = paste0("Contamination model failed convergence checks (R-hat <= 1.05 and bulk/tail ESS >= ",
                         100L * chains, " required). ",
                         "Inspect condition$model_diagnostics and increase iterations/burnin or revise the model."),
        call = NULL, model_diagnostics = model_diagnostics),
        class = c("catmech_convergence_error", "error", "condition"))
      stop(condition)
    }
    draws <- sampled$draws
    probability <- sampled$probability
    posterior_target <- sampled$posterior_target
    parameters <- catmech_contamination_parameter_summary(sampled$traces)
  }
  exact <- reference_sd == 0
  spec <- structure(list(version = "contamination1", engine = "reference_contamination",
    value = value, predictors = predictors, reference = reference,
    corruption_strata = corruption_strata, error_strata = error_strata,
    agreement_tol = agreement_tol, reference_sd = reference_sd,
    selection_delta = selection_delta, prior = resolved,
    burnin = as.integer(burnin), iterations = as.integer(iterations),
    max_iterations = as.integer(max_iterations), chains = as.integer(chains),
    target_model = "normal_target_given_complete_predictors",
    target_prior = if (exact) "p(beta,sigma_squared) proportional to 1/sigma_squared" else
      "flat beta; scaled-inverse-chi-squared(1, reference residual variance) for sigma_squared",
    recording_model = "exact_recording_or_additive_normal_error",
    reference_policy = if (exact) "observed_gold_preserved_exactly" else "declared_additive_normal_reference_error",
    reference_selection_assumption = "ignorable_given_predictors_recorded_value_and_strata",
    recorded_missingness_assumption = "ignorable_given_predictors_and_strata",
    error_assumption = "nondifferential_given_strata"),
    class = "catmech_contamination_spec")
  audits <- lapply(seq_len(m), function(j) {
    retained <- recorded & !is.na(draws[, j]) & draws[, j] == W
    exact_known <- known & exact
    reason <- ifelse(exact_known, "gold_reference", ifelse(!recorded,
      ifelse(known, "noisy_reference_model", "reference_missing"),
      ifelse(retained, "recorded_retained", "recorded_replaced")))
    row <- catmech_reference_audit(W, draws[, j], seq_len(nrow(data)), j, reason,
      reconstructed = !exact_known & !retained)
    row$reference_observed <- known
    row
  })
  diagnostics <- do.call(rbind, lapply(seq_len(m), function(j) {
    reason <- audits[[j]]$reason
    data.frame(imputation = j, requested = sum(!(known & exact)),
      reconstructed = sum(audits[[j]]$reconstructed), changed = sum(audits[[j]]$changed, na.rm = TRUE),
      reference_observed = sum(known), recorded_retained = sum(reason == "recorded_retained"),
      recorded_replaced = sum(reason == "recorded_replaced"),
      reference_missing = sum(reason == "reference_missing"), fallback_draws = 0L)
  }))
  result <- catmech_reference_assemble(data, value, predictors, if (exact) R else rep(NA_real_, nrow(data)),
    draws, audits, diagnostics, spec, m, seed, "conditional_contamination_model", model_diagnostics)
  result$manifest$reference_rows <- which(known)
  result$manifest$reference_column <- reference
  result$manifest$posterior_version <- as.character(utils::packageVersion("posterior"))
  result$manifest$iterations_used <- iterations_used
  result$contamination <- list(
    corruption = data.frame(row = seq_len(nrow(data)), recorded = W, reference = R,
      corruption_stratum = corruption$levels[corruption$index],
      evidence = ifelse(!recorded, "no_recorded_value", ifelse(!known | !exact, "posterior_model",
        ifelse(agree, "reference_agrees", "reference_disagrees"))),
      posterior_corruption = probability, posterior_target = posterior_target,
      stringsAsFactors = FALSE),
    parameters = parameters, recording_audit = recording_audit)
  result$mids$call <- match.call()
  result
}

catmech_contamination_reference_column <- function(data, value, reference) {
  if (!is.character(reference) || length(reference) != 1L || is.na(reference) ||
      !nzchar(reference) || !reference %in% names(data) || reference == value)
    stop("reference must name an existing column different from the target.", call. = FALSE)
  gold <- data[[reference]]
  if (!is.numeric(gold) || any(!is.na(gold) & !is.finite(gold)))
    stop("reference must contain independently measured finite numeric values or NA.", call. = FALSE)
  invisible(TRUE)
}

catmech_contamination_agreement <- function(W, R, tol) {
  abs(W - R) <= tol + 1e-8 * pmax(1, abs(R))
}

catmech_contamination_strata <- function(data, columns, rows, label) {
  if (is.null(columns)) return(list(index = rep(1L, nrow(data)), levels = "all"))
  if (!is.character(columns) || !length(columns) || anyNA(columns) ||
      anyDuplicated(columns) || any(!columns %in% names(data)))
    stop(label, " must name unique existing columns.", call. = FALSE)
  parts <- lapply(data[columns], as.character)
  if (any(vapply(parts, function(z) anyNA(z[rows]), logical(1))))
    stop(label, " columns must be complete on recorded rows.", call. = FALSE)
  key <- do.call(paste, c(unname(parts), sep = " | "))
  key[!rows] <- NA_character_
  labels <- sort(unique(key[rows]))
  list(index = match(key, labels), levels = labels)
}

catmech_contamination_prior <- function(prior, scale) {
  defaults <- list(pi_a = 0.5, pi_b = 0.5, mu0 = 0, kappa0 = 0.1, nu0 = 2, tau0 = scale)
  if (!is.list(prior) || (length(prior) && (is.null(names(prior)) || !all(names(prior) %in% names(defaults)))))
    stop("prior must be a named list using: ", paste(names(defaults), collapse = ", "), call. = FALSE)
  defaults[names(prior)] <- prior
  valid <- vapply(defaults, function(z) is.numeric(z) && length(z) == 1L && is.finite(z), logical(1))
  if (!all(valid) || any(unlist(defaults[c("pi_a", "pi_b", "kappa0", "nu0", "tau0")]) <= 0))
    stop("Prior values must be finite, with positive pi_a, pi_b, kappa0, nu0 and tau0.", call. = FALSE)
  defaults
}

# One univariate slice-sampling update with stepping out (Neal, 2003).
catmech_slice <- function(x, log_density, width, max_steps = 16L) {
  level <- log_density(x) - stats::rexp(1L)
  if (!is.finite(level)) return(x)
  left <- x - stats::runif(1L) * width
  right <- left + width
  steps_left <- floor(stats::runif(1L) * max_steps)
  steps_right <- max_steps - 1L - steps_left
  while (steps_left > 0L && isTRUE(log_density(left) > level)) {
    left <- left - width
    steps_left <- steps_left - 1L
  }
  while (steps_right > 0L && isTRUE(log_density(right) > level)) {
    right <- right + width
    steps_right <- steps_right - 1L
  }
  for (attempt in seq_len(100L)) {
    proposal <- stats::runif(1L, left, right)
    if (isTRUE(log_density(proposal) > level)) return(proposal)
    if (proposal < x) left <- proposal else right <- proposal
  }
  x
}

# Densities of latent recording rows with the target integrated out. Under A = 0
# a row has density exp(base). Under A = 1 it has density
# exp(add) * Normal(W; mean1 + mu, var1 + tau^2).
catmech_contamination_row_terms <- function(w, r, centre, sigma2, reference_variance) {
  if (reference_variance == 0 || !length(r)) {
    return(list(base = stats::dnorm(w, centre, sqrt(sigma2), log = TRUE), add = rep(0, length(w)),
                mean1 = centre, var1 = rep(sigma2, length(w))))
  }
  referenced <- !is.na(r)
  gain <- sigma2 / (sigma2 + reference_variance)
  base <- stats::dnorm(w, centre, sqrt(sigma2), log = TRUE)
  add <- rep(0, length(w))
  mean1 <- centre
  var1 <- rep(sigma2, length(w))
  base[referenced] <- base[referenced] + stats::dnorm(r[referenced], w[referenced], sqrt(reference_variance), log = TRUE)
  add[referenced] <- stats::dnorm(r[referenced], centre[referenced], sqrt(sigma2 + reference_variance), log = TRUE)
  mean1[referenced] <- centre[referenced] + gain * (r[referenced] - centre[referenced])
  var1[referenced] <- sigma2 - gain * sigma2
  list(base = base, add = add, mean1 = mean1, var1 = var1)
}

catmech_contamination_sampler <- function(s, prior, m, burnin, iterations, chains) {
  n <- length(s$W)
  exact <- s$reference_sd == 0
  reference_variance <- s$reference_sd^2
  latent <- which(s$recorded & (!s$known | !exact))
  free <- which(!s$recorded & !s$known)
  reference_only <- if (exact) integer() else which(!s$recorded & s$known)
  scored <- which(s$recorded)
  imputed <- which(!s$known | !exact)
  w_latent <- s$W[latent]
  r_latent <- s$R[latent]
  s_latent <- s$corruption$index[latent]
  g_latent <- s$errors$index[latent]
  unreferenced <- !s$known[latent]
  delta_latent <- s$delta * unreferenced
  extra_precision <- if (exact) rep(0, length(latent)) else (!unreferenced) / reference_variance
  extra_information <- if (exact) rep(0, length(latent)) else (!unreferenced) * ifelse(is.na(r_latent), 0, r_latent) / reference_variance
  delta_all <- s$delta * !s$known
  projector <- backsolve(s$full_design$R, t(qr.Q(s$full_design$qr)))
  pivot <- s$full_design$qr$pivot
  n_corruption <- length(s$corruption$levels)
  n_errors <- length(s$errors$levels)
  rows_s <- lapply(seq_len(n_corruption), function(j) which(s_latent == j))
  rows_g <- lapply(seq_len(n_errors), function(j) which(g_latent == j))
  # Exact references reveal the recording status and error on calibration rows.
  # Those rows enter the recording posterior only through these fixed summaries.
  target_start <- s$R
  known_scored <- which(s$known & s$recorded)
  known_k <- known_n <- numeric(n_corruption)
  known_count <- known_mean <- known_ss <- numeric(n_errors)
  if (exact) {
    known_n <- tabulate(s$corruption$index[known_scored], n_corruption)
    wrong <- known_scored[!s$agree[known_scored]]
    known_k <- tabulate(s$corruption$index[wrong], n_corruption)
    for (g in seq_len(n_errors)) {
      rows <- wrong[s$errors$index[wrong] == g]
      e <- s$W[rows] - s$R[rows]
      known_count[g] <- length(rows)
      known_mean[g] <- if (length(rows)) mean(e) else 0
      known_ss[g] <- sum((e - known_mean[g])^2)
    }
  }
  full_design <- s$full_design
  # A declared noisy reference leaves no row with exactly known target, so the
  # regression variance needs a proper prior to rule out a spike at zero.
  sigma_nu <- if (exact) 0 else 1
  sigma_scale <- if (exact) 0 else prior$tau0^2
  full_design$df <- full_design$df + sigma_nu
  log_prob_sum <- function(a, b) pmax(a, b) + log1p(exp(-abs(a - b)))
  # Updates corruption rates and error parameters with recording status and target
  # integrated out. The later joint draw of status and target completes a valid
  # two-block Gibbs scan, and avoids the slow mixing of conditioning on them here.
  update_recording <- function(recording, terms) {
    for (j in seq_len(n_corruption)) {
      rows <- rows_s[[j]]
      base <- terms$base[rows]
      wrong <- terms$add[rows] + stats::dnorm(w_latent[rows],
        terms$mean1[rows] + recording$mu[g_latent[rows]],
        sqrt(terms$var1[rows] + recording$tau2[g_latent[rows]]), log = TRUE)
      a_pi <- prior$pi_a + known_k[j]
      b_pi <- prior$pi_b + known_n[j] - known_k[j]
      recording$logit_pi[j] <- catmech_slice(recording$logit_pi[j], function(logit) {
        log_p <- stats::plogis(logit, log.p = TRUE)
        log_q <- stats::plogis(-logit, log.p = TRUE)
        a_pi * log_p + b_pi * log_q + sum(log_prob_sum(log_q + base, log_p + wrong))
      }, 1)
    }
    log_p_row <- stats::plogis(recording$logit_pi[s_latent], log.p = TRUE)
    log_q_row <- stats::plogis(-recording$logit_pi[s_latent], log.p = TRUE)
    for (g in seq_len(n_errors)) {
      rows <- rows_g[[g]]
      left <- log_q_row[rows] + terms$base[rows]
      right <- log_p_row[rows] + terms$add[rows]
      w_rows <- w_latent[rows]
      mean_rows <- terms$mean1[rows]
      var_rows <- terms$var1[rows]
      count <- known_count[g]
      density <- function(mu, tau2) {
        spread <- prior$nu0 * prior$tau0^2 + prior$kappa0 * (mu - prior$mu0)^2 +
          known_ss[g] + count * (known_mean[g] - mu)^2
        -(prior$nu0 + 1 + count) / 2 * log(tau2) - spread / (2 * tau2) +
          sum(log_prob_sum(left, right + stats::dnorm(w_rows, mean_rows + mu, sqrt(var_rows + tau2), log = TRUE)))
      }
      tau2 <- recording$tau2[g]
      recording$mu[g] <- catmech_slice(recording$mu[g], function(z) density(z, tau2), sqrt(tau2))
      mu <- recording$mu[g]
      tau2 <- exp(2 * catmech_slice(log(tau2) / 2, function(z) density(mu, exp(2 * z)), 1))
      # With few observed errors the prior ties the mean's scale to tau, a funnel that
      # coordinate updates cross slowly. This moves both along that ridge; the added z
      # is the Jacobian of the scale transformation (Liu and Sabatti, 2000).
      if (count < 5L) {
        offset <- mu - prior$mu0
        z <- catmech_slice(0, function(z) density(prior$mu0 + offset * exp(z), tau2 * exp(2 * z)) + z, 1)
        mu <- prior$mu0 + offset * exp(z)
        tau2 <- tau2 * exp(2 * z)
      }
      recording$mu[g] <- mu
      recording$tau2[g] <- tau2
    }
    if (any(!is.finite(c(recording$logit_pi, recording$mu, recording$tau2))) || any(recording$tau2 <= 0))
      stop("Nonfinite recording-model parameter draw.", call. = FALSE)
    recording
  }
  draw_latent <- function(fitted, sigma2, recording, terms) {
    target <- target_start
    probability <- numeric(n)
    sigma <- sqrt(sigma2)
    if (length(free)) target[free] <- stats::rnorm(length(free), fitted[free] + s$delta, sigma)
    if (length(reference_only)) {
      variance <- 1 / (1 / sigma2 + 1 / reference_variance)
      location <- variance * (fitted[reference_only] / sigma2 + s$R[reference_only] / reference_variance)
      target[reference_only] <- stats::rnorm(length(reference_only), location, sqrt(variance))
    }
    if (length(latent)) {
      mu <- recording$mu[g_latent]
      tau2 <- recording$tau2[g_latent]
      p <- stats::plogis(recording$logit_pi[s_latent] + terms$add +
        stats::dnorm(w_latent, terms$mean1 + mu, sqrt(terms$var1 + tau2), log = TRUE) - terms$base)
      probability[latent] <- p
      wrong <- stats::runif(length(latent)) < p
      target[latent] <- w_latent
      if (any(wrong)) {
        centre <- fitted[latent] + delta_latent
        variance <- 1 / (1 / sigma2 + 1 / tau2[wrong] + extra_precision[wrong])
        location <- variance * (centre[wrong] / sigma2 + (w_latent[wrong] - mu[wrong]) / tau2[wrong] +
                                  extra_information[wrong])
        target[latent[wrong]] <- stats::rnorm(sum(wrong), location, sqrt(variance))
      }
    }
    list(target = target, probability = probability)
  }
  start_recording <- function() {
    # Overdispersed start from the calibration rows alone.
    tau2 <- mu <- numeric(n_errors)
    for (g in seq_len(n_errors)) {
      count <- known_count[g]
      kappa <- prior$kappa0 + count
      location <- (prior$kappa0 * prior$mu0 + count * known_mean[g]) / kappa
      spread <- prior$nu0 * prior$tau0^2 + known_ss[g] + prior$kappa0 * count / kappa * (known_mean[g] - prior$mu0)^2
      tau2[g] <- spread / stats::rchisq(1L, prior$nu0 + count)
      mu[g] <- stats::rnorm(1L, location, sqrt(tau2[g] / kappa))
    }
    pi <- stats::rbeta(n_corruption, prior$pi_a + known_k, prior$pi_b + known_n - known_k)
    list(logit_pi = stats::qlogis(pmin(pmax(pi, 1e-9), 1 - 1e-9)), mu = mu, tau2 = tau2)
  }
  trace_names <- c(paste0("beta:", colnames(s$X)), "sigma", paste0("pi:", s$corruption$levels),
                   paste0("mu:", s$errors$levels), paste0("tau:", s$errors$levels), "latent_mean")
  traces <- array(NA_real_, dim = c(iterations, chains, length(trace_names)),
    dimnames = list(iteration = as.character(seq_len(iterations)),
      chain = as.character(seq_len(chains)), parameter = trace_names))
  allocation <- rep(m %/% chains, chains)
  if (m %% chains) allocation[seq_len(m %% chains)] <- allocation[seq_len(m %% chains)] + 1L
  completed <- matrix(NA_real_, n, m)
  saved <- 0L
  probability_total <- target_total <- numeric(n)
  for (chain in seq_len(chains)) {
    start <- catmech_model_parameter_draw(s$reference_design, s$reference_fit)
    recording <- start_recording()
    beta <- start$beta
    sigma2 <- start$sigma^2
    keep <- if (allocation[chain] == 0L) integer() else if (allocation[chain] == 1L) {
      as.integer(ceiling(iterations / 2))
    } else as.integer(round(seq(1, iterations, length.out = allocation[chain])))
    for (step in seq_len(burnin + iterations)) {
      fitted <- as.vector(s$X %*% beta)
      terms <- catmech_contamination_row_terms(w_latent, if (exact) numeric() else r_latent,
                                               fitted[latent] + delta_latent, sigma2, reference_variance)
      recording <- update_recording(recording, terms)
      state <- draw_latent(fitted, sigma2, recording, terms)
      adjusted <- state$target - delta_all
      coefficients <- numeric(s$full_design$p)
      coefficients[pivot] <- as.vector(projector %*% adjusted)
      sse <- sum((adjusted - as.vector(s$X %*% coefficients))^2)
      if (!is.finite(sse) || sse <= .Machine$double.eps^2 * max(1, sum(adjusted^2)) * n)
        stop("Trusted Gaussian regression has zero residual variance; no posterior draw is supported.", call. = FALSE)
      draw <- catmech_model_parameter_draw(full_design, list(beta = coefficients, sse = sse + sigma_nu * sigma_scale))
      beta <- draw$beta
      sigma2 <- draw$sigma^2
      if (step > burnin) {
        post <- step - burnin
        traces[post, chain, ] <- c(beta, draw$sigma, stats::plogis(recording$logit_pi), recording$mu,
                                   sqrt(recording$tau2), mean(state$target[imputed]))
        probability_total <- probability_total + state$probability
        target_total <- target_total + state$target
        if (post %in% keep) {
          saved <- saved + 1L
          completed[, saved] <- state$target
        }
      }
    }
  }
  if (saved != m || any(!is.finite(completed)))
    stop("Internal error: incomplete Gibbs output.", call. = FALSE)
  retained <- iterations * chains
  probability <- probability_total / retained
  probability[!s$recorded] <- NA_real_
  if (exact) probability[known_scored] <- as.numeric(!s$agree[known_scored])
  list(draws = completed, traces = traces, probability = probability,
       posterior_target = target_total / retained)
}

catmech_contamination_parameter_summary <- function(traces) {
  parameter <- dimnames(traces)[[3L]]
  do.call(rbind, lapply(seq_along(parameter), function(j) {
    z <- as.vector(traces[, , j])
    data.frame(parameter = parameter[j], mean = mean(z), sd = stats::sd(z),
      lower = unname(stats::quantile(z, .025)), upper = unname(stats::quantile(z, .975)))
  }))
}

#' Summarize what reference measurements reveal about recording errors
#'
#' A descriptive first step for any data with both recorded and reference values.
#' It fits no model, so it can be read before a modeling assumption is chosen.
#'
#' @param data A data frame.
#' @param value Name of the numeric recorded measurement.
#' @param reference Name of the numeric reference measurement.
#' @param by Optional categorical columns to compare, such as flag or site.
#' @param agreement_tol Recorded values within this distance of the reference
#'   count as correct. See [catmech_contamination()].
#' @return A list with `status`, `overall`, `by`, `heterogeneity`, `selection` and
#'   `notes`. Rates use Jeffreys intervals. `heterogeneity` holds exact tests of
#'   whether the disagreement rate differs across each `by` column; it is a screen
#'   for unmodeled structure, not a confirmatory test. `selection` compares recorded
#'   values on referenced and unreferenced rows; a large standardized difference
#'   shows the reference sample is selective, and agreement cannot show the reverse.
#' @export
catmech_recording_audit <- function(data, value, reference, by = NULL, agreement_tol = 0) {
  if (!is.data.frame(data)) stop("data must be a data frame.", call. = FALSE)
  if (!is.character(value) || length(value) != 1L || is.na(value) || !value %in% names(data) ||
      !is.numeric(data[[value]]) || any(!is.na(data[[value]]) & !is.finite(data[[value]])))
    stop("value must name a numeric recorded column with finite values or NA.", call. = FALSE)
  catmech_contamination_reference_column(data, value, reference)
  if (!is.numeric(agreement_tol) || length(agreement_tol) != 1L || !is.finite(agreement_tol) || agreement_tol < 0)
    stop("agreement_tol must be one nonnegative finite number.", call. = FALSE)
  if (!is.null(by) && (!is.character(by) || anyNA(by) || anyDuplicated(by) || any(!by %in% names(data))))
    stop("by must name unique existing columns.", call. = FALSE)
  W <- data[[value]]
  R <- data[[reference]]
  known <- !is.na(R)
  recorded <- !is.na(W)
  calibration <- known & recorded
  agree <- catmech_contamination_agreement(W, R, agreement_tol)
  disagree <- calibration & !agree
  summarize <- function(rows) {
    cal <- sum(calibration & rows)
    k <- sum(disagree & rows)
    e <- (W - R)[disagree & rows]
    data.frame(rows = sum(rows), reference_rows = sum(known & rows),
      reference_rate = sum(known & rows) / max(1L, sum(rows)), calibration_rows = cal,
      disagreements = k, disagreement_rate = if (cal) k / cal else NA_real_,
      lower = if (cal) stats::qbeta(.025, k + .5, cal - k + .5) else NA_real_,
      upper = if (cal) stats::qbeta(.975, k + .5, cal - k + .5) else NA_real_,
      error_mean = if (k) mean(e) else NA_real_,
      error_sd = if (k > 1L) stats::sd(e) else NA_real_)
  }
  overall <- summarize(rep(TRUE, nrow(data)))
  table_by <- heterogeneity <- NULL
  if (length(by)) {
    table_by <- do.call(rbind, lapply(by, function(column) {
      labels <- as.character(data[[column]])
      do.call(rbind, lapply(sort(unique(labels[!is.na(labels)])), function(level) {
        cbind(data.frame(column = column, level = level, stringsAsFactors = FALSE),
              summarize(!is.na(labels) & labels == level))
      }))
    }))
    heterogeneity <- do.call(rbind, lapply(by, function(column) {
      labels <- as.character(data[[column]])[calibration]
      use <- !is.na(labels)
      counts <- table(labels[use], disagree[calibration][use])
      p <- if (nrow(counts) < 2L || ncol(counts) < 2L) NA_real_ else
        catmech_contamination_exact_test(counts)
      data.frame(column = column, levels = nrow(counts), p_value = p, stringsAsFactors = FALSE)
    }))
  }
  smd <- NA_real_
  if (sum(recorded & known) > 1L && sum(recorded & !known) > 1L) {
    a <- W[recorded & known]
    b <- W[recorded & !known]
    spread <- sqrt((stats::var(a) + stats::var(b)) / 2)
    if (is.finite(spread) && spread > 0) smd <- (mean(a) - mean(b)) / spread
  }
  status <- if (!overall$calibration_rows) "no_calibration_rows" else
    if (!overall$disagreements) "no_errors_observed" else
      if (overall$disagreement_rate > .9) "mostly_disagreeing" else
        if (overall$disagreements < 5L) "few_errors_observed" else "errors_observed"
  notes <- c(
    if (status == "no_calibration_rows") "No row has both a recorded and a reference value; a recording model cannot be calibrated.",
    if (status == "no_errors_observed") "No disagreement was observed, so the error size distribution would come from the prior alone.",
    if (status == "few_errors_observed") "Fewer than five disagreements were observed, so the error distribution is weakly determined.",
    if (status == "mostly_disagreeing") "More than 90% of calibration rows disagree. Check units and tolerance, and whether the reference itself has measurement error.",
    if (!is.null(heterogeneity) && any(heterogeneity$p_value < .05, na.rm = TRUE))
      paste0("Disagreement rates differ across: ", paste(heterogeneity$column[which(heterogeneity$p_value < .05)], collapse = ", "),
             ". Consider declaring these as strata."),
    if (is.finite(smd) && abs(smd) > .2)
      "Recorded values differ between referenced and unreferenced rows, so reference selection is not simple random sampling.")
  list(status = status, overall = overall, by = table_by, heterogeneity = heterogeneity,
       selection = data.frame(measure = "recorded_value_standardized_difference", value = smd,
                              stringsAsFactors = FALSE),
       notes = as.character(notes))
}

catmech_contamination_exact_test <- function(counts) {
  exact <- tryCatch(stats::fisher.test(counts, workspace = 2e7)$p.value, error = function(e) NULL)
  if (!is.null(exact)) return(exact)
  restore <- preserve_rng()
  on.exit(restore(), add = TRUE)
  set.seed(1L)
  stats::fisher.test(counts, simulate.p.value = TRUE, B = 4000L)$p.value
}

#' Score an imputation model on reference rows it did not see
#'
#' Hides folds of reference values, refits the same workflow and compares its
#' predictions with the hidden measurements. This is the one check of imputation
#' calibration available on real data. It uses no simulation truth.
#'
#' @param fit A `catmech_result` from [catmech_contamination()], or from
#'   [catmech_reference()] with `method = "normal"`, with or without `use_recorded`.
#' @param folds Number of folds, at least two.
#' @param m Completions per refit, at least 20, used for predictive intervals.
#' @param level Nominal central predictive interval level.
#' @param seed Nonnegative integer seed, with caller RNG restored.
#' @return A list with `status`, `summary` (one row of metrics), `rows` (one per
#'   held-out reference), `flags`, `failures`, `advice` and `interpretation`. `advice` is
#'   nonempty when predictions are no better than the recorded value or the stratum error
#'   rate. On two real datasets that happened while `status` stayed clean, and simpler
#'   methods did as well or better. It is a reading from that study, not a tested rule.
#' @details
#' Folds are stratified by whether the reference agreed with the recorded value.
#' Reported measures include predictive interval coverage, randomized probability
#' integral transforms, root mean squared error of the model and of the recorded value, and
#' for the contamination model the Brier score against the stratum error rate as
#' well as expected against observed errors. `status` is `miscalibration_detected`
#' when coverage is significantly below nominal, the expected and observed error
#' counts differ at the 1% level, or a transform test rejects uniformity at the 1%
#' level. A clean result is not proof of validity: held-out rows share the reference
#' selection process, the check has limited power with few errors, and it does not
#' address nonignorable selection or a reference that is itself wrong.
#' @export
catmech_holdout_check <- function(fit, folds = 5L, m = 40L, level = 0.9, seed = 1L) {
  if (!inherits(fit, "catmech_result")) stop("Expected a catmech_result.", call. = FALSE)
  engine <- fit$manifest$engine
  if (!engine %in% c("reference_contamination", "reference_normal", "reference_normal_recorded"))
    stop("Hold-out checks support the contamination and normal reference workflows.", call. = FALSE)
  check_count(folds, "folds", 2L)
  check_count(m, "m", 20L)
  check_count(seed, "seed", 0L)
  if (!is.numeric(level) || length(level) != 1L || !is.finite(level) || level <= 0 || level >= 1)
    stop("level must be strictly between zero and one.", call. = FALSE)
  spec <- fit$spec
  data <- fit$original
  W <- data[[spec$value]]
  R <- data[[spec$reference]]
  contamination <- identical(engine, "reference_contamination")
  tol <- if (contamination) spec$agreement_tol else 0
  noise <- if (contamination) spec$reference_sd else 0
  rows <- which(!is.na(R))
  if (length(rows) < 2L * folds) stop("Too few reference rows for this many folds.", call. = FALSE)
  agree <- catmech_contamination_agreement(W, R, tol)
  group <- ifelse(is.na(W[rows]), "no_recorded_value", ifelse(agree[rows], "agrees", "disagrees"))
  restore <- preserve_rng()
  on.exit(restore(), add = TRUE)
  set.seed(as.integer(seed))
  assignment <- integer(length(rows))
  for (g in unique(group)) {
    members <- which(group == g)
    assignment[members[sample.int(length(members))]] <- sample(rep_len(sample.int(folds), length(members)))
  }
  refit <- function(train, fold_seed) {
    if (contamination) {
      catmech_contamination(train, spec$value, spec$reference, spec$predictors,
        corruption_strata = spec$corruption_strata, error_strata = spec$error_strata,
        agreement_tol = spec$agreement_tol, reference_sd = spec$reference_sd,
        selection_delta = spec$selection_delta, prior = if (is.null(spec$prior)) list() else spec$prior,
        m = m, seed = fold_seed,
        burnin = spec$burnin, iterations = spec$iterations, max_iterations = spec$max_iterations,
        chains = spec$chains)
    } else {
      catmech_reference(train, spec$value, spec$reference, spec$predictors, method = "normal",
                        m = m, seed = fold_seed, use_recorded = identical(engine, "reference_normal_recorded"),
                        cluster = spec[["cluster"]])
    }
  }
  results <- list()
  failures <- character()
  alpha <- (1 - level) / 2
  for (fold in seq_len(folds)) {
    hidden <- rows[assignment == fold]
    train <- data
    train[[spec$reference]][hidden] <- NA_real_
    fitted <- tryCatch(suppressWarnings(refit(train, as.integer(seed) + fold)), error = identity)
    if (inherits(fitted, "error")) {
      failures <- c(failures, paste0("fold ", fold, ": ", conditionMessage(fitted)))
      next
    }
    draws <- vapply(fitted$completed, function(z) z[[spec$value]][hidden], numeric(length(hidden)))
    draws <- matrix(draws, nrow = length(hidden))
    if (noise > 0) draws <- draws + matrix(stats::rnorm(length(draws), 0, noise), nrow = nrow(draws))
    # An agreeing record counts as a tie with its recorded value, so a
    # rounding tolerance cannot turn a correct record into a miss.
    truth <- ifelse(!is.na(W[hidden]) & agree[hidden] & noise == 0, W[hidden], R[hidden])
    below <- rowSums(draws < truth)
    ties <- rowSums(draws == truth)
    u <- stats::runif(length(hidden))
    corruption <- if (contamination) fitted$contamination$corruption$posterior_corruption[hidden] else rep(NA_real_, length(hidden))
    rate <- rep(NA_real_, length(hidden))
    if (contamination) {
      strata <- catmech_contamination_strata(train, spec$corruption_strata, !is.na(W), "corruption_strata")
      calibration <- !is.na(train[[spec$reference]]) & !is.na(W)
      k <- tabulate(strata$index[calibration & !agree], length(strata$levels))
      total <- tabulate(strata$index[calibration], length(strata$levels))
      rate <- ((k + .5) / (total + 1))[strata$index[hidden]]
    }
    results[[length(results) + 1L]] <- data.frame(fold = fold, row = hidden, reference = R[hidden],
      recorded = W[hidden], disagrees = ifelse(is.na(W[hidden]), NA, !agree[hidden]),
      predictive_mean = rowMeans(draws),
      lower = apply(draws, 1L, stats::quantile, probs = alpha, names = FALSE),
      upper = apply(draws, 1L, stats::quantile, probs = 1 - alpha, names = FALSE),
      covered = NA, pit = (below + u * (ties + 1)) / (ncol(draws) + 1),
      posterior_corruption = corruption, stratum_rate = rate)
    last <- length(results)
    results[[last]]$covered <- results[[last]]$lower <= truth & truth <= results[[last]]$upper
  }
  scored <- if (length(results)) do.call(rbind, results) else NULL
  flags <- character()
  if (is.null(scored)) {
    return(list(status = "unresolved_failed_folds", summary = NULL, rows = NULL,
                flags = flags, failures = failures, advice = character(),
                interpretation = "No fold produced a fit."))
  }
  n <- nrow(scored)
  covered <- sum(scored$covered)
  z <- stats::qnorm(.975)
  phat <- covered / n
  wilson_upper <- (phat + z^2 / (2 * n) + z * sqrt(phat * (1 - phat) / n + z^2 / (4 * n^2))) / (1 + z^2 / n)
  ks <- suppressWarnings(stats::ks.test(scored$pit, "punif")$p.value)
  has_w <- !is.na(scored$recorded)
  errors <- sum(scored$disagrees[has_w])
  expected <- brier <- brier_rate <- calibration_z <- NA_real_
  if (contamination && any(has_w)) {
    p <- scored$posterior_corruption[has_w]
    y <- as.numeric(scored$disagrees[has_w])
    expected <- sum(p)
    brier <- mean((p - y)^2)
    brier_rate <- mean((scored$stratum_rate[has_w] - y)^2)
    spread <- sum(p * (1 - p))
    if (spread > 0) calibration_z <- (errors - expected) / sqrt(spread)
  }
  if (wilson_upper < level) flags <- c(flags, "predictive_undercoverage")
  if (is.finite(ks) && ks < .01) flags <- c(flags, "predictive_distribution_rejected")
  if (is.finite(calibration_z) && abs(calibration_z) > stats::qnorm(.995))
    flags <- c(flags, "corruption_probability_miscalibrated")
  summary <- data.frame(engine = engine, holdout_rows = n, folds_completed = length(results),
    with_recorded = sum(has_w), errors_observed = errors, errors_expected = expected,
    calibration_z = calibration_z, brier_model = brier, brier_stratum_rate = brier_rate,
    nominal_level = level, predictive_coverage = phat, ks_p_value = ks,
    rmse_model = sqrt(mean((scored$predictive_mean - scored$reference)^2)),
    rmse_recorded = if (any(has_w)) sqrt(mean((scored$recorded[has_w] - scored$reference[has_w])^2)) else NA_real_,
    mean_interval_width = mean(scored$upper - scored$lower), stringsAsFactors = FALSE)
  status <- if (length(failures)) "unresolved_failed_folds" else
    if (length(flags)) "miscalibration_detected" else "no_miscalibration_detected"
  advice <- c(
    if (is.finite(summary$rmse_recorded) && summary$rmse_model >= summary$rmse_recorded)
      "Predictions are no closer to held-out references than the recorded values are.",
    if (is.finite(brier) && brier >= brier_rate)
      "Error probabilities are no better than the stratum error rate, so do not read them as per-record labels.")
  if (length(advice)) advice <- c(advice, "Consider catmech_reference(use_recorded = TRUE) or generalized raking.")
  list(status = status, summary = summary, rows = scored, flags = flags, failures = failures, advice = advice,
       interpretation = paste("Scores imputation against held-out references from the same selection process.",
         "A clean result is not proof of validity; see the function documentation."))
}
