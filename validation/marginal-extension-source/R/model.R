#' Gaussian completion with an explicit category observation model
#'
#' Replaces all flagged measurements and missing values using a conditional
#' Gaussian model. Trusted numeric measurements are kept exactly. The default
#' ignores categories; optional exact or externally specified noisy categories
#' enter through a likelihood, not by forcing trusted values to change.
#'
#' @param data A data frame. Original columns, values and category labels are
#'   preserved separately and in every completion except for the target cells.
#' @param value Name of the numeric target column.
#' @param predictors Explicit names of complete covariates, including the
#'   downstream outcome when appropriate. Use `character(0)` for an intercept.
#' @param flag Optional quality-flag column. All `flag_yes` measurements are
#'   treated as unavailable, including false-positive flags.
#' @param flag_yes,flag_no Labels for suspect and non-suspect measurements.
#'   Missing flags are treated as non-suspect, an explicit trust assumption.
#' @param category Optional category column. Ignored without validation when
#'   `category_mode = "ignore"`. Missing categories contribute no likelihood.
#' @param thresholds Numeric category endpoints. For informative categories,
#'   the first must be `-Inf` and the last `Inf`.
#' @param levels Ordered category labels with length one less than thresholds.
#' @param category_mode `"ignore"` (default), `"exact"`, or `"noisy"`.
#' @param confusion For noisy categories, a numeric matrix with named rows for
#'   true categories and named columns for observed categories, both matching
#'   `levels`. Entries are nonnegative and each row sums to one. These fixed
#'   probabilities must come from external validation or declared sensitivity
#'   assumptions; this function does not estimate them or their uncertainty.
#' @param m Number of completions, at least two.
#' @param seed Nonnegative integer seed. The caller's RNG state is restored.
#' @param burnin Discarded Gibbs iterations per chain.
#' @param iterations Retained diagnostic iterations per chain, at least 100
#'   when informative categories require Gibbs sampling.
#' @param chains Number of independent Gibbs chains, at least four when used.
#'
#' @details
#' The model is `Y | X, beta, sigma ~ Normal(X beta, sigma^2)`, with the reference
#' prior `p(beta, sigma^2) proportional to 1 / sigma^2`. The trusted design must
#' have full column rank, positive residual degrees of freedom and positive
#' residual variance. Globally aliased predictors must be removed explicitly.
#' Flag removal and missingness must be ignorable given the selected predictors.
#' This is not a model of the observed corrupted measurement or an estimate of
#' its probability of being wrong.
#'
#' Ignored categories use independent posterior predictive draws from trusted
#' data. Informative categories use Gibbs data augmentation: sample latent bins
#' with probabilities proportional to category likelihood times normal interval
#' mass; draw values from those truncated normals; and update regression
#' parameters using completed data. Equal rows in a confusion matrix contain no
#' information about the target and reduce to the independent-draw engine.
#'
#' Gibbs diagnostics use `posterior` rank-normalized split/folded R-hat and bulk
#' and tail effective sample size. An R-hat above 1.05, bulk or tail ESS below
#' `100 * chains`, or a nonfinite
#' diagnostic produces a `catmech_convergence_error`. Its `model_diagnostics`
#' field retains all traces and the diagnostic summary. These checks do not
#' prove convergence or calibration. Completions are spaced across retained
#' iterations and balanced across chains; they are not claimed independent.
#'
#' Exact-category conflicts with trusted values are errors. With noisy categories,
#' a trusted conflict is allowed only when its supplied likelihood is positive.
#' Recorded categories are never overwritten to make a completion agree with them.
#' Rubin pooling remains conditional on model adequacy and compatibility with
#' the intended analysis, especially outside a linear-Gaussian analysis setting.
#'
#' @return A `catmech_result` with `original`, `completed`, `mids`, `audit`,
#'   `diagnostics`, `spec`, `manifest`, and `model_diagnostics`. The `mids` object
#'   contains the target and selected predictors, plus a constant auxiliary
#'   column for an intercept-only model. Its completions are
#'   externally generated; it must not be continued using `mice.mids`.
#' @export
#' @examples
#' d <- simulate_catmech_lab(100, seed = 8)
#' fit <- catmech_model(d, "lab_value", c("age", "inflammation", "outcome"),
#'   flag = "outlier_flag", m = 3, seed = 12)
#' fit$manifest$engine
catmech_model <- function(data, value, predictors, flag = NULL,
                           flag_yes = "Yes", flag_no = "No",
                           category = NULL, thresholds = NULL, levels = NULL,
                           category_mode = c("ignore", "exact", "noisy"),
                           confusion = NULL, m = 20L, seed = 1L,
                           burnin = 200L, iterations = 400L, chains = 4L) {
  category_mode <- match.arg(category_mode)
  if (category_mode == "ignore" && !is.null(confusion))
    stop("confusion is unused when category_mode='ignore'; remove it or choose category_mode='noisy'.", call. = FALSE)
  scalar_name <- function(z) is.character(z) && length(z) == 1L &&
    !is.na(z) && nzchar(z)
  if (!is.data.frame(data) || nrow(data) < 2L || anyNA(names(data)) ||
      anyDuplicated(names(data)) || any(!nzchar(names(data)))) {
    stop("data must have at least two rows and unique nonempty column names.", call. = FALSE)
  }
  if (any(c(".imp", ".id") %in% names(data)))
    stop(".imp and .id are reserved column names.", call. = FALSE)
  if (any(!vapply(data, function(z) is.atomic(z) && is.null(dim(z)), logical(1))))
    stop("List and matrix columns are not supported.", call. = FALSE)
  if (!scalar_name(value) || (!is.null(flag) && !scalar_name(flag)) ||
      (!is.null(category) && !scalar_name(category)))
    stop("value, flag and category must be single nonempty names when supplied.", call. = FALSE)
  if (!scalar_name(flag_yes) || !scalar_name(flag_no) || flag_yes == flag_no)
    stop("flag_yes and flag_no must be distinct nonempty labels.", call. = FALSE)
  if (!is.character(predictors) || anyNA(predictors) || anyDuplicated(predictors) ||
      any(!nzchar(predictors)))
    stop("predictors must be an explicit character vector of unique column names.", call. = FALSE)
  special <- c(value, flag, category)
  if (anyDuplicated(special) || any(predictors %in% special))
    stop("Target, flag, category and predictor roles must use distinct columns.", call. = FALSE)
  needed <- c(value, flag, predictors, if (category_mode != "ignore") category)
  if (length(setdiff(needed, names(data))))
    stop("Missing data columns: ", paste(setdiff(needed, names(data)), collapse = ", "), call. = FALSE)
  check_count(m, "m", 2L)
  check_count(seed, "seed", 0L)
  check_count(burnin, "burnin", 0L)
  check_count(iterations, "iterations", 1L)
  check_count(chains, "chains", 1L)
  n <- nrow(data)
  y <- data[[value]]
  if (!is.numeric(y) || any(!is.na(y) & !is.finite(y)))
    stop("The target must contain finite numeric values or NA.", call. = FALSE)
  flag_values <- if (is.null(flag)) rep(flag_no, n) else as.character(data[[flag]])
  if (any(!is.na(flag_values) & !flag_values %in% c(flag_no, flag_yes)))
    stop("Unknown flag labels in the data.", call. = FALSE)
  suspect <- catmech_suspect_mask(flag_values, flag_yes)
  masked <- is.na(y) | suspect
  trusted <- !masked
  X <- catmech_model_matrix(data, predictors)
  trusted_design <- catmech_model_design(X[trusted, , drop = FALSE])
  # Checking the trusted posterior before starting makes its propriety requirement
  # explicit; missing/category-only rows cannot substitute for that anchor.
  trusted_fit <- catmech_model_fit(trusted_design, y[trusted])
  category_index <- rep(NA_integer_, n)
  H <- NULL
  if (category_mode != "ignore") {
    if (is.null(category)) stop("An informative category mode requires a category column.", call. = FALSE)
    if (!is.character(levels) || !length(levels) || anyNA(levels) ||
        anyDuplicated(levels) || any(!nzchar(levels)))
      stop("levels must be unique nonempty category labels.", call. = FALSE)
    thresholds <- normalize_thresholds(thresholds)
    make_interval_bounds(thresholds, levels)
    if (thresholds[1L] != -Inf || utils::tail(thresholds, 1L) != Inf)
      stop("The Gaussian category model requires outer endpoints -Inf and Inf.", call. = FALSE)
    raw_category <- as.character(data[[category]])
    if (any(!is.na(raw_category) & !raw_category %in% levels))
      stop("Unknown category labels in the data.", call. = FALSE)
    category_index <- match(raw_category, levels)
    if (category_mode == "exact") {
      if (!is.null(confusion)) stop("Supply confusion only for category_mode='noisy'.", call. = FALSE)
      H <- diag(length(levels))
      dimnames(H) <- list(levels, levels)
    } else H <- catmech_model_confusion(confusion, levels)
    known <- which(trusted & !is.na(category_index))
    if (length(known)) {
      true_bin <- as.integer(catmech_derive_category(y[known], thresholds, levels))
      likelihood <- H[cbind(true_bin, category_index[known])]
      if (any(likelihood <= 0))
        stop("A trusted value/category conflict has zero probability under the supplied category model.", call. = FALSE)
    }
  } else {
    thresholds <- NULL
    levels <- NULL
  }
  informative <- if (is.null(H)) rep(FALSE, n) else vapply(category_index, function(j) {
    if (is.na(j)) return(FALSE)
    any(H[, j] != H[1L, j])
  }, logical(1))
  # An observed category with zero probability for every true bin is impossible,
  # even when its column happens to be constant.
  if (!is.null(H) && any(vapply(category_index, function(j) {
    !is.na(j) && all(H[, j] == 0)
  }, logical(1)))) stop("An observed category has zero probability in every true category.", call. = FALSE)
  use_gibbs <- any(masked & informative)
  if (use_gibbs && (iterations < 100L || chains < 4L))
    stop("Informative-category Gibbs sampling requires iterations >= 100 and chains >= 4.", call. = FALSE)
  if (use_gibbs && m > iterations * chains)
    stop("m cannot exceed the number of retained Gibbs draws.", call. = FALSE)
  restore <- preserve_rng()
  on.exit(restore(), add = TRUE)
  set.seed(as.integer(seed))
  completed_y <- matrix(NA_real_, nrow = n, ncol = m)
  parameter_names <- c(paste0("beta:", colnames(X)), "sigma")
  if (!use_gibbs) {
    direct_traces <- matrix(NA_real_, m, length(parameter_names),
                            dimnames = list(NULL, parameter_names))
    for (j in seq_len(m)) {
      draw <- catmech_model_parameter_draw(trusted_design, trusted_fit)
      out <- y
      out[masked] <- stats::rnorm(sum(masked), as.vector(X[masked, , drop = FALSE] %*% draw$beta), draw$sigma)
      completed_y[, j] <- out
      direct_traces[j, ] <- c(draw$beta, draw$sigma)
    }
    model_diagnostics <- list(engine = "independent_normal", status = "not_applicable",
      convergence_pass = TRUE, summary = data.frame(parameter = parameter_names,
        rhat = NA_real_, ess_bulk = rep(as.numeric(m), length(parameter_names)),
        ess_tail = rep(as.numeric(m), length(parameter_names))),
      traces = direct_traces,
      reason = "Independent posterior-predictive draws; no informative category likelihood on masked rows.")
  } else {
    full_design <- catmech_model_design(X)
    trace_names <- c(parameter_names, "latent_mean")
    traces <- array(NA_real_, dim = c(iterations, chains, length(trace_names)),
      dimnames = list(iteration = as.character(seq_len(iterations)),
        chain = as.character(seq_len(chains)), parameter = trace_names))
    allocation <- rep(m %/% chains, chains)
    if (m %% chains) allocation[seq_len(m %% chains)] <- allocation[seq_len(m %% chains)] + 1L
    saved_index <- 0L
    for (chain in seq_len(chains)) {
      start <- catmech_model_parameter_draw(trusted_design, trusted_fit)
      current <- y
      current[masked] <- catmech_model_latent(
        as.vector(X[masked, , drop = FALSE] %*% start$beta), start$sigma,
        category_index[masked], thresholds, H)
      keep <- if (allocation[chain] == 0L) integer() else if (allocation[chain] == 1L) {
        as.integer(ceiling(iterations / 2))
      } else as.integer(round(seq(1, iterations, length.out = allocation[chain])))
      for (step in seq_len(burnin + iterations)) {
        fit <- catmech_model_fit(full_design, current)
        draw <- catmech_model_parameter_draw(full_design, fit)
        current[masked] <- catmech_model_latent(
          as.vector(X[masked, , drop = FALSE] %*% draw$beta), draw$sigma,
          category_index[masked], thresholds, H)
        if (step > burnin) {
          post <- step - burnin
          traces[post, chain, ] <- c(draw$beta, draw$sigma, mean(current[masked]))
          if (post %in% keep) {
            saved_index <- saved_index + 1L
            completed_y[, saved_index] <- current
          }
        }
      }
    }
    if (saved_index != m || any(!is.finite(completed_y)))
      stop("Internal error: incomplete Gibbs output.", call. = FALSE)
    model_diagnostics <- catmech_model_chain_diagnostics(traces)
    if (!model_diagnostics$convergence_pass) {
      condition <- structure(list(
        message = paste0("Gaussian category model failed convergence checks (R-hat <= 1.05 and bulk/tail ESS >= ",
                         100L * chains, " required). ",
                         "Inspect condition$model_diagnostics and increase iterations/burnin or revise the model."),
        call = NULL, model_diagnostics = model_diagnostics),
        class = c("catmech_convergence_error", "error", "condition"))
      stop(condition)
    }
  }
  if (any(!is.finite(completed_y))) stop("Model generated nonfinite target values.", call. = FALSE)
  completed <- lapply(seq_len(m), function(j) {
    out <- data
    out[[value]] <- completed_y[, j]
    out
  })
  selected <- which(masked)
  audits <- lapply(seq_len(m), function(j) data.frame(
    row = selected, reason = ifelse(is.na(y[selected]), "missing", "flag_masked"),
    reconstructed = rep(TRUE, length(selected)),
    changed = ifelse(is.na(y[selected]), NA, completed_y[selected, j] != y[selected]),
    donor_row = rep(NA_integer_, length(selected)), used_fallback = rep(FALSE, length(selected)),
    eligible_donors = rep(NA_integer_, length(selected)),
    imputation = rep(j, length(selected)), category_derived = rep(FALSE, length(selected)),
    stringsAsFactors = FALSE))
  audit <- do.call(rbind, audits)
  rownames(audit) <- NULL
  diagnostics <- do.call(rbind, lapply(seq_len(m), function(j) {
    known <- which(!is.na(category_index))
    disagreement <- if (is.null(H)) NA_integer_ else sum(as.integer(
      catmech_derive_category(completed_y[known, j], thresholds, levels)) != category_index[known])
    if (category_mode == "exact" && disagreement != 0L)
      stop("Internal error: exact category likelihood was violated.", call. = FALSE)
    data.frame(imputation = j, requested = length(selected), reconstructed = length(selected),
      changed = sum(audits[[j]]$changed, na.rm = TRUE), fallback_draws = 0L,
      category_disagreements = disagreement,
      threshold_violations = if (category_mode == "exact") disagreement else NA_integer_)
  }))
  adapter_names <- unique(c(predictors, value))
  # The mids constructor requires two columns even for a single target.
  # A complete constant auxiliary column is retained only in that adapter.
  base <- data[adapter_names]
  if (ncol(base) == 1L) {
    auxiliary <- ".catmech_intercept"
    while (auxiliary %in% names(base)) auxiliary <- paste0(auxiliary, "_")
    base[[auxiliary]] <- 1
  }
  where <- matrix(FALSE, n, ncol(base), dimnames = list(NULL, names(base)))
  where[, value] <- masked
  empty_method <- stats::setNames(rep("", ncol(base)), names(base))
  empty_predictors <- matrix(0L, ncol(base), ncol(base),
                             dimnames = list(names(base), names(base)))
  # Build a standard, non-iterating analysis adapter without asking mice to refit
  # or screen predictors that were already validated by this model.
  mids <- mice::mice(base, m = m, where = where, method = empty_method,
    predictorMatrix = empty_predictors, maxit = 0L, remove.constant = FALSE,
    remove.collinear = FALSE, printFlag = FALSE)
  target_imputations <- as.data.frame(completed_y[masked, , drop = FALSE])
  names(target_imputations) <- as.character(seq_len(m))
  rownames(target_imputations) <- rownames(base)[masked]
  mids$imp[[value]] <- target_imputations
  mids$method[value] <- "catmech_gaussian_external_draws"
  mids$call <- match.call()
  confusion_list <- if (is.null(H)) NULL else stats::setNames(lapply(seq_len(nrow(H)), function(k) {
    stats::setNames(as.list(as.numeric(H[k, ])), colnames(H))
  }), rownames(H))
  spec <- structure(list(version = "model1", engine = "gaussian_model", value = value,
    flag = flag, flag_yes = flag_yes, flag_no = flag_no, flag_policy = "replace",
    missing_flag_policy = "trust_observed", category = category,
    category_mode = category_mode, thresholds = thresholds, levels = levels,
    confusion = confusion_list, likelihood = "normal_target_given_complete_predictors",
    prior = "p(beta,sigma_squared) proportional to 1/sigma_squared",
    missingness_assumption = "ignorable_selection_given_predictors",
    category_parameters = "fixed_external_or_sensitivity_assumptions"), class = "catmech_model_spec")
  manifest <- list(schema_version = "model1", engine = "gaussian_model",
    assumption_hash = catmech_hash(spec),
    data_hash = digest::digest(serialize(data, NULL, version = 2), algo = "sha256", serialize = FALSE),
    data_hash_format = "R serialization version 2; environment-dependent data fingerprint",
    predictors = predictors, m = as.integer(m), seed = as.integer(seed), rng_kind = RNGkind(),
    burnin = if (use_gibbs) burnin else 0L, iterations = if (use_gibbs) iterations else m,
    chains = if (use_gibbs) chains else 0L,
    package_version = as.character(utils::packageVersion("catmech")),
    r_version = as.character(getRversion()), mice_version = as.character(utils::packageVersion("mice")),
    posterior_version = as.character(utils::packageVersion("posterior")),
    claim_level = "conditional_model_based", identification = "assumptions_not_individual_truth",
    assumptions = unclass(spec))
  structure(list(original = data, completed = completed, mids = mids, audit = audit,
    diagnostics = diagnostics, spec = spec, manifest = manifest,
    model_diagnostics = model_diagnostics), class = "catmech_result")
}

catmech_model_matrix <- function(data, predictors) {
  if (!length(predictors)) return(matrix(1, nrow(data), 1L,
                                       dimnames = list(NULL, "(Intercept)")))
  px <- data[predictors]
  if (anyNA(px) || any(!vapply(px, function(z) {
    is.numeric(z) || is.factor(z) || is.character(z) || is.logical(z)
  }, logical(1)))) stop("Predictors must be complete numeric, logical, character or factor columns.", call. = FALSE)
  for (nm in names(px)) {
    if (is.numeric(px[[nm]]) && any(!is.finite(px[[nm]])))
      stop("Predictors must be finite.", call. = FALSE)
    if (is.character(px[[nm]]) || is.factor(px[[nm]])) px[[nm]] <- droplevels(factor(px[[nm]]))
    if (length(unique(px[[nm]])) < 2L)
      stop("Remove constant predictor: ", nm, call. = FALSE)
  }
  stats::model.matrix(~ ., data = px)
}

catmech_model_confusion <- function(H, levels) {
  k <- length(levels)
  if (!is.matrix(H) || !is.numeric(H) || !identical(dim(H), c(k, k)) ||
      is.null(rownames(H)) || is.null(colnames(H)) ||
      anyDuplicated(rownames(H)) || anyDuplicated(colnames(H)) ||
      !setequal(rownames(H), levels) || !setequal(colnames(H), levels) ||
      any(!is.finite(H)) || any(H < 0) || any(abs(rowSums(H) - 1) > 1e-10))
    stop("confusion must be a nonnegative named square matrix with true rows and observed columns matching levels, each row summing to one.", call. = FALSE)
  H[levels, levels, drop = FALSE]
}

catmech_model_design <- function(X) {
  if (nrow(X) <= ncol(X))
    stop("The trusted regression needs more observations than coefficients.", call. = FALSE)
  q <- qr(X, tol = 1e-7)
  if (q$rank != ncol(X))
    stop("The trusted regression design must have full column rank; remove aliased or unsupported predictors.", call. = FALSE)
  list(X = X, qr = q, R = qr.R(q), df = nrow(X) - ncol(X), p = ncol(X))
}

catmech_model_fit <- function(design, y) {
  beta <- as.vector(qr.coef(design$qr, y))
  sse <- sum((y - as.vector(design$X %*% beta))^2)
  floor <- .Machine$double.eps^2 * max(1, sum(y^2)) * length(y)
  if (!is.finite(sse) || sse <= floor)
    stop("Trusted Gaussian regression has zero residual variance; no posterior draw is supported.", call. = FALSE)
  list(beta = beta, sse = sse)
}

catmech_model_parameter_draw <- function(design, fit) {
  sigma <- sqrt(fit$sse / stats::rchisq(1L, design$df))
  noise <- numeric(design$p)
  noise[design$qr$pivot] <- backsolve(design$R, stats::rnorm(design$p))
  beta <- fit$beta + sigma * noise
  if (!is.finite(sigma) || any(!is.finite(beta)))
    stop("Nonfinite Gaussian posterior parameter draw.", call. = FALSE)
  list(beta = beta, sigma = sigma)
}

catmech_model_log_interval <- function(lower, upper) {
  # Work in the nearer tail to prevent subtracting two probabilities near one.
  right <- lower > 0
  large <- stats::pnorm(upper, log.p = TRUE)
  small <- stats::pnorm(lower, log.p = TRUE)
  large[right] <- stats::pnorm(lower[right], lower.tail = FALSE, log.p = TRUE)
  small[right] <- stats::pnorm(upper[right], lower.tail = FALSE, log.p = TRUE)
  large + log(-expm1(small - large))
}

catmech_model_latent <- function(mu, sigma, category_index, thresholds, H) {
  n <- length(mu)
  if (!n) return(numeric())
  out <- numeric(n)
  free <- is.na(category_index)
  if (any(free)) out[free] <- stats::rnorm(sum(free), mu[free], sigma)
  rows <- which(!free)
  if (!length(rows)) return(out)
  k <- nrow(H)
  means <- mu[rows]
  log_mass <- matrix(vapply(seq_len(k), function(j) {
    catmech_model_log_interval((thresholds[j] - means) / sigma,
                              (thresholds[j + 1L] - means) / sigma)
  }, numeric(length(rows))), nrow = length(rows), ncol = k)
  log_weight <- log_mass + log(t(H[, category_index[rows], drop = FALSE]))
  maximum <- apply(log_weight, 1L, max)
  if (any(!is.finite(maximum)))
    stop("Category likelihood has no representable Gaussian interval probability.", call. = FALSE)
  weight <- exp(log_weight - maximum)
  probability <- weight / rowSums(weight)
  cumulative <- probability
  if (k > 1L) for (j in 2:k) cumulative[, j] <- cumulative[, j - 1L] + probability[, j]
  u <- stats::runif(length(rows))
  bin <- if (k == 1L) rep(1L, length(rows)) else
    1L + rowSums(u >= cumulative[, -k, drop = FALSE])
  for (j in unique(bin)) {
    recipients <- rows[bin == j]
    for (attempt in seq_len(100L)) {
      out[recipients] <- truncnorm::rtruncnorm(length(recipients),
        a = thresholds[j], b = thresholds[j + 1L], mean = mu[recipients], sd = sigma)
      good <- is.finite(out[recipients]) & out[recipients] > thresholds[j] &
        out[recipients] <= thresholds[j + 1L]
      recipients <- recipients[!good]
      if (!length(recipients)) break
    }
    if (length(recipients)) stop("No representable value in the sampled latent interval.", call. = FALSE)
  }
  out
}

catmech_model_chain_diagnostics <- function(traces) {
  parameter <- dimnames(traces)[[3L]]
  summary <- data.frame(parameter = parameter,
    rhat = vapply(seq_along(parameter), function(j) posterior::rhat(traces[, , j]), numeric(1)),
    ess_bulk = vapply(seq_along(parameter), function(j) posterior::ess_bulk(traces[, , j]), numeric(1)),
    ess_tail = vapply(seq_along(parameter), function(j) posterior::ess_tail(traces[, , j]), numeric(1)))
  minimum_ess <- 100L * dim(traces)[2L]
  pass <- all(is.finite(summary$rhat) & is.finite(summary$ess_bulk) &
                is.finite(summary$ess_tail) & summary$rhat <= 1.05 &
                summary$ess_bulk >= minimum_ess & summary$ess_tail >= minimum_ess)
  list(engine = "gibbs_category_likelihood", status = if (pass) "passed" else "failed",
    convergence_pass = pass, summary = summary, traces = traces,
    thresholds = list(max_rhat = 1.05, min_bulk_ess = minimum_ess,
                      min_tail_ess = minimum_ess))
}
