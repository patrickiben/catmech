# Synthetic one-analyte, one-visit validation fixture. Not exported package API.
# Conceptual transfer from the author's LLM_CDISC_Eval synthetic LB generator,
# frozen-truth protocol, and observable-defect checks. No private data are read.
# This implementation is independent of catmech's categorization/imputation code.
# Source this file to use the functions; Rscript tools/lab-fixture.R runs checks.

lab_reference_category <- function(x, lower = 0.60, upper = 1.30) {
  if (!is.numeric(x) || !is.numeric(lower) || !is.numeric(upper) ||
      length(lower) != 1L || length(upper) != 1L ||
      !is.finite(lower) || !is.finite(upper) || lower >= upper)
    stop("Supply numeric values and finite ordered reference limits.")
  factor(ifelse(is.na(x), NA_character_,
    ifelse(x < lower, "LOW", ifelse(x > upper, "HIGH", "NORMAL"))),
    levels = c("LOW", "NORMAL", "HIGH"))
}

lab_boundary_controls <- function() {
  data.frame(
    case = c("below_lln", "at_lln", "above_lln", "below_uln", "at_uln",
             "above_uln", "rounded_to_lln_from_below", "rounded_to_uln_from_above"),
    true_value = c(0.59, 0.60, 0.61, 1.29, 1.30, 1.31, 0.596, 1.304),
    recorded_2dp = c(0.59, 0.60, 0.61, 1.29, 1.30, 1.31, 0.60, 1.30),
    expected_reference = c("LOW", "NORMAL", "NORMAL", "NORMAL", "NORMAL",
                           "HIGH", "LOW", "HIGH"),
    expected_from_rounded = c("LOW", "NORMAL", "NORMAL", "NORMAL", "NORMAL",
                              "HIGH", "NORMAL", "NORMAL"),
    stringsAsFactors = FALSE
  )
}

lab_fixture_rng_restore <- function() {
  had <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  old <- if (had) get(".Random.seed", envir = .GlobalEnv, inherits = FALSE) else NULL
  kind <- RNGkind()
  function() {
    do.call(RNGkind, as.list(kind))
    if (had) assign(".Random.seed", old, envir = .GlobalEnv)
    else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
      rm(".Random.seed", envir = .GlobalEnv)
  }
}

make_lab_fixture <- function(n = 600L, replicate = 1L,
    phase = c("development", "calibration", "heldout"), seed = 20261031L,
    scenario = c("clean", "flagged_errors", "missed_flags", "site_shift",
                 "rounded_records")) {
  phase <- match.arg(phase)
  scenario <- match.arg(scenario)
  valid_integer <- function(x, minimum) is.numeric(x) && length(x) == 1L &&
    is.finite(x) && x >= minimum && x <= .Machine$integer.max && x == floor(x)
  if (!valid_integer(n, 30L) || !valid_integer(replicate, 1L) ||
      !valid_integer(seed, 0L)) stop("Invalid n, replicate or seed.")
  restore <- lab_fixture_rng_restore()
  on.exit(restore(), add = TRUE)
  RNGkind("L'Ecuyer-CMRG", "Inversion", "Rejection")
  set.seed(seed)
  state <- get(".Random.seed", envir = .GlobalEnv)
  # Independent full streams by phase, then independent substreams by repetition.
  stream_index <- match(phase, c("calibration", "development", "heldout"))
  for (i in seq_len(stream_index)) state <- parallel::nextRNGStream(state)
  if (replicate > 1L) for (i in seq_len(replicate - 1L))
    state <- parallel::nextRNGSubStream(state)
  assign(".Random.seed", state, envir = .GlobalEnv)
  initial_state <- state

  site <- factor(sample(rep(c("001", "002", "003"), length.out = n)),
                 levels = c("001", "002", "003"))
  age <- stats::rnorm(n, 50, 10)
  sex <- factor(sample(c("F", "M"), n, replace = TRUE), levels = c("F", "M"))
  inflammation <- stats::rnorm(n, 2, 1)
  biological_site_shift <- c(-0.025, 0, 0.025)[as.integer(site)]
  mean_lab <- 0.95 + 0.001 * (age - 50) + 0.025 * (sex == "M") +
    0.025 * (inflammation - 2) + biological_site_shift
  lab_true <- stats::rnorm(n, mean_lab, 0.20)
  # Positive laboratory values; the exact conditional outcome slope below is
  # unaffected by this rejection. The target model is otherwise approximately normal.
  while (any(lab_true <= 0)) {
    idx <- which(lab_true <= 0)
    lab_true[idx] <- stats::rnorm(length(idx), mean_lab[idx], 0.20)
  }
  outcome_site_shift <- c(-0.10, 0, 0.10)[as.integer(site)]
  outcome_noise <- stats::rnorm(n, 0, 0.20)
  outcome <- 0.8 * lab_true + 0.002 * (age - 50) + 0.10 * (sex == "M") +
    0.05 * (inflammation - 2) + outcome_site_shift + outcome_noise
  category <- lab_reference_category(lab_true)

  # Draw all latent error/flag uniforms before applying a scenario. Reusing a
  # phase/repetition across scenarios therefore preserves truth and missingness.
  hidden <- stats::runif(n) < 0.20
  error_u <- stats::runif(n)
  flag_u <- stats::runif(n)
  additive_error <- stats::rnorm(n, 0.35, 0.07)
  error_scenario <- scenario %in% c("flagged_errors", "missed_flags", "site_shift")
  random_error <- !hidden & error_u < if (error_scenario) 0.08 else 0
  site_error <- !hidden & site == "003" & scenario == "site_shift"
  value_shift <- ifelse(random_error, additive_error, 0) + ifelse(site_error, 0.20, 0)
  corrupted <- random_error | site_error
  sensitivity <- switch(scenario, flagged_errors = 1, missed_flags = 0.60,
                        site_shift = 0.20, 0)
  false_positive_rate <- if (error_scenario) 0.05 else 0
  flagged <- !hidden & flag_u < ifelse(corrupted, sensitivity, false_positive_rate)
  before_rounding <- lab_true + value_shift
  recorded <- if (scenario == "rounded_records") round(before_rounding, 2) else before_rounding
  rounding_error <- recorded - before_rounding
  recorded[hidden] <- NA_real_
  ids <- sprintf("SYN-CATMECH-LAB-%06d", seq_len(n))
  data <- data.frame(
    STUDYID = "SYN-CATMECH-LAB", USUBJID = ids, SITEID = site,
    VISIT = "BASELINE", VISITNUM = 1L, LBTESTCD = "CREAT",
    LBSTRESN = recorded, LBSTRESU = "mg/dL", LBSTNRLO = 0.60, LBSTNRHI = 1.30,
    LBNRIND = category, QCFLAG = factor(ifelse(flagged, "Yes", "No"),
                                      levels = c("No", "Yes")),
    age = age, sex = sex, inflammation = inflammation, outcome = outcome,
    stringsAsFactors = FALSE
  )
  truth <- data.frame(
    USUBJID = ids, LAB_TRUE = lab_true, CATEGORY_TRUE = category,
    MISSING = hidden, IS_CORRUPTED = corrupted, RANDOM_ERROR = random_error,
    SITE_ERROR = site_error, FLAGGED = flagged, VALUE_SHIFT = value_shift,
    ROUNDING_ERROR = rounding_error, OUTCOME_ERROR = outcome_noise,
    stringsAsFactors = FALSE
  )
  truth_hash <- digest::digest(truth[c("USUBJID", "LAB_TRUE", "CATEGORY_TRUE", "OUTCOME_ERROR")],
                               algo = "sha256", serializeVersion = 2)
  manifest <- list(
    fixture_version = "1", data_type = "synthetic_only", phase = phase,
    replicate = as.integer(replicate), n = as.integer(n), master_seed = as.integer(seed),
    rng_kind = RNGkind(), initial_rng_state = initial_state, scenario = scenario,
    truth_hash = truth_hash, truth_hash_format = "R serialization v2 via digest; not a portable semantic hash",
    known_slope = 0.8, slope_units = "synthetic outcome units per mg/dL",
    reference_category_source = "independent synthetic reference measurement before all injected errors",
    category_convention = "LOW: x < LLN; NORMAL: LLN <= x <= ULN; HIGH: x > ULN",
    lower = 0.60, upper = 1.30, target_mean_intercept = 0.95,
    target_residual_sd = 0.20, outcome_residual_sd = 0.20,
    missing_rate = 0.20, random_error_rate_observed = if (error_scenario) 0.08 else 0,
    flag_sensitivity = sensitivity, false_positive_rate = false_positive_rate,
    recording_shift_site003 = if (scenario == "site_shift") 0.20 else 0,
    record_rounding_digits = if (scenario == "rounded_records") 2L else NULL,
    analysis_formula = "outcome ~ LBSTRESN + age + sex + inflammation + SITEID",
    predictors = c("age", "sex", "inflammation", "SITEID", "outcome"),
    oracle_only_fields = names(truth)[-1L],
    clinical_closure_required = TRUE,
    source_transfer = list(
      project = "LLM_CDISC_Eval",
      generator_sha256 = "2f9aa33333f55bf3e7d34843ff8e5928f3c5785cde16a295eb4c5f9c54155124",
      freeze_sha256 = "684f248c6ed2f8522b25c2128fc487ff35a7197c9057191d656aec597f34921f",
      defects_sha256 = "9a9e7b326711561ee8ef0f4cdb52ea6d877486bbfc65794e97121093b4cb8dfa",
      scope = "conceptual controls and synthetic CREAT scale/limits; no private observations copied"
    )
  )
  structure(list(data = data, truth = truth, boundaries = lab_boundary_controls(),
                 manifest = manifest), class = "catmech_lab_fixture")
}

check_lab_fixture <- function(fixture) {
  d <- fixture$data
  t <- fixture$truth
  stopifnot(inherits(fixture, "catmech_lab_fixture"), nrow(d) == nrow(t),
            nrow(d) >= 30L, !anyDuplicated(d$USUBJID), identical(d$USUBJID, t$USUBJID),
            all(d$LBTESTCD == "CREAT"), all(d$VISITNUM == 1L),
            !any(names(t)[-1L] %in% names(d)))
  # Independent expression, not the generator's categorization function.
  expected <- rep("NORMAL", nrow(t))
  expected[t$LAB_TRUE < d$LBSTNRLO] <- "LOW"
  expected[t$LAB_TRUE > d$LBSTNRHI] <- "HIGH"
  stopifnot(identical(as.character(d$LBNRIND), expected),
            identical(is.na(d$LBSTRESN), t$MISSING))
  observed <- !t$MISSING
  stopifnot(isTRUE(all.equal(d$LBSTRESN[observed],
    (t$LAB_TRUE + t$VALUE_SHIFT + t$ROUNDING_ERROR)[observed], tolerance = 1e-13)))
  b <- fixture$boundaries
  stopifnot(identical(as.character(lab_reference_category(b$true_value)), b$expected_reference),
            identical(as.character(lab_reference_category(b$recorded_2dp)), b$expected_from_rounded))
  x <- d
  x$LBSTRESN <- t$LAB_TRUE
  x$outcome <- x$outcome - t$OUTCOME_ERROR
  exact <- stats::coef(stats::lm(outcome ~ LBSTRESN + age + sex + inflammation + SITEID, data = x))
  stopifnot(abs(exact[["LBSTRESN"]] - fixture$manifest$known_slope) < 1e-10)
  invisible(TRUE)
}

write_lab_fixture <- function(fixture, directory) {
  check_lab_fixture(fixture)
  if (dir.exists(directory) && length(list.files(directory, all.files = TRUE, no.. = TRUE)))
    stop("Fixture output directory must be new or empty.")
  if (!dir.exists(directory) && !dir.create(directory, recursive = TRUE))
    stop("Could not create fixture output directory.")
  saveRDS(fixture, file.path(directory, "fixture.rds"), version = 2)
  utils::write.csv(fixture$data, file.path(directory, "observed.csv"), row.names = FALSE)
  utils::write.csv(fixture$truth, file.path(directory, "oracle-truth.csv"), row.names = FALSE)
  utils::write.csv(fixture$boundaries, file.path(directory, "boundary-controls.csv"), row.names = FALSE)
  manifest <- fixture$manifest
  files <- c("fixture.rds", "observed.csv", "oracle-truth.csv", "boundary-controls.csv")
  manifest$files <- lapply(stats::setNames(files, files), function(name) list(
    sha256 = digest::digest(file = file.path(directory, name), algo = "sha256"),
    bytes = unname(file.info(file.path(directory, name))$size)))
  manifest$content_checks <- list(observed_rows = nrow(fixture$data),
    observed_columns = ncol(fixture$data), truth_rows = nrow(fixture$truth),
    unique_participants = length(unique(fixture$data$USUBJID)),
    truth_fields_absent_from_analysis_data = TRUE, independent_value_checks = "passed")
  jsonlite::write_json(manifest, file.path(directory, "manifest.json"),
                       auto_unbox = TRUE, pretty = TRUE, null = "null", digits = NA)
  invisible(normalizePath(directory))
}

lab_fixture_self_test <- function() {
  restore <- lab_fixture_rng_restore()
  on.exit(restore(), add = TRUE)
  set.seed(77)
  before <- get(".Random.seed", envir = .GlobalEnv)
  scenarios <- c("clean", "flagged_errors", "missed_flags", "site_shift", "rounded_records")
  fixtures <- lapply(scenarios, function(s) make_lab_fixture(scenario = s))
  names(fixtures) <- scenarios
  invisible(lapply(fixtures, check_lab_fixture))
  stopifnot(identical(before, get(".Random.seed", envir = .GlobalEnv)),
    identical(fixtures$clean, make_lab_fixture(scenario = "clean")),
    length(unique(vapply(fixtures, function(x) x$manifest$truth_hash, character(1)))) == 1L)
  # Demonstrate the injected mechanism in observable data, not only in labels.
  clean <- fixtures$clean$data$LBSTRESN
  for (name in scenarios[-1L]) {
    d <- fixtures[[name]]$data$LBSTRESN
    stopifnot(any(abs(d - clean) > 0, na.rm = TRUE))
  }
  t <- fixtures$missed_flags$truth
  stopifnot(any(t$IS_CORRUPTED & !t$FLAGGED), any(!t$IS_CORRUPTED & t$FLAGGED),
    all(fixtures$site_shift$truth$SITE_ERROR ==
      (!fixtures$site_shift$truth$MISSING & fixtures$site_shift$data$SITEID == "003")))
  different <- make_lab_fixture(phase = "heldout")
  stopifnot(!identical(different$manifest$truth_hash, fixtures$clean$manifest$truth_hash))
  # Known-wrong endpoint and rounding implementations must fail the boundary oracle.
  b <- lab_boundary_controls()
  wrong_closure <- as.character(cut(b$true_value, c(-Inf, .60, 1.30, Inf),
                                    labels = c("LOW", "NORMAL", "HIGH"), right = TRUE))
  stopifnot(any(wrong_closure != b$expected_reference),
            any(b$expected_from_rounded != b$expected_reference))
  invisible(data.frame(scenario = scenarios, checks = "pass", stringsAsFactors = FALSE))
}

if (sys.nframe() == 0L) {
  print(lab_fixture_self_test(), row.names = FALSE)
}
