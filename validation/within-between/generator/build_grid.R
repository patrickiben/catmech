# build_grid.R -- defines every cell and writes grid.csv, including each target.
#
# A target is stored from a closed form only where the closed form is exact (see
# closed_form below). Every other target comes from a fixed large Monte Carlo run
# of full-truth OLS, on seeds separate from those in check.R. Closed-form cells
# also get the same run, stored as target_mc, as an extra test of the algebra.
#
# Usage: Rscript build_grid.R [cohorts_closed] [cohorts_mc]

args <- commandArgs(FALSE)
script <- sub("^--file=", "", args[grep("^--file=", args)])
if (length(script)) setwd(dirname(normalizePath(script)))
source("generator.R")
library(parallel)

targs <- commandArgs(TRUE)
N_CLOSED <- if (length(targs) >= 1) as.integer(targs[1]) else 10000L
N_MC     <- if (length(targs) >= 2) as.integer(targs[2]) else 40000L
CORES    <- max(1L, detectCores() - 2L)

defaults <- list(
  # identity and design
  label = NA, role = NA, mechanism = NA, why = NA,
  K = 60, size_type = "balanced", size_mean = 10, size_cv = 0.5,
  n_strata = 0, strata_skew = 0, size_shift = 0, audit_tilt = 1,
  # covariate
  x_type = "normal", x_icc = 0.3, x_center = 0, x_size = 0, x_strat = 0,
  x_mean = 0, x_sd = 1,
  # truth
  truth_type = "linear", b0 = 12, b_x = 0.5, b_ctx = 0, b_x2 = 0, hfx = 0, hslope = 0,
  s_slope = 0, sd_u = 0.5, tail_u = "normal", u_size = 0, nom_sd = 0,
  sd_e = 1, tail_e = "normal", e_het = 0,
  # outcome (predictor role only)
  y_link = "identity", y_thr = 0, gamma = 0, delta = 0, c0 = 0, h_y = 0, g_slope = 0,
  ctx_truth = 0, sd_v = 0, tail_v = "normal", v_size = 0, sd_eps = 1, tail_eps = "normal",
  # recording errors (see draw_recorded in generator.R; 0 means off)
  berkson = 0, err_value = 0, prop_mean = 0, prop_sd = 0, mult_sd = 0,
  clu_mean = 0, clu_sd = 0, tail_clu = "normal", err_size = 0, err_v = 0,
  err_x = 0, err_x_clu = 0, err_diff = 0, drift_sd = 0,
  noise_sd = 0, tail_noise = "normal", noise_het_value = 0, noise_clu_sd = 0,
  err_s1_mean = 0, err_s1_sd = 0,
  compl_mean = 0, compl_sd = 0, inflate_share = 0, inflate = 0, inflate_x = 0,
  gaming_share = 0, gaming_thr = 0, gaming_band = 0, gaming_prob = 0,
  lod = 0, ceil = 0, fab_share = 0, unit_share = 0, unit_factor = 1,
  gross_rate = 0, gross_disp = 0, gross_type = "scale", gross_factor = 10,
  round_to = 0, round_share = 1, copy_rate = 0
)

mk <- function(label, mechanism, why, ...) {
  over <- list(...)
  bad <- setdiff(names(over), names(defaults))
  if (length(bad)) stop(label, ": unknown parameter ", paste(bad, collapse = ", "))
  x <- defaults
  x[names(over)] <- over
  x$label <- label
  x$mechanism <- mechanism
  x$why <- why
  x$role <- if (startsWith(label, "o_")) "outcome" else "predictor"
  if (x$role == "predictor") {
    # Default dose response unless the cell says otherwise.
    if (!"gamma" %in% names(over)) x$gamma <- 0.5
    if (!"delta" %in% names(over)) x$delta <- 0.3
    if (!"c0" %in% names(over)) x$c0 <- 2
    if (!"sd_v" %in% names(over)) x$sd_v <- 0.5
  }
  as.data.frame(x, stringsAsFactors = FALSE)
}

# Reusable scales. Hemoglobin g/dL is the default (b0 12). Others:
glu <- list(b0 = 100, b_x = 8, sd_u = 6, sd_e = 20)      # glucose mg/dL
bp  <- list(b0 = 130, b_x = 5, sd_u = 4, sd_e = 15)      # systolic BP mmHg
vl  <- list(b0 = 3, b_x = -0.5, sd_u = 0.4, sd_e = 1)    # log10 viral load
cnt <- list(truth_type = "count", b0 = 3, b_x = 0.3, sd_u = 0.5, sd_e = 0.2)  # monthly counts
cost <- list(truth_type = "lognormal", b0 = 6.5, b_x = 0.4, sd_u = 0.5, sd_e = 1.1)  # claim dollars
with_scale <- function(s, ...) {
  a <- list(...)                               # the cell's own values win over the scale's
  do.call(mk, c(a, s[setdiff(names(s), names(a))]))
}

cells <- list(
  # ---- OUTCOME ROLE -------------------------------------------------------
  mk("o_control_errorfree", "ERROR-FREE control; recorded equals true",
     "A negative control: audit finds nothing. Clean electronic import from the analyser."),
  mk("o_noise_noclusters", "Row assay noise; no cluster structure at all",
     "Imprecision of a single central lab: sites carry no information about value or error.",
     sd_u = 0, x_icc = 0, noise_sd = 0.5),
  mk("o_easy_large", "150 sites, mild noise, small offsets",
     "A well-run network with modest imprecision. Nothing should struggle here.",
     K = 150, sd_u = 0.3, x_icc = 0.2, noise_sd = 0.25, clu_sd = 0.1),
  mk("o_pairs_k150", "Two rows per cluster, site offsets and noise",
     "Paired records such as mother-infant or two-visit audits: many clusters, each tiny.",
     K = 150, size_mean = 2, sd_u = 1, x_icc = 0.5, noise_sd = 0.5, clu_sd = 0.5),
  mk("o_offset_x_within", "Large site offsets, covariate varies only within site",
     "Within-site allocation (a shift, a ward, a sex-balanced sample): site bias is large but orthogonal to the covariate.",
     K = 40, size_mean = 15, x_icc = 0, x_center = 1, sd_u = 1.5, clu_sd = 1.5, noise_sd = 0.2),
  mk("o_offset_x_site", "Binary site-level covariate, site offset dominates",
     "Facility ownership against a facility-calibrated value. Clustering dominates; effective sample is the 40 sites.",
     K = 40, size_mean = 8, x_type = "binary", x_icc = 1, b_x = 0.6, sd_u = 1, sd_e = 0.6, clu_sd = 1, noise_sd = 0.1),
  mk("o_offset_few_unbal", "20 sites, hospitals and posts, heavy-tailed offsets",
     "A small district network: a few big hospitals hold most rows, and one odd lab can sit far off.",
     K = 20, size_type = "two_tier", x_icc = 0.5, sd_u = 0.8, tail_u = "t4",
     clu_sd = 1, tail_clu = "t4", noise_sd = 0.3),
  with_scale(glu, "o_calib_slope", "Site calibration slope, error grows with the value",
     "Glucose analysers drift in percent, so high readings carry large site-shared error.",
     K = 40, size_type = "gamma", size_mean = 12, x_icc = 0.2, prop_mean = 0.03, prop_sd = 0.06, noise_sd = 2),
  mk("o_regress_to_mean", "Row-level compression toward the mean",
     "Self-report and recall: people and staff understate extremes, so error depends on the value.",
     K = 80, size_mean = 6, sd_u = 0.4, sd_e = 1.2, err_value = -0.4, noise_sd = 0.4),
  mk("o_bias_by_x", "Common bias depending on the covariate, three strata",
     "Busier or older-patient clinics under-record systematically, directly bending the slope.",
     K = 50, size_type = "gamma", size_mean = 12, n_strata = 3, hfx = 0.8, x_strat = 0.5,
     x_icc = 0.2, err_x = -0.3, noise_sd = 0.3),
  mk("o_site_bias_by_x", "Site-specific bias slope on the covariate",
     "One clinic misreads a subgroup and another the opposite: the error slope varies by site, zero on average.",
     K = 35, size_type = "gamma", size_mean = 12, x_icc = 0.2, err_x_clu = 0.5, noise_sd = 0.2),
  with_scale(glu, "o_hetero_value", "Noise grows with distance from the mean, sloppy sites",
     "Point-of-care meters lose precision at extreme readings, and some sites handle them worse.",
     K = 65, size_type = "gamma", size_mean = 10, x_type = "skew", x_mean = 50, x_sd = 12,
     noise_sd = 5, noise_het_value = 1, noise_clu_sd = 0.5),
  with_scale(bp, "o_digit_preference", "Rounding to 10 at most sites, unbalanced sizes",
     "Terminal-digit preference in manual blood pressure is widely seen; some sites use automated cuffs.",
     K = 60, size_type = "lognormal", size_cv = 0.7, round_to = 10, round_share = 0.6),
  mk("o_decimal_slip", "Sparse 10-fold slips, clumped by site",
     "Decimal point errors in keyed lab values: rare, enormous, and concentrated with certain clerks.",
     K = 80, size_type = "gamma", size_mean = 12, gross_rate = 0.02, gross_disp = 0.15,
     gross_factor = 10, noise_sd = 0.1),
  with_scale(glu, "o_unit_mixup_sites", "Whole sites report in another unit",
     "A site reports mmol/L in an mg/dL field. Every row at that site is off by 18-fold, and few sites do it.",
     K = 45, size_mean = 12, unit_share = 0.1, unit_factor = 1 / 18, noise_sd = 2),
  with_scale(vl, "o_lod_floor", "Limit-of-detection floor, binary covariate",
     "Viral load below the assay limit is reported at the limit: value-dependent, one-sided error.",
     K = 50, size_type = "gamma", size_mean = 10, x_type = "binary", x_icc = 0.2, lod = 2, noise_sd = 0.1),
  with_scale(cnt, "o_counts_underreport", "Register counts thinned by site completeness",
     "Monthly totals: events go missing from the summary form at a rate specific to the facility.",
     K = 60, size_type = "lognormal", size_cv = 0.8, compl_mean = 0.8, compl_sd = 0.5),
  with_scale(cnt, "o_counts_inflate_pay", "Sites paid by volume inflate, more at high covariate",
     "Performance-based financing rewards reported volume, and some facilities respond.",
     K = 70, size_type = "gamma", size_mean = 10, b0 = 3.2, inflate_share = 0.3, inflate = 0.25, inflate_x = 0.5),
  with_scale(cost, "o_claims_lognormal", "Heavy-tailed costs with percent errors by biller",
     "Claim amounts are long-tailed and billers vary in inflating them. Hard: a few rows drive everything.",
     K = 100, size_type = "lognormal", size_mean = 8, size_cv = 0.9, prop_mean = 0.05, prop_sd = 0.1, mult_sd = 0.2),
  mk("o_heavy_t3_k30", "t3 true values and t3 errors at both levels, 30 sites",
     "Mixed populations plus occasional bad runs give rare extremes at both patient and site level.",
     K = 30, size_mean = 12, tail_u = "t3", tail_e = "t3", clu_sd = 0.6, tail_clu = "t3",
     noise_sd = 0.5, tail_noise = "t3"),
  mk("o_contextual", "Contextual effect of site covariate mean",
     "Case mix: a site's composition moves everyone there (referral pattern), so between and within slopes differ.",
     K = 45, x_icc = 0.5, b_x = 0.3, b_ctx = 1, clu_sd = 0.5, noise_sd = 0.3),
  mk("o_size_informative", "Size linked to true value, covariate and error",
     "Big hospitals have sicker patients, different values and better records than small posts.",
     K = 40, size_type = "lognormal", size_cv = 0.9, u_size = 0.8, x_size = 0.5, err_size = -0.4,
     noise_sd = 0.3, noise_clu_sd = 0.5),
  mk("o_random_slope_k25", "Covariate effect varies by site, 25 sites",
     "A within-site allocation whose effect differs by local practice, so the row-level slope is a site average.",
     K = 25, size_mean = 16, x_type = "binary", x_icc = 0, x_center = 1, s_slope = 0.5, sd_u = 0.6, clu_sd = 0.4, noise_sd = 0.3),
  mk("o_strata_one_bad", "One small stratum bad and oversampled in the audit",
     "Private facilities outside the support system err more, and auditors deliberately sample more of them.",
     K = 90, size_type = "gamma", size_mean = 10, n_strata = 3, strata_skew = 1, hfx = 1, size_shift = 0.3,
     audit_tilt = 2, noise_sd = 0.1, err_s1_mean = 1, err_s1_sd = 1),
  mk("o_strata_slopes", "Covariate slope differs by stratum",
     "Urban and rural facilities differ in the covariate's effect, so the pooled slope is a weighted average.",
     K = 60, n_strata = 4, hslope = 0.4, x_strat = 0.6, hfx = 0.5, clu_sd = 0.4, noise_sd = 0.3),
  mk("o_copy_forward", "Copy-forward of the line above, 30 rows per site",
     "Register fatigue: a clerk repeats the previous entry, so errors are bounded in size but look real.",
     K = 50, size_mean = 30, x_icc = 0.2, copy_rate = 0.3, noise_sd = 0.2),
  mk("o_fabricated_sites", "A few sites fabricate believable values",
     "Desk-filled forms at remote sites: right range, no link to the patient.",
     K = 120, size_type = "two_tier", fab_share = 0.1, noise_sd = 0.1),
  with_scale(bp, "o_threshold_gaming", "Values nudged above a payment threshold at some sites",
     "Pay-for-performance targets: readings just below the cut are pushed over it, creating a heap.",
     K = 70, size_type = "gamma", size_mean = 10, b0 = 125, gaming_share = 0.4, gaming_thr = 130,
     gaming_band = 10, gaming_prob = 0.6),
  mk("o_zero_for_blank", "Blanks keyed as zero, rate clumped by site",
     "Registers with a mandatory field: an untested patient becomes 0. Common in aggregate and lab data.",
     K = 55, size_type = "gamma", size_mean = 10, gross_rate = 0.05, gross_disp = 0.2, gross_type = "zero"),
  mk("o_drift_within_run", "Site-specific drift across rows of a run",
     "Reagent lots and calibrations age through a run, at site-specific rates, unrelated to who is tested.",
     K = 40, size_mean = 25, x_icc = 0, drift_sd = 1, noise_sd = 0.2),
  mk("o_tiny_k20_n3", "20 sites of 3 rows, skewed covariate",
     "A small spot-check study. Few rows mean little to learn from and much to fear from any one site.",
     K = 20, size_mean = 3, x_type = "skew", x_mean = 45, x_sd = 12, b_x = 0.6, noise_sd = 0.5, clu_sd = 0.4),
  mk("o_quadratic_truth", "True value curved in a skewed covariate",
     "Biomarkers often move nonlinearly with age or load, so OLS's slope is only a projection.",
     K = 70, size_type = "gamma", size_mean = 10, x_type = "skew", b_x = 0.4, b_x2 = 0.15, noise_sd = 0.4, clu_sd = 0.3),
  mk("o_two_tier_k150", "150 sites, many posts and few hospitals, site offsets",
     "A national network: many clusters overall, yet most rows sit in a small share of them.",
     K = 150, size_type = "two_tier", clu_sd = 0.6, noise_sd = 0.3),

  # ---- PREDICTOR ROLE -----------------------------------------------------
  mk("p_control_errorfree", "ERROR-FREE control; recorded equals true",
     "A negative control for the predictor role. Clean electronic import."),
  mk("p_noise_noclusters", "Classical row noise; no cluster structure at all",
     "One central lab, independent patients, and a modest assay imprecision.",
     sd_u = 0, sd_v = 0, x_icc = 0, noise_sd = 0.5),
  mk("p_easy_large", "150 sites, mild noise",
     "A large, well-run network with reliable data. Nothing should struggle here.",
     K = 150, sd_u = 0.3, sd_v = 0.3, noise_sd = 0.2),
  mk("p_low_reliability", "Row noise as large as the true variation, unbalanced sizes",
     "Single spot measures of a variable quantity, such as one blood pressure reading for a habitual level.",
     K = 100, size_type = "lognormal", size_mean = 8, size_cv = 0.8, sd_u = 0.3, sd_v = 0.3, noise_sd = 1),
  mk("p_offset_mixed", "Site offset plus row noise at similar scale",
     "A multi-lab network: each lab has its own offset and its own imprecision.",
     K = 50, size_type = "gamma", size_mean = 12, x_icc = 0.4, sd_u = 0.7, clu_mean = 0.5, clu_sd = 0.7, noise_sd = 0.3),
  mk("p_site_dominant", "True value and error both nearly all between sites",
     "Value set by site (water source, altitude, local protocol) and the error too. Clustering dominates.",
     K = 30, size_mean = 8, x_icc = 0.8, sd_u = 1.5, sd_e = 0.3, clu_sd = 1, noise_sd = 0.05, sd_v = 1, sd_eps = 0.5),
  mk("p_offset_few_unbal", "20 sites, hospitals and posts, heavy-tailed site effects",
     "A small network where one big hospital sets most of the data and the audit may miss it.",
     K = 20, size_type = "two_tier", x_icc = 0.5, clu_sd = 0.8, tail_clu = "t4", noise_sd = 0.4,
     sd_v = 0.6, tail_v = "t4"),
  with_scale(glu, "p_calib_slope", "Site calibration slope on the predictor",
     "Percent-off analysers feeding a regression: the size of the error scales with the value.",
     K = 45, size_type = "gamma", size_mean = 12, gamma = 0.02, delta = 0.3, sd_v = 0.5, c0 = 2,
     prop_sd = 0.08, noise_sd = 3),
  mk("p_berkson_nominal", "Berkson: the register holds a nominal assigned value",
     "Prescribed dose or a kit's lot value is recorded and the actual value scatters around it.",
     K = 60, size_type = "gamma", size_mean = 10, nom_sd = 0.9, sd_e = 0.8, berkson = 1),
  mk("p_differential_outcome", "Error depends on the outcome",
     "Unblinded recording: a patient who did badly gets a value written with that in mind.",
     K = 70, size_type = "gamma", size_mean = 12, err_diff = 0.4, noise_sd = 0.3),
  mk("p_shared_site_quality", "Site quality moves both outcome and recording bias",
     "Well-run facilities get better outcomes and keep different records, so error links to the site's outcome residual.",
     K = 40, size_mean = 12, sd_v = 1, err_v = 0.8, clu_sd = 0.3, noise_sd = 0.2),
  mk("p_bias_by_x", "Bias depending on the covariate, two strata",
     "Older or sicker patients are measured differently, so the error is tied to another regressor.",
     K = 50, size_type = "gamma", size_mean = 10, n_strata = 2, hfx = 0.6, h_y = 0.5, x_icc = 0.3,
     err_x = 0.4, noise_sd = 0.3),
  mk("p_decimal_slip", "Sparse 10-fold slips on the predictor",
     "A few enormous keyed values have huge leverage when the value is a regressor.",
     K = 80, size_type = "gamma", size_mean = 12, gross_rate = 0.02, gross_disp = 0.15, noise_sd = 0.1),
  with_scale(vl, "p_lod_floor", "Limit-of-detection floor on the predictor",
     "A censored viral load is used as an exposure, a common practice with a value-dependent error.",
     K = 50, size_type = "gamma", size_mean = 10, x_type = "binary", x_icc = 0.2, gamma = -0.4,
     delta = 0.3, sd_v = 0.5, c0 = 2, lod = 2, noise_sd = 0.1),
  with_scale(bp, "p_digit_preference", "Rounding to 10 at most sites",
     "Hand-rounded blood pressure used as an exposure. Errors are bounded but tied to the value.",
     K = 60, size_type = "lognormal", size_cv = 0.7, gamma = 0.03, delta = 0.3, sd_v = 0.5, c0 = 2,
     round_to = 10, round_share = 0.6),
  mk("p_contextual_truth", "Outcome depends on site mean of the true value",
     "Peer or environmental effect: site average matters beyond the individual's value, so OLS mixes both.",
     K = 45, gamma = 0.3, ctx_truth = 0.8, sd_u = 0.7, clu_sd = 0.4, noise_sd = 0.4),
  mk("p_log_dose", "Lognormal predictor, outcome linear in its log",
     "Biomarkers such as ferritin act on a log scale but are analysed raw, so the target is a projection.",
     K = 60, size_type = "gamma", size_mean = 10, truth_type = "lognormal", b0 = 3, b_x = 0.3, sd_u = 0.3,
     sd_e = 0.6, y_link = "log", gamma = 0.8, prop_sd = 0.1, mult_sd = 0.15),
  with_scale(bp, "p_threshold_gaming", "Outcome steps at a threshold that recording games",
     "A clinical cut drives care and hence outcomes, and staff push values over that same cut.",
     K = 70, size_type = "gamma", size_mean = 10, b0 = 125, y_link = "step", y_thr = 130, gamma = 1,
     delta = 0.3, sd_v = 0.5, c0 = 2, gaming_share = 0.4, gaming_thr = 130, gaming_band = 10, gaming_prob = 0.6),
  mk("p_heavy_t3_k30", "t3 throughout: truth, outcome and errors, 30 sites",
     "Outliers everywhere, in values and in outcomes, with few sites to average them away.",
     K = 30, size_mean = 12, tail_u = "t3", tail_e = "t3", tail_v = "t3", tail_eps = "t3",
     clu_sd = 0.5, tail_clu = "t3", noise_sd = 0.5, tail_noise = "t3"),
  mk("p_size_informative", "Size linked to true value, outcome, covariate and error",
     "Large sites differ in everything. Which sites get audited changes what the audit looks like.",
     K = 40, size_type = "lognormal", size_cv = 0.9, u_size = 0.6, v_size = 0.6, x_size = 0.4, err_size = -0.4,
     noise_sd = 0.3),
  mk("p_strata_one_bad", "One small stratum bad and oversampled in the audit",
     "Private facilities err more, and auditors deliberately sample more of them.",
     K = 90, size_type = "gamma", size_mean = 10, n_strata = 3, strata_skew = 1, hfx = 1, h_y = 0.5,
     size_shift = 0.3, audit_tilt = 2, noise_sd = 0.1, err_s1_mean = 1, err_s1_sd = 1),
  mk("p_random_slope_k25", "Dose response varies by site, 25 sites",
     "Local care changes how the value relates to outcome, so the row-level slope is a site average.",
     K = 25, size_mean = 16, g_slope = 0.3, clu_sd = 0.3, noise_sd = 0.4),
  mk("p_copy_forward", "Copy-forward of the line above, 30 rows per site",
     "Repeat entries unlink value from patient, while the outcome stays with the patient.",
     K = 50, size_mean = 30, x_icc = 0.2, copy_rate = 0.3, noise_sd = 0.2),
  mk("p_fabricated_sites", "A few sites fabricate believable values",
     "Plausible desk-filled values, unrelated to the patient, while outcomes come from another system.",
     K = 60, size_type = "gamma", size_mean = 8, fab_share = 0.1, noise_sd = 0.1),
  with_scale(cnt, "p_counts_underreport", "Counts thinned by site completeness",
     "Facility volume from a register predicts an outcome, with completeness varying by site.",
     K = 120, size_type = "lognormal", size_cv = 0.8, gamma = 0.05, delta = 0.3, sd_v = 0.5, c0 = 2,
     compl_mean = 0.8, compl_sd = 0.5),
  mk("p_x_site_level", "Site-level covariate strongly tied to the true value",
     "A facility attribute predicts the value, so its effect and the value's are hard to separate.",
     K = 25, size_mean = 12, x_icc = 1, b_x = 1, sd_u = 0.5, clu_sd = 0.3, noise_sd = 0.4),
  mk("p_zero_for_blank", "Blanks keyed as zero, rate clumped by site",
     "An untested patient recorded as 0 becomes a huge outlier when the value is a regressor.",
     K = 55, size_type = "gamma", size_mean = 10, gross_rate = 0.05, gross_disp = 0.2, gross_type = "zero"),
  mk("p_two_tier_k150", "150 sites, many posts and few hospitals, site offsets",
     "A national network: many clusters overall, yet most rows sit in a small share of them.",
     K = 150, size_type = "two_tier", clu_sd = 0.6, noise_sd = 0.3),
  mk("p_tiny_k20_n3", "20 sites of 3 rows, skewed covariate",
     "A tiny spot-check. The audit sees 15 rows at the lower fraction.",
     K = 20, size_mean = 3, x_type = "skew", x_mean = 45, x_sd = 12, noise_sd = 0.5, clu_sd = 0.4)
)
grid <- do.call(rbind, cells)
if (anyDuplicated(grid$label)) stop("duplicate labels")

# Exact closed forms. OLS is conditionally unbiased when the residual has mean
# zero given the regressors, whatever its shape or dependence on the regressor's
# size. Anything that breaks linearity or ties a residual to size or a site mean
# goes to Monte Carlo.
closed_form <- function(p) {
  if (p$role == "outcome") {
    if (p$truth_type == "linear" && p$b_ctx == 0 && p$b_x2 == 0 && p$hslope == 0 && p$u_size == 0)
      return(if (p$x_type == "binary") p$b_x else p$b_x / p$x_sd)
    return(NA_real_)
  }
  if (p$y_link == "identity" && p$ctx_truth == 0 && p$v_size == 0) return(p$gamma)
  NA_real_
}

full_truth_coef <- function(sim) {
  d <- sim$input
  d$lab_value <- sim$truth
  unname(coef(lm(sim$formula, data = d))[sim$term])
}

SEED_BASE <- 100000000L
SEED_STEP <- 1000000L
run_target <- function(i, n) {
  cell <- as.list(grid[i, ])
  cell$target <- NA_real_
  splits <- parallel::splitIndices(n, CORES)
  res <- mclapply(splits, function(idx) {
    vapply(idx, function(r) full_truth_coef(
      simulate_cell(cell, SEED_BASE + i * SEED_STEP + r, 0.5)), numeric(1))
  }, mc.cores = CORES)
  b <- unlist(res)
  c(mean(b), sd(b) / sqrt(length(b)), sum(is.na(b)))
}

t0 <- Sys.time()
grid$target <- NA_real_
grid$target_method <- NA_character_
grid$target_mc <- NA_real_
grid$target_mc_se <- NA_real_
grid$target_mc_n <- NA_integer_
for (i in seq_len(nrow(grid))) {
  closed <- closed_form(as.list(grid[i, ]))
  n <- if (is.na(closed)) N_MC else N_CLOSED
  r <- run_target(i, n)
  if (r[3] > 0) stop(grid$label[i], ": rank-deficient cohorts in target run")
  grid$target[i] <- if (is.na(closed)) r[1] else closed
  grid$target_method[i] <- if (is.na(closed)) "monte_carlo" else "closed_form"
  grid$target_mc[i] <- r[1]
  grid$target_mc_se[i] <- r[2]
  grid$target_mc_n[i] <- n
  cat(sprintf("%-26s %-11s target %9.5f  run %9.5f se %.5f n %d\n", grid$label[i],
              grid$target_method[i], grid$target[i], r[1], r[2], n))
}
write.csv(grid, "grid.csv", row.names = FALSE)
writeLines(c(paste("generator_md5", unname(tools::md5sum("generator.R"))),
             paste("build_md5", unname(tools::md5sum("build_grid.R"))),
             paste("n_closed", N_CLOSED), paste("n_mc", N_MC),
             paste("seed_base", SEED_BASE), paste("seed_step", SEED_STEP),
             paste("elapsed_sec", round(as.numeric(difftime(Sys.time(), t0, units = "secs"))))),
           "target_run.txt")
cat("cells", nrow(grid), "predictor", sum(grid$role == "predictor"), "\n")
cat("elapsed", format(Sys.time() - t0), "\n")
