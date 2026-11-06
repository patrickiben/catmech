#' Simulate measurements, quality flags, and interval categories
#'
#' Synthetic data only. The coefficient of lab_true in the outcome model is
#' 0.8. Flags occur among observed values; a configurable fraction indicate
#' actual corruption. category_error_rate deliberately violates the package's
#' exact-category assumption. Truth fields must never enter imputation predictors.
#' @param n Number of rows, at least 20.
#' @param seed Nonnegative integer seed, with caller RNG restored.
#' @param corruption_probability Probability of corruption conditional on a flag.
#' @param missing_rate Probability the continuous value is missing (MCAR).
#' @param flag_rate Flag probability among observed measurements.
#' @param category_error_rate Probability of deliberately miscoding a category.
#' @param error_distribution "normal" or variance-standardized "heavy_tailed"
#'   Student t residuals with 3 degrees of freedom for the latent measurement.
#' @param thresholds Four interval endpoints for Low, Medium, High.
#' @return A data.frame with lab_true and corrupted for validation only.
#' @examples
#' dat <- simulate_catmech_lab(n = 50, seed = 123)
#' head(dat)
#' @export
simulate_catmech_lab <- function(n = 300L, seed = 1L, corruption_probability = 1,
                                missing_rate = 0.2, flag_rate = 0.1,
                                category_error_rate = 0,
                                error_distribution = c("normal", "heavy_tailed"),
                                thresholds = c(-Inf, 5.5, 7, Inf)) {
  check_count(n, "n", 20L)
  check_count(seed, "seed", 0L)
  check_probability(corruption_probability, "corruption_probability")
  check_probability(missing_rate, "missing_rate")
  check_probability(flag_rate, "flag_rate")
  check_probability(category_error_rate, "category_error_rate")
  error_distribution <- match.arg(error_distribution)
  make_interval_bounds(thresholds, c("Low", "Medium", "High"))
  restore <- preserve_rng()
  on.exit(restore(), add=TRUE)
  set.seed(seed)
  age <- stats::rnorm(n, 50, 10)
  sex <- factor(sample(c("F", "M"), n, replace=TRUE), levels=c("F", "M"))
  inflammation <- stats::rnorm(n, 2, 1)
  e <- if (error_distribution == "normal") stats::rnorm(n) else stats::rt(n, 3)/sqrt(3)
  lab_true <- 3 + 0.05*age + 0.5*inflammation + e
  category <- catmech_derive_category(lab_true, thresholds, c("Low", "Medium", "High"))
  wrong <- which(stats::runif(n) < category_error_rate)
  for (i in wrong) category[i] <- sample(setdiff(levels(category), as.character(category[i])), 1L)
  hidden <- stats::runif(n) < missing_rate
  flagged <- !hidden & stats::runif(n) < flag_rate
  corrupted <- flagged & stats::runif(n) < corruption_probability
  lab_value <- lab_true
  lab_value[hidden] <- NA_real_
  lab_value[corrupted] <- lab_value[corrupted] + stats::rnorm(sum(corrupted), 10, 2)
  outcome <- 0.8*lab_true + stats::rnorm(n)
  data.frame(patient_id=seq_len(n), age=age, sex=sex, inflammation=inflammation,
    lab_value=lab_value, lab_category=category,
    outlier_flag=factor(ifelse(flagged, "Yes", "No"), levels=c("No","Yes")),
    outcome=outcome, lab_true=lab_true, corrupted=corrupted)
}
