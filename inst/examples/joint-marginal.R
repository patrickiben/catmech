# Correlated site effects in a continuous predictor and the true measurement.
# Synthetic Gaussian demonstration; this is not a coverage study.
set.seed(822301)
marginal_site <- rep(seq_len(24), each = 10)
marginal_context <- stats::rnorm(24)
marginal_x <- marginal_context[marginal_site] + stats::rnorm(240)
marginal_site_error <- matrix(stats::rnorm(48), 24, 2) %*%
  chol(matrix(c(.5, .45, .45, .6), 2))
marginal_row_error <- matrix(stats::rnorm(480), 240, 2) %*%
  chol(matrix(c(.5, .55, .55, .8), 2))
marginal_truth <- .5 + .7 * marginal_x + .4 * marginal_context[marginal_site] +
  marginal_site_error[marginal_site, 1] + marginal_row_error[, 1]
marginal_data <- data.frame(
  recorded = .8 + 1.05 * marginal_x + marginal_site_error[marginal_site, 2] +
    marginal_row_error[, 2],
  reference = marginal_truth, x = marginal_x, site = as.character(marginal_site)
)

# Exact references on 60 rows in 18 sites; six sites have no references.
# Recorded values are also missing on some rows within each site.
marginal_data$reference[-sample(which(marginal_site <= 18), 60)] <- NA_real_
marginal_data$recorded[seq(1, 240, 3)] <- NA_real_
rm(marginal_truth, marginal_context, marginal_site_error, marginal_row_error,
   marginal_site, marginal_x)

# Jointly model (true, recorded, x), allowing their site effects to correlate.
# The predictors must be complete continuous numeric variables. This model
# assumes independent Gaussian sites/rows and noninformative site sizes.
marginal_fit <- catmech::catmech_reference(
  marginal_data, "recorded", "reference", "x", use_recorded = TRUE,
  cluster = "site", cluster_model = "joint_marginal", m = 40, seed = 822302
)
print(marginal_fit$model_diagnostics$attempts)

# Population marginal association, with conditional model-based credible
# intervals. In this generating model the marginal slope is 0.9, while the
# within-site structural coefficient is 0.7; these are different estimands.
# The summary uses all retained posterior draws independently of m.
marginal_posterior <- catmech::catmech_marginal_summary(marginal_fit)
print(marginal_posterior)

# Completed data retain every observed reference exactly and can be accessed
# through the existing mice adapter. This does not make other analyses inherit
# the marginal posterior summary's interpretation or validation.
marginal_completed <- mice::complete(marginal_fit$mids, 1)
