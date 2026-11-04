# Mixed recorded-present/recorded-missing rows within sites.
# Synthetic Gaussian demonstration; this is not a coverage study.
set.seed(801)
joint_data <- data.frame(site = rep(seq_len(40), each = 6),
                         visit = rep(seq_len(6), 40), x = stats::rnorm(240))
joint_truth <- 1 + 0.6 * joint_data$x +
  stats::rnorm(40, sd = 1.2)[joint_data$site] + stats::rnorm(240)
joint_data$recorded <- 0.7 * joint_truth +
  stats::rnorm(40, sd = 0.7)[joint_data$site] + stats::rnorm(240)
joint_data$reference <- ifelse(joint_data$site <= 20 & joint_data$visit <= 3,
                               joint_truth, NA_real_)
joint_data$recorded[joint_data$visit %in% c(3, 5)] <- NA_real_
# A wholly unobserved site can also be predicted.
joint_data$recorded[joint_data$site == 40] <- NA_real_
rm(joint_truth)

joint_fit <- catmech::catmech_reference(
  joint_data, "recorded", "reference", "x", use_recorded = TRUE,
  cluster = "site", cluster_model = "joint", m = 40, seed = 802
)
print(joint_fit$model_diagnostics$summary)
joint_analysis <- catmech::catmech_pool_lm(joint_fit, recorded ~ x)
print(joint_analysis)

# Population coefficients of the fitted joint model, with model-based credible
# intervals. This uses all posterior draws and the original imputation design.
joint_posterior <- catmech::catmech_joint_summary(joint_fit)
print(joint_posterior)
