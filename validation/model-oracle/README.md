# Development posterior check (superseded source)

This development run predates final input-validation edits. For the release source,
use [the final oracle record](../model-oracle-final/README.md).

Twelve trusted observations and eight interval-only observations were used in an intercept-only model.
The reference prior was p(mu,sigma²) proportional to 1/sigma². Six interval-only values were positive and two nonpositive.

The oracle integrates the trusted conjugate posterior over independent standard-normal and chi-square coordinates, weighted by the exact category likelihood. Base R Gaussian quadrature uses no package likelihood or Gibbs-update code.
Quadrature orders 64, 128 and 256 were compared, and the zero-category case was checked against closed-form posterior moments.

Four Gibbs chains used 1,500 burn-in and 5,000 retained iterations each. All four posterior moment checks passed: TRUE.
The known-wrong posterior omitting the category likelihood was detected: TRUE.
The acceptance tolerance is five Monte Carlo standard errors plus twice the last quadrature-order change, fixed before observing the chain results.

This checks the sampler's posterior moments for one known example. It does not establish frequentist confidence-interval coverage, robustness to misspecification, or posterior calibration across applications.

Reproduce from the repository root: `Rscript tools/check-model-oracle.R validation/new-model-oracle`.
See posterior-moment-checks.csv, quadrature-convergence.csv, conjugate-oracle-check.csv, negative-control.csv and source-hashes.csv.
