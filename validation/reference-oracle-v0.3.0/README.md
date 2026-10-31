# Historical v0.3.0 reference-model checks

This record refers to the frozen v0.3.0 sources. For the released guards, see
[the current oracle](../reference-oracle/README.md).

Overall result: **PASS**.

The normal path used 12000 independent completions of a fixed small dataset.
Direct linear algebra supplies the multivariate Student posterior predictive mean and covariance.
Each moment tolerance was prespecified as five estimated Monte Carlo standard errors plus 1e-12.
Covariances are estimated by products centered at the analytic mean, including cross-record terms.
No catmech fitting, posterior, or interval helper is used to calculate the analytic answer.

Data-contract checks test unused-recorded-value invariance, exact references, preservation,
PMM reference-donor membership, and the fully referenced ordinary-lm endpoint.
Source hashes were identical before and after the run.

These checks establish selected implementation properties. They do not establish repeated-sampling
coverage, correct reference selection, true gold-standard status, or clinical validity.
