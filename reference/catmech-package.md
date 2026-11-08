# catmech: Auditable Model-Based Imputation and Category Sensitivity

Auditable completion of a continuous measurement using Gaussian models,
predictive mean matching, an explicitly measured internal reference
subset, or a reference-calibrated contamination model of the recorded
values. Supports flag masking, optional exact or fallible interval
categories, posterior draws and per-record corruption probabilities,
convergence and donor-support diagnostics, reference hold-out checks,
finite-scenario sensitivity comparisons and 'mice' analysis adapters.
Retains a separate retain-or-replace sensitivity engine. Implements
multiple-imputation combining rules and finite-sample degrees of freedom
following Rubin (1987)
[doi:10.1002/9780470316696](https://doi.org/10.1002/9780470316696) and
Barnard and Rubin (1999)
[doi:10.1093/biomet/86.4.948](https://doi.org/10.1093/biomet/86.4.948) ,
and an unvalidated audit-size planning aid based on two-phase sampling
theory (Neyman 1938)
[doi:10.1080/01621459.1938.10503378](https://doi.org/10.1080/01621459.1938.10503378)
. Does not identify measurement error without assumptions or reference
information, or guarantee coverage under misspecification.

## See also

Useful links:

- <https://github.com/patrickiben/catmech>

- Report bugs at <https://github.com/patrickiben/catmech/issues>

## Author

**Maintainer**: Patrick Iben <patrickiben@gmail.com>
([ORCID](https://orcid.org/0009-0005-2443-7973))

Authors:

- Patrick Iben <patrickiben@gmail.com>
  ([ORCID](https://orcid.org/0009-0005-2443-7973))
