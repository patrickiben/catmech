# Validation record: catmech 0.10.1

Date: 2026-11-04. This is a numerical-budget update to the experimental marginal
joint model. The [v0.10.0 record](v0.10.0-VALIDATION.md), failed availability
screen and original source/evidence remain preserved.

The marginal default maximum increases to 64,000 retained draws per chain.
Every monitor, R-hat <= 1.05, bulk/tail ESS >= 400, four chains and the original
burn-in/doubling/restart rule remain fixed. Successful earlier attempts stop
as before; priors, likelihood, interval rules and conditional-joint defaults
are unchanged. Both original failures pass at 64k.

The complete dependent replay returned all 12,000 MCMC fits and 27,600 intervals
from 3600 cohorts. All four screens pass. A strict comparison finds all 27,596
previously successful analyses unchanged, with only the four formerly refused
analyses recovered. [Full report and reproducibility](MARGINAL_EXTENSION_SYNTHESIS.md).
This is not new independent confirmation or a guarantee for future inputs.

The initial forked replay is preserved as incomplete following a worker
segfault. Affected cohort 34 passed in independently started R workers; the
entire final grid was then run from scratch using PSOCK workers with explicit
missing-result checks. No interrupted records were reused.

Independent reporting/grid/diagnostic/arithmetic verification and deliberate
corruption checks are recorded beside the complete replay. All 1,526 assertions
pass on Ubuntu, Windows and macOS CI, including the resumed macOS check on a
pinned macos-15 image. Local source/archive/install check receipts are under
`v0.10.1-checks/`. The original model assumptions, one-predictor validation
scope, residual biases and unsupported data restrictions remain in force.


## Final source-package verification

The built v0.10.1 archive passed `R CMD check --as-cran --no-manual` with zero
errors, zero warnings and one routine “New submission” note. All 1,526 assertions
passed without skips. A fresh installation passed all seven bundled examples,
namespace-only native sampling, adapters and exports. Archive integrity verified
all 69 members and exact identity of all 20 scientific source files with the
completed numerical replay; no source-equivalence waiver is required.

| Check | Evidence |
|---|---|
| Source check and assertions | [check log](v0.10.1-checks/00check.log), [tests](v0.10.1-checks/tests.log) |
| Fresh installation and examples | [install](v0.10.1-checks/install.log), [smoke check](v0.10.1-checks/installed-smoke.log) |
| Archive identity and checksum | [verification](v0.10.1-checks/archive-verification.json), [SHA-256](v0.10.1-checks/source-package-sha256.txt) |
| macOS CI and all three jobs | [macOS log](v0.10.1-checks/macos-CI.log), [CI receipt](v0.10.1-checks/CI-jobs.json) |
| README and historical case claims | [fact check](v0.10.1-checks/fact-check.log) |
