# Numerical extension revalidation, v0.10.1

Frozen before the full replay on 2026-11-04. The v0.10.0 study, failed availability
screen and source archives remain unchanged. Both original refusals were bulk
ESS misses at32k, not evidence for changing priors or interval rules.

Change only the marginal engine default maximum from32k to64k retained draws per
chain. Keep4chains,1000burn-in,4000initial draws, exact seeded doubling/restarts,
R-hat<=1.05 and bulk/tail ESS>=400 for every existing monitor. Stop at the first
passing attempt. The likelihood, priors, starting states, diagnostics and
prediction law remain fixed. Existing conditional joint default stays32k.

First replay the two selected failures with original inputs/seeds and64k cap.
Verify their first four attempts match the recorded32k attempts, and export
only numerical diagnostics. This is selected engineering evidence, not coverage.

Then rerun ALL3600 original confirmation cohorts,12000 MCMC fits and27600
analyses, including all800 prior sets, with the same frozen seed mapping and
bootstrap seed. The separate run-marginal-extension.R and verifier are copied
from the original, changing candidate cap and provenance paths only. All
original prespecified screens and targets apply. This is a dependent numerical
revalidation of the existing confirmation grid, NOT a new independent scientific
confirmation. Do not replace original evidence or retry any new failures.

Record source hashes before/after; no source changes during replay. Independently
recompute all summaries, diagnostics, bootstrap contrasts and gates, including
deliberate corruptions. Compare old/new successful-fit estimates, intervals,
spread and diagnostics with tight numeric tolerance, excluding timing, source
fingerprints and control metadata. Any discrepancy needs investigation.

The method/DGM/prior ranges/targets and limitations remain as specified in the
original JOINT_MARGINAL_PLAN.md. No new scientific superiority claim arises
from using more draws only after a failed numerical screen.
