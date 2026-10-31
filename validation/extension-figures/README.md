# Measured references and remaining assumptions

A: a new random 25% accurate reference sample supports normal MI under its model when flag masking misses corruption. Shared reference results are shown once, alongside two recording-error scenarios from the same 1000 cohorts. The compound stress combines a site shift, random errors and 20% flag sensitivity; SITEID is included, so it does not isolate a site effect.

B: holding predictor information fixed and excluding categories, PMM improves this heavy-tail example but loses coverage under strong MAR. Its horizontal coverage axis is expanded; the other coverage panels start at zero.

C: more random references improve precision. Width is not a validity measure: experimental reference PMM has narrower intervals but confirmed undercoverage (74.8%, 85.1%, 90.8% at the tested 10%, 25%, 50% fractions). The 10% and 50% rows use 250 cohorts; the 25% row uses 1000. Subsets are nested within common cohorts.

D: normal reference MI fails when validation selection is nonignorable or supplied references are noisy. The three cases use a 25% expected reference fraction. Counts retain the 500-replicate noisy-assay design.

Coverage error bars are 95% Wilson intervals for Monte Carlo uncertainty, not individual scientific confidence intervals. Mean-width bars are mean +/- 1.96 MCSE. Exact MCSEs, counts, labels and every plotted source row are preserved in plotted-data.csv. Dashed vertical lines show nominal 95% coverage.

All results are synthetic, conditional on the stated designs, and use additional measured information for reference methods. Partial-reference PMM is experimental and its v0.3.1 pooled analysis requires explicit sensitivity-only opt-in. The figure does not certify general coverage or identify error without assumptions.

The script accepts completed bundles only and verifies all 12 frozen numerical source/plan files against validation/v0.3.0-source. figure-provenance.json records SHA-256 hashes of inputs, script and generated figure/data files.

Reproduce with `Rscript tools/plot-extension.R validation/new-extension-figures` from the repository root. Existing output directories are never overwritten. See ../EXTENSION_SYNTHESIS.md and ../EXTENSION_PLAN.md for the full evidence and prospective design.

Artifacts: extension-evidence.png, extension-evidence.pdf, plotted-data.csv, figure-provenance.json.
