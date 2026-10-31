# Coverage comparison figure

The figure compares legacy sensitivity PMM, the Gaussian model ignoring categories, and independently implemented trusted normal MI on the same held-out synthetic cohorts. The clinical-format panels expose remaining failures when contaminated measurements are missed by flags.

The site panel is a compound stress: site 003 has a fixed 0.20 mg/dL shift, 8% of observed records have random additive errors, and flag sensitivity is only 20%. The missed-flags panel has 60% sensitivity. SITEID is included in the analysis, where a constant site shift alone can be absorbed; this comparison does not isolate the effect of site bias.

Points show coverage among returned intervals. Error bars are 95% Wilson Monte Carlo intervals, which remain informative at zero or unit coverage; the accompanying CSV also gives each binomial MCSE and the returned/requested denominator. The dashed line is the nominal 95% target.

The legacy constraint method was not run on clinical-format fixtures because their endpoint convention differs. The figure is descriptive evidence for the specific simulated settings, not a clinical validation or universal coverage claim.

Artifacts: coverage-comparison.png, coverage-comparison.pdf, plotted-data.csv, figure-provenance.json.

Reproduce with `Rscript tools/plot-coverage.R validation/coverage-confirmatory validation/coverage-clinical validation/new-figures` from the repository root. Existing figure directories are never overwritten.
