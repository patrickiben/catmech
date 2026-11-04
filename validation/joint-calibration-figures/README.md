# Joint calibration confirmation figure

Generated only after the confirmation COMPLETE marker. The PDF is vector output; the PNG is 2800 by 2000 pixels.

Panels A and B show coverage among returned intervals, with Wilson Monte Carlo intervals and a 95% reference line. B is a deliberately misspecified contextual-X stress, separate from the original scenarios.

Panel C compares the mean absolute log width ratio across the prespecified prior scales on cohorts where both analyses returned all three intervals. Error bars are normal Monte Carlo intervals for each mean; the plotted-data file records the common denominator. These individual error bars do not replace the paired sensitivity analysis.

Panel D shows model-fit refusals over every request, with Wilson intervals; it is not the count of analysis-stage refusals. Full analysis refusal information remains in the study records.

Each plotted scalar is independently recomputed from saved raw records and checked against its saved summary before output. plotted-data.csv and figure-provenance.json retain numerical and source provenance.

The direct intervals are conditional model-based Bayesian credible intervals. CV3 intervals use cluster jackknife variance and Rubin pooling. Neither the figure nor diagnostic success establishes universal validity.
