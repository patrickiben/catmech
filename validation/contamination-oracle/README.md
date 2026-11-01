# Independent contamination-model oracle

A random-walk Metropolis sampler on the posterior with status and target integrated out, written in `tools/check-contamination-oracle.R` from the model statement, is compared with the package's slice-within-Gibbs sampler. Three synthetic cases (n = 70) cover one stratum, two strata with a selection shift, and a declared noisy reference. Run time 95 seconds.

- Posterior means and SDs: 38/38 within five combined MCSEs; largest |z| 1.96.
- Row-level corruption probabilities and target means: 36/36 within 0.015; largest difference 0.0053.
- Known-wrong likelihoods detected: 3/3; smallest largest-|z| 12.1.
- Data-contract controls: 15/15.

This is implementation evidence for these examples. It says nothing about repeated-sampling coverage, model adequacy or real data.
