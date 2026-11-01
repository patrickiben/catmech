# Prior-art and source checklist for audit targeting

Working notes from 2026-11-01, not a literature review. Nothing here supports a novelty claim. Verify each item at its source before citing it.

## Where the idea stands

The [go/no-go simulation](../validation/merged-audit-pilot/README.md) returned NO-GO for audits that follow the direction of spillover. A gap in the literature would not change that. That test used independent errors whose effects spread through exposure. Errors that spread themselves, through shared monitors, clerks, labs or systems, are a different idea that remains untested. Injecting such errors into a simulation and then recovering them would be circular. A real test needs record-level audit findings linked to shared-resource metadata, compared against random, severity-based, site-level and adaptive-cluster targeting, with a symmetrized comparison if direction is claimed.

**2026-11-01 result.** A [sponsor clustering screen](../validation/SPONSOR_CLUSTERING_SYNTHESIS.md) on the frozen ClinicalTrials.gov pool returned WEAK under a [plan fixed in advance](../validation/SPONSOR_CLUSTERING_PLAN.md). Differences between registered enrollment and participant-flow Started clustered within sponsors beyond chance (observed-to-null ratio 2.01, p = 0.0005), and a sponsor's revealed trials helped predict its other trials (AUC gain 0.088, 95% interval 0.049 to 0.126). Adaptive auditing did not beat risk ranking (yield ratio 0.97, 0.85 to 1.05, against a required 1.2). This shows clustering in two sponsor-entered registry fields only. It does not show errors spreading, because shared definitions, reporting habits, CROs, SOPs or training can produce the same pattern, and it does not cover clinical source-data errors. The test described above still needs record-level audit data.

## Sources

| Item | Status | Before use |
|---|---|---|
| Andersen et al. 2015: 0.45% overall error rate in three phase 3 trials | Found in a search summary | Confirm the denominator and definition |
| Giganti, Shepherd et al. 2019: 17.1% overall discrepancy rate across audited variables, n = 250, HIV sites | Found in a search summary | Confirm the definition. It is an overall rate, not a peak |
| Comparing those two rates | Not on one scale | Different settings and definitions show variation only |
| Uganda MIS 2018-19: paired RDT and microscopy, "97% of 7,311 eligible children" | Found in a search summary of the DHS report | Read the report. Check eligibility, pairing and variable coding |
| Paired-test analysis of that survey, false positives tied to regional transmission intensity | Found in a search summary | Confirm author, year, DOI and sample size. Do not use it for error contagion |
| Chao et al., arXiv 2310.02151: causal spillover under network membership misclassification, corrected with validation data | Abstract read | Close prior art. Read in full before claiming any gap |
| Shepherd et al. 2023, Biometrics: influence-function-based multiwave validation sampling | Found in a search summary | Read in full |
| Shepherd and Yu 2011: audit-based correction in linear regression | From memory only | Find and read |
| Venet et al. 2012: central statistical monitoring | Found in a search summary | Read in full |
| Thompson 1990: adaptive cluster sampling | From memory only | Find and read |
| Glickman et al. 2009 as a monitoring source | Unsafe. From memory, it is a NEJM paper on globalization of clinical research | Do not cite for monitoring unless read and confirmed |
| Houston et al. 2018: reported figures of 1.3%, about 40% and 50% conflict across notes | Unsettled | Do not use until checked at the source |
| The 0.1-1.5% and 3-10% error-rate bands in earlier notes | No source found | Removed |

## Questions for data managers and monitors

1. Do audit findings exist per record or field in a machine-readable form, or only as site-level summaries?
2. How are records chosen for audit: randomly, by risk, by extreme values, by site, by clerk or by another rule?
3. Can audited records be linked to monitor, site, clerk, lab, assay batch, EDC module or protocol version?
4. Only if direction is pursued: what is the lag between entry, monitoring, audit finding and correction?

Use public sources only.
