# Pooling guards and frozen-source numerical equivalence

Overall result: **PASS**.

54 supported scalar comparisons use identical completion objects in the frozen v0.3.0 pool and current pool.
They cover all five workflow engines, four complete-data endpoints, three supported scalar analysis forms, and 90%/95% confidence levels.
The comparison requires identical coefficient names and bit-identical estimate, standard error, degrees of freedom, limits, within/between variance, m, and n.
Partial reference PMM and the legacy engine receive explicit sensitivity opt-in; changed interpretation labels are recorded separately.

Guard checks reproduce the old silent-NULL multiresponse/zero-coefficient bugs and verify explicit rejection.
A single-column response and an intercept-only analysis remain supported. Partial-reference PMM is refused by default,
even with zero empirical between-imputation variance. Fully referenced PMM remains ordinary complete-data analysis.
The frozen source is from commit ea5ec6ae4688c5ff8a2f2086cc6f45e69cb5f21c; all current/frozen source hashes were unchanged during the check.

These are engineering and numerical-equivalence checks. They do not establish coverage or rehabilitate the undercoverage of reference PMM.
