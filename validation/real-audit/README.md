# Real-audit data

Two public source files, kept unchanged in `raw/`, with their hashes in `fetch-manifest.json`. Neither holds patient-level values.

**Rwanda facility recount.** `raw/pone.0235823.s001.dta` is S1 Dataset from Nshimyiryo A, Kirk CM, Sauer SM, Ntawuyirusha E, Muhire A, Sayinzoga F, Hedt-Gauthier B. Health management information system (HMIS) data verification: A case study in four districts in Rwanda. *PLOS ONE* 2020. [doi:10.1371/journal.pone.0235823](https://doi.org/10.1371/journal.pone.0235823). The article is distributed under the Creative Commons Attribution License (CC BY).

`rwanda-pool.csv` is that file converted to CSV, with value labels dropped and one column added, `facility_type`, decoded from `facility_cat`. No values were changed. The authors did not take part in this work and have not reviewed it.

**Effect-size audit.** `raw/codebook-primary-studies-complete.csv` comes from the OSF project for Maassen E, van Assen MALM, Nuijten MB, Olsson-Collentine A, Wicherts JM. Reproducibility of individual effect sizes in meta-analyses in psychology. *PLOS ONE* 2020. [doi:10.1371/journal.pone.0233107](https://doi.org/10.1371/journal.pone.0233107). [OSF project 7nsmd](https://osf.io/7nsmd/), under CC0 1.0 Universal.

`effect-size-pool.csv` selects nine columns from that file, renames some, and derives `family` from `efftype`. No values were changed.

`tools/real-audit-prepare.R` rebuilds both pools from `raw/`. What the tests do with them is in `../REAL_AUDIT_PLAN.md` and `../REAL_AUDIT_SYNTHESIS.md`.
