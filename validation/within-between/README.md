# Data for the test of within-cluster and between-cluster slopes

`snap-pool.csv` comes from `qc_pub_fy2020_per1.csv` in U.S. Department of Agriculture, Food and Nutrition Service, *Fiscal Year 2020 Supplemental Nutrition Assistance Program Quality Control Database*, Ag Data Commons, doi:10.15482/USDA.ADC/1528542, listed as U.S. Public Domain. `tools/within-between-prepare.py` builds it and checks the file's published MD5, d35bf400376f1774109c9f3b713458eb. The 27.2 MB source is not kept here.

It has one row per household review: `state`, `recorded` (`RAWBEN`, hundreds of dollars), `reference` (the benefit the reviewer determined, from `RAWBEN`, `AMTERR` and `STATUS`, hundreds of dollars), `size` (`CERTHHSZ`) and `gross` (`RAWGROSS`, hundreds of dollars).

The file's documentation says state reviewers re-interview households and that cases found ineligible, incomplete reviews and overissuances at least as large as the benefit are removed. So the pool leaves out some of the largest errors. It holds the state reviewer's finding, which is not a gold standard, and it is not a sample of all SNAP households in the sense of a weighted analysis, since the study's cohorts weight states equally and do not use the survey weights.
