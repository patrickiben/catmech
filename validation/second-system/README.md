# Second health system data

`raw/Zambia_data.xlsx` and `raw/Zambia_data_dictionary.docx` are unchanged copies from the Harvard Dataverse record for Grover D, Bauhoff S, Friedman J. Using supervised learning to select audit targets in performance-based financing in health: An example from Zambia. *PLOS ONE* 2019. [doi:10.1371/journal.pone.0211262](https://doi.org/10.1371/journal.pone.0211262). Data: [doi:10.7910/DVN/LHUIBO](https://doi.org/10.7910/DVN/LHUIBO), under CC0. They were downloaded on 2026-11-01, and the record's two README files were not read. The authors did not take part in this work.

`zambia-pool.csv` is built from the raw file by `tools/second-system-prepare.R`. It keeps the nine items' first-quarter values in long form with the clinic, district, province and catchment size, keeps rows with both values and a catchment size, and adds `both_zero` and `covariate`, the base-10 logarithm of catchment size. No values were changed. What the test does with it is in `../SECOND_SYSTEM_PLAN.md`.

The study's clinics were selected within districts with probability proportional to catchment size, so these rows are not a simple random sample of Zambian clinics.
