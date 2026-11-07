# catmech v1.0.0: Multimodal R-Package Educational Companion

Choose the learning modality and depth that fit your workflow:

1. **Concise 12-Lesson Study Guide (`DOCUMENT`):** [Read Markdown](STUDY_GUIDE.md) or [Download PDF](STUDY_GUIDE.pdf). Covers the clinical audit dilemma, two-phase contamination and linear reference MI models, Barnard-Rubin small-sample pooling across linear/GLM/Cox endpoints, synthetic vs real-world audit benchmarks, and regulatory provenance boundaries.
2. **Complete 30-Unit Reading Companion (`DOCUMENT`):** [Read Markdown](COMPANION.md) or [Download PDF](COMPANION.pdf). Covers 100% of [paper/manuscript.Rmd](../paper/manuscript.Rmd) (`lines 1-350`, SHA-256 `8152288b1c2af2a49103543db3c9324582177b1e842365b5c5ab83cdac8c1fd0`) and [paper/manuscript.pdf](../paper/manuscript.pdf) (`15` pages, SHA-256 `1d17480886d4f38c38e7cf5f900ea3c5f1ef48f3c4484fa801a72a83de05e0db`), with worked numerical reasoning, R one-liners, glossary, and a locked 12-item transfer assessment. See [COMPLETE_COVERAGE.md](COMPLETE_COVERAGE.md) for the line-by-line map.
3. **Interactive HTML5 Companion Studio & Browser Simulators (`INTERACTIVE`):**
   - [Interactive Companion Studio & 12-Item Practice Lab](companion_lab.html): Browse all 30 units (`C01`-`C30`), copy executable R snippets, and take the interactive 12-item Transfer & Boundary Practice Quiz with immediate attempt tracking, correctness feedback, answer explanations, and retry support.
   - [catmech Field Guide & 24-Scenario Method Picker](../docs/learn/index.html): Interactive decision tree generating verified R code snippets.
   - [Live Browser Linear-MI Audit Simulator](../docs/learn/simulator.html): Run 140 benchmark studies and custom audit fractions in the browser.
   - [Interactive Record Audit Planner](../docs/audit-planner/index.html): Explore two-phase precision ratios `rho(f)`.
4. **One-Page Visual Summaries (`ONE_PAGER`):** [Download Companion One-Pager PDF](ONE_PAGER.pdf) and the 1-page visual empirical summary [catmech-case-study.pdf](../docs/case-study/catmech-case-study.pdf).
5. **Executable R Lab Suite & Interactive R Console Tutor (`R`):**
   - Run all 16 live formula, workflow, and refusal labs in batch mode: `Rscript education/run_companion_labs.R`
   - Launch the interactive 12-lesson R console tutor: `Rscript -e "source('education/interactive_tutor.R'); catmech_tutor()"` (or `Rscript education/interactive_tutor.R --non-interactive`)

This educational companion explains the statistical definitions, benchmark archives, and R workflows in the frozen catmech manuscript; it does not establish universal interval coverage under unmodeled misspecification, validate per-record error probabilities for individual patient adjudication on diffuse real-world data, or certify FDA 21 CFR Part 11 or ICH E9 regulatory compliance.

To verify all G14 schema-v2 locks, PDF metadata hashes, readability ceilings, and live R labs:

```sh
python3 education/verify_educational_companion.py
```

See [ARTIFACTS.json](ARTIFACTS.json), [g14_companion_manifest.json](g14_companion_manifest.json), and [G14_HUMAN_SIGNOFF.md](G14_HUMAN_SIGNOFF.md) for the complete verification records.
