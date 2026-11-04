# Contributing

Report a reproducible example with the complete specification and package versions.
Use synthetic data. Bugs involving threshold boundaries, trust policies, missingness,
nonestimable predictions, or audit/provenance records should include a regression test.

Run `testthat::test_local()` and `R CMD check --no-manual` after changes. Documentation
is generated from roxygen comments. Benchmark changes require a dated plan and a new
output directory, including failures, paired baselines and Monte Carlo uncertainty.
Do not select a favorable q or data-generating process after inspecting results and
then present it as prespecified.

This initial revision used OpenAI Codex for code drafting, literature search,
adversarial review, test construction, and simulation orchestration. Runtime
validation uses ordinary R code and synthetic data, not an AI model. Human authors
retain responsibility for scientific claims, intended-use decisions and release.
Versions before 0.4 also used OpenAI Codex. Versions 0.4 and later used Anthropic Claude Code for the same kinds of
drafting, search, review and simulation orchestration, and it also wrote most later study plans, decision rules, scoring
and write-ups, at the owner's direction. OpenAI Codex reviewed version 0.7.1 and
implemented the three regression fixes and scientific wording corrections in
version 0.7.2, including focused tests and release verification. It implemented the
`pan`-backed joint clustered model and its validation in version 0.8.0, with separate
AI implementation, mathematical-review and simulation tasks. The write-ups say
which rounds, and no independent statistician has reviewed them.
Commit trailers record Claude Code's part, and no AI system is listed as a package or citation author. This is a software assessment and exploratory
validation record, not external peer review or a completed manuscript submission
certification.
