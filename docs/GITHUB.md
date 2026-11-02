# GitHub notes

The repository is `patrickiben/catmech`. Workflow results are recorded in `validation/VALIDATION.md`.

The history has been rewritten twice. Before the first push, one string was removed from a checking script. On 2026-11-02,
before any public release, one person's email address in a registry field of `validation/sponsor-clustering/sponsor-map.csv`
became a placeholder, and every recorded hash of that file and of the files that record it was updated to match. The same pass
reworded an old checking script, one checklist line and one commit message. No result changed. The order and dates of the
commits did not change, and [the commit map](commit-map.txt) gives old and new hashes, since older logs and CI records cite the
old ones. Commit order is local evidence, not third-party proof of timing.

The source includes an MIT license, Patrick Iben as author/maintainer and the account's GitHub no-reply email. Do not
add real participant data or generated analysis bundles from private data.

Local `R CMD check` evidence, workflow results and platform limitations are recorded in `validation/VALIDATION.md`.
A CRAN submission, DOI, independent scientific peer review and
clinical validation are separate future steps. Versions 0.2 and later include prespecified
held-out synthetic evaluations, and later versions add public real-data and real-audit tests, none of which is a clinical
source-data audit. The tests from version 0.5 on have not been reviewed by an independent statistician. The remaining failures must remain visible in any release description.
