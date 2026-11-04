# Namespace qualification after completed confirmation

The source package check found one unqualified `ave()` call in
`catmech_marginal_cluster_prepare()`. With the standard `stats` package attached,
that name resolves to `stats::ave`; with it detached, preparation fails. The
release qualifies that one call as `stats::ave`. No statistical specification,
native code, diagnostic threshold, interval rule or study result changes.

Confirmation completed first, from scientific source commit `295fea8`. Its
exact 35-file source tree was archived at `50b9ee6`, and the complete raw results,
including both convergence refusals, were committed at `83249c1` before the
qualification was applied. No confirmation fit was rerun or replaced.

The [equivalence receipt](release-equivalence.json) verifies all archived source
hashes and permits exactly this one textual substitution in the current source.
It regenerated all 3,600 public confirmation inputs from their frozen seeds,
matched their recorded hashes, and found every full prepared model object and
its serialized bytes identical. It also verified the old namespace failure and
corrected behavior with `stats` detached. The original name resolved to the
same `stats::ave` implementation when comparing preparations. This is an
engineering equivalence audit; it performs no fitting or accuracy evaluation.

The release source therefore differs from the confirmation source by this
explicit qualification. `tools/verify-release-archive.py` retains a strict
unchanged-source default and accepts this correction only with the completed
study, frozen sources and full equivalence receipt. The receipt cannot authorize
unrelated changes. Separate synthetic controls accepted the two intended cases
and rejected all seventeen deliberate corruptions. Final installed-package and
CRAN-style checks are recorded alongside this note.
