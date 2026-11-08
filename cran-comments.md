## Test environments

* macOS (aarch64-apple-darwin23.6.0, Apple Silicon), R 4.6.0 (local)
* macOS 15 (aarch64-apple-darwin23), R 4.6.2 (GitHub Actions)
* Ubuntu 24.04 (x86_64-pc-linux-gnu), R-devel (2026-11-07 r90850), R 4.6.2 and
  R 4.5.3 (GitHub Actions)
* Windows Server (x86_64-w64-mingw32), R-devel (2026-11-07 r90850 ucrt) and
  R 4.6.2 ucrt (GitHub Actions)

## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new submission.
* Words in DESCRIPTION that may be reported as possibly misspelled are the
  surnames of cited authors (Barnard, Neyman, Rubin) or correctly spelled terms
  (auditable, misspecification, unvalidated).

## Notes for the reviewer

* Seeds: ten exported functions take a `seed` argument with a fixed default
  (1, 42 or 123), so that an audit record can be reproduced without extra
  arguments. Each saves the user's random-number state on entry and restores it
  on exit, so calling it leaves the user's random stream unchanged. An internal
  Monte Carlo fallback for Fisher's exact test uses a fixed seed under the same
  save and restore.
* Global environment: the package writes to the global environment only to set
  and restore `.Random.seed`. The one `<<-`, in `R/audit-plan.R`, assigns to a
  local variable of the enclosing function.
* Parallel work: the experimental joint cluster sampler starts at most two
  worker processes at a time.

## Reverse dependencies

This is a new submission, so there are no reverse dependencies.
