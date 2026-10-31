# GitHub preparation

This directory is prepared as a local repository. Creating or publishing a remote
was outside the request to make the project ready for upload. The release archive
contains the package; the repository archive also includes assessment and validation.

After reviewing the README, assessment and limitations, a repository can be created
under `patrickiben/catmech`. Check name availability first; no remote availability
is asserted here. For example, from this directory with authenticated GitHub CLI:

```sh
gh repo create patrickiben/catmech --private --source=. --remote=origin --push
```

This example deliberately creates a private remote for an initial review. Choose
visibility explicitly when uploading. The source includes an MIT license, Patrick
Iben as author/maintainer and the account's GitHub no-reply email. Review those
metadata before public distribution. Do not add real participant data or generated
analysis bundles from private data.

The GitHub Actions workflow is included but cannot be reported as passing until it
runs remotely. Local `R CMD check` evidence and platform limitations are recorded
in `validation/VALIDATION.md`. A CRAN submission, DOI, scientific peer review and
confirmatory calibration study are separate future steps.
