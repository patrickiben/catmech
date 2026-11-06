# Building and checking the source package

Version 1.0.0 includes an original native sampler for the optional
`joint_marginal` engine. Building from source requires a C++17 compiler and the
R build toolchain for the installed R version, including its BLAS/LAPACK and
Fortran runtime link dependencies. `Rcpp` and `RcppArmadillo` are declared
dependencies; no separately installed Stan toolchain is used. The registered
native routines are built as part of the package installation.

Use the platform setup described in the official
[R Installation and Administration manual](https://cran.r-project.org/doc/manuals/r-release/R-admin.html).
The package uses the compiler and library settings supplied by R; it contains
no machine-specific compiler paths. A missing runtime library during linking
indicates an incomplete or mismatched local R toolchain, which should be
corrected in that environment.

From a repository checkout with the package dependencies installed:

```sh
R CMD build .
R CMD INSTALL --library=/path/to/isolated/R-library catmech_1.0.0.tar.gz
Rscript tools/smoke-installed.R /path/to/isolated/R-library
R CMD check --as-cran --no-manual catmech_1.0.0.tar.gz
```

Create the isolated library directory before installing. The installed smoke
check runs the bundled examples and verifies namespace-only calls, registered
native sampling, summaries, `mice` adapters and exported JSON/RDS records.
It requires the repository checkout because `tools/` and research validation
artifacts are deliberately excluded from the R source archive.

Validation runs also require the native library built from the same frozen
sources as their R code. Sourcing the R files alone does not supply the native
sampler. Source hashes and numerical checks are necessary alongside a successful
build; compilation does not establish statistical validity.
