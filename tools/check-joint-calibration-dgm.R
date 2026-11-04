#!/usr/bin/env Rscript
# Compare generators only; do not execute either study runner.
script <- sub('^--file=', '', grep('^--file=', commandArgs(), value = TRUE)[1L])
root <- normalizePath(file.path(dirname(script), '..'), mustWork = TRUE)
load_functions <- function(path, required, config_name) {
  e <- new.env(parent = globalenv())
  e[[config_name]] <- list(data_seed_base = 1530000000L, imputation_seed_base = 1540000000L)
  found <- character()
  for (z in parse(path)) if (is.call(z) && identical(z[[1L]], as.name('<-')) &&
    is.symbol(z[[2L]]) && as.character(z[[2L]]) %in% required) {
    eval(z, e); found <- c(found, as.character(z[[2L]]))
  }
  stopifnot(setequal(found, required));e
}
old <- load_functions(file.path(root,'validation/v0.8.0-source/validation/run-joint-cluster.R'),c('rmv','logistic_intercept','make_cohort'),'settings')
new <- load_functions(file.path(root,'validation/run-joint-calibration.R'),c('rmv','logit_intercept','make_cohort'),'config')
cells <- c('gaussian_balanced','gaussian_unbalanced','near_zero_icc','observed_mar','null_slope','student_t5')
for (i in seq_along(cells)) for(rep in c(1L,2L,30L)) {
  a <- old$make_cohort(cells[i],i,rep); b <- new$make_cohort(cells[i],i,rep)
  stopifnot(identical(a$input,b$input),identical(a$truth,b$truth),
    identical(serialize(a$input,NULL,version=2),serialize(b$input,NULL,version=2)),
    identical(serialize(a$truth,NULL,version=2),serialize(b$truth,NULL,version=2)),
    identical(a$meta$input_sha256,b$meta$input_sha256),
    identical(a$meta$data_seed,b$meta$data_seed),identical(a$meta$imputation_seed,b$meta$imputation_seed),
    identical(a$meta$truth,b$meta$target))
  cat(cells[i],rep,'input + private truth byte-identical; seeds/hash/target identical\n')
}
# Algebraic independent Gaussian check: (a,z,v_T,e_T) have diagonal covariance.
# x = a + z, T-.5 = (.7+.4)*a + .7*z + v_T + e_T.
variance <- diag(c(1,1,.5,.5)); x <- c(1,1,0,0); truth <- c(1.1,.7,1,1)
varx <- as.numeric(crossprod(x,variance%*%x))
covxt <- as.numeric(crossprod(x,variance%*%truth))
cat('contextual_x analytic Var(x)=',varx,' Cov(x,T)=',covxt,' marginal slope=',covxt/varx,'\n',sep='')
stopifnot(all.equal(covxt/varx,.9),identical(new$make_cohort('contextual_x',7L,1L)$meta$target,.9))
cat('Within-site structural coefficient .7; missingness uses ordinary random masks.\nNo fitting, coefficient accuracy, or coverage evaluation was performed.\n')
