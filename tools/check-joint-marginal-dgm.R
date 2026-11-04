#!/usr/bin/env Rscript
# Generator and estimand check only: never execute a runner or fit a model.
script <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1L])
root <- normalizePath(file.path(dirname(script), ".."), mustWork = TRUE)
load_generator <- function(path) {
  e <- new.env(parent = globalenv())
  e$config <- list(data_seed_base = 1670000000L, imputation_seed_base = 1680000000L)
  e$contextual_cells <- c("contextual_x", "contextual_unbalanced", "contextual_null")
  wanted <- c("rmv", "logit_intercept", "make_cohort"); found <- character()
  for (z in parse(path)) if (is.call(z) && identical(z[[1L]], as.name("<-")) &&
      is.symbol(z[[2L]]) && as.character(z[[2L]]) %in% wanted) {
    eval(z, e); found <- c(found, as.character(z[[2L]]))
  }
  stopifnot(setequal(wanted, found)); e
}
old <- load_generator(file.path(root, "validation/v0.9.0-source/validation/run-joint-calibration.R"))
new <- load_generator(file.path(root, "validation/run-joint-marginal.R"))
cells <- c("gaussian_balanced", "gaussian_unbalanced", "near_zero_icc", "observed_mar", "null_slope", "student_t5", "contextual_x")
for (i in seq_along(cells)) for (r in c(1L, 2L, 30L)) {
  a <- old$make_cohort(cells[i], i, r); b <- new$make_cohort(cells[i], i, r)
  stopifnot(identical(a, b), identical(serialize(a, NULL, version = 2), serialize(b, NULL, version = 2)))
  cat(cells[i], r, "public inputs, hidden truth, metadata, seeds and target byte-identical\n")
}
# Independent linear transformation of latent (a,z,v_T,v_W,e_T,e_W).
latent <- matrix(0, 6L, 6L); diag(latent)[1:2] <- 1
latent[3:4, 3:4] <- matrix(c(.5, .45, .45, .6), 2L)
latent[5:6, 5:6] <- matrix(c(.5, .55, .55, .8), 2L)
for (i in 7:9) {
  name <- c("contextual_x", "contextual_unbalanced", "contextual_null")[i - 6L]
  structural <- if (name == "contextual_null") -.2 else .7
  L <- rbind(T = c(structural + .4, structural, 1, 0, 1, 0),
             W = c(1.05, 1.05, 0, 1, 0, 1), X = c(1, 1, 0, 0, 0, 0))
  V <- L %*% latent %*% t(L)
  slope <- V["X", "T"] / V["X", "X"]
  p <- new$make_cohort(name, i, 1L)
  stopifnot(isTRUE(all.equal(slope, p$meta$target)), p$meta$n == 240L,
    p$meta$sites == 24L, p$meta$gold_n == 60L,
    p$meta$missing_recorded_n == if (name == "contextual_unbalanced") 84L else 72L,
    all(is.na(p$input$gold[is.na(p$input$gold)])))
  cat(name, ": independent Var(X)=", V["X", "X"], " Cov(X,T)=", V["X", "T"],
      " marginal slope=", slope, " structural slope=", structural, "\n", sep = "")
}
cat("21 same-seed v0.9 cohort objects agree byte-for-byte; 3 contextual targets agree with independent covariance algebra.\n")
cat("No fitted coefficient estimates or coverage were inspected.\n")
