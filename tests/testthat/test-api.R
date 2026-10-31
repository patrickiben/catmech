fixture <- function() {
  d <- simulate_catmech_lab(200, seed=99)
  d <- d[setdiff(names(d),c("lab_true","corrupted"))]
  s <- catmech_spec("lab_value","lab_category",c("Low","Medium","High"),
    c(-Inf,5.5,7,Inf), flag="outlier_flag", corruption_probability=1, conflict_policy="reconstruct")
  list(d=d,s=s,p=c("age","sex","inflammation","outcome"))
}

test_that("high-level completions preserve originals, audit and mids values", {
  a <- fixture()
  old <- a$d
  fit <- catmech_impute(a$d,a$s,a$p,m=3,seed=2)
  expect_identical(a$d,old)
  expect_identical(fit$original,old)
  expect_equal(fit$diagnostics$threshold_violations,rep(0,3))
  expect_false(anyNA(fit$diagnostics))
  for (j in 1:3) {
    x <- fit$completed[[j]]
    expect_false(anyNA(x$lab_value))
    expect_identical(as.character(catmech_derive_category(x$lab_value,a$s$thresholds,a$s$levels)),as.character(x$lab_category))
    expect_equal(mice::complete(fit$mids,j)$lab_value,x$lab_value)
    unaffected <- !is.na(old$lab_value) & as.character(old$outlier_flag)=="No"
    expect_equal(x$lab_value[unaffected],old$lab_value[unaffected])
  }
  expect_equal(mice::complete(fit$mids,0)$lab_value,old$lab_value)
  fit2 <- catmech_impute(a$d,a$s,a$p,m=3,seed=2)
  expect_identical(fit$completed,fit2$completed)
  expect_identical(fit$audit,fit2$audit)
})

test_that("seeded functions leave caller random state unchanged", {
  set.seed(999)
  before <- .Random.seed
  a <- fixture()
  expect_identical(.Random.seed,before)
  catmech_impute(a$d,a$s,a$p,m=2,seed=2)
  expect_identical(.Random.seed,before)
})

test_that("conflicts and unsupported missingness cannot pass silently", {
  a <- fixture()
  s <- a$s
  s$conflict_policy <- "error"
  expect_error(catmech_impute(a$d,s,a$p,m=2),"conflict")
  a$d$age[1] <- NA
  expect_error(catmech_impute(a$d,a$s,a$p,m=2),"complete")
  a <- fixture()
  a$d$lab_category[which(a$d$outlier_flag=="Yes")[1]] <- NA
  expect_error(catmech_impute(a$d,a$s,a$p,m=2),"flagged")
  a <- fixture()
  a$d$lab_category <- as.character(a$d$lab_category)
  a$d$lab_category[1] <- "unexpected"
  expect_error(catmech_impute(a$d,a$s,a$p,m=2),"Unknown category")
  a <- fixture()
  expect_error(catmech_impute(a$d,a$s,c(a$p,"lab_category"),m=2),"exclude")
})

test_that("only trusted numeric observations can supply a missing category", {
  a <- fixture()
  i <- which(!is.na(a$d$lab_value) & a$d$outlier_flag=="No")[1]
  expected <- as.character(a$d$lab_category[i])
  a$d$lab_category[i] <- NA
  fit <- catmech_impute(a$d,a$s,a$p,m=2)
  expect_true(is.na(fit$original$lab_category[i]))
  expect_equal(as.character(fit$completed[[1]]$lab_category[i]),expected)
  expect_equal(as.character(mice::complete(fit$mids,1)$lab_category[i]),expected)
  expect_equal(fit$manifest$category_derived_rows,i)
})

test_that("output bundles preserve signed infinity and refuse overwrite", {
  a <- fixture()
  fit <- catmech_impute(a$d,a$s,a$p,m=2)
  directory <- tempfile()
  write_catmech(fit,directory)
  expect_identical(readRDS(file.path(directory,"result.rds"))$original,a$d)
  manifest <- jsonlite::read_json(file.path(directory,"manifest.json"))
  expect_equal(manifest$assumptions$thresholds[[1]],"-Inf")
  expect_equal(manifest$assumptions$thresholds[[4]],"Inf")
  canonical <- paste(readLines(file.path(directory,"assumptions.canonical.json")),collapse="\n")
  expect_equal(digest::digest(canonical,algo="sha256",serialize=FALSE),fit$manifest$assumption_hash)
  expect_error(write_catmech(fit,directory),"new or empty")
})

test_that("positive unit transformations preserve donor selections", {
  a <- fixture()
  base <- catmech_impute(a$d,a$s,a$p,m=2,seed=12)
  # The target and thresholds have the same units; transform both together.
  d <- a$d
  d$lab_value <- 100*d$lab_value + 50
  s <- a$s
  s$thresholds <- 100*s$thresholds + 50
  scaled <- catmech_impute(d,s,a$p,m=2,seed=12)
  expect_equal(scaled$completed[[1]]$lab_value,100*base$completed[[1]]$lab_value+50,tolerance=1e-9)
  expect_equal(scaled$audit$donor_row,base$audit$donor_row)
  # Known-wrong control: keep the old thresholds after changing target units.
  s$thresholds <- a$s$thresholds
  s$conflict_policy <- "error"
  expect_error(catmech_impute(d,s,a$p,m=2,seed=12),"conflict")
})

test_that("exported low-level method works inside mice with immutable originals", {
  a <- fixture()
  d <- a$d[c(a$p,"lab_value")]
  original <- d$lab_value
  where <- is.na(d)
  suspect <- a$d$outlier_flag == "Yes"
  where[,"lab_value"] <- where[,"lab_value"] | suspect
  method <- stats::setNames(rep("",ncol(d)),names(d))
  method["lab_value"] <- "catmech_threshold_corruption"
  pred <- matrix(0,ncol(d),ncol(d),dimnames=list(names(d),names(d)))
  pred["lab_value",a$p] <- 1
  blots <- list(lab_value=list(category=a$d$lab_category,
    bounds=make_interval_bounds(a$s$thresholds,a$s$levels),
    flag=a$d$outlier_flag, corruption_probability=1, original_y=original,
    trusted=!is.na(original) & !suspect))
  imp <- mice::mice(d,m=2,maxit=2,where=where,method=method,
    predictorMatrix=pred,blots=blots,seed=22,printFlag=FALSE)
  for (j in 1:2) {
    v <- mice::complete(imp,j)$lab_value
    expect_identical(as.character(catmech_derive_category(v,a$s$thresholds,a$s$levels)),as.character(a$d$lab_category))
  }
})
