#!/usr/bin/env Rscript
# Selected numerical challenge, not independent coverage evidence.
root <- normalizePath(".")
pkgload::load_all(root, quiet=TRUE)
s <- jsonlite::read_json("validation/joint-marginal-confirmation/settings.json", simplifyVector=TRUE)
e <- new.env(parent=globalenv()); e$config <- s; e$contextual_cells <- s$contextual_cells
for (z in parse("validation/joint-marginal-confirmation-source/validation/run-joint-marginal.R")) {
 if(is.call(z) && identical(z[[1]],as.name("<-")) && is.symbol(z[[2]]) && as.character(z[[2]]) %in% c("rmv","logit_intercept","make_cohort")) eval(z,e)
}
out <- list()
for(case in list(c(2L,12L),c(9L,143L))) {
 p <- e$make_cohort(s$cells[case[1]],case[1],case[2])
 old <- read.csv("validation/joint-marginal-confirmation/fit-records.csv")
 row <- old[old$cell==s$cells[case[1]] & old$replicate==case[2] & old$fit_id=="marginal",]
 stopifnot(nrow(row)==1, row$input_sha256==p$meta$input_sha256)
 fit <- catmech_reference(p$input,"W","gold","x",m=40,seed=p$meta$imputation_seed,use_recorded=TRUE,cluster="site",cluster_model="joint_marginal",joint_control=list(max_iterations=64000L))
 d <- fit$model_diagnostics
 attempts <- jsonlite::fromJSON(row$attempts)
 stopifnot(identical(as.character(jsonlite::toJSON(d$attempts[1:4,],auto_unbox=TRUE)),row$attempts),
 isTRUE(all.equal(d$attempts$max_rhat[4],row$max_rhat,tolerance=1e-12)),
 isTRUE(all.equal(d$attempts$min_bulk_ess[4],row$min_ess_bulk,tolerance=1e-12)),d$convergence_pass)

 out[[length(out)+1L]] <- list(cell=s$cells[case[1]],replicate=case[2],input_sha256=p$meta$input_sha256,iterations=d$iterations_used,attempts=d$attempts,diagnostics=d$summary)
 cat(s$cells[case[1]],case[2],"PASS at",d$iterations_used,"retained draws per chain\n")
}
dir.create("validation/marginal-extension-challenge",showWarnings=FALSE)
jsonlite::write_json(list(scope="selected numerical challenge; no accuracy evaluated",cases=out),"validation/marginal-extension-challenge/results.json",pretty=TRUE,auto_unbox=TRUE,digits=NA)
