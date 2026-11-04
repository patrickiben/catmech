#!/usr/bin/env Rscript
# Fabricated reporting controls only; this file never fits a statistical model.
script <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1L])
root <- normalizePath(file.path(dirname(script), ".."), mustWork = TRUE)
args <- commandArgs(trailingOnly = TRUE)
out <- if (length(args)) args[1L] else tempfile("catmech-marginal-reporting-controls-")
dir.create(out, showWarnings=FALSE)
e <- new.env(parent=globalenv())
wanted <- c('cells','prior_cells','contextual_cells','primary_cells','se_mean','wilson','summary_rows','paired_rows','prior_rows','calibration_screen','availability_screen')
for(z in parse(file.path(root,'validation/run-joint-marginal.R'))) if(is.call(z) && identical(z[[1]],as.name('<-')) && is.symbol(z[[2]]) && as.character(z[[2]]) %in% wanted) eval(z,e)
e$config <- list(prior_reps=30L,primary_reps=30L,bootstrap_seed=1770000000L,bootstrap_reps=2000L)
records <- fits <- list()
for(ci in seq_along(e$cells)) for(rep in 1:30) {
  cell <- e$cells[ci]
  ids <- c('baseline','marginal','oracle',if(cell %in% e$prior_cells)c('baseline_prior_quarter','baseline_prior_four','marginal_prior_quarter','marginal_prior_four','marginal_regression_half','marginal_regression_double'))
  for(id in ids) {
    status <- if((rep==1 && cell=='gaussian_balanced' && startsWith(id,'marginal')) || (cell=='student_t5' && id=='marginal') || (rep==2 && cell=='near_zero_icc' && id=='marginal_regression_half')) 'error' else 'ok'
    fits[[length(fits)+1]] <- data.frame(cell=cell,replicate=rep,fit_id=id,fit_status=status)
    for(analysis in if(id=='oracle')'cv3' else c('direct','cv3')) {
      family <- sub('_(prior|regression).*','',id)
      method <- sub(family,paste0(family,'_',analysis),id,fixed=TRUE)
      base <- if(family=='baseline')-.1 else .005
      shift <- if(endsWith(id,'quarter'))-.03 else if(endsWith(id,'four')).04 else if(endsWith(id,'half'))-.06 else if(endsWith(id,'double')).08 else 0
      if(family=='marginal')shift <- shift/3
      target <- if(cell %in% c('null_slope','contextual_null'))0 else if(cell %in%e$contextual_cells).9 else .7
      estimate <- target + base + shift + .12*sin(rep+ci)
      spread <- .1 + .02*(analysis=='cv3') + .01*abs(shift)
      lower <- estimate-1.96*spread;upper <- estimate+1.96*spread
      row <- data.frame(cell=cell,replicate=rep,method=method,input_sha256=paste(ci,rep),status=status,estimate=estimate,spread=spread,lower=lower,upper=upper,width=upper-lower,df=if(analysis=='cv3')20 else NA_real_,within_variance=if(analysis=='cv3').008 else NA_real_,between_variance=if(analysis=='cv3').003 else NA_real_,error=estimate-target,covered=as.integer(lower<=target&&upper>=target),interval_score=upper-lower+40*max(0,lower-target)+40*max(0,target-upper),mcse_mean=if(analysis=='direct').0005 else NA_real_,mcse_lower=if(analysis=='direct').001 else NA_real_,mcse_upper=if(analysis=='direct').001 else NA_real_)
      if(status!='ok') { row[setdiff(names(row),c('cell','replicate','method','input_sha256','status'))] <- NA_real_;row$covered<-0L }
      records[[length(records)+1]]<-row
    }
  }
}
records<-do.call(rbind,records);fits<-do.call(rbind,fits)
wr<-function(x,name)write.csv(x,file.path(out,name),row.names=FALSE,na='')
wr(records,'records.csv');wr(fits,'fit-records.csv');wr(e$summary_rows(records),'summary.csv');wr(e$paired_rows(records),'paired.csv')
p<-e$prior_rows(records);wr(p$triples,'prior-triples.csv');wr(p$summary,'prior-summary.csv')
c<-e$calibration_screen(records);wr(c$bootstrap,'calibration-bootstrap.csv');wr(c$contextual_pairs,'contextual-paired.csv')
a<-e$availability_screen(fits,records);wr(a$cells,'availability.csv')
jsonlite::write_json(list(contextual=c$contextual,ordinary=c$ordinary,conditional_coverage=c$conditional,bias=c$bias,prior=p$gate,availability=a$gate),file.path(out,'controls.json'),auto_unbox=TRUE,pretty=TRUE,na='null',digits=16)
cat('Synthetic reporting-control rows only; no catmech fits or actual study outcomes used.\n')
