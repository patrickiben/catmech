#!/usr/bin/env Rscript
# Evidence plots consume completed confirmation only, never partial or pilot data.
args <- commandArgs(trailingOnly = TRUE)
script <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1L])
root <- normalizePath(file.path(dirname(script), ".."), mustWork = TRUE)
input <- if (length(args)) args[1L] else file.path(root, "validation/joint-marginal-confirmation")
output <- if (length(args) > 1L) args[2L] else file.path(root, "validation/joint-marginal-figures")
if (!file.exists(file.path(input, "COMPLETE.json"))) stop("Confirmation is incomplete; no outcome files were read.")
complete <- jsonlite::read_json(file.path(input, "COMPLETE.json"), simplifyVector = TRUE)
if (!identical(complete$status, "complete") || !identical(complete$mode, "confirmation"))
  stop("Only a completed confirmation stage can produce this evidence figure.")
files <- c("COMPLETE.json", "settings.json", "records.csv", "fit-records.csv", "summary.csv",
  "prior-triples.csv", "prior-summary.csv", "availability.csv", "adoption.json",
  "source-hashes-before.csv", "source-hashes-after.csv")
if (!all(file.exists(file.path(input, files)))) stop("Required evidence is missing.")
read_csv <- function(name) utils::read.csv(file.path(input, name), stringsAsFactors = FALSE)
config <- jsonlite::read_json(file.path(input, "settings.json"), simplifyVector = TRUE)
records <- read_csv("records.csv"); fits <- read_csv("fit-records.csv")
summary <- read_csv("summary.csv"); triples <- read_csv("prior-triples.csv")
prior_summary <- read_csv("prior-summary.csv"); availability <- read_csv("availability.csv")
cells <- c("gaussian_balanced", "gaussian_unbalanced", "near_zero_icc", "observed_mar", "null_slope", "student_t5",
  "contextual_x", "contextual_unbalanced", "contextual_null")
prior_cells <- c("gaussian_balanced", "near_zero_icc", "observed_mar", "contextual_x")
methods <- c("baseline_direct", "marginal_direct")
stopifnot(identical(config$mode, "confirmation"), isTRUE(config$accuracy_exported),
  identical(config$cells, cells), identical(config$prior_cells, prior_cells),
  identical(read_csv("source-hashes-before.csv"), read_csv("source-hashes-after.csv")),
  nrow(records) == complete$analysis_requests, nrow(fits) == complete$fitted_workflow_requests,
  !anyDuplicated(records[c("cell", "replicate", "method")]),
  !anyDuplicated(fits[c("cell", "replicate", "fit_id")]),
  all(records$source_fingerprint == config$source_fingerprint),
  all(fits$source_fingerprint == config$source_fingerprint),
  identical(config$source_fingerprint, complete$source_fingerprint))
check <- function(actual, expected, label) {
  actual <- as.numeric(actual); expected <- as.numeric(expected)
  good <- length(actual) == length(expected) && identical(is.na(actual), is.na(expected))
  if (good) {
    observed <- !is.na(expected)
    good <- all(is.finite(actual[observed]) & is.finite(expected[observed]) &
      abs(actual[observed] - expected[observed]) <= 1e-12 + 1e-11 * abs(expected[observed]))
  }
  if (!good) stop("Plotted-data check failed: ", label)
}
wilson <- function(k, n) {
  if (!n) return(c(NA_real_, NA_real_))
  z <- stats::qnorm(.975); p <- k/n; den <- 1 + z^2/n
  center <- (p + z^2/(2*n))/den
  half <- z*sqrt(p*(1-p)/n + z^2/(4*n^2))/den
  c(max(0, center-half), min(1, center+half))
}
se_mean <- function(x) if (length(x) > 1L) stats::sd(x)/sqrt(length(x)) else NA_real_
mean_or_na <- function(x) if (length(x)) mean(x) else NA_real_
plot_data <- list()
add <- function(panel, cell, method, estimate, interval, n, events = NA_integer_, uncertainty) {
  plot_data[[length(plot_data)+1L]] <<- data.frame(panel=panel, cell=cell, method=method,
    estimate=estimate, lower=interval[1L], upper=interval[2L], n=n, events=events,
    uncertainty=uncertainty, stringsAsFactors=FALSE)
}
for (cell in cells) for (method in methods) {
  r <- records[records$cell == cell & records$method == method, ]
  stopifnot(nrow(r) == config$primary_reps, setequal(r$replicate, seq_len(config$primary_reps)),
    all(r$status %in% c("ok", "error")), all(r$covered %in% c(0,1)), all(r$covered[r$status != "ok"] == 0))
  expected_target <- if (cell %in% c("null_slope", "contextual_null")) 0 else
    if (cell %in% c("contextual_x", "contextual_unbalanced")) .9 else .7
  check(r$target, rep(expected_target, nrow(r)), paste(cell, method, "target"))
  good <- r[r$status == "ok", ]; n <- nrow(good); k <- sum(r$covered)
  coverage <- if (n) k/n else NA_real_; ci <- wilson(k,n)
  saved <- summary[summary$cell == cell & summary$method == method, ]; stopifnot(nrow(saved) == 1L)
  check(unlist(saved[c("requests", "intervals", "failures", "coverage_conditional", "coverage_conditional_lower", "coverage_conditional_upper")]),
    c(nrow(r),n,nrow(r)-n,coverage,ci), paste(cell,method,"coverage"))
  add("A",cell,method,coverage,ci,n,k,"95% Wilson Monte Carlo interval; returned intervals")
  error <- good$estimate-expected_target
  check(good$error,error,paste(cell,method,"error"))
  bias <- mean_or_na(error); bias_se <- se_mean(error)
  check(unlist(saved[c("bias","bias_mcse")]),c(bias,bias_se),paste(cell,method,"bias"))
  add("B",cell,method,bias,bias+c(-1,1)*stats::qnorm(.975)*bias_se,n,
    uncertainty="95% normal Monte Carlo interval for bias; returned intervals")
}
prior_noise <- list()
for (cell in prior_cells) {
  saved_triples <- triples[triples$cell == cell & triples$analysis == "direct", ]
  stopifnot(nrow(saved_triples) == config$prior_reps, !anyDuplicated(saved_triples$replicate),
    setequal(saved_triples$replicate,seq_len(config$prior_reps)))
  sensitivity <- matrix(NA_real_,config$prior_reps,2L,dimnames=list(NULL,methods))
  mean_mcse <- matrix(NA_real_,config$prior_reps,2L,dimnames=list(NULL,methods))
  for (replicate in seq_len(config$prior_reps)) {
    r <- records[records$cell == cell & records$replicate == replicate, ]
    saved <- saved_triples[saved_triples$replicate == replicate, ]
    for (method in methods) {
      suffixes <- c("","_prior_quarter","_prior_four",if(method=="marginal_direct")c("_regression_half","_regression_double"))
      selected <- r[match(paste0(method,suffixes),r$method), ]
      stopifnot(nrow(selected)==length(suffixes),!anyNA(selected$method))
      good <- all(selected$status=="ok"); family <- sub("_direct$","",method)
      stopifnot(identical(good,saved[[paste0(family,"_complete")]]))
      value <- NA_real_
      if(good) {
        hw_shift <- mean(abs(selected$estimate[2:3]-selected$estimate[1L]))
        value <- if(method=="marginal_direct")max(hw_shift,mean(abs(selected$estimate[4:5]-selected$estimate[1L])))else hw_shift
      }
      check(saved[[paste0(family,"_coefficient_sensitivity")]],value,paste(cell,replicate,method,"coefficient S"))
      sensitivity[replicate,method] <- value
      if(good)mean_mcse[replicate,method] <- mean(selected$mcse_mean)
    }
    check(saved$coefficient_difference,sensitivity[replicate,2L]-sensitivity[replicate,1L],paste(cell,replicate,"paired coefficient S"))
  }
  common <- apply(is.finite(sensitivity),1L,all)
  values <- sensitivity[common,,drop=FALSE]
  saved <- prior_summary[prior_summary$cell == cell & prior_summary$analysis == "direct", ];stopifnot(nrow(saved)==1L)
  delta <- values[,2L]-values[,1L]
  check(unlist(saved[c("complete_pairs","mean_baseline_coefficient_sensitivity","mean_marginal_coefficient_sensitivity","coefficient_difference","coefficient_mcse")]),
    c(nrow(values),mean_or_na(values[,1L]),mean_or_na(values[,2L]),mean_or_na(delta),se_mean(delta)),paste(cell,"prior summary"))
  for(method in methods) {
    value <- mean_or_na(values[,method]); se <- se_mean(values[,method])
    add("C",cell,method,value,value+c(-1,1)*stats::qnorm(.975)*se,nrow(values),
      uncertainty="95% normal Monte Carlo interval for mean absolute coefficient shift; common complete old-three/new-five sets")
    prior_noise[[length(prior_noise)+1L]] <- data.frame(cell=cell,method=method,common_sets=nrow(values),
      mean_absolute_coefficient_shift=value,mean_single_summary_mcse=mean_or_na(mean_mcse[common,method]),
      note="Finite-chain posterior means add simulation noise to absolute shifts; MCSE is not subtracted.")
  }
}
for(cell in cells) for(family in c("baseline","marginal")) {
  r <- fits[fits$cell == cell & fits$fit_id == family, ]
  stopifnot(nrow(r)==config$primary_reps,setequal(r$replicate,seq_len(config$primary_reps)),all(r$fit_status %in% c("ok","error")))
  n<-nrow(r); k<-sum(r$fit_status!="ok")
  saved<-availability[availability$cell==cell, ];stopifnot(nrow(saved)==1L)
  check(saved$requests,n,paste(cell,"fit requests"))
  check(saved[[paste0(family,"_failures")]],k,paste(cell,family,"fit refusals"))
  check(saved[[paste0(family,"_failure_rate")]],k/n,paste(cell,family,"fit refusal rate"))
  add("D",cell,paste0(family,"_direct"),k/n,wilson(k,n),n,k,"95% Wilson Monte Carlo interval; all requested fits")
}
plot_data <- do.call(rbind,plot_data)
prior_noise <- do.call(rbind,prior_noise)
labels <- c(gaussian_balanced="Gaussian balanced",gaussian_unbalanced="Gaussian unbalanced",near_zero_icc="Near-zero ICC",
  observed_mar="Observed-data MAR",null_slope="Null slope",student_t5="Student-t5",
  contextual_x="Contextual X",contextual_unbalanced="Contextual unbalanced",contextual_null="Contextual null")
colors <- c(baseline_direct="#737C8C",marginal_direct="#007E87")
offsets <- c(baseline_direct=.12,marginal_direct=-.12)
points <- function(z,y,multiplier=1) {
  color <- unname(colors[z$method]);good<-is.finite(z$estimate);ci<-good&is.finite(z$lower)&is.finite(z$upper)
  graphics::segments(multiplier*z$lower[ci],y[ci],multiplier*z$upper[ci],y[ci],col=color[ci],lwd=1.5)
  for(end in c("lower","upper"))graphics::segments(multiplier*z[[end]][ci],y[ci]-.04,multiplier*z[[end]][ci],y[ci]+.04,col=color[ci])
  graphics::points(multiplier*z$estimate[good],y[good],pch=16,cex=1.1,col=color[good])
}
axes <- function(xlim,cells,title,subtitle,xlab,reference=NULL) {
  n<-length(cells)
  graphics::plot(NA_real_,NA_real_,type="n",xlim=xlim,ylim=c(.45,n+.55),axes=FALSE,xlab=xlab,ylab="",xaxs="i",yaxs="i")
  ticks<-pretty(xlim);ticks<-ticks[ticks>=xlim[1L]&ticks<=xlim[2L]]
  graphics::abline(v=ticks,col="#E8EBEF",lwd=.8)
  if(!is.null(reference))graphics::abline(v=reference,col="#536171",lty=2,lwd=1.2)
  if(n==9L)graphics::abline(h=3.5,col="#CCD2DA",lty=3)
  graphics::axis(1,at=ticks,col="#9BA4B0",col.axis="#273142",cex.axis=.85)
  graphics::axis(2,at=n:1,labels=labels[cells],las=1,tick=FALSE,col.axis="#273142",cex.axis=.85)
  graphics::title(main=title,adj=0,line=2.0,cex.main=1.06,col.main="#172B45")
  graphics::mtext(subtitle,side=3,adj=0,line=.45,cex=.77,col="#526073")
}
range_finite <- function(z,zero=FALSE,padding=.08) {
  r<-range(c(z$lower[is.finite(z$lower)],z$upper[is.finite(z$upper)],if(zero)0),na.rm=TRUE)
  if(!all(is.finite(r)))return(c(-.01,.01))
  span<-max(diff(r),.001);r+c(-1,1)*padding*span
}
render <- function() {
  graphics::par(mfrow=c(2,2),mar=c(4.4,10,3.8,1.5),oma=c(4.2,0,5.2,0),family="sans",cex=1,col.lab="#273142",mgp=c(2.7,.65,0))
  z<-plot_data[plot_data$panel=="A", ]
  coverage_lower<-max(0,min(c(.9,z$lower[is.finite(z$lower)]))*100-2)
  axes(c(coverage_lower,101),cells,"A  Coverage of the population marginal slope",
    sprintf("%d cohorts per scenario; denominator = returned intervals",config$primary_reps),"Coverage (%)",95)
  points(z,10-match(z$cell,cells)+offsets[z$method],100)
  z<-plot_data[plot_data$panel=="B", ]
  axes(range(pretty(range_finite(z,TRUE))),cells,"B  Bias for the same population target",
    "Contextual targets: 0.9, 0.9 and 0; no within-site target substitution","Mean estimate minus target",0)
  points(z,10-match(z$cell,cells)+offsets[z$method])
  z<-plot_data[plot_data$panel=="C", ]
  limits<-range(pretty(range_finite(z,TRUE)));limits[1L]<-min(0,min(z$lower[is.finite(z$lower)]))
  axes(limits,prior_cells,"C  Sensitivity of coefficient estimates",
    "Old IW axis; new maximum of HW and regression-scale axes","Mean sensitivity, strongest declared axis (original units)",0)
  points(z,5-match(z$cell,prior_cells)+offsets[z$method])
  z<-plot_data[plot_data$panel=="D", ]
  upper<-max(c(1,z$upper[is.finite(z$upper)]*100))*1.08
  axes(c(-.025*upper,upper),cells,"D  Numerical fit refusals",
    sprintf("All %d default-prior fit requests per scenario",config$primary_reps),"Refused fits (%)",0)
  points(z,10-match(z$cell,cells)+offsets[z$method],100)
  graphics::mtext(if(startsWith(config$source_fingerprint,"SYNTHETIC"))
    "Layout check: fabricated numbers, not study evidence" else
    "Joint marginal model: fresh confirmation evidence",side=3,outer=TRUE,line=3.3,font=2,cex=1.35,col="#172B45")
  graphics::mtext("Prespecified synthetic scenarios; posterior credible intervals remain conditional on each model and prior",side=3,outer=TRUE,line=1.9,cex=.9,col="#526073")
  graphics::mtext("Error bars show 95% Monte Carlo uncertainty: Wilson intervals for proportions; normal intervals for mean bias and sensitivity.",side=1,outer=TRUE,line=1.2,cex=.78,col="#526073")
  graphics::mtext("Sensitivity uses common old-three/new-five sets and the larger new sensitivity axis. Paired inference and all analysis refusals remain in the numerical evidence files.",side=1,outer=TRUE,line=2.4,cex=.78,col="#526073")
  graphics::par(fig=c(0,1,.914,.945),mar=rep(0,4),oma=rep(0,4),new=TRUE)
  graphics::plot.new()
  graphics::legend("center",legend=c("v0.9 conditional joint model","New marginal joint model + hierarchical prior"),
    col=colors,pch=16,horiz=TRUE,bty="n",cex=.92,x.intersp=.8)
}
if(dir.exists(output)&&length(list.files(output,all.files=TRUE,no..=TRUE)))stop("Use a new or empty figure output directory.")
dir.create(output,recursive=TRUE,showWarnings=FALSE)
utils::write.csv(plot_data,file.path(output,"plotted-data.csv"),row.names=FALSE,na="")
utils::write.csv(prior_noise,file.path(output,"prior-shift-mcse.csv"),row.names=FALSE,na="")
grDevices::pdf(file.path(output,"joint-marginal-evidence.pdf"),width=16,height=12,family="sans",useDingbats=FALSE)
render();invisible(grDevices::dev.off())
grDevices::png(file.path(output,"joint-marginal-evidence.png"),width=3200,height=2400,res=200,type=if(capabilities("cairo"))"cairo"else"quartz")
render();invisible(grDevices::dev.off())
jsonlite::write_json(list(stage="confirmation",source_fingerprint=complete$source_fingerprint,
  script_sha256=digest::digest(file=script,algo="sha256"),
  input_files=setNames(lapply(file.path(input,files),function(f)digest::digest(file=f,algo="sha256")),files),
  points=nrow(plot_data),self_checks="All plotted values recomputed from raw records and checked against saved summaries; prior comparisons use common complete old-three/new-five sets.",
  uncertainty="95% Monte Carlo intervals; not posterior coefficient intervals for any individual cohort."),
  file.path(output,"figure-provenance.json"),auto_unbox=TRUE,pretty=TRUE)
writeLines(c("# Joint marginal confirmation figure","",
  "This figure consumes a completed confirmation stage. The PDF is vector output; PNG resolution is 3200 by 2400.","",
  "A: coverage among returned intervals, with Wilson Monte Carlo intervals. B: bias against the same marginal population target for both methods, with normal Monte Carlo intervals. Contextual targets stay 0.9, 0.9 and zero.","",
  "C: the old model uses mean absolute shifts at inverse-Wishart scales 0.25/four versus one. The new model uses the larger of that HW-axis response and its regression-prior response at scales five/20 versus ten. A cohort enters only when all three old and five new fits return. Individual mean error bars do not replace the paired sensitivity test. Prior families differ; this is a declared range comparison, not proof of prior independence. prior-shift-mcse.csv shows the finite-chain mean MCSE alongside shift magnitudes.","",
  "D: fit refusals among every requested default-prior fit. Summary-only failures are separate in the study records. Zero observed refusals do not imply zero future failure probability.","",
  "The first six scenarios and three contextual scenarios are separated visually, but all nine are reported. Near-zero ICC and Student-t5 are stress results; the complete prospective plan defines each adoption screen.","",
  "All plotted numbers and error bars are recomputed from raw records, checked against saved summaries and retained in plotted-data.csv. No visual or numerical result establishes universal calibration."),file.path(output,"README.md"))
cat(sprintf("Generated %d self-checked plotted points and %d prior-shift/MCSE summaries from completed confirmation.\n",nrow(plot_data),nrow(prior_noise)))
