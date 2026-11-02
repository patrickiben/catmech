# Usage: Rscript check_interface.R setups/your_file.R
# Checks the interface and timing on the 12 fixture cohorts. It does not judge the method.
args <- commandArgs(TRUE); source(args[1])
fx <- readRDS("fixtures.rds"); res <- lapply(seq_along(fx), function(i) { p <- fx[[i]]; p$truth <- NULL; p$target <- NULL
  t0 <- Sys.time(); got <- tryCatch(fit_method(p, seed = 1000L + i), error = function(e) e); s <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  ok <- !inherits(got, "error") && all(c("estimate","se","df") %in% names(got)) && all(is.finite(unlist(got[c("estimate","se","df")]))) && got$se > 0 && got$df > 0
  data.frame(cell = fx[[i]]$cell, ok = ok, seconds = round(s, 3), note = if (inherits(got,"error")) substr(conditionMessage(got),1,80) else "") })
print(do.call(rbind, res))
