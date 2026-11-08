# Runs the illustrative workflow printed in paper/manuscript.Rmd and paper/paper.md, exactly
# as printed, and checks two things: that it raises no warning (in particular none of the
# pooling functions' outcome-role warnings), and that print(pool_survival) at a console
# width of 100 gives the Cox table shown after the code. Needs the installed catmech
# version the paper describes. Run from the repository root:
#   Rscript paper/check_workflow_example.R

suppressPackageStartupMessages({
  library(catmech)
  library(survival)
})

paper_dir <- "paper"
if (!dir.exists(paper_dir)) stop("Run from the repository root.", call. = FALSE)

fenced_block <- function(lines, from, opener) {
  start <- which(lines == opener & seq_along(lines) > from)[1L]
  if (is.na(start)) stop("No ", opener, " block after line ", from, call. = FALSE)
  end <- which(lines == "```" & seq_along(lines) > start)[1L]
  list(text = lines[(start + 1L):(end - 1L)], end = end)
}

check_paper <- function(file) {
  lines <- readLines(file.path(paper_dir, file), warn = FALSE)
  heading <- grep("^# Illustrative clinical (trial|audit) workflow", lines)
  if (length(heading) != 1L) stop(file, ": workflow heading not found once", call. = FALSE)
  code <- fenced_block(lines, heading, "```r")
  printed <- fenced_block(lines, code$end, "```text")$text

  workdir <- tempfile("catmech-workflow-")
  dir.create(workdir)
  old_wd <- setwd(workdir)
  on.exit(setwd(old_wd), add = TRUE)
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)

  env <- new.env(parent = globalenv())
  warnings_seen <- character()
  withCallingHandlers(
    eval(parse(text = code$text), envir = env),
    warning = function(w) {
      warnings_seen <<- c(warnings_seen, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  old_width <- options(width = 100)
  on.exit(options(old_width), add = TRUE)
  shown <- utils::capture.output(print(env$pool_survival))

  ok_warn <- length(warnings_seen) == 0L
  ok_table <- identical(sub("\\s+$", "", shown), sub("\\s+$", "", printed))
  cat(sprintf("%s: %d warning(s); printed Cox table %s\n", file, length(warnings_seen),
              if (ok_table) "matches the run" else "DIFFERS from the run"))
  for (w in unique(warnings_seen)) cat("  warning:", w, "\n")
  if (!ok_table) {
    cat("  printed in the paper:\n", paste0("    ", printed, "\n"), sep = "")
    cat("  produced by the run:\n", paste0("    ", shown, "\n"), sep = "")
  }
  ok_warn && ok_table
}

cat("catmech", as.character(utils::packageVersion("catmech")), "\n")
results <- vapply(c("manuscript.Rmd", "paper.md"), check_paper, logical(1))
if (!all(results)) quit(status = 1)
cat("Workflow example check: PASS\n")
