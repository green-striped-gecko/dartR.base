## Runs every help-page example of a dartRverse package twice, in memory and
## file-backed (environment variable dartR_fbm=true, which makes the examples
## convert their data with gl.gen2fbm), and compares the outcomes.
##
## Usage (from any folder; the package must be installed):
##   Rscript run_examples_fbm.R <package source folder> <output csv>
##           [timeout per example in s, default 600] [topic regex, default all]
##
## Exit status 1 when any example fails file-backed but not in memory (a
## regression of the FBM code; used by the fbm-examples GitHub workflow);
## examples failing in both modes are reported but do not fail the run.
##
## Each example runs in its own R session (Rscript), with \donttest sections
## included and \dontrun sections excluded, in a temporary working folder,
## with plots sent to a null device. An example that fails in both modes
## points to the environment (a missing package or external program); one
## that fails only file-backed points to the FBM code.

args <- commandArgs(TRUE)
pkg_dir <- normalizePath(args[1], winslash = "/")
out_csv <- args[2]
timeout <- if (length(args) >= 3) as.numeric(args[3]) else 600
topics <- if (length(args) >= 4) args[4] else "."
pkg <- read.dcf(file.path(pkg_dir, "DESCRIPTION"), fields = "Package")[1, 1]

rd_files <- list.files(file.path(pkg_dir, "man"), pattern = "\\.Rd$",
                       full.names = TRUE)
rd_files <- rd_files[grepl(topics, sub("\\.Rd$", "", basename(rd_files)))]
work <- file.path(tempdir(), "fbm_examples")
dir.create(work, showWarnings = FALSE)

run_one <- function(ex_file, fbm) {
  log <- tempfile(fileext = ".log")
  runner <- tempfile(fileext = ".R")
  # The variable is set inside the script, before the package is loaded
  # (dartR.base reads it in .onLoad); system2(env = ) does not work on
  # Windows.
  writeLines(c(
    sprintf("Sys.setenv(dartR_fbm = '%s')", if (fbm) "true" else "false"),
    sprintf("suppressPackageStartupMessages(library(%s))", pkg),
    "suppressPackageStartupMessages(library(dartR.data))",
    "pdf(NULL)",
    sprintf("setwd(%s)", deparse(work)),
    "t0 <- proc.time()[['elapsed']]",
    sprintf("r <- tryCatch({ source(%s, echo = FALSE, local = new.env()); 'OK' },",
            deparse(normalizePath(ex_file, winslash = "/"))),
    "  error = function(e) paste('ERROR:', conditionMessage(e)))",
    "cat('\\n@@RESULT@@', gsub('[\\r\\n]+', ' ', r), '\\n')",
    "cat('@@TIME@@', proc.time()[['elapsed']] - t0, '\\n')"
  ), runner)
  status <- suppressWarnings(system2(
    file.path(R.home("bin"), "Rscript"), runner, stdout = log, stderr = log,
    timeout = timeout))
  txt <- readLines(log, warn = FALSE)
  res <- grep("@@RESULT@@", txt, value = TRUE)
  tim <- grep("@@TIME@@", txt, value = TRUE)
  list(
    result = if (length(res)) trimws(sub(".*@@RESULT@@", "", res[1])) else
      if (identical(status, 124L)) "TIMEOUT" else
        paste("CRASH:", paste(tail(txt, 3), collapse = " ")),
    time = if (length(tim)) as.numeric(sub(".*@@TIME@@", "", tim[1])) else NA
  )
}

rows <- list()
for (rd in rd_files) {
  ex <- tempfile(fileext = ".R")
  tools::Rd2ex(rd, ex, commentDonttest = FALSE, commentDontrun = TRUE)
  if (!file.exists(ex)) next                  # no examples
  name <- sub("\\.Rd$", "", basename(rd))
  uses_fbm <- any(grepl("dartR_fbm", readLines(ex, warn = FALSE)))
  g <- run_one(ex, fbm = FALSE)
  f <- run_one(ex, fbm = TRUE)
  rows[[name]] <- data.frame(topic = name, uses_fbm_switch = uses_fbm,
                             gen = g$result, fbm = f$result,
                             gen_s = g$time, fbm_s = f$time)
  cat(sprintf("%-40s gen: %-6s fbm: %-6s\n", name,
              substr(g$result, 1, 6), substr(f$result, 1, 6)))
}
tab <- do.call(rbind, rows)
utils::write.csv(tab, out_csv, row.names = FALSE)

ok <- function(r) r == "OK"
cat("\n", pkg, ":", nrow(tab), "help pages with examples\n")
cat("  both OK:        ", sum(ok(tab$gen) & ok(tab$fbm)), "\n")
cat("  fail only FBM:  ", sum(ok(tab$gen) & !ok(tab$fbm)), "\n")
cat("  fail only gen:  ", sum(!ok(tab$gen) & ok(tab$fbm)), "\n")
cat("  fail both:      ", sum(!ok(tab$gen) & !ok(tab$fbm)), "\n")
cat("  examples with the dartR_fbm switch:", sum(tab$uses_fbm_switch), "\n")

only_fbm <- tab$topic[ok(tab$gen) & !ok(tab$fbm)]
both <- tab$topic[!ok(tab$gen) & !ok(tab$fbm)]
if (length(both)) {
  cat("\nFailing in both modes (environment, e.g. a missing package):\n")
  print(tab[tab$topic %in% both, c("topic", "gen")], row.names = FALSE)
}
if (length(only_fbm)) {
  cat("\nFailing only file-backed (FBM regression):\n")
  print(tab[tab$topic %in% only_fbm, c("topic", "fbm")], row.names = FALSE)
  quit(status = 1)
}
