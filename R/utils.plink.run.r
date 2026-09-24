#' @name utils.plink.run
#' @title Runs PLINK from within R
#'
#' @description Runs PLINK from within R. 
#' @param dir.in The path where the data files are
#' @param plink.cmd The 'name' to call plink. This will depend on the file name 
#' (without the extension '.exe' if on windows) or the name of the PATH variable
#' @param plink.path The path where the executable is. If plink is listed in
#'   the PATH then there is no need for this. This is what the option "path"
#'   means
#' @param out The root of the output file name
#' @param syntax the flags to pass to plink call; any file names in it must
#'   be quoted by the caller (e.g. with \code{shQuote()})
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' progress log; 3, progress and results summary; 5, full report [default NULL].
#' @return A character vector with the command used for PLINK. Stops with
#'   an error, quoting PLINK's output, when PLINK exits with an error.
#' @details
#' PLINK needs to be installed on the 
#'   machine and syntax used need to be appropriate for the version installed.
#' @references
#' Purcell, Shaun, et al. 'PLINK: a tool set for whole-genome association and
#' population-based linkage analyses.' The American journal of human genetics
#' 81.3 (2007): 559-575.
#' @keywords internal
#' @export
#' @author Custodian: Carlo Pacioni and Luis Mijangos (Post to
#'  \url{https://groups.google.com/d/forum/dartr})

utils.plink.run <- function(dir.in,  
                            plink.cmd = "plink",  
                            plink.path = "path",  
                            out = "hapmap1",
                            syntax, 
                            verbose = NULL) {
  
  # SET VERBOSITY
  verbose <- gl.check.verbosity(verbose)
  
  # FLAG SCRIPT START
  funname <- match.call()[[1]]
  utils.flag.start(func = funname, verbose = verbose)
  
  # DO THE JOB
  old.wd <- getwd()
  on.exit(setwd(old.wd))
  # I set wd to where the user wants the output file and copy there the 
  # executable because providing the path to PLINK doesn't work for me on 
  # Windows and I couldn't find a better solution
  
  setwd(dir.in)

  # plink.path == "path" means the executable is found via the system PATH
  if (plink.path == "path") {
    exe <- plink.cmd
  } else {
    exe <- file.path(plink.path, plink.cmd)
  }
  # The executable and output name are quoted, so paths containing spaces
  # reach PLINK as single arguments; paths inside syntax must be quoted by
  # the caller
  cmd <- paste(shQuote(exe), syntax, "--out", shQuote(out))
  # Capture PLINK's output, stderr included (2>&1), and print it only at
  # verbose >= 3; it bypasses sink() and would print even at verbose = 0
  # (gl.read.PLINK review, F5). A non-zero exit stops with PLINK's own
  # message rather than letting the caller fail later on a missing file.
  # system() itself errors when the command cannot be started at all
  # (e.g. a wrong plink.path); report that the same way
  plink.out <- tryCatch(
    suppressWarnings(system(paste(cmd, "2>&1"), intern = TRUE)),
    error = function(e) structure(conditionMessage(e), status = 127L)
  )
  status <- attr(plink.out, "status")
  if (verbose >= 3) {
    cat(report(paste(plink.out, collapse = "\n"), "\n"))
  }
  if (!is.null(status) && status != 0) {
    stop(error(
      "Fatal Error: PLINK exited with status", status, "running\n", cmd,
      "\nLast lines of the PLINK output:\n",
      paste(utils::tail(plink.out, 5), collapse = "\n"), "\n"
    ))
  }
  
  # FLAG SCRIPT END
  
  if (verbose > 0) {
    cat(report("Completed:", funname, "\n"))
  }
  
  return(cmd)
}
