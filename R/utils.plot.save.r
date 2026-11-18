#' @name utils.plot.save
#' @title An internal function to save a ggplot object to disk in RDS binary format 
#' @family utilities

#' @description
#' WARNING: UTILITY SCRIPTS ARE FOR INTERNAL USE ONLY AND SHOULD NOT BE USED BY END USERS AS THEIR USE OUT OF CONTEXT COULD LEAD TO UNPREDICTABLE OUTCOMES.

#' @param x Name of the ggplot object.
#' @param dir Name of the directory to save the file.
#' @param file Name of the file to save the plot to (omit file extension)
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' progress log; 3, progress and results summary; 5, full report
#'  [default NULL, unless specified using gl.set.verbosity]
#' @param ... Ignored; retained for call compatibility.
#'  
#' @details
#' An internal function to save a ggplot object to disk in RDS binary format.
#' Uses saveRDS() to save the file with an .RDS extension; can be reloaded with gl.load().

#' @author Custodian: Arthur Georges (Post to \url{https://groups.google.com/d/forum/dartr})

#' @keywords internal
#' @export
#' @return NULL, invisibly


utils.plot.save <- function(
  x,
  dir=NULL,
  file=NULL,
  verbose=NULL,
  ...
){

  # SET VERBOSITY
  verbose <- gl.check.verbosity(verbose)

  if(is.null(file)){
    if (verbose >= 2) {
      cat(report("  No plot saved\n"))
    }
  } else {

    if(is.null(dir)){
      dir <- paste0(getwd(),"/")
    }

    if(!dir.exists(dir)){
      dir <- tempdir()
      if(verbose >= 2){cat(warn("Specified directory does not exist, ggplot saved to", dir, "\n"))}
    }
    
    filespec <- file.path(dir, file)
    filespec <- paste0(filespec,".RDS")
    if(verbose >= 2){cat(report("ggplot object saved as RDS binary to",filespec,"using saveRDS()\n"))}
    saveRDS(x, filespec)

  }

  return(invisible(NULL))

}
