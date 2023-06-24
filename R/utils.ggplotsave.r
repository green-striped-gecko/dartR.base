#' @name utils.ggplotsave
#' @title An internal function to save a ggplot file/object to disk. 
#' @family utilities

#' @description
#' WARNING: UTILITY SCRIPTS ARE FOR INTERNAL USE ONLY AND SHOULD NOT BE USED BY END USERS AS THEIR USE OUT OF CONTEXT COULD LEAD TO UNPREDICTABLE OUTCOMES.

#' @param x Name of the ggplot object.
#' @param type Type of file to save 
#' @param dir Name of the directory to save the file.
#' @param file Name of the file to save the plot to (omit file extension)
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' progress log; 3, progress and results summary; 5, full report
#'  [default NULL, unless specified using gl.set.verbosity]
#' @param ... Parameters passed to function \link[ggplot2]{ggsave}, 
#'  such as width and height, when the ggplot is to be saved.
#'  
#'  @details
#'  #' Uses saveRDS() and ggsave().
#' Additional details ..... options for saving are specified by the parameter
#' type, which can be one of 
#' "RDS", eps", "ps", "tex" (pictex), "pdf", "jpeg", "tiff", "png", 
#' "bmp", "svg" or "wmf" (windows only). If type is specified, Whether or not "RDS", 
#' the function also saves the ggplot object as an RDS binary file using gl.save(); 
#' can be reloaded with gl.load().
#' 
#' @author Maintainer: Bernd Gruber (Post to \url{https://groups.google.com/d/forum/dartr})

#' @export
#' @return NULL

utils.ggplotsave <- function(
  x,
  type=NULL,
  dir=NULL,
  file=NULL,
  verbose=NULL,
  ...
){
  
  errorflag <- 0
  
  # If saving ------------------
  if(!is.null(type)){
    
    typelist <- c("RDS","eps", "ps", "tex", "pdf", "jpeg", "tiff", "png", "bmp", "svg", "wmf")
    if(type=="jpg"){type<-"jpeg"}
    if(type=="tif"){type<-"tiff"}
    if(type=="rds"){type<-"RDS"}
    if(!(type %in% typelist)){
      errorflag <- 1
      if(verbose >= 2){cat(error("type specified as",type,"but not included in list of acceptable types\n"))}
      if(verbose >= 2){cat(error("  Should be one of",typelist,"\n"))}
    }
    
    if (!is.ggplot(x) & errorflag==0) {
      if(verbose >= 2){cat(error(deparse(substitute(x))," is not a ggplot object, no object saved\n"))}
      errorflag <- 1
    }
    
    if(!is.null(dir) & errorflag==0){
      if(!file.exists(dir)){
        if(verbose >= 2){cat(error("Directory to receive the saved plot file does not exist or is misspecified. Defaulting to working directory\n"))}
        dir <- getwd()
      }
    } else {
      dir <- getwd()
    }
    
    if(is.null(file) & errorflag==0){
      if(verbose >= 2){cat(error("No file name provided for the plot file.\n"))}
      file <- NULL
      errorflag <- 1
    }
    
    if(errorflag==0){
      filespec <- file.path(dir, file)
      filespec.ext <- paste0(filespec,".",type)
      filespec.RDS <- paste0(filespec,".RDS")
      if(type=="RDS"){
        if(verbose >= 2){cat(report("ggplot object will be saved as RDS to",filespec.RDS,"using saveRDS()\n"))}
        saveRDS(x, filespec.RDS)
      } else {
        if(verbose >= 2){cat(report("ggplot file will be saved as",type,"to",filespec.ext,"using ggsave()\n"))}
        if(verbose >= 2){cat(report("ggplot object will also be saved as RDS binary to",filespec.RDS,"using saveRDS()\n"))}
        ggsave(x,filename=filespec.ext,device=type,dpi=600)
        saveRDS(x, filespec.RDS)
      }
    } 
  } 
  
  if(is.null(type) | errorflag==1){
    cat(error("    No plot saved\n"))
  }

return(NULL)
  
}
