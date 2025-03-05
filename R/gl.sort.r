#'@name gl.sort
#'@title Sorts genlight objects
#' @family data manipulation
#'@description 
#' This function provides the ability to sort genotypes in a genlight object by  
#' individual name or population name. 
#' 
#'@param x Genlight object containing SNP/Silicodart genotypes [required].
#'@param sort.by Specify to sort the genotypes by either 'ind', "pop" [default 'pop'].
#'@param order.by Vector used to order genotypes [default NULL]
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' progress log; 3, progress and results summary; 5, full report
#'  [default NULL, unless specified using gl.set.verbosity]
#'
#'@details 
#'This is a function to sort genotypes in a genlight object by individual name
#'or population name. This will be useful if you want to visualise the structure
#'across populations (bands) in a gl.smearplot; order of genotypes is important.
#'
#' The order.by parameter needs to be a vector upon which to effect the sort, of length of 
#' nPop(gl) if sort.by is 'pop' or  nInd(gl) if sort.by is 'ind'.  For sort.by='ind' 
#' order.by can be a vector such as a variable in gl@other$ind.metrics.
#' 
#' If not specified by nominating a vector with order.by, alphabetical order of populations or individuals is used. 
#'  
#'@author Bernd Gruber (Post to \url{https://groups.google.com/d/forum/dartr})
#'
#'@examples 
#'#sort by populations
#'bc <- gl.sort(bandicoot.gl)
#'#sort from West to East
#'bc2 <- gl.sort(bandicoot.gl, sort.by="pop" ,
#'order.by=c("WA", "SA", "VIC", "NSW", "QLD"))
#'#sort by missing values
#'miss <- rowSums(is.na(as.matrix(bandicoot.gl)))
#'bc3 <- gl.sort(bandicoot.gl, sort.by="ind", order.by=miss)
#'gl.smearplot(bc3)
#'@export 
#'
#'@return Returns a reordered genlight object. Sorts also the ind/loc.metrics 
#'and coordinates accordingly

gl.sort <- function(x,
                  sort.by = "pop",
                  order.by = NULL,
                  verbose = NULL) {
  
  # SET VERBOSITY
  verbose <- gl.check.verbosity(verbose)
  
  # FLAG SCRIPT START
  funname <- match.call()[[1]]
  utils.flag.start(func=funname,build="v.2023.2",verbose=verbose)
  
  # CHECK DATATYPE
  datatype <- utils.check.datatype(x, verbose=verbose)
  
  if (!is(x, "dartR")) {
    class(x) <- "dartR"  
    if (verbose>2) {
      cat(warn("Warning: Standard adegenet genlight object encountered. Converted to compatible dartR genlight object\n"))
      cat(warn("                    Should you wish to convert it back to an adegenet genlight object for later use outside dartR, 
                 please use function dartR2gl\n"))
    }
  }
  
  # FUNCTION SPECIFIC ERROR CHECKING
  
  if (is.na(pmatch(sort.by,c("pop","ind")))){
    stop(error("sort.by is not either 'pop', 'ind'. Please specify one of these options."))
  }
  
  if (!is.null(order.by)) {
    if (sort.by=="pop") len <- nPop(x)
    if (sort.by=="ind") len <- nInd(x)
    #if (sort.by=="other ") len <- nLoc(x)
    if (len!=length(order.by)) stop(error("Length of vector order.by does not match length of",sort.by, ". Check your input parameters and make sure length of your sort.by selection matches the length of the vector provided by order.by"))
  }
  
  # DO THE JOB
  
  if (sort.by=="pop") {
    if (is.null(order.by)) index <- order(pop(x)) else {
      if (!(sum(pop(x) %in% order.by))==nInd(x)) stop(error("order.by does not contain all levels of sort.by."))
      index  <- order(factor(pop(x), levels=order.by))
    } 
    xx <- x[index,]
    xx@other$latlon <- x@other$latlon[index,]
    xx@other$ind.metrics <- x@other$ind.metrics[index,]
    if (!is.null(order.by)) pop(xx)<- factor(pop(xx), levels=order.by)
    }
    if (sort.by=="ind") {
    if (is.null(order.by)) index <- order(indNames(x)) else {
      if (!(length(order.by)==nInd(x))) stop(error("order.by does not contain all levels of sort.by."))
      index  <- order(order.by)
    } 
    
  # Apply sorting
    xx <- x[index,]
    xx@other$latlon <- x@other$latlon[index,]
    xx@other$ind.metrics <- x@other$ind.metrics[index,]
  }
  
  # ADD TO HISTORY
  if (is(xx,"genlight")) {
    nh <- length(xx@other$history)
    xx@other$history[[nh + 1]] <- c(match.call())
  } 
  
  return(xx)
}
    
