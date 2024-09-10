#' @name gl.add.indmetrics
#' @title Adds metadata into a genlight object
#' @description
#' This function adds the metadata information to the slot ind.metrics and
#' populates population and coordinates information slots if the they are
#' found in the metadata.
#' @param x Name of the genlight object containing the SNP data, or the genind
#'  object containing the SilocoDArT data [required].
#' @param ind.metafile path and name of CSV file containing the metadata
#' information for each individual (see details for explanation) [required].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' progress log ; 3, progress and results summary; 5, full report
#' [default 2, unless specified using gl.set.verbosity].
#'
#' @details
#' The ind.metadata file needs to have very specific headings. First a column
#' with a heading named 'id'. Here the ids must match the ids in the genlight
#'  object, e.g. \code{indNames(your_genlight)}. The following column headings
#'  are optional:
#'  \itemize{
#'  \item 'pop' - specifies the population membership of each individual.
#'  \item 'lat' - latitude coordinates (in decimal degrees WGS1984 format).
#'  \item 'lon' - longitude coordinates (in decimal degrees WGS1984 format).
#'  }
#'
#'  Additional columns with individual metadata can be imported (e.g. age,
#'  sex, etc).
#' @author Custodian: Luis Mijangos -- Post to
#' \url{https://groups.google.com/d/forum/dartr}
#' @examples
#' dartfile <- system.file('extdata','testset_SNPs_2Row.csv', package='dartR.data')
#' metadata <- system.file('extdata','testset_metadata.csv', package='dartR.data')
#' gl <- gl.read.dart(dartfile, probar=TRUE)
#' gl <- gl.add.indmetrics(gl, ind.metafile = metadata)
#' @export
#' @return A genlight object with metadata information for each individual.

gl.add.indmetrics <- function(x,
                              ind.metafile,
                              verbose = NULL) {
  # SET VERBOSITY
  verbose <- gl.check.verbosity(verbose)
  
  # FLAG SCRIPT START
  funname <- match.call()[[1]]
  utils.flag.start(func = funname,
                   build = "v.2023.2",
                   verbose = verbose)
  
  # CHECK DATATYPE
  datatype <- utils.check.datatype(x, verbose = verbose)
  
  # DO THE JOB
  if (!is.null(ind.metafile)) {
    if (verbose >= 2) {
      cat(report(paste(
        "Adding individual metrics:", ind.metafile, ".\n"
      )))
    }
    ###### population and individual file to link AAnumbers to populations...
    ind.cov <-
      read.csv(ind.metafile,
               header = T,
               stringsAsFactors = T)
    # is there an entry for every individual
    
    id.col <- match("id", names(ind.cov))
    
    if (is.na(id.col)) {
      stop(error("Fatal Error: There is no id column\n"))
    } else {
      ind.cov[, id.col] <-
        trimws(ind.cov[, id.col], which = "both")  #trim spaces
      
      if (length(ind.cov[, id.col]) != length(unique(ind.cov[, id.col]))) {
        cat(error(
          "Individual names are not unique. You need to change them!\n"
        ))
        stop()
      }
      
      # reorder
      if (length(ind.cov[, id.col]) != length(indNames(x))) {
        cat(
          warn(
            "Ids for individual metadata does not match the number of ids in the SNP data file. Maybe this is fine if a subset matches.\n"
          )
        )
        nam.indmeta <- ind.cov[, id.col]
        nam.dart <- indNames(x)
        
        nm.indmeta <-
          nam.indmeta[!nam.indmeta %in% nam.dart]
        nm.inddart <-
          nam.dart[!nam.dart %in% nam.indmeta]
        if (length(nm.indmeta) > 0) {
          cat(warn("ind.metafile ids not matched were:\n"))
          print(nm.indmeta)
        }
        if (length(nm.inddart) > 0) {
          cat(warn("DArT file ids not matched were:\n"))
          print(nm.inddart)
        }
      }
      
      ord <- match(indNames(x), ind.cov[, id.col])
      ord <- ord[!is.na(ord)]
      
      if (length(ord) > 1 & length(ord) <= nInd(x)) {
        if (verbose >= 2) {
          cat(report(
            paste(
              "  Ids for individual metadata (at least a subset of) are matching!\n"
            )
          ))
          cat(report(
            paste(
              "  Found ",
              length(ord == nInd(x)),
              "matching ids out of",
              nrow(ind.cov),
              "ids provided in the ind.metadata file.\n "
            )
          ))
        }
        ord2 <-
          match(ind.cov[ord, id.col], indNames(x))
        x <- x[ord2,]
      } else {
        stop(error("Fatal Error: Individual ids are not matching!!!!\n"))
      }
    }
    
    pop.col <- match("pop", names(ind.cov))
    
    if (is.na(pop.col)) {
      if (verbose >= 1) {
        cat(
          warn(
            "Warning: There is no pop column, created one with all pop1 as default for all individuals\n"
          )
        )
      }
      pop(x) <- factor(rep("pop1", nInd(x)))
    } else {
      pop(x) <- as.factor(ind.cov[ord, pop.col])
      if (verbose >= 2) {
        cat(report(" Added population assignments.\n"))
      }
    }
    
    lat.col <- match("lat", names(ind.cov))
    lon.col <- match("lon", names(ind.cov))
    if (verbose >= 2) {
      if (is.na(lat.col)) {
        cat(warn(
          "Warning: Individual metrics do not include a latitude (lat) column\n"
        ))
      }
      if (is.na(lon.col)) {
        cat(warn(
          "Warning: Individual metrics do not include a longitude (lon) column\n"
        ))
      }
    }
    if (!is.na(lat.col) & !is.na(lon.col)) {
      x@other$latlon <- ind.cov[ord, c(lat.col, lon.col)]
      rownames(x@other$latlon) <- ind.cov[ord, id.col]
      if (verbose >= 2) {
        cat(report("  Added latlon data.\n"))
      }
    }
    
    other.col <- names(ind.cov)
    if (length(other.col) > 0) {
      # conserving previous ind.metrics
      x@other$ind.metrics <- as.data.frame(cbind(x@other$ind.metrics,ind.cov[ord, other.col, drop = FALSE]))
      rownames(x@other$ind.metrics) <- ind.cov[ord, id.col]
      if (verbose >= 2) {
        cat(report(
          paste(" Added ",
                other.col,
                " to the other$ind.metrics slot.\n")
        ))
      }
    }
  }
  
  # ADD TO HISTORY
  nh <- length(x@other$history)
  x@other$history[[nh + 1]] <- match.call()
  
  # FLAG SCRIPT END
  
  if (verbose > 0) {
    cat(report("Completed:", funname, "\n"))
  }
  
  return(x)
  
}