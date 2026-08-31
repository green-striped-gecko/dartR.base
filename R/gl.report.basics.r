#' @name gl.report.basics
#' @title Basic statistics for a genlight object
#' @family unmatched report

#' @description
#' Calculates basic statistics for a genlight object. 
#' 
#' @param x Name of the genlight object  [required].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' progress log; 3, progress and results summary; 5, full report
#' [default 2, unless specified using gl.set.verbosity].
#' 
#' @author Custodian: Arthur Georges (Post to
#' \url{https://groups.google.com/d/forum/dartr})
#' 
#' @examples
#' if (isTRUE(getOption("dartR_fbm"))) platypus.gl <- gl.gen2fbm(platypus.gl)
#' gl.report.basics(platypus.gl)
#' 
#' @return Returns NULL invisibly; the results are printed to the console.
#' @export

gl.report.basics <- function(x,
                           verbose = NULL) {
    # SET VERBOSITY
    verbose <- gl.check.verbosity(verbose)
    
    # FLAG SCRIPT START
    funname <- match.call()[[1]]
    utils.flag.start(func = funname,
                     build = "v.2023.3",
                     verbose = verbose)
    
    # CHECK DATATYPE
    datatype <- utils.check.datatype(x, verbose = verbose)

    if(verbose >= 1){
      cat(report("\nSUMMARY STATISTICS\n\n"))
      
      # Report the datatype
      cat(report("Datatype: "))
          cat(datatype,"\n")
      
      # Report dimensions
      cat(report("Loci: "))
      cat(nLoc(x),"\n")
      cat(report("Individuals: "))
      cat(nInd(x),"\n")
      cat(report("Populations: "))
      cat(nPop(x),"\n")
      cat("\n")
      
      # Report average read depth
      cat(report("Average Read Depth: "))
      if (is.null(x@other$loc.metrics$rdepth)) {
        cat("not available")
      } else {
        cat(round(mean(x@other$loc.metrics$rdepth, na.rm = TRUE), 2))
      }
      cat("\n")

      # Report composition (tabulated over explicit levels so the table
      # has a fixed shape for both datatypes and for data in which a
      # genotype class is entirely absent)
      cat(report("Values: "))
      cat(nLoc(x)*nInd(x),"\n")
      mat <- as.matrix(x)
      if (datatype == "SilicoDArT") {
        lv <- c("0", "1")
      } else {
        lv <- c("0", "1", "2")
      }
      tmp <- table(factor(mat, levels = lv), useNA = "always")
      tmp <- round(tmp * 100 / (nLoc(x) * nInd(x)), 1)
      tmp <- as.data.frame(t(as.matrix(tmp)))
      colnames(tmp) <- c(lv, "NA")
      rownames(tmp) <- "percent"
      print(tmp)
      cat("\n")

      # Report Monomorphic loci (all-NA loci count as monomorphic,
      # consistent with gl.filter.monomorphs)
      non.na <- colSums(!is.na(mat))
      if (datatype == "SilicoDArT") {
        mono.count <- sum(colSums(mat == 0, na.rm = TRUE) == non.na |
                            colSums(mat == 1, na.rm = TRUE) == non.na)
      } else {
        mono.count <- sum(colSums(mat == 0, na.rm = TRUE) == non.na |
                            colSums(mat == 2, na.rm = TRUE) == non.na)
      }
      cat(report("Monomorphic Loci: "))
      cat(mono.count,"\n")

      # Report allNA
      # Loci
      cat(report("Loci all NA: "))
      cat(sum(non.na == 0),"\n")
      # Individuals
      ind.all.na <- indNames(x)[rowSums(!is.na(mat)) == 0]
      cat(report("Individuals all NA: "))
      cat(length(ind.all.na),"\n")
      if(length(ind.all.na) > 0){cat(ind.all.na,"\n")}
      cat("\n")

      # Report sample sizes
      cat(report("Sample Sizes:\n"))
      tmp <- table(pop(x))
      print(tmp)
      cat("\n")

      # Report loci and individuals allNA by pop (single pass)
      loc.na.pop <- array(NA,nPop(x))
      ind.na.pop <- array(NA,nPop(x))
      for (i in 1:nPop(x)){
        tmpop <- mat[as.character(pop(x)) == popNames(x)[i], , drop = FALSE]
        loc.na.pop[i] <- sum(colSums(!is.na(tmpop)) == 0)
        ind.na.pop[i] <- sum(rowSums(!is.na(tmpop)) == 0)
      }
      tmp <- t(as.data.frame(loc.na.pop))
      colnames(tmp) <- popNames(x)
      rownames(tmp) <- ""
      cat(report("Loci all NA across individuals by Population\n"))
      print(tmp)
      cat("\n")

      tmp <- t(as.data.frame(ind.na.pop))
      colnames(tmp) <- popNames(x)
      rownames(tmp) <- ""
      cat(report("Individuals all NA across loci by Population\n"))
      print(tmp)
      cat("\n")
    }
    
    if(verbose >= 3){
      
      # List locus metrics
      cat(report("Locus Metrics:\n"))
      cat(names(x@other$loc.metrics),"\n",sep=", ")
      print(tibble::tibble(x@other$loc.metrics))
      cat("\n\n")
      
      # List individual metrics
      cat(report("Individual Metrics\n"))
      cat(names(x@other$ind.metrics),"\n",sep=", ")
      print(tibble::tibble(x@other$ind.metrics))
      cat("\n\n")
      
      # Report history
      cat(report("History:\n")) 
      print(x@other$history)
      
    }
    
   # FLAG SCRIPT END
    
    if (verbose >= 1) {
        cat(report("Completed:", funname, "\n"))
    }
    
    # RETURN
    
    invisible(NULL)
}
