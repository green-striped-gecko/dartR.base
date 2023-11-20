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
#' gl.report.basics(platypus.gl)
#' 
# @import tibble
#' @export
#' @return NULL

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
    datatype <- utils.check.datatype(x, verbose = 0)

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
      cat(mean(x@other$loc.metrics$rdepth))
      cat("\n")
      
      # Report composition
      cat(report("Values: "))
      cat(nLoc(x)*nInd(x),"\n")
      tmp <- table(as.matrix(x),useNA='always')
      tmp <- round(tmp*100/(nLoc(x)*nInd(x)),1)
      tmp <- t(as.data.frame(tmp))
      colnames(tmp) <- c("0","1","2","NA")
      tmp <- tmp[2,]
      tmp <- as.data.frame(t(tmp))
      rownames(tmp) <- "percent"
      print(tmp)
      cat("\n")
      
      # Report Monomorphic loci
      if (datatype == "SilicoDArT") {
        mat <- as.matrix(x)
        mono.count <- 0
        for (i in 1:nLoc(x)) {
          row <- mat[, i]  # Row for each locus
          if (all(row == 0, na.rm = TRUE) |
              all(row == 1, na.rm = TRUE)) {
            mono.count <- mono.count + 1
          }
        }
      }
      
      # SNP data
      if (datatype == "SNP") {
        mat <- as.matrix(x)
        mono.count <- 0
        for (i in 1:nLoc(x)) {
          row <- mat[, i]  # Row for each locus
          if (all(row == 0, na.rm = TRUE) |
              all(row == 2, na.rm = TRUE)) {
            mono.count <- mono.count + 1
          }
        }
      }
      cat(report("Monomorphic Loci: "))
      cat(mono.count,"\n")
      
      # Report allNA
      # Loci
      na.counter <- 0
      matrix <- as.matrix(x)
      for (i in 1:nLoc(x)) {
        row <- matrix[, i]  # Row for each locus
        if (all(is.na(row))) {
          na.counter <-na.counter + 1
        }
      }
      cat(report("Loci all NA: "))
      cat(na.counter,"\n")
      # Individuals
      na.counter <- 0
      ind.list <- array(NA, nInd(x))
      matrix <- as.matrix(x)
      i.names <- indNames(x)
      for (i in 1:nInd(x)) {
        col <- matrix[i, ]  # Row for each individual
        if (all(is.na(col))) {
          ind.list[i] <- i.names[i]
          na.counter <-na.counter + 1
        }
      }
      cat(report("Individuals all NA: "))
      cat(na.counter,"\n")
      if(na.counter > 0){cat(ind.list,"\n")}
      cat("\n")
      
      # Report sample sizes
      cat(report("Sample Sizes:\n"))
      tmp <- table(pop(x))
      print(tmp)
      cat("\n")
      
      # Report loci allNA by pop
      tmp <- array(NA,nPop(x))
      for (i in 1:nPop(x)){
        tmpop <- as.matrix(gl.keep.pop(x,popNames(x)[i],verbose=0))
        tmpsums <- apply(tmpop,2,function(x){all(is.na(x))})
        tmp[i] <- sum(tmpsums)
      }
      tmp <- t(as.data.frame(tmp))
      colnames(tmp) <- popNames(x)
      rownames(tmp) <- ""
      cat(report("Loci all NA across individuals by Population\n"))
      print(tmp)
      cat("\n")
      
      # Report individuals allNA by pop
      tmp <- array(NA,nPop(x))
      for (i in 1:nPop(x)){
        tmpop <- as.matrix(gl.keep.pop(x,popNames(x)[i],verbose=0))
        tmpsums <- apply(tmpop,1,function(x){all(is.na(x))})
        tmp[i] <- sum(tmpsums)
      }
      tmp <- t(as.data.frame(tmp))
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
