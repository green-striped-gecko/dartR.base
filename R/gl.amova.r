#' @name gl.amova
#' @title Performs AMOVA using genlight data
#' @family basic statistics

#' @description 
#' This script performs an AMOVA based on the genetic distance matrix from
#' stamppNeisD() [package StAMPP] using the amova() function from the package
#' PEGAS for exploring within and between population variation. For detailed
#' information use their help pages: ?pegas::amova, ?StAMPP::stamppAmova. Be
#' aware due to a conflict of the amova functions from various packages I had
#' to 'hack' StAMPP::stamppAmova to avoid a namespace conflict.

#' @param x Name of the genlight containing the SNP genotypes, with
#' population information [required].
#' @param distance Distance matrix between individuals (if not provided NeisD
#' from StAMPP::stamppNeisD is calculated) [default NULL].
#' @param permutations Number of permutations to perform for hypothesis
#' testing [default 100]. Please note should be set to 1000 for analysis.
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' progress log ; 3, progress and results summary; 5, full report
#' [default 2, unless specified using gl.set.verbosity].
#' 
#' @author Bernd Gruber (bugs? Post to 
#' \url{https://groups.google.com/d/forum/dartr})
#' 
#' @examples
#' #permutations should be higher, here set to 1 because of speed
#' out <- gl.amova(bandicoot.gl, permutations=1)
#' 
#' @export
#' @return An object of class 'amova' which is a list with a table of sums of
#' square deviations (SSD), mean square deviations (MSD), and the number of
#' degrees of freedom, and a vector of variance components.


# Function to perform AMOVA (Analysis of Molecular Variance) using genlight data
gl.amova <- function(x,
                     distance = NULL,
                     permutations = 100,
                     verbose = NULL) {
  
  # SET VERBOSITY: Check and set the verbosity level
  verbose <- gl.check.verbosity(verbose)
  
  # FLAG SCRIPT START: Start recording function execution for profiling
  funname <- match.call()[[1]]
  utils.flag.start(func = funname,
                   build = "v.2023.3",
                   verbose = verbose)
  
  # CHECK DATATYPE: Ensure the input data x is of the right datatype
  datatype <- utils.check.datatype(x, verbose = verbose)
  
  # CHECK IF PACKAGES ARE INSTALLED: Ensure the required package 'pegas' is installed
  pkg <- "pegas"
  if (!(requireNamespace(pkg, quietly = TRUE))) {
    cat(error(
      "Package",
      pkg,
      " needed for this function to work. Please install it.\n"
    ))
    return(-1)
  }
  
  # CALCULATE DISTANCES: Calculate distances if not provided as an argument
  if (is.null(distance)) {
    class(x)<- "genlight"  # Needs to be 'genlight' due to stampp
    dd <- StAMPP::stamppNeisD(x, FALSE)
  } else {
    dd <- distance
  }
  
  # PROCESS 'genlight' OBJECT: Convert and process 'genlight' object if present
  if (is(x,"genlight")) {
    geno2 <- x
    geno <- as.matrix(geno2)
    sample <- row.names(geno)
    pop.names <- pop(geno2)
    ploidy <- ploidy(geno2)
    geno <- geno * (1 / ploidy)
    geno[is.na(geno)] <- NaN
    format <- vector(length = length(geno[, 1]))
    format[1:length(geno[, 1])] <- "genlight"
    pops <- unique(pop.names)
    pop.num <- vector(length = length(geno[, 1]))
    for (i in 1:length(geno[, 1])) {
      pop.num[i] <- which(pop.names[i] == pops)
    }
    genoLHS <-
      as.data.frame(cbind(sample, pop.names, pop.num, ploidy, format))
    geno <- cbind(genoLHS, geno)
    geno[, 2] <- as.character(pop.names)
    geno[, 4] <- as.numeric(as.character(geno[, 4]))
    row.names(geno) <- NULL
  }
  
  # PREPARE DATA FOR AMOVA: Assign relevant variables to temporary environment
  pop.names <- geno[, 2]
  pop.names <- factor(pop.names)
  temp <- new.env()
  assign("distance", dd, envir = temp)
  assign("pop.names", pop.names, envir = temp)
  assign("permutations", permutations, envir = temp)
  
  # PERFORM AMOVA: Execute AMOVA using pegas::amova function
  res <- with(temp, pegas::amova(distance ~ pop.names, nperm = permutations))
  rm(pop.names, permutations, temp)
  
  # FLAG SCRIPT END: End recording function execution for profiling
  
  if (verbose >= 1) {
    cat(report("Completed:", funname, "\n"))
  }
  
  # RETURN: Return the result of AMOVA
  return(res)
}

