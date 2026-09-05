#' @name gl2demerelate
#' @title Creates a dataframe suitable for input to package \{Demerelate\} from a
#' genlight \{adegenet\} object
#' @family linker
#'
#' @description
#' Converts a genlight object containing SNP data into a dataframe in the
#' format expected by package \{Demerelate\} for estimating pairwise
#' relatedness: one row per individual, with the sample identifier and
#' population followed by two adjacent allele columns per locus.
#'
#' @details
#' Each SNP genotype is recoded as a pair of alleles labelled 1 (reference)
#' and 2 (alternate): 0 becomes 1/1, 1 becomes 1/2, 2 becomes 2/2, and
#' missing genotypes are retained as NA in both allele columns, which
#' Demerelate reads as missing. The two allele columns for each locus are
#' suffixed _1 and _2 and kept adjacent. Locus names are sanitized for
#' Demerelate: '-' and '|' are replaced with '_' and '/' is removed.
#'
#' @param x Name of the genlight object containing the SNP data [required].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' progress log; 3, progress and results summary; 5, full report
#' [default NULL, adopting the global verbosity set by gl.set.verbosity(),
#' or 2 if no global is set].
#'
#' @return A dataframe suitable as input to package \{Demerelate\}
#'
#' @author Author(s): Luis Mijangos. Custodian: Luis Mijangos -- Post to
#' \url{https://groups.google.com/d/forum/dartr}
#'
#' @examples
#' if (isTRUE(getOption("dartR_fbm"))) testset.gl <- gl.gen2fbm(testset.gl)
#' df <- gl2demerelate(testset.gl)
#'
#' @export

gl2demerelate <- function(x,
                          verbose = NULL) {
    # SET VERBOSITY
    verbose <- gl.check.verbosity(verbose)
    
    # FLAG SCRIPT START
    funname <- match.call()[[1]]
    utils.flag.start(func = funname,
                     build = "v.2023.2",
                     verbose = verbose)
    
    # CHECK DATATYPE
    datatype <- utils.check.datatype(x, accept = "SNP", verbose = verbose)
    
    # DO THE JOB
    
    # Convert the genlight data to a form expected by demerelate Strategy is to
    # create two identical genlight objects, converted to
    # matrices, then to have one hold one allelic state for each locus and the
    #other to hold the alternate allelic state for each locus.
    # The locus names are relabeled _1 and _2 in each matrix and then the two
    #matrices are concatenated side by side. The resultant matrix
    # is ordered on locus to bring the two allelic states for a locus back to
    #adjacency. Format tidying up and bob's your uncle.
    
    x1 <- as.matrix(x) # x is a genlight object
    x2 <- as.matrix(x)
    
    x1[x1 == 2] <- 2 # homozygote alternate
    x1[x1 == 1] <- 1 # heterozygote
    x1[x1 == 0] <- 1 # homozygote reference
    colnames(x1) <- gsub(" ", "", paste(colnames(x1), "_1"))
    
    x2[x2 == 2] <- 2 # homozygote alternate
    x2[x2 == 1] <- 2 # heterozygote
    x2[x2 == 0] <- 1 # homozygote reference
    colnames(x2) <- gsub(" ", "", paste(colnames(x2), "_2"))
    
    x_temp <- cbind(x1, x2)
    x_temp <- x_temp[, order(colnames(x_temp))]
    #x_temp[is.na(x_temp)] <- 0 # Related uses zero as missing
    
    # Tidy up the locus names
    colnames(x_temp) <- gsub("-", "_", colnames(x_temp), fixed = TRUE)
    colnames(x_temp) <- gsub("/", "", colnames(x_temp), fixed = TRUE)
    colnames(x_temp) <- gsub("|", "_", colnames(x_temp), fixed = TRUE)
    
    #Add the individual names
    
    df <-
        cbind.data.frame(row.names(x_temp), factor(pop(x)), x_temp, stringsAsFactors = FALSE)
    df[, 1] <- as.character(df[, 1])
    
    # Convert to a dataframe suitable for input to package demerelate
    #df <- data.frame(x_temp, stringsAsFactors = FALSE)
    names(df)[1] <- "Sample-ID"
    names(df)[2] <- "Population"
    #  df[,2] <- factor(df[,2])
    
    # FLAG SCRIPT END
    
    if (verbose > 0) {
        cat(report("Completed:", funname, "\n"))
    }
    
    return(df)
}
