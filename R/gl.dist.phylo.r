#' @name gl.dist.phylo
#' @title Generates a distance matrix taking into account a substitution model
#' 
#' @family phylogeny
#' 
#' @description
#' Generates a distance matrix for individuals or populations in a genlight object
#'  using one of a selection of substitution models.
#'  
#' @param xx Name of the genlight object containing the SNP data [required].
#' 
# Parameters that govern creation of the sequence file
#' @param min.tag.len Minimum tag length of sequence tags to be used in the analysis [default NULL]
#' 
# Parameters that govern the creation of the distance matrix
#' @param subst.model The evolutionary model of nucleotide substitutions to employ in calculating genetic
#' distance between individuals [default "F81"]
#' @param by.pop If TRUE, the distance matrix is based on comparing populations; if FALSE, on individuals [default TRUE].
#' @param pairwise.missing How to handle missing sequences [default TRUE]
#' 
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' progress log; 3, progress and results summary; 5, full report
#' [default 2, unless specified using gl.set.verbosity].
#' 
#' @details
#' The script takes a genlight object as input, creates a set of sequences from the trimmed sequence
#' tags for each individual, calculates distances between the individuals and then optionally averages those distances between
#' the populations defined in the genlight object (typically OTUs). 
#' 
#'      min.tag.length : Sequence tags can vary considerably in length, which results in large numbers of Ns in
#'      alignments. This can have an impact of distance measures depending on how missing values are managed.
#'      To minimize this effect, you might elect to filter on tag length using this parameter.
#'      
#'      subst.model : Use this parameter to specify the substitution model, selecting from the list used by
#'      {ape}
#'      
#'      \itemize{
#'      \item
#'       raw: This is simply the proportion or the number of sites that differ between each pair of sequences. This may be useful to draw “saturation plots”. The options variance and gamma have no effect, but pairwise.deletion can.
#'      \item
#'       TS, TV: These are the numbers of transitions and transversions, respectively.
#'      \item
#'       JC69: This model was developed by Jukes and Cantor (1969). It assumes that all substitutions (i.e. a change of a base by another one) have the same probability. This probability is the same for all sites along the DNA sequence. This last assumption can be relaxed by assuming that the substition rate varies among site following a gamma distribution which parameter must be given by the user. By default, no gamma correction is applied. Another assumption is that the base frequencies are balanced and thus equal to 0.25.
#'      \item
#'       K80: The distance derived by Kimura (1980), sometimes referred to as “Kimura's 2-parameters distance”, has the same underlying assumptions than the Jukes–Cantor distance except that two kinds of substitutions are considered: transitions (A <-> G, C <-> T), and transversions (A <-> C, A <-> T, C <-> G, G <-> T). They are assumed to have different probabilities. A transition is the substitution of a purine (C, T) by another one, or the substitution of a pyrimidine (A, G) by another one. A transversion is the substitution of a purine by a pyrimidine, or vice-versa. Both transition and transversion rates are the same for all sites along the DNA sequence. Jin and Nei (1990) modified the Kimura model to allow for variation among sites following a gamma distribution. Like for the Jukes–Cantor model, the gamma parameter must be given by the user. By default, no gamma correction is applied.
#'      \item
#'       F81: Felsenstein (1981) generalized the Jukes–Cantor model by relaxing the assumption of equal base frequencies. The formulae used in this function were taken from McGuire et al. (1999).
#'      \item
#'       K81: Kimura (1981) generalized his model (Kimura 1980) by assuming different rates for two kinds of transversions: A <-> C and G <-> T on one side, and A <-> T and C <-> G on the other. This is what Kimura called his “three substitution types model” (3ST), and is sometimes referred to as “Kimura's 3-parameters distance”.
#'      \item
#'       F84: This model generalizes K80 by relaxing the assumption of equal base frequencies. It was first introduced by Felsenstein in 1984 in Phylip, and is fully described by Felsenstein and Churchill (1996). The formulae used in this function were taken from McGuire et al. (1999).
#'      \item
#'       BH87: Barry and Hartigan (1987) developed a distance based on the observed proportions of changes among the four bases. This distance is not symmetric.
#'      \item
#'       T92: Tamura (1992) generalized the Kimura model by relaxing the assumption of equal base frequencies. This is done by taking into account the bias in G+C content in the sequences. The substitution rates are assumed to be the same for all sites along the DNA sequence.
#'      \item
#'       TN93: Tamura and Nei (1993) developed a model which assumes distinct rates for both kinds of transition (A <-> G versus C <-> T), and transversions. The base frequencies are not assumed to be equal and are estimated from the data. A gamma correction of the inter-site variation in substitution rates is possible.
#'      \item
#'       GG95: Galtier and Gouy (1995) introduced a model where the G+C content may change through time. Different rates are assumed for transitons and transversions.
#'      \item
#'       logdet: The Log-Det distance, developed by Lockhart et al. (1994), is related to BH87. However, this distance is symmetric. Formulae from Gu and Li (1996) are used. dist.logdet in phangorn uses a different implementation that gives substantially different distances for low-diverging sequences.
#'      \item
#'       paralin: Lake (1994) developed the paralinear distance which can be viewed as another variant of the Barry–Hartigan distance.
#'      }
#'      
#'      pairwise.missing : If TRUE, then missing values in the sequence (NNNs) will be accommodated in the calculations pair of taxa at a time; otherwise, the 
#'      deletion of data at positions in the sequence will be global (deleted if any missing data at the position in any individual).
#'      
#' @author Custodian: Arthur Georges -- Post to 
#' \url{https://groups.google.com/d/forum/dartr}
#' 
#' @examples
#' 
#' \dontrun{
#' tmp <- gl.filter.monomorphs(testset.gl)
#' gl.dist.phylo(xx=tmp,subst.model="F80")
#' }
#' 
#' @import ape
#' 
#' @export
#' @return The distance matrix as an object of class dist.

gl.dist.phylo <- function(xx,
                      subst.model="F81",
                      min.tag.len=NULL,
                      pairwise.missing=TRUE,
                      by.pop=TRUE,
                      verbose = NULL) {
  
  # SET VERBOSITY
  verbose <- gl.check.verbosity(verbose)
  
  # FLAG SCRIPT START
  funname <- match.call()[[1]]
  utils.flag.start(func = funname,
                   build = "v2023.3",
                   verbose = verbose)
  
  # STANDARD ERROR CHECKING
  # Check datatype
  datatype <- utils.check.datatype(xx, verbose = verbose)
  
  # Check for monomorphic loci
  tmp <- gl.filter.monomorphs(xx, verbose = 0)
  if ((nLoc(tmp) < nLoc(xx))) {
    if(verbose >= 2){cat(warn("  Warning: genlight object contains monomorphic loci\n"))}
  }
  
 # DEFINE FUNCTIONS
  
  avg.dist <- function(
    gl,
    dist){
    # A function to collapse a distance matrix to populations (by averaging)
    mat <- as.matrix(dist)
    # Number of populations or OTUs
    n_pop = nPop(gl)
    map <- pop(gl)
    names(map)<- indNames(gl)
    # Initialize an empty matrix to store aggregated distances
    aggregated_matrix = matrix(NA, nrow = n_pop, ncol = n_pop)
    rownames(aggregated_matrix) = colnames(aggregated_matrix) = popNames(gl)
    # Iterate over each pair of populations/OTUs
    for (i in 1:(n_pop - 1)) {
      for (j in (i + 1):n_pop) {
        # Indices of individuals in each population/OTU
        inds_i = which(map == rownames(aggregated_matrix)[i])
        inds_j = which(map == colnames(aggregated_matrix)[j])
        # Mean distance between populations/OTUs
        aggregated_matrix[i, j] = aggregated_matrix[j, i] = mean(mat[inds_i, inds_j], na.rm = TRUE)
      }
    }
    dist <- as.dist(aggregated_matrix)
    # dist <- as.matrix(dist)
    # rownames(dist) <- NULL
    # colnames(dist) <- NULL
    # names <- sprintf("%-10s", popNames(gl))
    # dist <- cbind(names,dist)
    return(dist)
  }
  
  # DO THE JOB

  if(!is.null(min.tag.len)){
    if(verbose >= 2){cat(report("  Filtering sequence tags with length less than",min.tag.len,"\n"))}
    xx <- gl.filter.taglength(xx,lower=min.tag.len,verbose=0)
  }
  
  # Create the sequences in a form amenable to analysis by ape
  hold <- getwd()
  setwd(tempdir())
  if(verbose >= 2){cat(report("  Converting sequence tags to input format for package {ape}\n"))}
  gl2fasta(xx,outfile="tmp.fas",outpath = tempdir(),verbose=0)
    sequences <- ape::read.dna("tmp.fas", format = "fasta")
  setwd(hold)  
  
  # Calculate distances between individuals
  if(verbose >= 2){
    cat(report("  Calculating distances between individuals\n"))
    cat(report("    Substitution model:",subst.model,"\n"))
    cat(report("    Pairwise missing value deletion:",pairwise.missing,"\n"))
  }
  D <- ape::dist.dna(sequences,model=subst.model,pairwise.deletion = pairwise.missing)
  
  if(by.pop){
    #Calculate average distances for pairwise populations
    if(verbose >= 2){cat(report("  Calculating average distances between populations\n"))}
    D <- avg.dist(gl=xx,dist=D)
    
  #   hold <- getwd()
  #   setwd(tempdir())
  #     write(nPop(xx),file="infile")
  #     write.table(D,file="infile",row.names=FALSE,col.names=FALSE,quote=FALSE,append=TRUE)
  #   setwd(hold)
  # } else {
  #   hold <- getwd()
  #   setwd(tempdir())
  #     write(nInd(xx),file="infile")
  #     write.table(D,file="infile",row.names=FALSE,col.names=FALSE,quote=FALSE,append=TRUE)
  #   setwd(hold)
  }
  
  # Convert to a dist object
  # D <- D[,-1]
  # rownames(D) <- popNames(xx)
  # colnames <- rownames
  # D <- as.dist(D)
  
  # FLAG SCRIPT END
  
  if (verbose > 0) {
    cat(report("Completed:", funname, "\n"))
  }

  return(D)
}
