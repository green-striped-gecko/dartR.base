#'@name gl.pcoa
#'@title Ordination applied to genotypes in a genlight object (PCA), in an fd
#'object, or to a distance matrix (PCoA)
#'@description
#' This function takes the genotypes for individuals and undertakes a Pearson
#' Principal Component analysis (PCA) on SNP or Sequence tag P/A (SilicoDArT) data; it
#' undertakes a Gower Principal Coordinate analysis (PCoA) if supplied with a
#' distance matrix. 

#' @param x Name of the genlight object or fd object containing the SNP data, or
#'  a distance matrix of type dist [required].
#' @param nfactors Number of axes to retain in the output of factor scores
#' [default 5].
#' @param pc.select Method for identifying substantial PC axes. One of Kaiser-Guttman,
#' broken-stick, Tracy-Widom [default broken-stick]
#' @param correction Method applied to correct for negative eigenvalues, either
#'  'lingoes' or 'cailliez' [Default NULL].
#' @param mono.rm If TRUE, remove monomorphic loci [default TRUE].
#' @param parallel TRUE if parallel processing is required (does fail under
#' Windows) [default FALSE].
#' @param n.cores Number of cores to use if parallel processing is requested
#' [default 16].
#' @param plot.out If TRUE, a diagnostic plot is displayed showing a scree plot
#' for the "informative" axes and a histogram of eigenvalues of the remaining 
#' "noise" axes [Default TRUE].
#' @param plot.theme Theme for the plot. See Details for options
#'  [default theme_dartR()].
#' @param plot.colors List of two color names for the borders and fill of the
#'  plot [default gl.colors(2)].
#' @param plot.dir Directory to save the plot RDS files [default as specified 
#' by the global working directory or tempdir()]
#' @param plot.file Name for the RDS binary file to save (base name only, exclude extension) [default NULL]
#' @param verbose verbose= 0, silent or fatal errors; 1, begin and end; 2,
#' progress log; 3, progress and results summary; 5, full report
#' [default 2 or as specified using gl.set.verbosity].

#' @details
#'The function is essentially a wrapper for glPca \{adegenet\} or pcoa \{ape\}
#'with default settings apart from those specified as parameters in this
#'function.
#'\strong{ Sources of stress in the visual representation }
#' While, technically, any distance matrix can be represented in an ordinated
#' space, the representation will not typically be exact.There are three major
#' sources of stress in a reduced-representation of distances or dissimilarities
#' among entities using PCA or PCoA. By far the greatest source comes from the
#' decision to select only the top two or three axes from the ordinated set of
#' axes derived from the PCA or PCoA. The representation of the entities such a
#' heavily reduced space will not faithfully represent the distances in the
#' input distance matrix simply because of the loss of information in deeper
#' informative dimensions. For this reason, it is not sensible to be too
#' precious about managing the other two sources of stress in the visual
#' representation.

#' The measure of distance between entities in a PCA is the Pearson Correlation
#' Coefficient, essentially a standardized Euclidean distance. This is both a
#' metric distance and a Euclidean distance and so the distances in the final 
#' ordination are faithful to those between entities in
#' the dataset. Note that missing values are imputed in PCA, and that this can be
#' a source of disparity between the distances between the entities in the dataset
#' and the distances in the ordinated space.
#' 
#' In PCoA, the second source of stress is the choice of distance 
#' measure or dissimilarity measure. While pretty much any distance or dissimilarity matrix 
#' can be represented in an ordinated space, the distances between entities can 
#' be faithfully represented in that space (that is, without stress) only if the 
#' distances are metric. Furthermore, for distances between entities to be 
#' faithfully represented in a rigid Cartesian space, the distance measure 
#' needs to be Euclidean. If this is not the case, the distances between the 
#' entities in the ordinated visualized space will not #' exactly represent 
#' the distances in the input matrix (stress will be non-zero).
#' This source of stress will be evident as negative eigenvalues in the deeper
#' dimensions. 
#' 
#' A third source of stress arises from having a sparse dataset, one with
#' missing values. This affects both PCA and PCoA. If the original data matrix
#' is not fully populated, that is, if there are missing values, then even a
#' Euclidean distance matrix will not necessarily be 'positive definite'. It
#' follows that some of the eigenvalues may be negative, even though the
#' distance metric is Euclidean. This issue is exacerbated when the number of
#' loci greatly exceeds the number of individuals, as is typically the case when
#' working with SNP data. The impact of missing values can be minimized by
#' stringently filtering on Call Rate, albeit with loss of data. An alternative
#' is given in a paper 'Honey, I shrunk the sample covariance matrix' and more
#' recently by Ledoit and Wolf (2018), but their approach has not been
#'  implemented here.
#'  
#'  Options for imputing missing values while minimizing distortion are provided in 
#'  the function gl.impute().

#' The good news is that, unless the sum of the negative eigenvalues, arising
#' from a non-Euclidean distance measure or from missing values, approaches
#' those of the final PCA or PCoA axes to be displayed, the distortion is
#' probably of no practical consequence and certainly not comparable to the
#' stress arising from selecting only two or three final dimensions out of
#' several informative dimensions for the visual representation.

#'\strong{ Function's output }

#' Two diagnostic plots are produced. The first is a Scree Plot, showing the
#' percentage variation explained by each of the PCA or PCoA axes, for those
#' axes are considered informative. The scree plot informs a decision on the number of
#' dimensions to be retained in the visual summaries. Various approaches are available
#' for identifying which axes are informative (in terms of containing biologically
#' significant variation) and which are noise axes. The simplest method is to consider
#' only those axes that explain more variance than the original variables on average as
#' being informative (pc.select="Kaiser-Guttman"). A second method (the default) is
#' the broken-stick method (pc.select="broken-stick"). A third method is the Tracy-Widom
#' statistical approach (pc.select="Tracy-Widom").
#' 
#' Once you have the informative axes identified, a judgement call is made as to how
#' many dimensions to retain and present as results. This requires a decision on how 
#' much information on structure in the data is to be discarded. Retaining at least those axes 
#' that explain 10% or more of variation is used by some as a rule of thumb.

#' The second graph is for diagnostic purposes only. It shows the distribution of 
#' eigenvalues for the remaining uninformative (noise) axes, including those with 
#' negative eigenvalues. If a plot.file is given, the ggplot arising from this function 
#' is saved as an binary file using gl.save(); can be reloaded with gl.load(). 
#' A file name must be  specified for the plot to be saved.
#' If a plot directory (plot.dir) is specified, the ggplot binary is saved to that
#' directory; otherwise to the tempdir(). 
#' Action is recommended (verbose >= 2) if the negative eigenvalues are
#' dominant, their sum approaching in magnitude the eigenvalues for axes
#' selected for the final visual solution.

#' Output is a glPca object conforming to adegenet::glPca but with only the
#' following retained.
#'\itemize{
#'\item  $call - The call that generated the PCA/PCoA
#'\item  $eig - Eigenvalues -- All eigenvalues (positive, null, negative).
#'\item  $scores - Scores (coefficients) for each individual
#'\item  $loadings - Loadings of each SNP for each principal component
#'    }
#'
#'
#'  Examples of other themes that can be used can be consulted in \itemize{
#'  \item \url{https://ggplot2.tidyverse.org/reference/ggtheme.html} and \item
#'  \url{https://yutannihilation.github.io/allYourFigureAreBelongToUs/ggthemes/}
#'  }
#'
#' PCA was developed by Pearson (1901) and Hotelling (1933), whilst the best
#' modern reference is Jolliffe (2002). PCoA was developed by Gower (1966) while
#' the best modern reference is Legendre & Legendre (1998).

#'@return An object of class pcoa containing the eigenvalues and factor scores
#'@author Author(s): Arthur Georges and Jesus Castrejon. Custodian: Arthur Georges (Post to
#'\url{https://groups.google.com/d/forum/dartr})
#'@examples
#' # PCA (using SNP genlight object)
#' gl <- possums.gl
#' pca <- gl.pcoa(possums.gl[1:50,],verbose=2)
#' gl.pcoa.plot(pca,gl)
#' \donttest{
#' gs <- testset.gs
#' levels(pop(gs))<-c(rep('Coast',5),rep('Cooper',3),rep('Coast',5),
#' rep('MDB',8),rep('Coast',6),'Em.subglobosa','Em.victoriae')
#' # PCA (using SilicoDArT genlight object)
#' pca <- gl.pcoa(gs)
#' gl.pcoa.plot(pca,gs)
#' # Using a distance matrix
#' D <- gl.dist.ind(testset.gs, method='jaccard')
#' pcoa <- gl.pcoa(D,correction="cailliez")
#' gl.pcoa.plot(pcoa,gs)
#' }
#'@references
#'\itemize{
#'\item Cailliez, F. (1983) The analytical solution of the additive constant
#'problem. Psychometrika, 48, 305-308.
#'\item Gower, J. C. (1966) Some distance properties of latent root and vector
#'methods used in multivariate analysis. Biometrika, 53, 325-338.
#'\item Hotelling, H., 1933. Analysis of a complex of statistical variables into
#' Principal Components. Journal of Educational Psychology 24:417-441, 498-520.
#'\item Jolliffe, I. (2002) Principal Component Analysis. 2nd Edition, Springer,
#' New York.
#'\item Ledoit, O. and Wolf, M. (2018). Analytical nonlinear shrinkage of
#'large-dimensional covariance matrices. University of Zurich, Department of
#'Economics, Working Paper No. 264, Revised version. Available at SSRN:
#'https://ssrn.com/abstract=3047302 or http://dx.doi.org/10.2139/ssrn.3047302
#'\item Legendre, P. and Legendre, L. (1998). Numerical Ecology, Volume 24, 2nd
#'Edition. Elsevier Science, NY.
#'\item Lingoes, J. C. (1971) Some boundary conditions for a monotone analysis
#'of symmetric matrices. Psychometrika, 36, 195-203.
#'\item Pearson, K. (1901). On lines and planes of closest fit to systems of
#'points in space. Philosophical Magazine. Series 6, vol. 2, no. 11, pp.
#'559-572.
#' }
#' @seealso \code{\link{gl.pcoa.plot}}
#' @family data exploration functions
#' @importFrom ape pcoa
#' @export

gl.pcoa <- function(x,
                    nfactors = 5,
                    pc.select="broken-stick",
                    correction = NULL,
                    mono.rm = TRUE,
                    parallel = FALSE,
                    n.cores = 1,
                    plot.out = TRUE,
                    plot.theme = theme_dartR(),
                    plot.colors = gl.colors(2),
                    plot.file=NULL,
                    plot.dir=NULL,
                    verbose = NULL) {
  
    # SET VERBOSITY
    verbose <- gl.check.verbosity(verbose)
    
    # FLAG SCRIPT START
    funname <- match.call()[[1]]
    utils.flag.start(func = funname,
                     build = "2024V1",
                     verbose = verbose)
    
    # SET WORKING DIRECTORY
    plot.dir <- gl.check.wd(plot.dir,verbose=0)
    
    # CHECK DATATYPE
    datatype <-
        utils.check.datatype(x,
                             accept = c("SNP", "SilicoDArT", "dist", "fd"),
                             verbose = verbose)
    
    if (datatype == "SNP" | datatype == "SilicoDArT") {
        x <- gl.filter.allna(x, verbose = 0)
    }
    
    if (datatype == "fd") {
        datatype <- utils.check.datatype(x$gl, verbose = verbose)
        x <- x$gl
        x <- gl.filter.allna(x, verbose = 0)
    }
    
    # SCRIPT SPECIFIC ERROR CHECKING
    
    if (mono.rm == TRUE &
        (datatype == "SNP" | datatype == "SilicoDArT")) {
        x <- gl.filter.monomorphs(x, verbose = 0)
    }
    
    if((datatype == "SNP" | datatype == "SilicoDArT")){
        if (nInd(x) < 2) {
        stop(
            error(
                "Fatal Error: Only one individual or no individuals present in the dataset\n"
            )
        )
        }
   
        if (nLoc(x) < nInd(x)) {
        if(verbose >=2){cat(
            warn(
                "  Warning: Number of loci is less than the number of individuals to be represented\n"
            )
        )}
        }
    }
    
    if (is.null(correction)) {
        correction <- "none"
    } else {
        correction <- tolower(correction)
        if (correction != "lingoes" &&
            correction != "cailliez") {
            if (verbose >= 2) {
                cat(
                    warn(
                        "  Warning: Correction if specified needs to be lingoes or cailliez, set to the default 'None'\n"
                    )
                )
            }
            correction <- "none"
        }
    }
    
    # if(pc.select=="Tracy-Widom"){
    #   cat(warn("  Note: The Tracy-Widom Criterion for determining informative axes is not yet implemented. Selecting 'broken-stick criterion'\n"))
    #   pc.select <- "broken-stick"
    # }
    
    # FUNCTIONS
    
    twtest <- function(x){
      # Returns the p-value of the Tracy Widom statistics 
      # (righ-side cumulative function). The function uses linear interpolation 
      # over the tables reported originally by Patterson et al (2006).
      
      # The package RMTstat provides the same values (1-ptw(x)) with more precision but for smaller values
      
      # twtest(0.6)
      
      # Author: Jesus Castrejon
      
      # This array holds a pre-calculated table for the Tracy-Widom distribution.
      # Each group of three numbers corresponds to a value `x` from the distribution,
      # the density at `x`, and the cumulative distribution up to `x`.
      twtable <- c(-8.000000000,1.000000000,0.000000000,-7.900000000,1.000000000,0.000000000,-7.800000000,1.000000000,0.000000000,-7.700000000,1.000000000,0.000000000,-7.600000000,1.000000000,0.000000000,
                   -7.500000000,1.000000000,0.000000001,-7.400000000,1.000000000,0.000000002,-7.300000000,0.999999999,0.000000005,-7.200000000,0.999999999,0.000000010,-7.100000000,0.999999997,0.000000019,
                   -7.000000000,0.999999995,0.000000039,-6.900000000,0.999999989,0.000000076,-6.800000000,0.999999978,0.000000146,-6.700000000,0.999999958,0.000000276,-6.600000000,0.999999920,0.000000511,
                   -6.500000000,0.999999849,0.000000932,-6.400000000,0.999999723,0.000001670,-6.300000000,0.999999498,0.000002942,-6.200000000,0.999999105,0.000005097,-6.100000000,0.999998431,0.000008683,
                   -6.000000000,0.999997293,0.000014554,-5.900000000,0.999995401,0.000024005,-5.800000000,0.999992309,0.000038969,-5.700000000,0.999987331,0.000062279,-5.600000000,0.999979441,0.000098012,
                   -5.500000000,0.999967125,0.000151923,-5.400000000,0.999948187,0.000231995,-5.300000000,0.999919496,0.000349097,-5.200000000,0.999876655,0.000517756,-5.100000000,0.999813597,0.000757035,
                   -5.000000000,0.999722082,0.001091485,-4.900000000,0.999591101,0.001552137,-4.800000000,0.999406175,0.002177466,-4.700000000,0.999148569,0.003014256,-4.600000000,0.998794427,0.004118267,
                   -4.500000000,0.998313849,0.005554591,-4.400000000,0.997669962,0.007397591,-4.300000000,0.996818016,0.009730295,-4.200000000,0.995704571,0.012643159,-4.100000000,0.994266851,0.016232112,
                   -4.000000000,0.992432322,0.020595851,-3.900000000,0.990118582,0.025832397,-3.800000000,0.987233631,0.032034971,-3.700000000,0.983676579,0.039287325,-3.600000000,0.979338843,0.047658716,
                   -3.500000000,0.974105853,0.057198759,-3.400000000,0.967859270,0.067932445,-3.300000000,0.960479677,0.079855636,-3.200000000,0.951849687,0.092931337,-3.100000000,0.941857369,0.107087044,
                   -3.000000000,0.930399881,0.122213418,-2.900000000,0.917387157,0.138164458,-2.800000000,0.902745495,0.154759279,-2.700000000,0.886420892,0.171785501,-2.600000000,0.868381957,0.189004169,
                   -2.500000000,0.848622271,0.206156009,-2.400000000,0.827162053,0.222968755,-2.300000000,0.804049066,0.239165233,-2.200000000,0.779358684,0.254471803,-2.100000000,0.753193114,0.268626779,
                   -2.000000000,0.725679802,0.281388431,-1.900000000,0.696969061,0.292542221,-1.800000000,0.667231036,0.301906945,-1.700000000,0.636652122,0.309339558,-1.600000000,0.605430961,0.314738516,
                   -1.500000000,0.573774198,0.318045543,-1.400000000,0.541892124,0.319245849,-1.300000000,0.509994383,0.318366852,-1.200000000,0.478285870,0.315475570,-1.100000000,0.446962951,0.310674866,
                   -1.000000000,0.416210105,0.304098784,-0.900000000,0.386197065,0.295907232,-0.800000000,0.357076521,0.286280263,-0.700000000,0.328982392,0.275412215,-0.600000000,0.302028689,0.263505933,
                   -0.500000000,0.276308949,0.250767272,-0.400000000,0.251896179,0.237400053,-0.300000000,0.228843301,0.223601597,-0.200000000,0.207183986,0.209558915,-0.100000000,0.186933854,0.195445624,
                   0.000000000,0.168091934,0.181419571,0.100000000,0.150642330,0.167621190,0.200000000,0.134556018,0.154172511,0.300000000,0.119792709,0.141176787,0.400000000,0.106302721,0.128718659,
                   0.500000000,0.094028817,0.116864772,0.600000000,0.082907953,0.105664756,0.700000000,0.072872924,0.095152500,0.800000000,0.063853860,0.085347620,0.900000000,0.055779577,0.076257058,
                   1.000000000,0.048578763,0.067876743,1.100000000,0.042180992,0.060193257,1.200000000,0.036517582,0.053185457,1.300000000,0.031522284,0.046826015,1.400000000,0.027131832,0.041082856,
                   1.500000000,0.023286351,0.035920459,1.600000000,0.019929640,0.031301023,1.700000000,0.017009350,0.027185487,1.800000000,0.014477062,0.023534398,1.900000000,0.012288293,0.020308645,
                   2.000000000,0.010402429,0.017470054,2.100000000,0.008782605,0.014981856,2.200000000,0.007395547,0.012809046,2.300000000,0.006211384,0.010918644,2.400000000,0.005203434,0.009279861,
                   2.500000000,0.004347977,0.007864200,2.600000000,0.003624031,0.006645482,2.700000000,0.003013114,0.005599836,2.800000000,0.002499018,0.004705636,2.900000000,0.002067590,0.003943413,
                   3.000000000,0.001706520,0.003295741,3.100000000,0.001405143,0.002747112,3.200000000,0.001154255,0.002283795,3.300000000,0.000945945,0.001893694,3.400000000,0.000773431,0.001566204,
                   3.500000000,0.000630927,0.001292071,3.600000000,0.000513508,0.001063253,3.700000000,0.000416999,0.000872795,3.800000000,0.000337871,0.000714702,3.900000000,0.000273152,0.000583831,
                   4.000000000,0.000220344,0.000475784,4.100000000,0.000177359,0.000386816,4.200000000,0.000142452,0.000313749,4.300000000,0.000114170,0.000253894,4.400000000,0.000091308,0.000204987,
                   4.500000000,0.000072871,0.000165125,4.600000000,0.000058035,0.000132716,4.700000000,0.000046124,0.000106431,4.800000000,0.000036582,0.000085163,4.900000000,0.000028955,0.000067996,
                   5.000000000,0.000022872,0.000054172,5.100000000,0.000018030,0.000043066,5.200000000,0.000014185,0.000034164,5.300000000,0.000011138,0.000027045,5.400000000,0.000008728,0.000021365,
                   5.500000000,0.000006826,0.000016843,5.600000000,0.000005328,0.000013250,5.700000000,0.000004151,0.000010403,5.800000000,0.000003228,0.000008151,5.900000000,0.000002505,0.000006374,
                   6.000000000,0.000001941,0.000004974,6.100000000,0.000001501,0.000003874,6.200000000,0.000001158,0.000003011,6.300000000,0.000000892,0.000002336,6.400000000,0.000000686,0.000001809,
                   6.500000000,0.000000527,0.000001398,6.600000000,0.000000403,0.000001078,6.700000000,0.000000308,0.000000830,6.800000000,0.000000235,0.000000638,6.900000000,0.000000179,0.000000489,
                   7.000000000,0.000000136,0.000000375,7.100000000,0.000000104,0.000000286,7.200000000,0.000000079,0.000000218,7.300000000,0.000000059,0.000000166,7.400000000,0.000000045,0.000000126,
                   7.500000000,0.000000034,0.000000096,7.600000000,0.000000025,0.000000073,7.700000000,0.000000019,0.000000055,7.800000000,0.000000014,0.000000041,7.900000000,0.000000011,0.000000031,
                   8.000000000,0.000000008,0.000000023)
      # If the input `x` is 8 or greater, the function returns 0 because the table 
      # does not include values at or beyond this point.
      if(x >= 8){
        return(0)
      }
      
      i <- 0
      # Search through twtable to find the right index where the value of `x` 
      # fits between entries in the table.
      while(i < 162 & x >= twtable[(i)*3+1]){
        i <- i+1
      }
      # lower boundary 
      z0 = twtable[(i-1)*3 + 2]
      # upper boundary
      z1 = twtable[(i)*3 + 2]
      # If at the end of the table, return the density at the lower boundary.
      if(i == 162){
        return(z0)
        # x is less than the smallest x in the table, return the density at 
        # the upper boundary
      } else if(i == 0){
        return(z1)
        # x is within the range of the table, interpolate
      } else{
        # Linear interpolation to find the density at x
        result <- z0 + (z1- z0)*(x - twtable[(i-1)*3+1])/(twtable[(i)*3+1] - twtable[(i-1)*3+1])
      }
      return(result)
    }

    tw.statistics <- function(eigenvalues,alpha=0.05,plot=TRUE) {
      
      # Returns two vectors, one holding the eigenvalues of "informative" dimensions,
      # the other holding the eigenvalues of the "noise" dimensions. Criterion is
      # statistical significance under the Tracy-Widom distribution (implemented by
      #Patterson, 2006.
      
      # Author: Jesus Castrejon
      
      # Reference: Patterson, N., Price, A. L., & Reich, D. (2006). Population 
      # structure and eigenanalysis. PLoS Genetics, 2, e190. 
      # https://doi.org/10.1371/journal.pgen.0020190
      
      # Takes M as the number of eigenvalues
      M <- length(eigenvalues)
      
      # Sum eigenvalues and squared eigenvalues
      s  <- sum(eigenvalues)
      s2 <- sum(eigenvalues**2.)
      
      # Initialise lists to store values
      index <- 1:length(eigenvalues)
      p.value <- c()
      tw.stat <- c()
      
      # Checks TW statistics for each eigenvalue
      for(i in index){
        Mp <- M-i+1
        nhat <- (Mp+2)*s**2./(Mp*s2 - s**2.)
        lambda <- eigenvalues[i]*Mp/s
        # Fix effective n for very small eigenvalues
        if(nhat<1.){
          nhat <- 1.
        }
        mu <- (sqrt(nhat-1.) + sqrt(Mp))**2./nhat
        sigma <- (sqrt(nhat-1.) + sqrt(Mp))*(1./sqrt(nhat-1) + 1./sqrt(Mp))**(1./3)/nhat
        # TW statistic
        x <- (lambda-mu)/sigma
        # Gets pvalue of the Tracy-Widom statistics
        pvalue <- twtest(x)
        # Updates sum of eigenvalues
        s <- s - eigenvalues[i]
        s2 <- s2 - eigenvalues[i]**2.
        p.value[i] <- pvalue
        tw.stat[i] <- x
      }
      
      # Adjust pvalues
      p.value <- stats::p.adjust(p.value)
      
      # Creates dataframe with eigenvalues and filter into noisy/structure base on pvalue>alpha
      df <- data.frame(index,eigenvalues,tw.stat,p.value)
      idx <- which(p.value>alpha)[1]
      struc <- df[1:(idx-1),]
      noise <- df[(idx):length(eigenvalues),]
      struc$structure <- 'structured'
      noise$structure <- 'noisy'
      
      # Plots eigenvalues
      if(plot){
        p1 <- ggplot() +  
          geom_point(data=struc, aes(x=index,y = eigenvalues,color = structure)) +
          geom_point(data=noise, aes(x=index,y = eigenvalues,color = structure)) +
          scale_y_continuous(trans='log10') + 
          scale_color_manual(values=c("blue","red", "green"))+
          theme(legend.title=element_blank())
        print(p1)
      }
      return(list('struct'=struc,'noise'=noise))
    }
    
    bs.statistics <- function(eigenvalues,plot=FALSE){
      # Returns two vectors, one holding the eigenvalues of "informative" dimensions,
      # the other holding the eigenvalues of the "noise" dimensions. Criterion is
      # based on the broken-stick algorithm (Macarthur, 1957).
      
      # Author: Jesus Castrejon
      
      # Reference: Macarthur, R. (1957). On the relative abundance of bird species. 
      # Proceedings of the National Academy of Sciences USA, 43, 293–295. 
      # https://doi.org/10.1073/pnas.43.3.293
      
      # Set number of eigenvalues
      n <- length(eigenvalues)
      index<-1:n
      
      # Calculate expected lengths under Broken stick approach
      broken.sticks <- sum(eigenvalues)*rev(cumsum(1/n:1))/n
      
      # Creates dataframe with eigenvalues and filter noise/structure based on BS lengths
      df <- data.frame(index,eigenvalues,broken.sticks)
      idx <- which(eigenvalues<broken.sticks)[1]
      struc <- df[1:idx-1,]
      noise <- df[idx:length(eigenvalues),]
      struc$structure <- 'structured'
      noise$structure <- 'noisy'
      
      # Creates Plot of eigenvalues
      if(plot){
        p1 <- ggplot() +  
          geom_point(data=struc, aes(x=index,y = eigenvalues,color = structure)) +
          geom_point(data=noise, aes(x=index,y = eigenvalues,color = structure)) +
          geom_line(data=df, aes(x=index,y = broken.sticks,color = 'broken sticks'),linetype="dashed") +
          scale_y_continuous(trans='log10') + 
          scale_color_manual(values=c("black","blue", "red")) +
          theme(legend.title=element_blank())
        print(p1)
      }
      
      return(list('struct'=struc,'noise'=noise))
    }
    
    # DO THE JOB
    
    ######## DISTANCE ANALYSIS
    
    if (datatype == "dist")
    {
        D <- x
        
        # Calculate the pcoa
        if (verbose >= 2) {
            if (correction == "none") {
                cat(
                    report(
                        "  Performing a PCoA, individuals as entities, no correction applied\n"
                    )
                )
                title <-
                    paste0("PCoA on Distance Matrix (no correction)\nScree Plot\n (informative axes only -- ",pc.select,"criterion)")
            } else {
                cat(
                    report(
                        "  Performing a PCoA, individuals as entities,",
                        correction,
                        "correction for stress (-ive eigenvalues) applied\n"
                    )
                )
                title <-
                    paste0(
                        "PCoA on Distance Matrix(",
                        correction,
                        ")\nScree Plot (informative axes only -- ",pc.select,"criterion)")
            }
        }
        pco <-ape::pcoa(D, correction = correction, rn = labels(D))
        
        # Extract relevant variables
        
        if (correction == "none") {
            eig.raw <- pco$values$Eigenvalues
        } else {
            eig.raw <- pco$values$Corr_eig
        }
        # M=length(eig.raw)+1
        
        # Identify the number of axes with explanatory value
        
        if(pc.select=="Kaiser-Guttman"){
          eig.raw.pos <- eig.raw[eig.raw >= 0]
          eig.raw.pos.pc <- eig.raw.pos * 100 / sum(eig.raw.pos)
          eig.top <- eig.raw.pos[eig.raw.pos > mean(eig.raw.pos)] 
          eig.top.pc <- round(eig.top * 100 / sum(eig.raw.pos), 1)
          eig.raw.noise <- eig.raw[eig.raw <= mean(eig.raw)]
        } else if(pc.select=="broken-stick"){
          eig.raw.pos <- eig.raw[eig.raw >= 0]
          tmp <- bs.statistics(eig.raw.pos)
          eig.top <- tmp$struct$eigenvalues
          eig.top.pc <- round(eig.top * 100 / sum(eig.raw.pos), 1)
          eig.raw.noise <- tmp$noise$eigenvalues
        } else if(pc.select=="Tracy-Widom"){
          # M <- nInd(x)
          eig.raw.pos <- eig.raw[eig.raw >= 0]
          tmp <- tw.statistics(eig.raw.pos)
          eig.top <- tmp$struct$eigenvalues
          eig.top.pc <- round(eig.top * 100 / sum(eig.raw.pos), 1)
          eig.raw.noise <- tmp$noise$eigenvalues
        } else {
          stop("Error: incorrect specification for model of significant eigenvalues. Set to Tracy-Widom\n")
          pc.select <- "Tracy-Widom"
        }  
        
        if (any(eig.raw < 0)) {
            if (verbose >= 2) {
                problem <- (-sum(eig.raw[eig.raw < 0]) / mean(eig.raw[1:3])) * 100
                cat(
                    warn(
                        "  Warning: Some eigenvalues negative -- sum to",
                        round(problem, 2),
                        "% of the mean eigenvalue for PCoA axes 1-3\n"
                    )
                )
                cat(
                    report(
                        "    Tolerable negative eigenvalues should sum to much less than the eigenvalues of displayed PCoA axes (say, less than 20%)\n"
                    )
                )
                if (problem > 20) {
                    cat(
                        report(
                            "    If the stress (negative eigenvalues) is considered a problem, and you might reasonably choose to ignore it, you have the following options:\n"
                        )
                    )
                    cat(
                        report(
                            "    (a) Apply more stringent filtering on Call Rate before generating the distance matrix; or\n"
                        )
                    )
                    cat(
                        report(
                            "    (b) Select an alternate distance measure, preferably a metric distance or better still, a Euclidean distance, if you have not already; or\n"
                        )
                    )
                    cat(
                        report(
                            "    (c) Apply a transformation (correction) to eliminate the negative eigenvalues. If this was already done, try another correction; or\n"
                        )
                    )
                    cat(
                        report(
                            "    (d) Interpret the visual representation of the ordination with caution, seeking corroborating evidence.\n"
                        )
                    )
                }
            }
        }
        
        # Provide a summary
        if (verbose >= 3) {
            if (correction == "lingoes" | correction == "cailliez") {
                cat(" Correction",
                    correction,
                    "applied to remove negative eigenvalues\n")
                cat(
                    paste(
                        "  Uncorrected ordination yielded",
                        length(eig.top),
                        "informative dimensions(",pc.select,"criterion) from",
                        nInd(x) - 1,
                        "original dimensions\n"
                    )
                )
            } else {
                cat(
                    paste(
                        "  Ordination yielded",
                        length(eig.top),
                        "informative dimensions(",pc.select,"criterion) from",
                        dim(as.matrix(D))[1],
                        "original dimensions\n"
                    )
                )
            }
            cat(paste(
                "    PCoA Axis 1 explains",
                round(eig.raw.pos.pc[1], 1),
                "% of the total variance\n"
            ))
            cat(
                paste(
                    "    PCoA Axis 1 and 2 combined explain",
                    round(eig.raw.pos.pc[1] + eig.raw.pos.pc[2], 1),
                    "% of the total variance\n"
                )
            )
            cat(paste(
                "    PCoA Axis 1-3 combined explain",
                round(
                    eig.raw.pos.pc[1] + eig.raw.pos.pc[2] + eig.raw.pos.pc[3],
                    1
                ),
                "% of the total variance\n"
            ))
        }
        # Construct a universal output file
        p.object <- list()
        p.object$scores <- pco$vectors[, 1:nfactors]
        p.object$eig <- pco$values$Eigenvalues
        p.object$loadings <- pco$vectors.cor[, 1:nfactors]
        p.object$call <- match.call()
        
    }  ######## END DISTANCE DATA
    
    ######## SNP or P/A DATA, PCA
    
    if (datatype == "SNP" || datatype == "SilicoDArT")
    {
        if (datatype == "SNP") {
            if (verbose >= 2) {cat(
                    report(
                        "  Performing a PCA, individuals as entities, loci as attributes, SNP genotype as state\n"
                    )
                )}
                title <-
                    paste0("PCA on SNP Genotypes\nScree Plot\n (informative axes only -- ",pc.select," criterion)")
            }
            if (datatype == "SilicoDArT") {
                if (verbose >= 2) {cat(
                    report(
                        "  Performing a PCA, individuals as entities, loci as attributes, Tag P/A as state\n"
                    )
                )}
                title <-
                    paste0("PCA on Tag P/A Data\nScree Plot\n (informative axes only -- ",pc.select," criterion)")
            }
            
        pca <-
            glPca(x,
                  nf = nfactors,
                  parallel = parallel,
                  n.cores = n.cores)
        
        # # Identify the number of axes with explanatory value greater than the original variables on average
        # 
        # 
        # eig.raw.pos <- eig.raw[eig.raw >= 0]
        # eig.raw.pos.pc <- eig.raw.pos * 100 / sum(eig.raw.pos)
        # eig.top <- eig.raw.pos[eig.raw.pos >= mean(eig.raw.pos)]
        # eig.top.pc <- round(eig.top * 100 / sum(eig.raw.pos), 1)
        # eig.raw.noise <- eig.raw[eig.raw <= mean(eig.raw)]
        
        # Extract relevant variables
        eig.raw <- pca$eig

        # Identify the number of axes with explanatory value
        # Identify the number of axes with explanatory value
        
        if(pc.select=="Kaiser-Guttman"){
          eig.raw.pos <- eig.raw[eig.raw >= 0]
          eig.raw.pos.pc <- eig.raw.pos * 100 / sum(eig.raw.pos)
          eig.top <- eig.raw.pos[eig.raw.pos > mean(eig.raw.pos)] 
          eig.top.pc <- round(eig.top * 100 / sum(eig.raw.pos), 1)
          eig.raw.noise <- eig.raw[eig.raw <= mean(eig.raw)]
        } else if(pc.select=="broken-stick"){
          eig.raw.pos <- eig.raw[eig.raw >= 0]
          tmp <- bs.statistics(eig.raw.pos)
          eig.top <- tmp$struct$eigenvalues
          eig.top.pc <- round(eig.top * 100 / sum(eig.raw.pos), 1)
          eig.raw.noise <- tmp$noise$eigenvalues
        } else if(pc.select=="Tracy-Widom"){
          # M <- nInd(x)
          eig.raw.pos <- eig.raw[eig.raw >= 0]
          tmp <- tw.statistics(eig.raw.pos)
          eig.top <- tmp$struct$eigenvalues
          eig.top.pc <- round(eig.top * 100 / sum(eig.raw.pos), 1)
          eig.raw.noise <- tmp$noise$eigenvalues
        } else {
          stop("Error: incorrect specification for model of significant eigenvalues. Set to Tracy-Widom\n")
          pc.select <- "Tracy-Widom"
        }  
        
        if (any(eig.raw < 0) && verbose >= 2) {
          problem <- (-sum(eig.raw[eig.raw < 0]) / mean(eig.raw[1:3])) * 100
          cat(warn("  Warning: Some eigenvalues negative -- they sum to", round(problem, 2), "% of the mean eigenvalue for PCA axes 1-3\n"))
          cat(report("    Tolerable negative eigenvalues should sum to much less than the eigenvalues of displayed PCA axes (say, less than 20%)\n"))
          if (problem > 20) {
            cat(report("    If the stress (negative eigenvalues) is considered a problem, and you might reasonably choose to ignore it, you have the following options:\n"))
            cat(report("    (a) Apply more stringent filtering on Call Rate and repeat the PCA; or\n"))
            cat(report("    (b) Undertake a PCoA with an appropriate distance measure and a transformation (correction) to eliminate the negative eigenvalues; or\n"))
            cat(report("    (c) Interpret the visual representation of the ordination with caution, seeking corroborating evidence.\n"))
          }
        }
        
        e <- pca$eig[pca$eig > sum(pca$eig / length(pca$eig))]
        e <- eig.top
        e <- round(e * 100 / sum(pca$eig), 1)
        if (verbose >= 3) {
          cat(paste("  Ordination yielded", length(e), "informative dimensions(", pc.select, "criterion) from", nInd(x) - 1, "original dimensions\n"))
          cat(paste("    PCA Axis 1 explains", e[1], "% of the total variance\n"))
          cat(paste("    PCA Axis 1 and 2 combined explain", e[1] + e[2], "% of the total variance\n"))
          cat(paste("    PCA Axis 1-3 combined explain", e[1] + e[2] + e[3], "% of the total variance\n"))
        }
        
        # Construct a universal output file
        p.object <- list()
        p.object$scores <- pca$scores
        p.object$eig <- pca$eig
        p.object$loadings <- pca$loadings
        p.object$call <- match.call()
        
    }  #### END SNP ANALYSIS
    
    # PLOT THE DIAGNOSTICS
    
    # Plot Scree plot avoid no visible binding probl
    eigenvalue <- percent <- NULL
    
    df <-
        data.frame(eigenvalue = seq(1:length(eig.top.pc)), percent = eig.top.pc)
    if (datatype == "SNP") {
        xlab <- paste("PCA Axis")
    } else {
        xlab <- paste("PCoA Axis")
    }
    ylab <- paste("Percentage Contribution")
    
    p1 <-
      ggplot(df, aes(x = eigenvalue, y = percent)) + 
      geom_line(color = plot.colors[2]) + 
      geom_point(color = plot.colors[1], size = 4) +
        geom_hline(yintercept = 10, color = "blue") + 
      plot.theme + xlab(xlab) + 
      ylab(ylab) + 
      ggtitle(title)
    
    if (any(eig.raw < 0)) {
        main <- "Noise Axes -- Warning: some eigenvalues < 0"
    } else {
        main <- "Noise Axes -- all eigenvalues positive"
    }
    
    p2 <-
      ggplot(as.data.frame(eig.raw.noise), aes(x = eig.raw.noise)) + 
      geom_histogram(bins = 50, color = plot.colors[1], fill = plot.colors[2]) +
        geom_vline(xintercept = 0, color = "blue") +
      plot.theme + 
      xlab("Eigenvalue") + 
      ylab("Count") + 
      ggtitle(main)
    
    if(any(eig.raw.noise==0)){
      if(verbose>=2){
        cat(warn("Warning: Some eigenvalues are zero, which suggests redundancy in the SNP loci\n"))
        cat(warn("  This may arise if there are duplicate individuals or clones in the dataset. Worth checking.\n"))    
      }
    }
    
    # printing outputs
    p3 <- (p1 / p2)
    if (verbose >= 1) {
        if(plot.out){print(p3)}
    }

    # Optionally save the plot ---------------------
    
    if(!is.null(plot.file)){
      tmp <- utils.plot.save(p3,
                             dir=plot.dir,
                             file=plot.file,
                             verbose=verbose)
    }

    # FLAG SCRIPT END
    
    if (verbose > 0) {
        cat(report("Completed:", funname, "\n"))
    }
    
    class(p.object) <- "glPca"
    invisible(p.object)
}
