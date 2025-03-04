#' @name gl2genepop
#' @title Converts a genlight object into genepop format (and file)
#' @family linker

#' @description
#' The genepop format is used by several external applications (for example
#' Neestimator2. Unfortunatelly, the software seems to be no longer easily 
#' available. To install use the \link[dartRverse]{gl.download.binary} function.
#' So the main idea is to create the genepop file and then run the other
#' software externally. As a feature, the genepop file is also returned as an
#' invisible data.frame by the function.
#' 
#' @param x Name of the genlight object containing the SNP data [required].
#' @param outfile File name of the output file [default 'genepop.gen'].
#' @param outpath Path where to save the output file [default global working 
#' directory or if not specified, tempdir()].
#' @param pop.order Order of the output populations either "alphabetic" or a 
#' vector of population names in the order required by the user (see examples)
#' [default "alphabetic"].
#' @param output.format Whether to use a 2-digit format ("2_digits") or 3-digits
#'  format ("3_digits") [default "2_digits"].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' progress log; 3, progress and results summary; 5, full report
#' [default 2, unless specified using gl.set.verbosity].
#' 
#' @author Custodian: Bernd Gruber (Post to
#' \url{https://groups.google.com/d/forum/dartr})
#' 
#' @examples
#' \donttest{
#' require("dartR.data")
#' # SNP data
#' geno <- gl2genepop(possums.gl[1:3,1:9], outpath = tempdir())
#' head(geno)
#' test <- gl.filter.callrate(platypus.gl,threshold = 1)
#' popNames(test)
#' gl2genepop(test, pop.order = c("TENTERFIELD","SEVERN_ABOVE","SEVERN_BELOW"),
#'            output.format="3_digits", outpath = tempdir())
#' }
#' 
#' @export
#' @return Invisible data frame in genepop format

gl2genepop <- function (x,
                        outfile = "genepop.gen",
                        outpath = NULL,
                        pop.order = "alphabetic",
                        output.format = "2_digits", 
                        verbose = NULL) {
  # SET VERBOSITY
  verbose <- gl.check.verbosity(verbose)
  
  # SET WORKING DIRECTORY
  outpath <- gl.check.wd(outpath,verbose=0)
  
  # FLAG SCRIPT START
  funname <- match.call()[[1]]
  utils.flag.start(func = funname,
                   build = "v.2023.2",
                   verbose = verbose)
  
  # CHECK DATATYPE
  datatype <- utils.check.datatype(x, verbose = verbose)
  
  # FUNCTION SPECIFIC ERROR CHECKING
  
  #works only with SNP data
  if (datatype != "SNP") {
    cat(error(
      "  Only SNPs (diploid data can be transformed into genepop format!\n"
    ))
    stop()
  }
  
  if (is.null(pop(x))) {
    cat(
      important(
        "Your genlight object does not have a population definition. Therefore the function assumes the whole genlight object has one population!\n"
      )
    )
    pop(x) <- rep("Pop1", nInd(x))
  }
  
  # DO THE JOB
  #ordering populations
  # filtering all loci with all NAs to avoid crashing the function
  x <- gl.filter.allna(x,verbose =0 )
  
  if(length(pop.order)==1){
    x <- x[order(pop(x)), ]
  }
  
  if(length(pop.order)>1){
    
    tmp <- seppop(x)
    tmp_2 <-  tmp[match(pop.order, names(tmp))]
    x <- Reduce(rbind,tmp_2)
  }

  #convert to genind
  x <- gl2gi(x, verbose = 0, probar = FALSE)
  data <- as.matrix(x)
  pop_names <- x@pop
  
  if (all(unlist(unique(x@all.names)) %in% c("A", "T",
                                             "C", "G"))) {
    m_type <- "snp"
    colnames(data) <- gsub(colnames(data),
                           pattern = "\\.A",
                           replacement = ".01")
    colnames(data) <- gsub(colnames(data),
                           pattern = "\\.T",
                           replacement = ".02")
    colnames(data) <- gsub(colnames(data),
                           pattern = "\\.C",
                           replacement = ".03")
    colnames(data) <- gsub(colnames(data),
                           pattern = "\\.G",
                           replacement = ".04")
  }
  
  loci_names_l <- x@loc.fac
  loc_all <- data.frame(col = colnames(data))
  if (all(stringr::str_count(colnames(data), "\\.") ==
          1) != TRUE) {
    stop(
      error(
        "The columns' names of x@tab must be of form 'locus.allele' with only 1\n'.' between locus and allele"
      )
    )
  }
  loc_all <- tidyr::separate(
    loc_all,
    col = 1,
    sep = "\\.",
    into = c("locus", "allele")
  )
  
  if(output.format == "3_digits"){
    loc_all$allele <- paste0("0",loc_all$allele)
  }
  
  loci_names <- as.character(loci_names_l[-which(duplicated(loci_names_l))])
  n.loci <- length(loci_names_l[-which(duplicated(loci_names_l))])
  data_gpop <- data.frame(id = paste(pop_names, "_",
                                     row.names(data), ",", sep = ""))
  for (i in 1:n.loci) {
    loc <- loci_names[i]
    a <- c()
    for (j in 1:nrow(data)) {
      col_loc <- which(loc_all[, "locus"] == loc)
      hom <- which(data[j, col_loc] == 2)
      het <- which(data[j, col_loc] == 1)
      if (length(hom) != 0) {
        a[j] <- paste(loc_all[col_loc[hom], "allele"],
                      loc_all[col_loc[hom], "allele"], sep = "")
      } else if (length(het) != 0) {
        # disabling the order of alleles
        # if (as.character(loc_all[col_loc[het[1]], "allele"]) <
        #     as.character(loc_all[col_loc[het[2]], "allele"])) {
          a[j] <- paste(loc_all[col_loc[het[1]], "allele"],
                        loc_all[col_loc[het[2]], "allele"],
                        sep = "")
        # } else {
        #   a[j] <- paste(loc_all[col_loc[het[2]], "allele"],
        #                 loc_all[col_loc[het[1]], "allele"],
        #                 sep = "")
        # }
      } else {
        if (nchar(loc_all[1, "allele"]) == 3) {
          a[j] <- "000000"
        } else {
          a[j] <- "0000"
        }
      }
    }
    data_gpop <- cbind(data_gpop, a)
  }
  colnames(data_gpop) <- c("ID", as.character(loci_names))
  data_gpop[,] <- apply(data_gpop, c(1, 2), as.character)
  
  dummy <- paste("Genepop output. Loci:", nLoc(x), "Populations:", nPop(x))
  dummy[2] <- paste(locNames(x), collapse = ",")
  cs <- c(cumsum(table(pop(x))))
  from <-c(1, (cs[-length(cs)] + 1))
  to <-cs
  
  for (i in 1:nPop(x)) {
    da <-
      apply(data_gpop[from[i]:to[i], ], 1, function(y)
        paste(y, collapse = " "))
    dummy <- c(dummy, "Pop", da)
  }
  
  data_gpop2 <- dummy
  utils::write.table(
    data_gpop2,
    file = file.path(outpath, outfile),
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE
  )
  cat(report(
    "  The genepop file is saved as: ",
    file.path(outpath, outfile, "\n")
  ))
  
  # FLAG SCRIPT END
  
  if (verbose >= 1) {
    cat(report("Completed:", funname, "\n"))
  }
  
  # RETURN
  return(invisible(data.frame(lines = data_gpop2)))
  
}
