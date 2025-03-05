
# rs snp identifier # contains the SNP identifier;
# alleles contains SNP alleles according to NCBI database dbSNP;
# chrom contains the chromosome that the SNP was mapped;
# pos contains the respective position of this SNP on chromosome;
# strand contains the orientation of the SNP in the DNA strand. Thus, SNPs could be in the forward (+) or in the reverse (-) orientation relative to the reference genome
# assembly contains the version of reference sequence assembly (from NCBI);
# center contains the name of genotyping center that produced the genotypes
# protLSID contains the identifier for HapMap protocol;
# assayLSID contain the identifier HapMap assay used for genotyping;
# panelLSID contains the identifier for panel of individuals genotyped;
# QCcode contains the quality control for all entries;


gl2hapmap <- function(x){
  x_mat <- as.matrix(x[, ])
  homs1 <- paste(substr(x@loc.all, 1, 1), "/", substr(x@loc.all, 1, 1), sep = "")
  hets <- x@loc.all
  homs2 <- paste(substr(x@loc.all, 3, 3), "/", substr(x@loc.all, 3, 3), sep = "")
  xx <- matrix(NA, ncol = ncol(x_mat), nrow = nrow(x_mat))
  for (i in 1:nrow(x_mat)) {
    for (ii in 1:ncol(x_mat)) {
      inp <- x_mat[i, ii]
      if (!is.na(inp)) {
        if (inp == 0)
          xx[i, ii] <- homs1[ii]
        else if (inp == 1)
          xx[i, ii] <- hets[ii]
        else if (inp == 2)
          xx[i, ii] <- homs2[ii]
      } else{
        xx[i, ii] <-"0/0"
      }
    }
  }
  xx <- gsub("/", "", xx)
  xx <- as.data.frame(xx)
  xx <- t(xx)
  colnames(xx) <- indNames(x)
  
  x$chromosome <- as.factor(as.numeric(x$chromosome))
  
  geno_tmp <- data.frame(rs = locNames(x),
                         alleles= x$loc.all,
                         chrom= x$chromosome,
                         pos= x$position,
                         strand="+",
                         assembly="Oilpalm",
                         center= NA,
                         protLSID= NA,
                         assayLSID= NA,
                         panel=NA,
                         QCcode=NA)
  
  res_output <- cbind(geno_tmp,xx)
  res_output <- as.matrix(res_output)
  res_output[] <- as.character(res_output)
  res_output <- as.data.frame(rbind(colnames(res_output),res_output))
  
  return(res_output)
}


