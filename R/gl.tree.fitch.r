#' @name gl.tree.fitch
#' @title Generates a distance phylogeny
#' 
#' @family phylogeny
#' 
#' @description
#' Generates a distance phylogeny from a distance object
#'  using the Fitch-Margoliash algorithm in Phylip.

#' @param D Name of the distance matrix for tree building [required]
#' @param x Name of the genlight object containing the SNP data [required for bootstrapping, default NULL].
#' @param phylip.path Path to the directory that holds the Phylip executables [required].
#' @param out.path Path to the directory to save files produced by the analysis [default tempdir()]

# Parameters that govern tree selection
#' @param tree.method Algorithm used for constructing trees and selecting the best tree [default "FM"]
#' @param outgroup Name of the outgroup taxon [default NULL, no outgroup, tree not rooted]
#' @param global.rearrange If TRUE, undertake global rearrangements when generating the tree [default FALSE].
#' @param randomize If TRUE, randomize the order of the input taxa [default FALSE].
#' @param n.jumble Number of randomizations of the input order, must be odd [default 9]
#' @param bstrap Number of bootstrap replicates [default 1000]
#' 
# Parameters that govern the appearance of the tree
#' @param plot.type One of 'phylogram','cladogram','unrooted','fan','tidy','radial' [default "phylogram"]
#' @param bstrap.threshold Threshold for bootstrap values to be displayed on the tree [default 0.8]
#' @param branch.width Width of the branches [default 2]
#' @param branch.color Colour of the branches [default "blue"]
#' @param node.label.color Colour of the node labels [default "red"]
#' @param terminal.label.cex Height of the taxon label text [default 0.8]
#' @param node.label.cex Height of the node label text [default 0.8]
#' @param offset Horizontal offset of the node labels from the node [default 1.8]
#' 
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' progress log; 3, progress and results summary; 5, full report
#' [default 2, unless specified using gl.set.verbosity].
#' 
#' @details
#' The script takes a distance object as input. This distance object is typically created with gl.dist.phylo().
#' The script then creates a file consistent with what is expected
#' by program fitch in the Phylip suite of executables. It then runs fitch to generate the "best" phylogenetic
#' tree. Program fitch is run again with bstrap replicates to generate bootstrap support for each node
#' in the tree and plots these on the tree.
#' 
#'      tree.method : Currently only Fitch-Margoliash is implemented.
#'      
#'      outgroup : Name the taxon to be used as outgroup. Must be among the names of the populations defined in the genlight object.
#' 
#' @author Custodian: Arthur Georges -- Post to 
#' \url{https://groups.google.com/d/forum/dartr}
#' 
#' @examples
#' 
#' \dontrun{
#' tmp <- gl.filter.monomorphs(testset.gl)
#' D <- gl.dist.phylo(testset.gl,subst.model="F81")
#' gl.phylip(D=D,x=tmp,phylip.path="D:/workspace/R/phylip-3.695/exe",plot.type="unrooted",
#' node.label.cex=0.5,terminal.label.cex=0.6,global.rearrange = FALSE, bstrap=10)
#' }
#' 
# Testing
# gl <- testset.gl
# gl <- gl.filter.callrate(gl)
# gl <- gl.filter.taglength(gl,lower=40)
# dd <- gl.dist.phylo(testset.gl,subst.model="F81")
# gl.tree.fitch(D=dd,phylip.path="D:/workspace/R/phylip-3.695/exe")
# gl.tree.fitch(D=dd,phylip.path="D:/workspace/R/phylip-3.695/exe",global.rearrange = TRUE)
# gl.tree.fitch(D=dd,phylip.path="D:/workspace/R/phylip-3.695/exe",outgroup="EmvicVictJasp")
# gl.tree.fitch(D=dd,x=gl,plot.type="tidy",bstrap=5,phylip.path="D:/workspace/R/phylip-3.695/exe",outgroup="EmvicVictJasp",verbose=3)
#' 
#' @import ape
#' 
#' @export
#' @return The tree file in newick format.

gl.tree.fitch <- function(D,
                      x=NULL,
                      phylip.path,
                      out.path=tempdir(),
                      tree.method="FM",
                      outgroup=NULL,
                      global.rearrange=FALSE,
                      randomize=FALSE,
                      n.jumble=9,
                      bstrap=1, # No bootstrapping unless set to > 1, e.g. 1000
                      plot.type="phylogram",
                      bstrap.threshold=0.8,
                      branch.width=2,
                      branch.color="blue",
                      node.label.color="red",
                      terminal.label.cex=0.8,
                      node.label.cex=0.8,
                      offset=1.2,
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
  
  datatype1 <- utils.check.datatype(D, accept=c('dist'), verbose = verbose)
  
  if(bstrap > 1){
    if(is.null(x)){
      cat(warn("  Bootstrapping requires access to the genlight object, not provided, bootstrapping disabled\n"))
      bstrap <- 1
    } else {
      datatype2 <- utils.check.datatype(x, verbose = verbose)
    }
  } 
  
  # n.jumble must be odd
    if(n.jumble %% 2 == 0){
      n.jumble <- n.jumble+1
      cat(warn("  Parameter n.jumble must be an odd integer, adding 1 to give",n.jumble,"\n"))}
  
  # # Check for monomorphic loci
  # tmp <- gl.filter.monomorphs(x, verbose = 0)
  # if ((nLoc(tmp) < nLoc(x))) {
  #   if(verbose >= 2){cat(warn("  Warning: genlight object contains monomorphic loci\n"))}
  # }
  
  tree.method <- toupper(tree.method)
  if(tree.method=="F-M"){tree.method <- "FM"}
  if(tree.method != "FM"){
    cat(warn(  "  Using Fitch-Margoliash method of tree construction\n"))
  }
  
  if(!is.null(outgroup)){
    if(!outgroup %in% attr(D,"Labels")){
      cat(warn("  Specified outgroup not in list of available taxa, set to unrooted\n"))
      outgroup <- NULL
    } else {
      outgroup <- which(attr(D,"Labels") == outgroup)
      outgroup <- outgroup[1]
    }
  }  
  
  # Specify the Phylip program fitch and command line to run it
  #  on different platforms
  
  if(Sys.info()["sysname"] == "Windows") {
    prog <- "fitch.exe"
    cmd.fm <- paste0("fitch.exe < fitch.cmd")
  }
  
  if (Sys.info()["sysname"] == "Linux") {
    prog <- "fitch"
    cmd.fm <- "./fitch < fitch.cmd"
  }
  
  if (Sys.info()["sysname"] == "Darwin") {
    prog <- "fitch"
    cmd.fm <- "./fitch < fitch.cmd"
  }
  
  # Transfer the Phylip executable fitch.exe to the tempdir
  
  if (file.exists(file.path(phylip.path, prog))) {
    file.copy(file.path(phylip.path, prog),
              to = tempdir(),
              overwrite = TRUE)
    if(verbose >= 2){cat(report("  Transferred fitch executables to tempdir\n"))}
  } else{
    cat(
      error(
        " Cannot find",
        prog,
        "in the specified folder :",
        phylip.path,
        "\n"
      )
    )
    stop()
  }
  
  # Specify the Phylip program consense and command line to run it
  #  on different platforms
  
  if(Sys.info()["sysname"] == "Windows") {
    prog <- "consense.exe"
    cmd.cns <- paste0("consense.exe < consense.cmd")
  }
  
  if (Sys.info()["sysname"] == "Linux") {
    prog <- "consense"
    cmd.cns <- "./consense < consense.cmd"
  }
  
  if (Sys.info()["sysname"] == "Darwin") {
    prog <- "consense"
    cmd.cns <- "./consense < consense.cmd"
  }
  
  # Transfer the Phylip executable consense.exe to the tempdir
  
  if (file.exists(file.path(phylip.path, prog))) {
    file.copy(file.path(phylip.path, prog),
              to = tempdir(),
              overwrite = TRUE)
    if(verbose > 1){cat(report("  Transferred consense executables to tempdir\n"))}
  } else{
    cat(
      error(
        " Cannot find",
        prog,
        "in the specified folder :",
        phylip.path,
        "\n"
      )
    )
    stop()
  }
  
  # DO THE JOB
  
  # # Calculate the distance matrix
  # D <- gl.dist.phylo(x)
  
  # Convert dist object to matrix of form expected by phylip
  mat <- as.matrix(D)
  rownames(mat) <- NULL
  colnames(mat) <- NULL
  names <- sapply(attr(D,"Labels"),function(x) substr(x, start = 1, stop = 10))
  names <- sprintf("%-10s",names)
  mat <- cbind(names,mat)
  
  # Output the distance matrix in phylip format to tempdir()
  hold <- getwd()
  setwd(tempdir())
    write(length(attr(D,"Labels")),file="infile")
    write.table(mat,file="infile",row.names=FALSE,col.names=FALSE,quote=FALSE,append=TRUE)
    #write.table(D,file="infile",quote=FALSE,append=TRUE)
    #cat("Hello Vietnam\n")
  setwd(hold)
  
  #### TREE BUILDING
  #### METHOD = FITCH

  # Delete previous versions of output and outtree 
  hold <- getwd()
  setwd(tempdir())
  if(file.exists(file.path(tempdir(),"outfile"))){
    shell(paste("del outfile"))
  }
  if(file.exists(file.path(tempdir(),"outtree"))){
    shell(paste("del outtree"))
  }
  setwd(hold)
  
  # Create the command file for fitch
  vector <- rep(NA,7)
  if(!is.null(outgroup)){
    vector[1] <- "O"
    vector[2] <- outgroup
  }
  if(global.rearrange){
    vector[3] <- "G"
  }
  if(randomize){
    vector[4] <- "J"
    vector[5] <- 12345
    vector[6] <- n.jumble
  }
  vector[7] <- "Y"
  vector <- na.omit(vector)
  
  fitch.cmd <- file.path(tempdir(), "fitch.cmd")
  sink(fitch.cmd)

  cat(vector,sep="\n")
  
  sink()
  
  # Execute the command file by passing it to fitch
  if(verbose >= 2){cat(report("  Running fitch from the Phylip executables\n"))}
  hold <- getwd()
  setwd(tempdir())
     shell(cmd.fm)
  setwd(hold)
  # This will produce an output file called outfile, and a newick file called outtree
  
  # Display the newick file for the phylogeny    
  if(verbose>=2){
    newick <- readLines(file.path(tempdir(),"outtree"))
    if(verbose > 2){
      cat(report("Newick format tree file:\n"))
      cat(newick, sep = "\n")
    }
  }
  
  # Save the 'best' tree
  tree <- read.tree(file.path(tempdir(),"outtree"))
  tree_1 <- tree
  
  # BOOTSTRAPS -- USING CONSENSE
  
  if(bstrap > 1){ 
    
    if (verbose >= 2) {
      cat(report("Writing bootstrap distance matricies to the Phylip input file\n"))
      cat(report(
        "Repeating calculations for",
        bstrap,
        "iterations\n"
      ))
    }
    
    # First the n=bstrap trees
    
    count <- 0
    hold <- getwd()
    setwd(tempdir())
    for(i in 1:bstrap){
      count <- count + 1
      if(verbose >= 2){cat("  Bootstrap replicate",count,"\n")}
      # Create the sequences
      tmp <- gl.subsample.loc(x,n = nLoc(x),replace=TRUE,verbose=0)
      
      dd <- gl.dist.phylo(tmp,by.pop=TRUE,verbose=verbose)
      mat <- as.matrix(dd)
      rownames(mat) <- NULL
      colnames(mat) <- NULL
      names <- sapply(attr(dd,"Labels"),function(x) substr(x, start = 1, stop = 10))
      names <- sprintf("%-10s", names)
      mat <- cbind(names,mat)
      # gl2fasta(x,outfile="tmp.fas",outpath = tempdir(),verbose=0)
      # sequences <- ape::read.dna("tmp.fas", format = "fasta")
      # # Calculate distances
      # D <- ape::dist.dna(sequences,model=subst.model,pairwise.deletion = pairwise.missing)
      # #Calculate average distances for pairwise populations
      # D <- avg.dist(x,D,verbose=0)
      if(count==1){
        write(length(attr(dd,"Labels")),file="infile")
      } else {
        write(length(attr(dd,"Labels")),file="infile",append=TRUE)
      }
      write.table(mat,file="infile",row.names=FALSE,col.names=FALSE,quote=FALSE,append=TRUE)
      #write.table(D,file="infile",quote=FALSE,append=TRUE)
    }
    setwd(hold) ### infile contains the multiple bootstrapped matricies
    
    
    #### Run fitch to generate the outtree file to pass to consense
    
    #Delete previous versions of output and outtree
    hold <- getwd()
    setwd(tempdir())
    if(file.exists(file.path(tempdir(),"outtree"))){
      cat("    Deleting old outtree file\n")
      shell(paste("del outtree"))
    }
    
    # Create the command file for fitch
    vector <- rep(NA,11)
    vector[1] <- "A"
    if(!is.null(outgroup)){
      vector[2] <- "O"
      vector[3] <- outgroup
    }
    # if(global.rearrange){
    #   vector[4] <- "G"
    # }
    # if(randomize){
    #   vector[5] <- "J"
    #   vector[6] <- 12345
    #   vector[7] <- n.jumble
    # }
    vector[8] <- "M"
    vector[9] <- bstrap
    vector[10] <- 331
    vector[11] <- "Y"
    vector <- as.vector(na.omit(vector))
    
    fitch.cmd <- file.path(tempdir(), "fitch.cmd")
    sink(fitch.cmd)
    cat(vector,sep="\n")
    sink()
    
    # Execute the command file by passing it to fitch
    
    shell(cmd.fm)
    
    # Copy the outtree file from fitch to the intree for consense
    shell(paste("cp outtree intree"))
    shell(paste("del outtree"))
    
    # Create the consensus command file
    vector <- rep(NA,4)
    vector[1] <- "A"
    if(!is.null(outgroup)){
      vector[2] <- "O"
      vector[3] <- outgroup
    }
    vector[4] <- "Y"
    vector <- na.omit(vector)
    
    consense.cmd <- file.path(tempdir(), "consense.cmd")
    sink(consense.cmd)
    cat(vector,sep="\n")
    sink()
    
    # Execute the command file by passing it to consense
    # hold <- getwd()
    #  setwd(tempdir())
    shell(cmd.cns)
    #setwd(hold)
    # This will append to an output file called outfile, and output a newick file called outtree
    
    tree.bstraps <- ape::read.tree(file.path(tempdir(),"outtree"))
    if(verbose==3){
      newick <- readLines(file.path(tempdir(),"outtree"))
      if(verbose > 2){
        cat(report("Newick format for the bootstrap tree file showing node values:\n"))
        cat(newick, sep = "\n")
      }
    }
    # Extracting branch lengths
    branch_lengths <- tree.bstraps$edge.length
    # The number of tips is the same as the number of leaf labels
    num_tips <- length(tree.bstraps$tip.label)
    # Branch lengths associated with internal nodes start after all the tips in the tree structure
    # Because in a fully bifurcating tree, the number of internal nodes is one less than the number of tips
    node_values <- branch_lengths[(num_tips + 1):length(branch_lengths)]/bstrap
    node_values[node_values < bstrap.threshold] <- 0
    node_values <- as.character(node_values)
    node_values[node_values=="0"] <- " "
    
    # cat(report(node_values))
    # cat("\n")
  }
  
  # Plot the tree
  
  plot(tree_1, type = plot.type, edge.width = branch.width, edge.color = branch.color, cex = terminal.label.cex)
  if(bstrap > 1){
    nodelabels(node_values, frame = "none", cex = node.label.cex, col = node.label.color,adj=offset)
  }
  if(verbose > 2){cat(report("Plotting a Fitch-Margoliash Distance Phylogeny\n"))}
  
  # FLAG SCRIPT END
  
  if (verbose > 0) {
    cat(report("Output files from Phylip written to", tempdir(),"\n"))
    cat(report("Completed:", funname, "\n"))
  }

  return(tree_1)
}
