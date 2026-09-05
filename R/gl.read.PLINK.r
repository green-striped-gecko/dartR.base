#' @name gl.read.PLINK
#' @title Reads PLINK data file into a genlight object
#' @family io
#' @description This function imports PLINK data into a genlight object and
#' append available metadata.
#'
#' @details This function handles .ped or .bed file (with the associate files -
#' e.g. .fam, .bim). However, if a .ped file is provided, PLINK needs to be
#' installed and it is used to convert the .ped into a .bed, which is then
#'  converted into a genlight. The conversion is run in a temporary
#' directory, so no files are written into the input directory.
#'
#' Genotypes are coded as the count of \code{allele.2} (column 6 of the
#' .bim file, the major allele under PLINK 1.x defaults): 0 = homozygous
#' \code{allele.1}, 2 = homozygous \code{allele.2}. Note that this is the
#' opposite orientation to \code{\link{gl.read.vcf}}, which counts the ALT
#' allele.
#'
#' Additional metadata can be included passing .csv files. These will be
#' appended to the existing metadata present in the PLINK files.
#'
#' The locus metadata needs to be in a csv file with headings, with a mandatory
#' column headed AlleleID corresponding exactly to the locus identity labels
#' provided with the SNP data.
#'
#' The individual metadata needs a mandatory column headed id whose entries
#' match the individual labels in the PLINK .fam file; rows are matched to
#' individuals by id, so their order in the file does not matter. If a pop
#' column is present it is used to assign populations.
#'
#' @param filename Fully qualified path to PLINK input file (without including
#' the extension) [required].
#' @param ind.metafile Name of the csv file containing the metrics for
#' individuals [default NULL].
#' @param loc.metafile Name of the csv file containing the metrics for
#' loci [default NULL].
#' @param fbm Logical. If TRUE, the genotypes of the returned object are held
#' in a file-backed matrix (fbm) in its @fbm slot rather than in @gen. This
#' is useful for very large datasets that do not fit into RAM. Note that
#' using fbm objects requires the package bigsnpr to be installed. To back
#' convert use \code{gl.fbm2gen()} [default FALSE].
#' @param plink.flags additional possible parameters passed on to plink
#' [default NULL].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' progress log; 3, progress and results summary; 5, full report
#' [default NULL, adopting the global verbosity set by gl.set.verbosity(),
#' or 2 if no global is set].
#' @inheritParams utils.plink.run
#' @return A genlight object with the SNP data and associated metadata included.
#' @author Author(s): Carlo Pacioni. Custodian: Carlo Pacioni -- Post to
#' \url{https://groups.google.com/d/forum/dartr}
#' @examples
#' \donttest{
#' # Create a small PLINK binary fileset (4 individuals x 4 SNPs) and read it
#' base <- file.path(tempdir(), "toy")
#' writeLines(c("1\tsnp1\t0\t100\tG\tA", "1\tsnp2\t0\t200\tT\tG",
#'              "2\tsnp3\t0\t150\tA\tC", "2\tsnp4\t0\t300\tC\tT"),
#'            paste0(base, ".bim"))
#' writeLines(c("FAM1 ind1 0 0 1 1", "FAM1 ind2 0 0 2 1",
#'              "FAM2 ind3 0 0 1 2", "FAM2 ind4 0 0 2 2"),
#'            paste0(base, ".fam"))
#' writeBin(as.raw(c(0x6c, 0x1b, 0x01, 0xcb, 0x2f, 0xab, 0xb7)),
#'          paste0(base, ".bed"))
#' gl <- gl.read.PLINK(base, verbose = 0)
#' as.matrix(gl)
#' }
#' @export
#' @rawNamespace import(data.table, except = c(melt,dcast))
#' @importFrom snpStats read.plink write.SnpMatrix row.summary

gl.read.PLINK <- function(filename,
                          ind.metafile = NULL,
                          loc.metafile = NULL,
                          plink.cmd = "plink",
                          plink.path = "path",
                          fbm=FALSE,
                          plink.flags = NULL,
                          verbose = NULL) {
  # SET VERBOSITY
  verbose <- gl.check.verbosity(verbose)
  
  # FLAG SCRIPT START
  funname <- match.call()[[1]]
  utils.flag.start(func = funname,
                   build = "Jackson",
                   verbose = verbose)
  
  # FUNCTION SPECIFIC ERROR CHECKING
  
  # check if packages are installed
  pkg <- "snpStats"
  if (!(requireNamespace(pkg, quietly = TRUE))) {
    stop(error(
      "Package",
      pkg,
      " needed for this function to work. Please install it."
    ))
  }
  
  if (is.null(loc.metafile) & verbose > 0) {
    cat(
      warn(
        "Warning: Locus metafile not provided, locus metrics will be
        calculated where this is possible\n"
      )
    )
  }
  
  if (is.null(ind.metafile) & verbose > 0) {
    cat(
      warn(
        "Warning: Individual metafile not provided, pop set to 'A' for all individuals\n"
      )
    )
  }
  
  # DO THE JOB
  dir.out <- tempdir() # use a tmp dir to handle transformation
  plink.fns <- list.files(
    path = dirname(filename),
    pattern = paste0("^", basename(filename)),
    full.names = TRUE,
    recursive = TRUE
  )
  plink.fn <- plink.fns[grep(".bed$", plink.fns)]
  if (length(plink.fn) > 1) {
    # If there is more than one .bed
    stop(error("Found more than one .bed file and don't know which one to use\n"))
  } else {
    if (length(plink.fn) == 0) {
      # If there is no .bed
      plink.fn <- plink.fns[grep(".ped$", plink.fns)] # look for .ped
      if (length(plink.fn) == 0 |
          length(plink.fn) > 1) {
        # If there is no or >1 .ped
        stop(
          error(
            "Found no .bed files and",
            length(plink.fn),
            ".ped file(s). This function needs either one .bed or .ped file\n"
          )
        )
      } else {
        # This should mean that there is one .ped
        if (length(plink.fn) == 1) {
          # But it is checked here for safety
          # Copy the fileset to the temporary directory so the conversion
          # does not write into the user's data directory (F5)
          file.copy(plink.fns, dir.out, overwrite = TRUE)
          utils.plink.run(
            dir.in = dir.out,
            plink.cmd = plink.cmd,
            plink.path = plink.path,
            out = basename(filename),
            syntax = paste(
              "--file",
              basename(filename),
              "--make-bed",
              plink.flags
            ),
            verbose = verbose
          )
          plink.fn <- file.path(dir.out, paste0(basename(filename), ".bed"))
          if (!file.exists(plink.fn)) {
            stop(error("PLINK did not produce a .bed file. Check that PLINK",
                       "is installed and the input files are valid\n"))
          }
        } else {
          # If we got it wrong stop the FUN
          stop(error("Couldn't find any .bed or .ped file\n"))
        }
        
      } # Close if there is 1 .ped
    } # close if there is no .bed
  } # Close the else L85 from here onwards there should be only 1 .bed in plink.fn
  
  snpMatrix <- snpStats::read.plink(bed = sub(x = plink.fn, pattern =
                                                ".bed$", ""))
  gen <- snpMatrix[[1]]
  fam <- snpMatrix[[2]]
  map <- snpMatrix[[3]]
  row.names(gen) <- fam$member
  snpStats::write.SnpMatrix(gen, file.path(dir.out, "snpMatrixComb.txt"))
  
  suppressWarnings(genCombdt <- data.table::fread(file = file.path(dir.out, "snpMatrixComb.txt")))
  setnames(genCombdt, "V1", "id")
  
  genCombNA <- as.matrix(genCombdt[, -1, with = FALSE])
  row.names(genCombNA) <- fam$member
  gl <- new(
    "genlight",
    gen = genCombNA,
    ind.names = fam$member,
    loc.names = map$snp.name,
    chromosome = map$chromosome,
    position = map$position,
    # Need to confirm that we want here
    #the chr position and that the SNP position on the read
    # goes only in the loc.metrics
    other = list(
      ind.metrics = cbind(
        id = fam$member,
        snpStats::row.summary(gen),
        # This adds Call.rate, Certain.calls, Heterozygosity
        Family = fam$pedigree,
        Father = fam$father,
        Mother = fam$mother,
        Sex = fam$sex
      ),
      loc.metrics = map
    )
  )
  
  if (length(unique(fam$member)) != length(fam$member)) {
    stop(error(
      "Fatal Error: Individual labels are not unique, check and edit your input file\n"
    ))
  }

  if (length(unique(map$snp.name)) != length(map$snp.name)) {
    stop(error(
      "Fatal Error: AlleleID not unique, check and edit your input file\n"
    ))
  }
  
  pop(gl) <- array("A", nInd(gl))
  names(gl@other$loc.metrics)[2] <- "AlleleID"
  # Fix order of cols
  gl@other$loc.metrics <-
    gl@other$loc.metrics[, c("AlleleID", names(gl@other$loc.metrics)[-2])]
  gl <- gl.compliance.check(gl, verbose = verbose)
  
  # NOW THE LOCUS METADATA
  
  if (!is.null(loc.metafile)) {
    loc.metrics <-
      read.csv(file = loc.metafile,
               header = TRUE,
               stringsAsFactors = TRUE)
    if (!("AlleleID" %in% names(loc.metrics))) {
      stop(error(
        "Fatal Error: mandatory AlleleID column absent from the locus metrics file\n"
      ))
    }
    
    if (nrow(loc.metrics) != nLoc(gl))
      stop(
        error(
          "Fatal Error: the locus metrics file does not have the same number of loci of the input data file\n"
        )
      )
    
    if (!all(loc.metrics[, "AlleleID"] %in% gl@other$loc.metrics$AlleleID))
      stop(
        error(
          "Fatal Error: AlleleID in the locus metrics file does not correspond with",
          "AlleleID in the input data file\n"
        )
      )
    
    row.names(loc.metrics) <- loc.metrics$AlleleID
    gl@other$loc.metrics <- cbind(gl@other$loc.metrics, loc.metrics[map$snp.name, ])
  }
  
  gl <- gl.recalc.metrics(gl, verbose = 0)
  
  if (verbose >= 2) {
    cat(report(
      paste(
        " Added or updated ",
        names(gl@other$loc.metrics),
        "to the other$ind.metrics slot.\n"
      )
    ))
  }
  
  # NOW THE INDIVIDUAL METADATA
  
  if (!is.null(ind.metafile)) {
    ind.metrics <-
      read.csv(
        file = ind.metafile,
        header = TRUE,
        stringsAsFactors = TRUE,
        fileEncoding = "UTF-8-BOM"
      )
    if (!("id" %in% names(ind.metrics))) {
      cat(
        error(
          "Fatal Error: mandatory id column absent from the individual metadata file\n"
        )
      )
      stop()
    }
    fam$member %in% ind.metrics[, "id"]
    if (sum(fam$member %in% ind.metrics[, "id"]) < length(fam$member))
      stop(error(
        paste(
          "Fatal Error: there are",
          sum(fam$member %in% ind.metrics[, "id"]),
          "individuals id that match the ones in the data, but",
          length(fam$member),
          "individuals genotyped\n"
        )
      ))
    
    
    if (!("pop" %in% names(ind.metrics))) {
      cat(
        warn(
          "  Warning: pop column absent from the individual metadata file, setting to 'A'\n"
        )
      )
      ind.metrics$pop <- array("A", nInd(gl))
    }

    # Match the metafile rows to the individuals in .fam order by id, so
    # the row order of the metafile does not matter (F2)
    ind.metrics <- ind.metrics[match(fam$member, ind.metrics$id), ]

    gl@other$ind.metrics <- cbind(gl@other$ind.metrics, ind.metrics)
    pop(gl) <- as.character(ind.metrics$pop)
    if (verbose >= 2) {
      cat(report(
        paste(
          " Added ",
          names(gl@other$ind.metrics),
          " to the other$ind.metrics slot.\n"
        )
      ))
    }
  }
  
  # MAKE COMPLIANT
  gl <- gl.compliance.check(gl, verbose = verbose)
  
  # ADD TO HISTORY (add the first entry)
  gl@other$history <- list()
  gl@other$history[[1]] <- match.call()
  
  # FLAG SCRIPT END
  
  if (verbose > 0) {
    cat(report("Completed:", funname, "\n"))
  }
  
  # convert to fbm only when requested; otherwise the object remains
  # @gen-backed (F1)
  if (fbm) {
    gl <- gl.gen2fbm(gl, verbose = verbose)
    if (verbose >= 2) {
      cat(report(" Created a file-backed matrix (fbm) dartR object\n"))
    }
  }

  return(gl)
  
}
