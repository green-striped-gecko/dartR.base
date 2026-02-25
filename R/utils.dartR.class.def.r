#' @import bigstatsr
#' @import bigsnpr
#' @importFrom methods callNextMethod slotNames
#' @import adegenet

#seppop needs to be imported to work for dartR
#also internal functions for "[" methods
seppop <- getFromNamespace("seppop", "adegenet")
.seppop_internal <- getFromNamespace(".seppop_internal", "adegenet")
.get_pop_inds <- getFromNamespace(".get_pop_inds", "adegenet")

## Allow the slot to be either an FBM or NULL
setClassUnion("FBMcode256_or_NULL", c("NULL", "FBM.code256"))

## returns NULL if the 'fbm' slot is missing OR is NULL
.fbm_or_null <- function(x) {
    if (methods::.hasSlot(x, "fbm")) {
    val <- methods::slot(x, "fbm")
    return(if (is.null(val)) NULL else val)
  }
  NULL
}

.has_fbm <- function(x) !is.null(.fbm_or_null(x))

fbm_or_gen <- function(x) {
  # early guard: must be an S4 object
  
  fbm_present <- methods::.hasSlot(x, "fbm")
  gen_present <- methods::.hasSlot(x, "gen")
  
  fbm_nonempty <- FALSE
  gen_nonempty <- FALSE
  
  if (fbm_present) {
    fbm_val <- tryCatch(methods::slot(x, "fbm"), error = function(e) NULL)
    fbm_nonempty <- !is.null(fbm_val)
  }
  
  if (gen_present) {
    gen_val <- tryCatch(methods::slot(x, "gen"), error = function(e) NULL)
    # gen is usually a list; consider it empty if length==0
    gen_nonempty <- !(is.list(gen_val) && length(gen_val) == 0L)
  }
  
  if (fbm_nonempty && !gen_nonempty) return("fbm")
  if (gen_nonempty && !fbm_nonempty) return("gen")
  
  # if both non-empty, that violates the XOR rule
  if (fbm_nonempty && gen_nonempty) {
    stop("Both @fbm and @gen are populated - invalid object.")
  }
  
  # if neither non-empty, report the missing situation clearly
  stop("Neither @fbm nor @gen slot found or populated in object.")
}

## (Re)define your subclass; it CONTAINS genlight and ADDS a slot
setClass("dartR",
         contains = "genlight",
         slots = c(fbm="ANY"
          
           #seqposition = "intOrNULL" 
           ),   # <-- snp position data set
         prototype = prototype(fbm = NULL), #, seqpostion = NULL),
         validity = function(object) {
           fbm     <- .fbm_or_null(object)
           has_fbm <- !is.null(fbm)
           has_gen <- length(object@gen) > 0L
           
           ## Enforce XOR: exactly one populated
           if ((has_fbm && has_gen) || (!has_fbm && !has_gen))
             return("Exactly one of @fbm or @gen must be populated (mutual exclusivity).")
           
           ## If using FBM, check dimensions against names/slots
           if (has_fbm) {
             if (!inherits(object@fbm, "FBM.code256"))
               return("@fbm must inherit from 'FBM.code256'.")
             dm <- dim(object@fbm)
             if (length(dm) != 2L) return("@fbm must be a 2D FBM.code256.")
             ## When FBM is the source of truth, ensure nInd/nLoc are consistent
             if (length(object@ind.names) > 0 && length(object@ind.names) != dm[1])
               return("length(@ind.names) must match nrow(@fbm).")
             if (length(object@loc.names) > 0 && length(object@loc.names) != dm[2])
               return("length(@loc.names) must match ncol(@fbm).")
           } 
           TRUE
         }
)

### adding new slots...
#setClass("dartR",slots=list(what="integer"), contains="genlight")
###

### adding new slots...
#setClass("dartR",slots=list(what="integer"), contains="genlight")
###

########################
## show dartR ##
########################
setMethod ("show", "dartR", function(object) {
  ## HEADER
  cat(" ********************\n")
  cat(" *** DARTR OBJECT ***\n")
  cat(" ********************")
  marker <- "mixed markers"
  if (all(!is.na(ploidy(object)))) {
    if (all(ploidy(object) == 2))
      marker <- "SNPs"
    if (all(ploidy(object) == 1))
      marker <- "silicoDarts (P/A) "
  }
  
  cat(
    "\n\n **",
    format(nInd(object), big.mark = ","),
    "genotypes, ",
    format(nLoc(object), big.mark = ","),
    marker,
    ", size:",
    format(object.size(object), units = "auto")
  )
  
  ## MISSING DATA
  if (!is.null(object@gen)) {
  temp <- sapply(object@gen, function(e)
    length(e@NA.posi))
  if (length(temp > 1)) {
    cat("\n\n    missing data: ",
        sum(temp),
        " (=",
        round((sum(temp) / (
          nInd(object) * nLoc(object)
        )) * 100, 2),
        " %) scored as NA",
        sep = "")
  }
  } 
  fbm <- .fbm_or_null(object)
  if (!is.null(fbm)) {
    fbm <- object@fbm
    n <- nInd(object); p <- nLoc(object)
    if (n > 0 && p > 0) {
      nNA <- sum(is.na(fbm[]))
      cat("\n\n    missing data: ",
          nNA,
          " (=",
          round((nNA / (n * p)) * 100, 2),
          " %) scored as NA",
          sep = "")
    }
  }
  
  ## BASIC CONTENT
  cat("\n\n ** Genetic data")
  if (!is.null(object@gen))  cat("\n   @gen: list of", length(object@gen), "SNPbin")
  if (!is.null(fbm))  cat("\n   @fbm: fbm.code256 object storing genotypes of ",nrow(object@fbm)," individuals and ",ncol(object@fbm)," loci")
  
  if (!is.null(object@ploidy)) {
    ploidytxt <-
      paste("(range: ", paste(range(object@ploidy), collapse = "-"), ")", sep =
              "")
    cat("\n   @ploidy: ploidy of each individual ", ploidytxt)
  }
  
  ## Additional data
  cat("\n\n ** Additional data")
  optional <- FALSE
  
  if (!is.null(object@ind.names)) {
    optional <- TRUE
    cat("\n   @ind.names: ",
        length(object@ind.names),
        "individual labels")
  } else
    cat("\n   @ind.names: ", "no individual labels")
  
  if (!is.null(object@loc.names)) {
    optional <- TRUE
    cat("\n   @loc.names: ", length(object@loc.names), "locus labels")
  } else
    cat("\n   @loc.names: ", "no locus labels")
  
  if (!is.null(object@loc.all)) {
    optional <- TRUE
    cat("\n   @loc.all: ", length(object@loc.all), "allele labels")
  } else
    cat("\n   @loc.all: ", " no allele labels")
  
  if (!is.null(object@chromosome)) {
    optional <- TRUE
    cat("\n   @chromosome: factor storing chromosomes of the", marker)
  }
  
  if (!is.null(object@position)) {
    optional <- TRUE
    cat("\n   @position: integer storing positions of the",
        marker,
        "[within 69 base sequence]")
  }
  if (!is.null(object@pop)) {
    optional <- TRUE
    poptxt <-
      paste("(group size range: ", paste(range(table(object@pop)), collapse =
                                           "-"), ")", sep = "")
    cat("\n   @pop:", paste("population of each individual", poptxt))
  } else
    cat("\n   @pop:", "no population lables for individuals")
  
  if (!is.null(object@strata)) {
    optional <- TRUE
    cat("\n   @strata: ")
    levs <- names(object@strata)
    if (length(levs) > 6) {
      levs <- paste(paste(head(levs), collapse = ", "), "...", sep = ", ")
    } else {
      levs <- paste(levs, collapse = ", ")
    }
    cat("a data frame with",
        length(object@strata),
        "columns (",
        levs,
        ")")
  }
  
  if (!is.null(object@hierarchy)) {
    optional <- TRUE
    cat("\n   @hierarchy:", paste(object@hierarchy, collapse = ""))
  }
  
  if (!is.null(object@other)) {
    optional <- TRUE
    cat("\n   @other: ")
    cat("a list containing: ")
    cat(ifelse(
      is.null(names(object@other)),
      "elements without names",
      paste(names(object@other), collapse = ", ")
    ), "\n")
  }
  
  if (!is.null(object@other$ind.metrics)) {
    optional <- TRUE
    cat("    @other$ind.metrics: ")
    cat(ifelse(
      is.null(names(object@other$ind.metrics)),
      "elements without names",
      paste(names(object@other$ind.metrics), collapse = ", ")
    ), "\n")
  }
  
  if (!is.null(object@other$ind.metrics)) {
    optional <- TRUE
    cat("    @other$loc.metrics: ")
    cat(ifelse(
      is.null(names(object@other$loc.metrics)),
      "elements without names",
      paste(names(object@other$loc.metrics), collapse = ", ")
    ), "\n")
  }
  
  if (!optional) {
    cat("\n   - empty -")
  }
  
  cat("   @other$latlon[g]:")
  if (!is.null(object@other$latlon)) {
    if (nrow(object@other$latlon) == nInd(object)) {
      cat(" coordinates for all individuals are attached")
    } else{
      cat(" number of coordinates does not match number of individuals")
    }
  } else {
    cat(" no coordinates attached")
  }
  cat("\n")
  
}) # end show method


.fbmsub_copy <- function(G, i, j, backingfile = tempfile("geno_"), col_block = 1024L) {
  stopifnot(inherits(G, "FBM.code256"))
  
  nr <- nrow(G); nc <- ncol(G)
  
  # Normalize i, j to integer indices
  if (is.logical(i)) i <- which(rep_len(i, nr))
  if (is.logical(j)) j <- which(rep_len(j, nc))
  if (anyNA(i) || anyNA(j)) stop("Indices i/j contain NAs after normalization.")
  
  G2 <- bigstatsr::FBM.code256(length(i), length(j), code = code)
  
  # Chunk over columns to limit memory footprint
  if (length(j) > 0L) {
    starts <- seq.int(1L, length(j), by = col_block)
    for (s in starts) {
      e <- min(s + col_block - 1L, length(j))
      j_chunk <- j[s:e]
      # pull (all rows, chunked columns), then subset rows
      chunk <- G[, j_chunk]
      if (!is.null(i)) chunk <- chunk[i, , drop = FALSE]
      G2[, s:e] <- chunk
    }
  }
  G2
}

#################
## subset dartR
#################

#' indexing dartR objects correctly...
#'
#' @param x dartR object
#' @param i index for individuals
#' @param j index for loci
#' @param ... other parameters
#' @param pop list of populations to be kept
#' @param treatOther elements in other (and ind.metrics & loci.metrics) as 
#' indexed as well. default: TRUE
#' @param quiet warnings are suppressed. default: TRUE
#' @param drop reduced to a vector if a single individual/loci is selected. 
#' default: FALSE [should never set to TRUE]
#' @return dartR object
#' 
## Helper: subset an FBM.code256 by (rows i, cols j) without materializing everything

## dartR: subset method with FBM support
setMethod("[", signature(x = "dartR", i = "ANY", j = "ANY", drop = "ANY"),
          function(x, i, j, ..., pop = NULL, treatOther = TRUE, quiet = TRUE, drop = FALSE) {
            
            if (drop) stop("drop=TRUE not supported for dartR objects.")
            
            ## Resolve missing i/j
            if (missing(i)) i <- TRUE
            if (missing(j)) j <- TRUE
            
            ## Current dims (works with your nInd/nLoc methods that already prefer FBM)
            ori.n <- nInd(x)
            ori.p <- nLoc(x)
            
            ## Recycle logical vectors
            if (is.logical(i)) i <- rep_len(i, ori.n)
            if (is.logical(j)) j <- rep_len(j, ori.p)
            
            ## Optionally subset by population names
            if (!is.null(pop) && !is.null(pop(x))) {
              i <- .get_pop_inds(x, pop)  # returns integer indices
            }
            
            ## Character j = match locus names
            if (is.character(j)) {
              j <- match(j, x@loc.names, nomatch = 0L)
            }
            
            ## Normalize i/j to integer indices now (keeps logicals if you prefer)
            ii <- if (is.logical(i)) which(i) else i
            jj <- if (is.logical(j)) which(j) else j
            
            if (length(ii)==0 ) {
              stop(error("Subsetting resulted in zero individuals."))
            
            }
            if (length(jj)==0 ) {
              stop(error("Subsetting resulted in zero loci."))
            
            }
            
            ## -------- FBM-backed branch --------
            fbm <- .fbm_or_null(x)
            if (!is.null(fbm)) {
              # Subset FBM first (keeps XOR: @gen remains empty)
              #x@fbm <- .fbmsub_copy(x@fbm, ii, jj, code = x@fbm$code256)
              
                dummyG <-bigstatsr::big_copy(x@fbm, ind.row = ii, ind.col = jj, backingfile = tempfile("geno_"))
                x@fbm <- dummyG
              #unlink temp file
              #bigstatsr::big_unlink(dummyG)
              
              ## Update names / locus/ind attributes
              if (!is.null(x@ind.names)) x@ind.names <- x@ind.names[ii]
              if (!is.null(x@loc.names)) x@loc.names <- x@loc.names[jj]
              if (!is.null(x@chromosome)) x@chromosome <- x@chromosome[jj]
              if (!is.null(x@position))   x@position   <- x@position[jj]
              if (!is.null(x@loc.all))    x@loc.all    <- x@loc.all[jj]
              if (!is.null(x@n.loc))   x@n.loc <- length(jj)
              
              ## Ploidy / pop / strata
              if (!is.null(x@ploidy))  ploidy(x) <- ploidy(x)[ii]
              if (!is.null(pop(x)))    pop(x)    <- factor(pop(x)[ii])
              if (!is.null(x@strata))  x@strata  <- x@strata[ii, , drop = FALSE]
              
              ## other: subset elements shaped by individuals; keep flags
              if (treatOther && length(x@other)) {
                namesOther <- names(x@other)
                flags_tmp  <- x$other$loc.metrics.flags
                counter <- 0L
                f1 <- function(obj) {
                  counter <<- counter + 1L
                  if (!is.null(dim(obj)) && nrow(obj) == ori.n) {
                    obj[ii, , drop = FALSE]
                  } else if (length(obj) == ori.n) {
                    obj <- obj[ii]
                    if (is.factor(obj)) factor(obj) else obj
                  } else {
                    if (!quiet) warning(paste("cannot treat the object", namesOther[counter]))
                    obj
                  }
                }
                x@other <- lapply(x@other, f1)
                x$other$loc.metrics.flags <- flags_tmp
              }
              
              ## Locus-wise metrics in @other
              if (!is.null(x@other$loc.metrics)) {
                x@other$loc.metrics <- x@other$loc.metrics[jj, , drop = FALSE]
              }
              
              ## Ensure XOR: keep @gen empty
              x@gen <- vector("list", 0L)
              
              return(x)
            }
            
            ## -------- SNPbin (genlight) branch (original logic) --------
            ## genotypes
            x@gen <- x@gen[i]
            
            ## ind names
            x@ind.names <- x@ind.names[ii]
            
            ## ploidy
            if (!is.null(x@ploidy)) {
              ploidy(x) <- ploidy(x)[ii]
            }
            
            ## pop
            if (!is.null(pop(x))) {
              pop(x) <- factor(pop(x)[ii])
            }
            
            ## strata
            if (!is.null(x@strata)) {
              x@strata <- x@strata[ii, , drop = FALSE]
            }
            
            ## HANDLE 'OTHER' SLOT ##
            nOther <- length(other(x))
            namesOther <- names(other(x))
            flags_tmp <- x$other$loc.metrics.flags
            counter <- 0L
            if (treatOther && !(is.logical(i) && all(i))) {
              f2 <- function(obj, n = ori.n) {
                counter <<- counter + 1L
                if (!is.null(dim(obj)) && nrow(obj) == ori.n) {
                  obj[ii, , drop = FALSE]
                } else if (length(obj) == ori.n) {
                  obj <- obj[ii]
                  if (is.factor(obj)) factor(obj) else obj
                } else {
                  if (!quiet) warning(paste("cannot treat the object", namesOther[counter]))
                  obj
                }
              }
              other(x) <- lapply(x@other, f2)
              x$other$loc.metrics.flags <- flags_tmp
            }
            
            ## SUBSET LOCI ##
            x@loc.names  <- x@loc.names[jj]
            x@chromosome <- chr(x)[jj]
            x@position   <- position(x)[jj]
            x@loc.all    <- alleles(x)[jj]
            x@gen        <- lapply(x@gen, function(e) e[jj])
            x@n.loc      <- x@gen[[1]]@n.loc
            
            if (!is.null(x@other$loc.metrics))
              x@other$loc.metrics <- x@other$loc.metrics[jj, , drop = FALSE]
            if (!.has_fbm(x)) x@fbm <- NULL
            
            x
          })

###############################################################
#' @name cbind.dartR
#' @title cbind for dartR objects
#' @description cbind is a bit lazy and does not take care for the metadata (so data in the
#' other slot is lost). You can get most of the loci metadata back using
#' gl.compliance.check.
#' @param ... list of dartR objects
#' @param backingfile prefix for the backing file of the resulting FBM
#' @param code code mapping to use for the resulting FBM=CODE_012; if NULL, inherits from the first FBM input
#' @param chunk number of columns to process in a block when copying from FBMs
#' @param quiet suppress warnings. default: TRUE
#' @examples
#' t1 <- platypus.gl
#' class(t1) <- "dartR"
#' t2 <- cbind(t1[,1:10],t1[,11:20])
#' @return A genlight object
#' @export

cbind.dartR <- function(...,
                        backingfile = tempfile("geno_"),
                        code = NULL,          # if NULL: inherit from first FBM or use CODE_DOSAGE
                        chunk = 2048L,        # columns per block to copy from FBMs
                        quiet = TRUE) {
  
  dots <- list(...)
  
  ## Accept both dartR and genlight; flatten a single list input
  objs <- dots[sapply(dots, function(x) inherits(x, "genlight"))]
  if (length(objs) == 1L && is.list(objs[[1L]])) objs <- objs[[1L]]
  if (!length(objs)) stop("No dartR/genlight objects supplied to cbind().")
  
  ## Remove empty objects (0 loci or 0 inds)
  objs <- objs[sapply(objs, nLoc) > 0 & sapply(objs, nInd) > 0]
  if (!length(objs)) {
    if (!quiet) message("All objects are empty; returning NULL.")
    return(NULL)
  }
  
  ## Basic checks
  n.ind <- unique(sapply(objs, nInd))
  if (length(n.ind) != 1L) stop("Objects have different numbers of individuals.")
  n.ind <- n.ind[[1]]
  
  ## Check ploidy consistency (keep ploidy of first object)
  ploidy_mat <- as.matrix(as.data.frame(lapply(objs, ploidy)))
  if (nrow(ploidy_mat) && any(apply(ploidy_mat, 1, function(r) length(unique(r))) > 1))
    stop("Non-consistent ploidy across datasets.")
  ploidy_out <- ploidy(objs[[1L]])
  
  ## Determine if any object has FBM present
  has_fbm_vec <- vapply(objs, function(o) !is.null(.fbm_or_null(o)), logical(1))
  any_fbm <- any(has_fbm_vec)
  
  ## Collect names/metadata to set on the result
  ind_names <- indNames(objs[[1L]])
  loc_names <- make.unique(unlist(lapply(objs, locNames)))
  alleles_out <- unlist(lapply(objs, alleles), use.names = FALSE)
  pop_out <- pop(objs[[1L]])          # keep pop from the first; warn if differs
  if (any(vapply(objs, function(o) !identical(pop(o), pop_out), logical(1)))) {
    if (!quiet) warning("Populations differ across inputs; using pop() from the first object.")
  }
  strata_out <- objs[[1L]]@strata
  
  ## --- PATH 1: FBM present in at least one input = build a new FBM target ----
  if (any_fbm) {
    total_p <- sum(sapply(objs, nLoc))
    
    ## Choose code
    if (is.null(code)) {
      first_fbm <- .fbm_or_null(objs[[which(has_fbm_vec)[1L]]])
      if (!is.null(first_fbm)) {
        code <- first_fbm$code256
        ## ensure all FBMs share the same code
        mismatch <- which(has_fbm_vec & !vapply(objs, function(o) identical(.fbm_or_null(o)$code256, code), logical(1)))
        if (length(mismatch)) stop("FBMs have different code mappings; specify a common 'code=' or recode.")
      } else {
        ## no FBM actually? (shouldn't happen) default to dosage
        if (!requireNamespace("bigsnpr", quietly = TRUE))
          stop("Package 'bigsnpr' required to set default code (CODE_DOSAGE).")
        code <- bigsnpr::CODE_012
      }
    }
    #    code <- bigsnpr::CODE_012
    ## Create target FBM
    Gnew <- bigstatsr::FBM.code256(n.ind, total_p, code = code, backingfile = backingfile)
    
    ## Fill columns from each object using big_apply
    col_off <- 0L
    for (k in seq_along(objs)) {
      ok    <- objs[[k]]
      p.k   <- nLoc(ok)
      fbm.k <- .fbm_or_null(ok)
      
      if (!is.null(fbm.k)) {
        ## FBM -> FBM copy in column blocks
        bigstatsr::big_apply(
          fbm.k,
          a.FUN = function(X, ind, ind.col) {
            blk <- X[ind, ind.col]                 # decoded chunk (n.ind x length(ind.col))
            # If you want to store the raw missing byte explicitly, uncomment:
            blk[is.na(blk)] <- 3
            Gnew[ind, col_off + ind.col] <- blk    # write into shifted columns of Gnew
            NULL
          },
          a.combine  = "c",
          ind        = bigstatsr::rows_along(fbm.k),
          ind.col    = bigstatsr::cols_along(fbm.k),
          block.size = max(1L, as.integer(chunk))
        )
      } else {
        ## Non-FBM input: materialize once, then write in blocks via big_apply on Gnew
        Xk <- as.matrix(ok)
        bigstatsr::big_apply(
          Gnew,
          a.FUN = function(Y, ind, ind.col) {
            blk <- Xk[ind, ind.col, drop = FALSE]
             blk[is.na(blk)] <- 3
            Y[ind, col_off + ind.col] <- blk
            NULL
          },
          a.combine  = "c",
          ind        = seq_len(n.ind),
          ind.col    = seq_len(p.k),
          block.size = max(1L, as.integer(chunk))
        )
      }
      
      col_off <- col_off + p.k
    }
    
    ## Build output dartR with XOR: fbm set, gen empty
    out <- new("dartR",
               gen        = vector("list", 0L),
               ind.names  = ind_names,
               loc.names  = loc_names,
               ploidy     = ploidy_out,
               pop        = pop_out,
               other      = objs[[1L]]@other,  # keep 'other' from first; adjust below
               strata     = strata_out,
               chromosome = if (!is.null(objs[[1L]]@chromosome)) do.call(c, lapply(objs, chr)) else NULL,
               position   = if (!is.null(objs[[1L]]@position))   do.call(c, lapply(objs, position)) else NULL,
               loc.all    = alleles_out,
               fbm        = Gnew)
    
    ## Merge @other$loc.metrics if present and compatible (row-bind by loci)
    lm_list <- lapply(objs, function(o) o@other$loc.metrics)
    if (all(vapply(lm_list, function(dm) is.data.frame(dm) && nrow(dm) > 0, logical(1)))) {
      ## try to rbind on common columns
      common_cols <- Reduce(intersect, lapply(lm_list, names))
      if (length(common_cols)) {
        out@other$loc.metrics <- do.call(rbind, lapply(lm_list, function(dm) dm[, common_cols, drop = FALSE]))
        rownames(out@other$loc.metrics) <- loc_names
      } else if (!quiet) {
        warning("loc.metrics columns differ; not merging.")
      }
    }
    
    return(out)
  }
  
  ## --- PATH 2: No FBM at all = fall back to original SNPbin-based cbind ----
  myList <- objs
  
  ## (This part mirrors your original implementation)
  if (length(myList) == 1 && is.list(myList[[1]])) myList <- myList[[1]]
  if (!all(sapply(myList, function(x) inherits(x, "genlight"))))
    stop("Some objects are not genlight/dartR objects")
  
  ## remove empty again (defensive)
  myList <- myList[sapply(myList, nLoc) > 0 & sapply(myList, nInd) > 0]
  if (length(myList) == 0) {
    if (!quiet) message("All objects are empty")
    return(NULL)
  }
  
  ## consistency already checked above
  ori.ploidy <- ploidy(myList[[1]])
  
  ## merge one individual at a time (SNPbin concatenation)
  res <- vector("list", length = n.ind)
  for (i in seq_len(n.ind)) {
    res[[i]] <- Reduce(function(a, b) { cbind(a, b, checkPloidy = FALSE) },
                       lapply(myList, function(e) e@gen[[i]]))
  }
  
  ## build output
  out <- new("dartR",
             gen        = res,
             ind.names  = indNames(myList[[1]]),
             loc.names  = make.unique(unlist(lapply(myList, locNames))),
             ploidy     = ori.ploidy,
             pop        = pop(myList[[1]]),
             other      = myList[[1]]@other,
             strata     = myList[[1]]@strata,
             chromosome = do.call(c, lapply(myList, chr)),
             position   = do.call(c, lapply(myList, position)),
             loc.all    = unlist(lapply(myList, alleles)),
             fbm        = NULL)
  
  ## (Optionally) merge loc.metrics similarly as FBM path
  if (!is.null(myList[[1]]@other$loc.metrics)) {
    lm_list <- lapply(myList, function(o) o@other$loc.metrics)
    if (all(vapply(lm_list, function(dm) is.data.frame(dm) && nrow(dm) > 0, logical(1)))) {
      common_cols <- Reduce(intersect, lapply(lm_list, names))
      if (length(common_cols)) {
        out@other$loc.metrics <- do.call(rbind, lapply(lm_list, function(dm) dm[, common_cols, drop = FALSE]))
        rownames(out@other$loc.metrics) <- locNames(out)
      }
    }
  }
  
  out
} # end cbind.dartR
###############################################################
##################################################################
#' @name rbind.dartR
#' @title rbind for dartR objects
#' @description rbind is a bit lazy and does not take care for the metadata (so data in the
#' other slot is lost). You can get most of the loci metadata back using
#' gl.compliance.check.
#' @param ... list of dartR objects
#' @param backingfile prefix for the backing file of the resulting FBM
#' @param code code mapping to use for the resulting FBM=CODE_012; if NULL, inherits from the first FBM input
#' @param chunk number of columns to process in a block when copying from FBMs
#' @param quiet suppress warnings. default: TRUE
#' @examples
#' t1 <- platypus.gl
#' class(t1) <- "dartR"
#' t2 <- rbind(t1[1:5,],t1[6:10,])
#' @return A genlight object 
#' @export
rbind.dartR <- function(...,
                        backingfile   = tempfile("geno_rbind_"),
                        code          = NULL,        # if NULL: inherit from first FBM
                        chunk         = 2048L,       # columns per block
                        quiet         = TRUE) {
  dots <- list(...)
  objs <- dots[sapply(dots, function(x) inherits(x, "genlight"))]
  if (length(objs) == 1L && is.list(objs[[1L]])) objs <- objs[[1L]]
  if (!length(objs)) stop("No dartR/genlight objects supplied to rbind().")
  
  ## Drop empties
  objs <- objs[sapply(objs, nLoc) > 0 & sapply(objs, nInd) > 0]
  if (!length(objs)) {
    if (!quiet) message("All objects are empty; returning NULL.")
    return(NULL)
  }
  
  ## Reference locus set (from first object)
  p.ref <- nLoc(objs[[1L]])
  ref_loc <- locNames(objs[[1L]])
  if (is.null(ref_loc)) ref_loc <- paste0("l", seq_len(p.ref))
  
  ## Build per-object column mapping to reference
  colmap <- vector("list", length(objs))
  for (k in seq_along(objs)) {
    ok <- objs[[k]]
    if (nLoc(ok) != p.ref) {
      stop("Objects have different numbers of loci (", nLoc(ok), " vs ", p.ref, ").")
    }
    lk <- locNames(ok)
    if (!is.null(lk) && !is.null(ref_loc)) {
      idx <- match(ref_loc, lk)
      if (anyNA(idx)) {
        stop("Locus names don't match between objects; cannot align by names.")
      }
      colmap[[k]] <- idx
    } else {
      ## Assume identical order
      colmap[[k]] <- seq_len(p.ref)
    }
  }
  
  ## Concatenate per-individual metadata
  ind_names_in <- unlist(lapply(objs, indNames), use.names = FALSE)
  ind_names    <- make.unique(ind_names_in, "__ind")
  ploidy_out   <- unlist(lapply(objs, ploidy), use.names = FALSE)
  
  ## Pop / strata (row-wise)
  pop_list   <- lapply(objs, pop)
  pop_out    <- if (all(vapply(pop_list, function(p) !is.null(p), logical(1)))) {
    factor(unlist(pop_list), levels = unique(unlist(lapply(pop_list, levels))))
  } else NULL
  
  strata_out <- {
    has_str <- vapply(objs, function(o) !is.null(o@strata), logical(1))
    if (any(has_str)) {
      # simple base rbind on common columns
      dfs <- lapply(objs[has_str], function(o) o@strata)
      common <- Reduce(intersect, lapply(dfs, names))
      if (length(common)) do.call(rbind, lapply(dfs, function(d) d[, common, drop = FALSE])) else NULL
    } else NULL
  }
  
  ## Locus-wise metadata: keep from reference (first object)
  loc_names <- ref_loc
  alleles_out  <- alleles(objs[[1L]])
  chrom_out    <- chr(objs[[1L]])
  pos_out      <- position(objs[[1L]])
  loc_metrics  <- objs[[1L]]@other$loc.metrics
  other_out    <- objs[[1L]]@other
  if (!is.null(loc_metrics) && nrow(loc_metrics) == length(loc_names)) {
    rownames(loc_metrics) <- loc_names
    other_out$loc.metrics <- loc_metrics
  }
  
  ## Any FBM present?
  has_fbm_vec <- vapply(objs, function(o) !is.null(.fbm_or_null(o)), logical(1))
  any_fbm <- any(has_fbm_vec)
  
  if (any_fbm) {
    total_n <- sum(sapply(objs, nInd))
    p       <- p.ref
    
    ## Choose code
    if (is.null(code)) {
      first_fbm <- .fbm_or_null(objs[[which(has_fbm_vec)[1L]]])
      code <- first_fbm$code256
      mismatch <- which(has_fbm_vec & !vapply(objs, function(o) identical(.fbm_or_null(o)$code256, code), logical(1)))
      if (length(mismatch)) stop("FBMs have different code mappings; supply a common 'code=' or recode inputs.")
    }
    
    ## Target FBM
    Gnew <- bigstatsr::FBM.code256(total_n, p, code = code)
    
    ## Copy each object into row slices of Gnew using big_apply
    row_off <- 0L
    for (k in seq_along(objs)) {
      ok     <- objs[[k]]
      n.k    <- nInd(ok)
      fbm.k  <- .fbm_or_null(ok)
      cmap   <- colmap[[k]]              # fbm columns to read (in ref order)
      
      if (!is.null(fbm.k)) {
        ## Read from source FBM in column blocks; write to destination rows with offset
        pos <- integer(ncol(fbm.k)); pos[cmap] <- seq_len(p)   # map source col index -> destination col
        bigstatsr::big_apply(
          fbm.k,
          a.FUN = function(X, ind, ind.col) {
            ## ind is 1:n.k; ind.col is a subset of cmap
            blk <- X[ind, ind.col]
            blk[is.na(blk)] <- 3
            dest_cols <- pos[ind.col]  # positions in 1:p
            Gnew[row_off + ind, dest_cols] <- blk
            integer(0)  # keep 'c' combiner happy
          },
          a.combine  = "c",
          ind        = bigstatsr::rows_along(fbm.k),
          ind.col    = cmap,
          block.size = max(1L, as.integer(chunk))
        )
      } else {
        ## Non-FBM: materialize matrix (reorder columns if needed), then write via big_apply on Gnew
        Xk <- as.matrix(ok)
        Xk[is.na(Xk)] <- 3
        if (!identical(colnames(Xk), ref_loc)) {
          Xk <- Xk[, cmap, drop = FALSE]
        }
        dest_rows <- (row_off + seq_len(n.k))
        bigstatsr::big_apply(
          Gnew,
          a.FUN = function(Y, ind, ind.col) {
            Y[dest_rows, ind.col] <- Xk[, ind.col, drop = FALSE]
            integer(0)
          },
          a.combine  = "c",
          ind        = dest_rows,
          ind.col    = seq_len(p),
          block.size = max(1L, as.integer(chunk))
        )
      }
      
      row_off <- row_off + n.k
    }
    
    ## Build output (XOR: fbm set, gen empty)
    out <- methods::new("dartR",
                        gen        = vector("list", 0L),
                        fbm        = Gnew,
                        ind.names  = ind_names,
                        loc.names  = loc_names,
                        ploidy     = ploidy_out,
                        pop        = pop_out,
                        other      = other_out,
                        strata     = strata_out,
                        chromosome = chrom_out,
                        position   = pos_out,
                        loc.all    = alleles_out
    )
    if (!.has_fbm(out)) out@fbm <- NULL
    methods::validObject(out)
    return(out)
  }
  
  ## ---------------- No FBM: SNPbin-based fallback (original logic) ----------------
  myList <- objs
  if (length(unique(sapply(myList, nLoc))) != 1L)
    stop("objects have different numbers of SNPs")
  
  dots <- list()
  dots$Class <- "dartR"
  dots$gen <- Reduce(c, lapply(myList, function(e) e@gen))
  out <- do.call(methods::new, dots)
  locNames(out) <- locNames(myList[[1L]])
  alleles(out)  <- alleles(myList[[1L]])
  indNames(out) <- ind_names
  pop(out)      <- if (!is.null(pop_out)) pop_out else pop(myList[[1L]])
  out@strata    <- strata_out
  out@chromosome <- chr(myList[[1L]])
  out@position   <- position(myList[[1L]])
  ploidy(out)    <- ploidy_out
  if (!.has_fbm(out)) out@fbm <- NULL
  methods::validObject(out)
  out
}
# end of rbind.dartR
###############################################################

#now set methods to use FBM when present.
as.matrix.dartR <- function(x, ...) {
  fbm <- .fbm_or_null(x)
  if (!is.null(fbm)) return(x@fbm[])
  # fall back to genlight
  methods::as(methods::as(x, "genlight"), "matrix")
}

if (!methods::isGeneric("nInd"))
  methods::setGeneric("nInd", function(x) standardGeneric("nInd"))
if (!methods::isGeneric("nLoc"))
  methods::setGeneric("nLoc", function(x) standardGeneric("nLoc"))

## Use FBM dims when present
setMethod("nInd", "dartR", function(x) {
  fbm <- .fbm_or_null(x)
  if (is.null(fbm)) callNextMethod() else  nrow(x@fbm) 
})
setMethod("nLoc", "dartR", function(x) {
  fbm <- .fbm_or_null(x)
  if (is.null(fbm)) callNextMethod() else ncol(x@fbm)
})

if (!methods::isGeneric("as.matrix")) {
  methods::setGeneric("as.matrix", useAsDefault = base::as.matrix)
}

setMethod("as.matrix", "dartR", function(x, ...) {
  fbm <- .fbm_or_null(x)
  if (is.null(fbm)) methods::as(methods::as(x, "genlight"), "matrix") else {
    
    dummy <- x@fbm[, drop=FALSE]
    colnames(dummy)<- locNames(x)
    rownames(dummy)<- indNames(x)
    return(dummy) 
    }
})
# 3) Also register an explicit coercion (used by some internals)
methods::setAs("dartR", "matrix", function(from) {
  if (!is.null(.fbm_or_null(from))) return(from@fbm[])
  methods::as(methods::as(from, "genlight"), "matrix")
})

############glSum############################################################################

## Ensure the generic exists (adegenet defines it; this is safe if already present)
#if (!methods::isGeneric("glSum")) {
#  methods::setGeneric("glSum", function(x, alleleAsUnit = TRUE, useC=FALSE) standardGeneric("glSum"))
#}



#' @name glSum
#' @title glSum for dartR objects
#' @description glSum is necessary as adegenet is using it internally and we need one for fbm projects
#' @param x a dartR object 
#' @param alleleAsUnit logical; if TRUE, the mean is calculated per allele,
#' if FALSE, per individual
#' @param useC FALSE, default if set to true not sure what happens ;-)
#' @return A numeric vector of sum of second allele per locus
#' @export
glSum <- function(x, alleleAsUnit = TRUE, useC=FALSE) {
  fbm <- .has_fbm(x)

  ## If no FBM (old object or genlight mode), delegate to next method (genlight)
  if (!fbm){
    class(x)<- "genlight"
    res <- adegenet::glSum(x, alleleAsUnit = alleleAsUnit, useC=useC); return(res)
  }


## ---- FBM-aware glSum for dartR ----

  
  if (alleleAsUnit) {
    res <- integer(nLoc(x))
    for (e in 1:nInd(x)) {
      temp <- as.integer(x@fbm[e,])
      temp[is.na(temp)] <- 0L
      res <- res + temp
    }
  }  else {
    res <- numeric(nLoc(x))
    myPloidy <- ploidy(x)
    for (i in 1:nInd(x)) {
      temp <- as.integer(x@fbm[i,])/myPloidy[i]
      temp[is.na(temp)] <- 0
      res <- res + temp
    }
  }

  names(res) <- locNames(x)
  return(res)
}

setMethod("glNA", signature(x = "dartR"), function(x, alleleAsUnit = TRUE)  {
  fbm <- .fbm_or_null(x)
  ## If no FBM (old object or genlight mode), delegate to next method (genlight)
  if (is.null(fbm)) return(callNextMethod())
  res <- integer(nLoc(x))
  temp <- NA.posi(x)
  if (alleleAsUnit) {
    for (i in 1:length(temp)) {
      if (length(temp[[i]]) > 0) {
        res[temp[[i]]] <- res[temp[[i]]] + ploidy(x)[i]
      }
    }
  }
  else {
    for (e in temp) {
      if (length(e) > 0) {
        res[e] <- res[e] + 1
      }
    }
  }
  names(res) <- locNames(x)
  return(res)
})

#if (!methods::isGeneric("glMean")) {
#  methods::setGeneric("glMean", function(x, ...) standardGeneric("glMean"))
#}
#setMethod("glMean", signature(x = "dartR"),  function(x, alleleAsUnit = TRUE) {
#' @name glMean
#' @title glMean for dartR objects
#' @description glMean is necessary as adegenet is using it internally and we need one for fbm projects
#' @param x a dartR object 
#' @param alleleAsUnit logical; if TRUE, the mean is calculated per allele,
#' if FALSE, per individual
#' @return A numeric vector of means per locus
#' @export
glMean <- function(x, alleleAsUnit = TRUE) {
  fbm <- .has_fbm(x)
  ## If no FBM (old object or genlight mode), delegate to next method (genlight)
  if (!fbm){
    class(x)<- "genlight"
    res <- adegenet::glMean(x, alleleAsUnit = alleleAsUnit); return(res)
    }
    
  
  if (alleleAsUnit) {
    
    N <- sum(ploidy(x)) - glNA(x, alleleAsUnit = TRUE)
    res <- glSum(x, alleleAsUnit = TRUE)/N
  } else {
    N <- nInd(x) - glNA(x, alleleAsUnit = FALSE)
    res <- glSum(x, alleleAsUnit = FALSE)/N
  }
  names(res) <- locNames(x)
  return(res)
}



#############NA.posi##########################################################
## Ensure the generic exists (adegenet defines it; safe if already present)
if (!methods::isGeneric("NA.posi")) {
  methods::setGeneric("NA.posi", function(x, ...) standardGeneric("NA.posi"))
}

setMethod("NA.posi", signature(x = "dartR"), function(x, ...) {
  fbm <- .fbm_or_null(x)
  if (is.null(fbm)) return(callNextMethod())
  
  n <- nInd(x)
  p <- nLoc(x)
  
  res <- vector("list", n)
  nm  <- if (!is.null(x@ind.names) && length(x@ind.names) == n) x@ind.names else NULL
  if (!is.null(nm)) names(res) <- nm
  
  ## optional chunk size via 'chunk=' in ...
  dots  <- list(...)
  chunk <- if (!is.null(dots$chunk)) as.integer(dots$chunk) else 4096L
  if (!is.finite(chunk) || chunk <= 0L) chunk <- 4096L
  
  for (s in seq.int(1L, p, by = chunk)) {
    e <- min(s + chunk - 1L, p)
    
    blk <- fbm[, s:e]  # may be vector if s==e, depending on FBM subsetting
    
    ## CRITICAL: enforce 2D shape (n x (e-s+1)) even when (e-s+1)==1
    if (is.null(dim(blk))) {
      blk <- matrix(blk, nrow = n, ncol = e - s + 1L)
    }
    
    NAblk <- is.na(blk)
    
    ## Collect missing positions only where they exist (more efficient than scanning all rows)
    pos <- which(NAblk, arr.ind = TRUE)  # columns: row, col (within block)
    if (nrow(pos) > 0L) {
      rows <- pos[, 1L]
      cols <- pos[, 2L] + (s - 1L)       # offset to locus index in full matrix
      by_row <- split(cols, rows)
      
      for (ri in names(by_row)) {
        i <- as.integer(ri)
        res[[i]] <- c(res[[i]], by_row[[ri]])
      }
    }
  }
  
  ## Ensure integer vectors (empty = integer())
  res <- lapply(res, function(v) if (length(v)) as.integer(v) else integer())
  res
})

## KEY PART: initialize method
setMethod("initialize", "dartR", function(.Object, ...) {
  args <- list(...)
  has_fbm <- !is.null(args$fbm)
  
  if (has_fbm) {
    ## Build a minimal genlight-like shell WITHOUT letting genlight's initialize reset things
    ## (i.e., do NOT call callNextMethod() here)
    slots <- slotNames(.Object)
    fill <- function(nm, default = NULL) if (!is.null(args[[nm]])) args[[nm]] else default
    
    ## Set parent slots explicitly
    if ("gen"        %in% slots) methods::slot(.Object, "gen")        <- vector("list", 0L)  # XOR: empty
    if ("ind.names"  %in% slots) methods::slot(.Object, "ind.names")  <- fill("ind.names")
    if ("loc.names"  %in% slots) methods::slot(.Object, "loc.names")  <- fill("loc.names")
    if ("ploidy"     %in% slots) methods::slot(.Object, "ploidy")     <- fill("ploidy")
    if ("pop"        %in% slots) methods::slot(.Object, "pop")        <- fill("pop")
    if ("other"      %in% slots) methods::slot(.Object, "other")      <- fill("other", list())
    if ("strata"     %in% slots) methods::slot(.Object, "strata")     <- fill("strata")
    if ("chromosome" %in% slots) methods::slot(.Object, "chromosome") <- fill("chromosome")
    if ("position"   %in% slots) methods::slot(.Object, "position")   <- fill("position")
    if ("loc.all"    %in% slots) methods::slot(.Object, "loc.all")    <- fill("loc.all")
    if ("n.loc"      %in% slots && !is.null(args$fbm)) methods::slot(.Object, "n.loc") <- ncol(args$fbm)
    
    ## Set our FBM slot last
    methods::slot(.Object, "fbm") <- args$fbm
    
    ## Done
    return(.Object)
  }
  
  ## No FBM provided -> normal genlight initialization
  .Object <- methods::callNextMethod()
  .Object
})







