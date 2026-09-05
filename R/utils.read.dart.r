#' @name utils.read.dart
#' @title Utility to import DarT data to R
#' @family io

#' @description 
#' WARNING: UTILITY SCRIPTS ARE FOR INTERNAL USE ONLY AND SHOULD NOT BE USED BY END USERS AS THEIR USE OUT OF CONTEXT COULD LEAD TO UNPREDICTABLE OUTCOMES.

#' @param filename Path to file (csv file only currently) [required].
#' @param nas A character specifying NAs [default '-'].
#' @param topskip A number specifying the number of rows to be skipped. If not
#' provided the number of rows to be skipped are 'guessed' by the number of rows
#' with '*' at the beginning [default NULL].
#' @param service.row The row number in which the information of the DArT
#' service is contained [default 1].
#' @param plate.row The row number in which the information of the plate
#' location is contained [default 3].
#' @param lastmetric Specifies the last non genetic column [default 'RepAvg'].
#' Be sure to check if that is true, otherwise the number of individuals will
#' not match. You can also specify the last column by a number.
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' progress log ; 3, progress and results summary; 5, full report [default NULL].
#' 
#' @details
#' Internal function called by gl.read.dart()
#' 
#' @author Custodian: Bernd Gruber (Post to \url{https://groups.google.com/d/forum/dartr})
#' 
# @export
#' @return A list of length 5. #dart format (one or two rows) #individuals,
#' #snps, #non genetic metrics, #genetic data (still two line format, rows=snps,
#'  columns=individuals)


utils.read.dart <- function(filename,
                            nas = "-",
                            topskip = NULL,
                            lastmetric = "RepAvg",
                            service.row = 1,
                            plate.row = 3,
                            verbose = NULL) {
    # SET VERBOSITY
    verbose <- gl.check.verbosity(verbose)
    
    # FLAG SCRIPT START
    funname <- match.call()[[1]]
    utils.flag.start(func = funname,
                     build = "v.2023.2",
                     verbose = verbose)
    
    # DO THE JOB
    
    if (is.null(topskip)) {
        if (verbose >= 2) {
            cat(report("  Topskip not provided.\n "))
        }
        tdummy <-
            read.csv(
                filename,
                na.strings = nas,
                check.names = FALSE,
                nrows = 20,
                header = FALSE,
                stringsAsFactors = TRUE
            )
        
        nskip <- sum(tdummy[, 1] == "*")
        if (nskip > 0) {
            topskip <- nskip
            if (verbose >= 2) {
                cat(report(paste(
                    "Setting topskip to", nskip, ".\n"
                )))
            }
        } else {
            stop(
                error(
                    "Could not determine the number of rows that need to be skipped. Please provide it manually by setting the topskip parameter.\n"
                )
            )
        }
    }
    
    if (verbose >= 2) {
        cat(report("  Reading in the SNP data\n"))
    }
    snpraw <-
        read.csv(
            filename,
            na.strings = nas,
            skip = topskip,
            check.names = FALSE,
            stringsAsFactors = TRUE
        )
    
    if (is.character(lastmetric)) {
        lmet <- which(lastmetric == names(snpraw))
        if (length(lmet) == 0) {
            stop(error(
                paste(
                    "Could not determine number of data columns based on",
                    lastmetric,
                    "!\n"
                )
            ))
        }
    } else {
        lmet <- lastmetric
    }
    
    service <- NA
    plate_location <- NA
    if(exists("tdummy")){
    # The header block holds topskip rows; indexing tdummy beyond it would
    # paste the column-header row and data rows into the service/plate
    # fields, so out-of-range rows yield NA instead.
    if (service.row <= topskip) {
    # extracting service information
    service <- tdummy[service.row, (lmet + 1):ncol(tdummy)]
    } else if (verbose >= 1) {
        cat(warn(
            "  Warning: service.row", service.row, "lies beyond the",
            topskip, "header rows; service set to NA.\n"
        ))
    }
    if ((plate.row + 2) <= topskip) {
    # extracting plate information
    plate <-
        unlist(unname(tdummy[plate.row, (lmet + 1):ncol(tdummy)]))
    plate.row_res <-
        unlist(unname(tdummy[(plate.row + 1), (lmet + 1):ncol(tdummy)]))
    plate_col_res <-
        unlist(unname(tdummy[(plate.row + 2), (lmet + 1):ncol(tdummy)]))
    plate_location <-
        paste0(plate, "-", plate.row_res, plate_col_res)
    } else if (verbose >= 1) {
        cat(warn(
            "  Warning: plate.row", plate.row, "plus its two following rows lie beyond the",
            topskip, "header rows; plate location set to NA.\n"
        ))
    }
    }
    
    ind.names <- colnames(snpraw)[(lmet + 1):ncol(snpraw)]
    ind.names <-
        trimws(ind.names, which = "both")  #trim for spaces
    if (length(ind.names) != length(unique(ind.names))) {
        if (verbose >= 1) {
            cat(
                warn(
                    "Warning: Individual names are not unique, adding '_n' to replicates (but not the first instance) to render them unique.\n"
                )
            )
        }
        ind.names <-
            make.unique(as.character(ind.names), sep = "_")
    }
    
    stdmetricscols <- 1:lmet
    # if (length(stdmetricscols) != length(stdmetrics)) { cat(paste('\nCould not find all standard metrics.\n',stdmetrics[!(stdmetrics
    # %in% names(snpraw) )] ,' is missing.\n Carefully check the spelling of your headers!\n')) stop() } if (!is.null(addmetrics)) {
    # addmetricscols <- which( names(snpraw) %in% addmetrics ) if (length(addmetricscols) != length(addmetrics)) { cat(paste('\nCould not
    # find all additional metrics.\n',addmetrics[!(addmetrics %in% names(snpraw) )] ,' is missing.\n Carefully check the spelling of your
    # headers! or set addmetrics to NULL\n')) stop() } stdmetricscols <- c(stdmetricscols, addmetricscols) }
    
    covmetrics <- snpraw[, stdmetricscols]
    
    ##### Various checks (first are there two rows per allele?  we do not need cloneid any more....  covmetrics$CloneID =
    ##### as.character(covmetrics$CloneID) check that there are two lines per locus... covmetrics = separate(covmetrics, CloneID, into =
    ##### c('clid','clrest'),sep = '\\|', extra='merge')
    
    # covmetrics$AlleleID = as.character(covmetrics$AlleleID)
    
    # check that there are two lines per locus... covmetrics = separate(covmetrics, AlleleID, into = c('allid','alrest'),sep = '\\|',
    # extra='merge')
    metrics_names <- colnames(covmetrics)
    
    if("AlleleID" %in% metrics_names){
      covmetrics$clone <- (sub("\\|.*", "", covmetrics$AlleleID, perl = T))
      spp <- (sub(".+-+(\\d{1,3}):.+", "\\1", covmetrics$AlleleID))
    }
    
    if("MarkerName" %in% metrics_names){
      covmetrics$clone <- (sub("\\|.*", "", covmetrics$MarkerName, perl = T))
      spp <- (sub(".+-+(\\d{1,3}):.+", "\\1", covmetrics$MarkerName))
    }
  
    #### find uid within allelid
    covmetrics$uid <- paste(covmetrics$clone, spp, sep = "-")
    ### there should be only twos (and maybe fours)
    tt <- table(table(covmetrics$uid))

    # Testing for SNPs with the same ID
    # Jason Carling e-mail
    # the marker discovery and calling pipeline which we use (called ds14) has a 
    # mechanism to prevent this situation from occurring. When the clusters are 
    # parsed, the software will not allow more than a single marker with the same 
    # CloneID, SNP position, and SNP variant. However, there are a number of 
    # scenarios I can think of which could occur after the initial marker outputs
    # are produced which could violate this rule, resulting in the situation you 
    # have described. One example could be the combination of markers from multiple
    # software outputs into a single new report file. This can occur either 
    # manually, or using SNP re-calling software which calls the markers based on an 
    # input definition file, rather than running the marker discovery de novo. 
    # Secondly, marker renaming is sometimes undertaken, when newly discovered 
    # sequences are added to our CloneID database, and these IDs are retrospectively
    # incorporated into the marker report to replace temporary IDs. I have not yet 
    # determined where the marker name clash was introduced in this case, but I can 
    # answer your question by saying that unfortunately, we cannot completely 
    # prevent this sort of problem for the reasons indicated above due to the
    # potential for processes which modify names or combine data sets after the 
    # initial marker discovery and output.
    
    if(length(tt)>1){
      # Mixed row-counts per uid. The modal count is the expected number of
      # rows per locus (2 in 2-row format, 1 in 1-row format); only uids
      # departing from it are malformed (e.g. a marker-name clash yielding
      # 4 rows, or an orphaned allele row yielding an odd count). Remove
      # those uids only - never the valid majority. (The previous logic
      # removed every uid with count > 1 whenever the second count class
      # was not 4, which wiped every locus of a 2-row file.)
      uid.counts <- table(covmetrics$uid)
      modal.count <- as.numeric(names(tt)[which.max(tt)])
      bad.uid <- names(uid.counts)[uid.counts != modal.count]
      bad.rows <- which(covmetrics$uid %in% bad.uid)
      # removing SNPs whose row count departs from the modal count
      covmetrics <- covmetrics[-bad.rows,]
      snpraw <- snpraw[-bad.rows,]

      if (verbose >= 1) {
        cat(warn(
          "  There were",length(bad.uid),"SNPs whose number of rows differed from the expected",modal.count,"per locus. These SNPs have been removed. SNP names are:",paste(bad.uid,collapse = " "),"\n"
        ))
      }

      # recompute the row-count table so downstream format checks and
      # reports reflect the cleaned data
      tt <- table(table(covmetrics$uid))
    }
      
    datas <- snpraw[, (lmet + 1):ncol(snpraw)]
    # Apply the trimmed, uniquified individual names to the genotype
    # columns, so the '_n' suffixes promised by the warning above are what
    # downstream code and ind.metafile matching actually see (previously
    # the make.unique result was computed but never applied).
    colnames(datas) <- ind.names

    nrows <-NULL
    if (is.null(nrows)) {
        gnrows <-3 - max(datas, na.rm = TRUE)  #if max(datas==1) then two row format, if two then one row format
        
        if (gnrows == 1 | gnrows == 2) {
            nrows <- gnrows
            if (verbose >= 2) {
                cat(report(paste(
                    "  Detected", nrows, "row format.\n"
                )))
            }
        } else {
            stop(
                error(
                    "The DArT format must be either 1row or 2row. This does not seem to be the case here.\n"
                )
            )
        }
        
    }
    
    # Cross-check the two independent format signals: the genotype range
    # (which set nrows above) and the number of rows per locus ID. A genuine
    # 1-row file whose called genotypes happen to be only 0/1 (no
    # heterozygotes) would otherwise be silently misread as 2-row, pairing
    # unrelated adjacent rows into fabricated genotypes.
    uid.rows <- as.numeric(names(tt)[which.max(tt)])
    if (nrows != uid.rows) {
      stop(error(
        paste(
          "The genotype codes suggest a", nrows, "row format, but each locus ID occurs",
          uid.rows, "time(s), which contradicts that format. The file would be misread.",
          "Check the file: a 1-row report with no heterozygous calls, or a corrupted",
          "2-row report, can produce this disagreement.\n"
        )
      ))
    }
    
    if (verbose >= 2) {
      cat(report(
        paste(
          "Number of rows per clone (should be only ",
          nrows,
          "s):",
          names(tt),
          "\n "
        )
      ))
    }
    
    if (verbose >= 2) {
        cat(report("Added the following locus metrics:\n"))
        cat(report(paste(
            paste(names(snpraw)[stdmetricscols], collapse = " "), ".\n"
        )))
    }
    
    nind <- ncol(datas)
    nsnp <- nrow(covmetrics) / nrows
    
    if (verbose >= 2) {
        cat(important(
            paste(
                "Recognised:",
                nind,
                "individuals and",
                nsnp,
                " SNPs in a",
                nrows,
                "row format using",
                filename,
                "\n"
            )
        ))
    }
    
    out <-
        list(
            nrows = nrows,
            nind = nind,
            nsnp = nsnp,
            covmetrics = covmetrics,
            gendata = datas,
            service = service,
            plate_location = plate_location
        )
    
    # FLAG SCRIPT END
    
    if (verbose >= 1) {
        cat(report("Completed:", funname, "\n"))
    }
    
    return(out)
    
}
