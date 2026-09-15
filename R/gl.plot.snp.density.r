#' @name gl.plot.snp.density
#' @title Plot SNP density along chromosomes (heat-map)
#' @family graphics
#'
#' @description
#' Generates a tiled heat-map of single-nucleotide polymorphisms (SNPs)
#' across chromosomes or scaffolds in a genlight object. SNPs are binned
#' into fixed-width windows and coloured by SNP count, optionally
#' annotating chromosomes with their SNP number and the position of their
#' last SNP (Mb).
#'
#' @details
#' Chromosome names are taken from \code{x@chromosome} and SNP positions
#' from \code{x@position}. Packaged datasets ship with these slots empty;
#' fill them from the locus metrics first, as in the example. Loci with a
#' missing chromosome or a missing or zero position are ignored.
#'
#' A chromosome is plotted when it carries at least \code{min.snps}
#' positioned SNPs and its last SNP lies at or beyond \code{min.length}
#' bp; the "length" used here and in the labels is the position of the
#' last SNP, not the assembly length. Chromosomes are ordered
#' alphabetically, numeric-aware (\code{chr2} before \code{chr10}), with
#' the first name at the top of the plot. Every bin from the start of the
#' chromosome to its last SNP is drawn; bins containing no SNPs are
#' rendered in the lowest colour of the palette. The function does not
#' modify the input genlight object.
#'
#' @param x A genlight object with chromosome names in \code{@chromosome}
#' and SNP positions in \code{@position} [required].
#' @param bin.size Width (bp) of the genomic bins used to count SNPs
#' [default 1e6].
#' @param min.snps Minimum number of positioned SNPs a chromosome must have
#' to be plotted [default 50].
#' @param min.length Minimum position (bp) of the last SNP for a chromosome
#' to be plotted [default 1e6].
#' @param color.palette A function returning a vector of colours, passed to
#' ggplot2; typically viridis::viridis [default viridis::viridis].
#' @param chr.info If TRUE, append (N SNPs, L Mb) to the chromosome labels
#' [default TRUE].
#' @param plot.title Main title for the plot [default NULL].
#' @param plot.theme User specified theme [default theme_dartR()].
#' @param plot.display Specify if plot is to be displayed in the graphics
#' window [default TRUE].
#' @param plot.file Filename (minus extension) for the RDS plot file
#' [Required for plot save].
#' @param plot.dir Directory to save the plot RDS file [default as specified
#' by the global working directory or tempdir()].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' brief progress messages; 3, progress and results summary; 5, full report
#' [default 2, unless specified using gl.set.verbosity].
#'
#' @return A ggplot object (invisibly) displaying the SNP-density
#' heat-map.
#'
#' @author Author(s): Luis Mijangos. Custodian: Luis Mijangos -- Post to
#' \url{https://groups.google.com/d/forum/dartr}
#'
#' @export
#' @examples
#' t1 <- platypus.gl
#' t1$chromosome <- t1$other$loc.metrics$Chrom_Platypus_Chrom_NCBIv1
#' t1$position   <- t1$other$loc.metrics$ChromPos_Platypus_Chrom_NCBIv1
#' gl.plot.snp.density(t1,
#'                     bin.size   = 5e6,
#'                     min.snps   = 10,
#'                     min.length = 2e6,
#'                     plot.title = "Platypus SNP density")

gl.plot.snp.density <- function(x,
                                bin.size      = 1e6,
                                min.snps      = 50,
                                min.length    = 1e6,
                                color.palette = viridis::viridis,
                                chr.info      = TRUE,
                                plot.title    = NULL,
                                plot.theme    = theme_dartR(),
                                plot.display  = TRUE,
                                plot.file     = NULL,
                                plot.dir      = NULL,
                                verbose       = NULL) {
  
  pos <- n_snps <- chr_size <- bin_center <- chr_label <- NULL
  
  # SET VERBOSITY
  verbose <- gl.check.verbosity(verbose)
  
  # FLAG SCRIPT START
  funname <- match.call()[[1]]
  utils.flag.start(func = funname, verbose = verbose)
  
  # CHECK DATATYPE
  datatype <- utils.check.datatype(x, accept = "SNP", verbose = verbose)
  
  # SET WORKING DIRECTORY
  plot.dir <- gl.check.wd(plot.dir, verbose = 0)
  
  # DEPENDENCY CHECKS
  needed_pkgs <- c("ggplot2", "dplyr", "viridis")
  for (pkg in needed_pkgs) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop(error("Package ", pkg,
                 " is required but not installed. Please install it.\n"))
    }
  }
  
  # FUNCTION-SPECIFIC CHECKS
  if (bin.size <= 0) {
    stop(error("Parameter bin.size must be > 0 bp.\n"))
  }
  if (min.snps < 1) {
    stop(error("Parameter min.snps must be >= 1.\n"))
  }
  if (min.length < 1) {
    stop(error("Parameter min.length must be >= 1 bp.\n"))
  }
  
  n_chr <- length(x@chromosome)
  n_pos <- length(x@position)
  if (n_chr != nLoc(x) || n_pos != nLoc(x)) {
    stop(error(
      "x needs a chromosome name and a position for every locus; found",
      n_chr, "chromosome entries and", n_pos, "position entries for",
      nLoc(x), "loci. Fill the slots from the locus metrics first, e.g.",
      "x$chromosome <- x$other$loc.metrics$<chromosome column>;",
      "x$position <- x$other$loc.metrics$<position column>\n"
    ))
  }
  
  # DO THE JOB
  
  # Extract valid chromosome / position pairs
  df_info <- data.frame(
    chr = as.character(x@chromosome),
    pos = x@position,
    stringsAsFactors = FALSE
  )
  
  # keep only rows where chr is non-NA/non-empty and pos > 0
  df_info <- df_info[
    !is.na(df_info$chr) &
      df_info$chr != ""    &
      !is.na(df_info$pos)  &
      df_info$pos > 0,
  ]
  
  if (nrow(df_info) == 0) {
    stop(error("No valid chromosome/position data found in x.\n"))
  }
  
  if (verbose >= 2) {
    cat(report("  Retained", nrow(df_info), "SNPs after initial filtering\n"))
  }
  
  # Summarise chromosomes & apply filters 
  chr_all <- df_info |>
    dplyr::group_by(chr) |>
    dplyr::summarise(chr_size = max(pos, na.rm = TRUE),
                     n_snps   = dplyr::n(),
                     .groups  = "drop")
  
  chr_stats <- chr_all |>
    dplyr::filter(n_snps   >= min.snps,
                  chr_size >= min.length)
  
  if (verbose >= 2) {
    n_few <- sum(chr_all$n_snps < min.snps)
    n_short <- sum(chr_all$n_snps >= min.snps & chr_all$chr_size < min.length)
    cat(report(
      "  Chromosomes with positioned SNPs:", nrow(chr_all), "; dropped",
      n_few, "with fewer than", min.snps, "SNPs and", n_short,
      "with the last SNP below",
      format(min.length, scientific = FALSE, big.mark = ","),
      "bp; retained", nrow(chr_stats), "\n"
    ))
  }
  
  if (nrow(chr_stats) == 0) {
    stop(error("No chromosomes meet the min.snps / min.length criteria.\n"))
  }
  
  # Order chromosomes alphabetically, numeric-aware; the first name is
  # drawn at the top, so the factor levels run in reverse
  chr_order <- rev(stringr::str_sort(chr_stats$chr, numeric = TRUE))
  chr_stats <- chr_stats[match(chr_order, chr_stats$chr), ]
  
  # Build y axis labels
  chr_stats <- chr_stats |>
    dplyr::mutate(
      chr_label = if (chr.info) {
        sprintf("%s (%.0f SNPs, %.1f Mb)", chr, n_snps, chr_size / 1e6)
      } else {
        chr
      }
    )
  
  # Bin SNPs; every bin up to the last SNP of each chromosome is kept, so
  # empty bins are drawn in the lowest colour
  counts <- df_info |>
    dplyr::filter(chr %in% chr_stats$chr) |>
    dplyr::mutate(bin_center = floor(pos / bin.size) * bin.size + bin.size / 2) |>
    dplyr::count(chr, bin_center, name = "n_snps")
  
  grid <- do.call(rbind, lapply(seq_len(nrow(chr_stats)), function(i) {
    last_center <- floor(chr_stats$chr_size[i] / bin.size) * bin.size + bin.size / 2
    data.frame(chr = chr_stats$chr[i],
               bin_center = seq(bin.size / 2, last_center, by = bin.size),
               stringsAsFactors = FALSE)
  }))
  
  plot_dat <- dplyr::left_join(grid, counts, by = c("chr", "bin_center"))
  plot_dat$n_snps[is.na(plot_dat$n_snps)] <- 0L
  plot_dat$chr_label <- factor(chr_stats$chr_label[match(plot_dat$chr, chr_stats$chr)],
                               levels = chr_stats$chr_label)
  plot_dat <- plot_dat[, c("chr_label", "bin_center", "n_snps")]
  
  if (verbose >= 3) {
    n_bins <- as.vector(table(factor(plot_dat$chr_label, levels = chr_stats$chr_label)))
    summary_tab <- data.frame(
      chromosome = chr_stats$chr,
      n_snps = chr_stats$n_snps,
      last_snp_Mb = round(chr_stats$chr_size / 1e6, 1),
      n_bins = n_bins,
      stringsAsFactors = FALSE
    )
    cat(report("  Chromosomes plotted (top to bottom):\n"))
    print(summary_tab[rev(seq_len(nrow(summary_tab))), ], row.names = FALSE)
  }
  
  xmax <- max(plot_dat$bin_center) + bin.size / 2
  
  # Draw heat-map
  p1 <- ggplot2::ggplot(plot_dat,
                        ggplot2::aes(x = bin_center,
                                     y = chr_label,
                                     fill = n_snps)) +
    ggplot2::geom_tile(width = bin.size, height = 0.9) +
    ggplot2::scale_x_continuous(
      limits = c(0, xmax),
      expand = ggplot2::expansion(mult = c(0, 0))
    ) +
    ggplot2::scale_fill_gradientn(
      colours = color.palette(255),
      name = paste0("SNPs per ",
                    format(bin.size, big.mark = ","),
                    " bp")
    ) +
    plot.theme +
    ggplot2::labs(
      x = "Genomic position (bp)",
      y = "Chromosome",
      title = plot.title
    ) +
    ggplot2::theme(
      panel.grid   = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_line()
    )
  
  # PRINTING OUTPUTS
  if (plot.display) {
    print(p1)
  }
  
  # Optionally save the plot
  if (!is.null(plot.file)) {
    tmp <- utils.plot.save(p1,
                           dir = plot.dir,
                           file = plot.file,
                           verbose = verbose)
  }
  
  # FLAG SCRIPT END
  if (verbose >= 1) {
    cat(report("Completed:", funname, "\n"))
  }
  
  # RETURN 
  invisible(p1)
}
