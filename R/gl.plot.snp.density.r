#' @name gl.plot.snp.density
#'
#' @title Plot SNP density along chromosomes (heat-map)
#'
#' @description
#' Generates a tiled heat-map of single-nucleotide polymorphisms (SNPs)
#' across chromosomes or scaffolds in a genlight object. SNPs are binned
#' into fixed-width windows and coloured by SNP count, optionally
#' annotating chromosomes with their length (Mb) and SNP number.
#'
#' @param x            A genlight object with chromosome information in
#'                     @chromosome and SNP positions in @position
#'                     [required].
#' @param bin.size     Width (bp) of the genomic bins used to count SNPs
#'                     [default 1 000 000].
#' @param min.snps     Minimum number of SNPs a chromosome must possess to
#'                     be plotted [default 50].
#' @param min.length   Minimum chromosome length (bp) to include
#'                     [default 1 000 000].
#' @param color.palette A function returning a vector of colours to be
#'                     passed to ggplot2; typically viridis::viridis
#'                     [default viridis::viridis].
#' @param chr.info     Logical; if TRUE append (N SNPs, L Mb) to
#'                     chromosome labels [default TRUE].
#' @param plot.title   Optional main title for the plot [default NULL].
#' @param plot.theme   ggplot2 theme applied to the plot
#'                     [default theme_dartR()].
#' @param save2tmp     Logical; save the ggplot object to tempdir()
#'                     for later retrieval with gl.print.reports()
#'                     [default FALSE].
#' @param verbose      Verbosity: 0 = silent; 1 = begin/end; 2 = progress;
#'                     3 = progress + summary; 5 = full report
#'                     [default 2 or as set by gl.set.verbosity()].
#'
#' @details
#' Chromosomes are ordered from longest (bottom) to shortest (top) so that
#' density patterns can be compared visually.  Bins containing no SNPs are
#' rendered in the lowest colour of the palette.  The function does not
#' modify the input genlight object.
#'
#' @return A ggplot object (invisibly) displaying the SNP-density
#'         heat-map.
#'
#' @author Custodian: Luis Mijangos -- Post to
#' \url{https://groups.google.com/d/forum/dartr}
#' 
#' @export
#' @examples
#' 
#'   t1 <- platypus.gl
#'   t1$chromosome <- t1$other$loc.metrics$Chrom_Platypus_Chrom_NCBIv1
#'   t1$position   <- t1$other$loc.metrics$ChromPos_Platypus_Chrom_NCBIv1
#'   gl.plot.snp.density(t1,
#'                       bin.size   = 5e6,
#'                       min.snps   = 10,
#'                       min.length = 2e6,
#'                       plot.title = "Platypus SNP density")


gl.plot.snp.density <- function(x,
                                bin.size      = 1e6,
                                min.snps      = 50,
                                min.length    = 1e6,
                                color.palette = viridis::viridis,
                                chr.info      = TRUE,
                                plot.title    = NULL,
                                plot.theme    = theme_dartR(),
                                save2tmp      = FALSE,
                                verbose       = NULL) {
  
  pos <- n_snps <- chr_size <- bin_start <- chr_label <- bin_center <- NULL
  
  # SET VERBOSITY
  verbose <- gl.check.verbosity(verbose)
  
  # FLAG SCRIPT START
  funname <- match.call()[[1]]
  utils.flag.start(func = funname, verbose = verbose)
  
  # CHECK DATATYPE
  utils.check.datatype(x, accept = "SNP", verbose = verbose)
  
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
    stop(error("Parameter min.snps must be > 1.\n"))
  }
  if (min.length < 1) {
    stop(error("Parameter min.length must be > 1 bp.\n"))
  }
  
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
  chr_stats <- df_info |>
    dplyr::group_by(chr) |>
    dplyr::summarise(chr_size = max(pos, na.rm = TRUE),
                     n_snps   = dplyr::n(),
                     .groups  = "drop") |>
    dplyr::filter(n_snps   >= min.snps,
                  chr_size >= min.length) |>
    dplyr::arrange(dplyr::desc(chr_size))
  
  if (nrow(chr_stats) == 0) {
    stop(error("No chromosomes meet the min.snps / min.length criteria.\n"))
  }
  
  # Build y axis labels & factor ordering
  chr_stats <- chr_stats |>
    dplyr::mutate(
      chr_label = if (chr.info) {
        sprintf("%s (%.0f SNPs, %.1f Mb)", chr, n_snps, chr_size / 1e6)
      } else {
        chr
      }
    )

  # Bin SNPs
  plot_dat <- df_info |>
    dplyr::inner_join(chr_stats, by = "chr") |>
    dplyr::mutate(
      bin_start  = floor(pos / bin.size) * bin.size,
      bin_center = bin_start + bin.size / 2
    ) |>
    as.data.frame() |>
    plyr::count(vars = c("chr_label", "bin_center")) |>
    plyr::rename(c("freq" = "n_snps"))
  
  plot_dat$chr_label <- factor(plot_dat$chr_label,
                               levels = sort(unique(plot_dat$chr_label),
                                             decreasing = T))
  
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
  
  # DISPLAY OUTPUT
  print(p1)
  
  # SAVE TO TEMPDIR (optional)
  if (save2tmp) {
    tmp_plot <- tempfile(pattern = paste0("dartR_plot_", funname, "_"),
                         fileext = ".rds")
    saveRDS(p1, file = tmp_plot)
    if (verbose >= 2) {
      cat(report("  Saved ggplot object to", tmp_plot, "\n"))
      cat(report("  Retrieve with gl.print.reports()\n"))
    }
  }
  
  # FLAG SCRIPT END
  if (verbose >= 1) {
    cat(report("Completed:", funname, "\n"))
  }
  
  # RETURN 
  invisible(p1)
}
