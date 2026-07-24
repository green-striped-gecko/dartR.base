#' @name gl.test.heterozygosity
#' @title Tests the difference in heterozygosity between populations taken
#'  pairwise
#' @description
#' Tests whether the unbiased expected heterozygosity (uHe) differs between
#' populations, for every pair of populations in a genlight object. The uHe of
#' each population is bootstrapped (resampling individuals or loci, see Details),
#' and for each pair the function returns the observed difference (pop1 - pop2),
#' a bootstrap confidence interval, a significance label and a two-sided p value
#' (raw and Benjamini-Hochberg adjusted across pairs). uHe is corrected for
#' sample size following equation 2 of Nei (1978).
#' 
#' @param x A genlight object containing the SNP genotypes [required].
#' @param nreps Number of bootstrap replicates used to build the sampling
#' distribution of each pairwise difference [default 1,000].
#' @param boot.method Resample across individuals ("ind") or across loci ("loc")
#' when bootstrapping the heterozygosity of each population [default "ind"].
#' @param alpha1 First significance level for comparison with diff=0 on plot
#' [default 0.05].
#' @param alpha2 Second significance level for comparison with diff=0 on plot
#' [default 0.01].
#' @param conf Confidence level for the bootstrap confidence interval of each
#' pairwise difference [default 0.95].
#' @param plot.out If TRUE, plots a sampling distribution of the differences for
#' each comparison [default TRUE].
#' @param max_plots Maximum number of plots to print per page [default 6].
#' @param plot.theme Theme for the plot. See Details for options
#'  [default theme_dartR()].
#' @param plot.colors List of two color names for the borders and fill of the
#'  plots [default gl.colors(2)].
#' @param plot.file Name for the RDS binary file to save (base name only, exclude extension) [default NULL]
#' @param plot.dir Directory to save the plot RDS files [default as specified 
#' by the global working directory or tempdir()]
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' progress log; 3, progress and results summary; 5, full report
#' [default NULL, unless specified using gl.set.verbosity].
#' @details
#' \strong{Choosing the resampling unit (boot.method)}
#'
#' The uHe of each population is bootstrapped with the same engine as
#' \code{\link{gl.report.heterozygosity}}, and \code{boot.method} sets what is
#' resampled. The two units answer different questions and can yield different
#' confidence intervals and p values for the same data:
#' \itemize{
#'  \item \code{"ind"} (default) resamples individuals. It treats the sampled
#'  individuals as a sample of each population and tests whether the two
#'  populations differ in heterozygosity. This is the usual question.
#'  \item \code{"loc"} resamples loci, with the individuals held fixed. It tests
#'  whether the genome-wide difference is consistent across loci rather than
#'  driven by a few of them, and does not generalise to the populations the
#'  individuals came from. Use it only when the individuals are effectively a
#'  census of each group (for example a fully sampled captive lineage or a
#'  complete cohort), where there is no individual-sampling variance to account
#'  for; with a sample of individuals from each population, use \code{"ind"}.
#'  }
#' For \code{"loc"} the loci are resampled independently within each population
#' and treated as independent, so the interval ignores linkage and is
#' approximate.
#'
#' \strong{Interpreting the plot}
#'
#' With \code{plot.out = TRUE} one histogram is drawn per pair of populations.
#' The histogram is the bootstrap distribution of the difference in uHe
#' (pop1 - pop2), with four vertical lines:
#' \itemize{
#'  \item green: the observed difference.
#'  \item blue: zero (no difference).
#'  \item two dark red lines: the significance limits at level \code{alpha1}
#'  (default 0.05).
#'  \item two bright red lines: the limits at level \code{alpha2} (default 0.01).
#'  }
#' Read significance from the blue zero line, not the green one. If zero lies
#' outside the red lines of a given level, the difference is significant at that
#' level; if it lies inside the \code{alpha1} lines, the difference is not
#' significant. The green observed line sits near the centre of its own bootstrap
#' distribution by construction, so its position is not the test. Each panel
#' subtitle gives the significance label and the p value.
#'
#' \strong{Saving the plot}
#'
#' If \code{plot.file} is given, the ggplot is saved as an RDS file with
#' \code{saveRDS()} and can be reloaded with \code{readRDS()}. It is written to
#' \code{plot.dir} if given, otherwise to \code{tempdir()}. Other ggplot themes
#' can be consulted at
#' \url{https://ggplot2.tidyverse.org/reference/ggtheme.html} and
#' \url{https://yutannihilation.github.io/allYourFigureAreBelongToUs/ggthemes/}.
#' @return A dataframe with one row per population pair containing the two
#'  population labels, the observed difference in unbiased expected
#'  heterozygosity (pop1 - pop2), the lower and upper bounds of its bootstrap
#'  confidence interval, the significance label and the two-sided bootstrap
#'  p value (raw and Benjamini-Hochberg adjusted).
#' @author Custodian: Luis Mijangos (Post to
#'  \url{https://groups.google.com/d/forum/dartr})
#' @references 
#' Nei, M. (1978). Estimation of average heterozygosity and genetic distance 
#' from a small number of individuals. Genetics, 89(3), 583-590.
#' 
#' @examples
#' if (isTRUE(getOption("dartR_fbm"))) platypus.gl <- gl.gen2fbm(platypus.gl)
#' out <- gl.test.heterozygosity(platypus.gl, nreps=100, verbose=3, plot.out=TRUE)
#' @family Genetic variation within populations
#' @import patchwork
#' @export

gl.test.heterozygosity <- function(x,
                                   nreps = 1000,
                                   boot.method = "ind",
                                   alpha1 = 0.05,
                                   alpha2 = 0.01,
                                   conf = 0.95,
                                   plot.out = TRUE,
                                   max_plots = 6,
                                   plot.theme = theme_dartR(),
                                   plot.colors = gl.select.colors(ncolors=2, verbose=0),
                                   plot.file=NULL,
                                   plot.dir=NULL,
                                   verbose = NULL) {
    # TRAP COMMAND
    
    funname <- match.call()[[1]]
    
    # SET VERBOSITY
    
    verbose <- gl.check.verbosity(verbose)
    
    # SET WORKING DIRECTORY
    plot.dir <- gl.check.wd(plot.dir,verbose=0)
    
    # CHECKS DATATYPE
    
    datatype <- utils.check.datatype(x, verbose = verbose)
    
    # SCRIPT SPECIFIC ERROR CHECKING
    
    # Upper and lower significance boundaries
    
    if (alpha1 < 0 || alpha1 > 1) {
        cat(warn(
            "Warning: First alpha value should be between 0 and 1, set to 0.05\n"
        ))
        alpha1 <- 0.05
    }
    
    if (alpha2 < 0 || alpha2 > 1) {
        cat(warn(
            "Warning: Second alpha value should be between 0 and 1, set to 0.01\n"
        ))
        alpha2 <- 0.01
    }

    # Bootstrap resampling unit
    if (!(boot.method == "ind" || boot.method == "loc")) {
        cat(warn(
            "  Warning: boot.method must be 'ind' or 'loc', set to 'ind'\n"
        ))
        boot.method <- "ind"
    }
    
    upper1 <- 1 - alpha1  # significance level #1
    upper2 <- 1 - alpha2  # significance level #2
    if (upper1 > upper2) {
        tmp <- upper2
        upper2 <- upper1
        upper1 <- tmp
    }
    lower1 <- 1 - upper1
    lower2 <- 1 - upper2
    
    # Set a population if none is specified (such as if the genlight object has been generated manually)
    if (is.null(pop(x)) |
        is.na(length(pop(x))) | length(pop(x)) <= 0) {
        if (verbose >= 2) {
            cat(
                warn(
                    "  No population assignments detected,
                             individuals assigned to a single population labelled 'pop1'\n"
                )
            )
        }
        pop(x) <- array("pop1", dim = nInd(x))
        pop(x) <- as.factor(pop(x))
    }
    
    # Check for monomorphic loci
    
    if (x@other$loc.metrics.flags$monomorphs == FALSE) {
        if (verbose >= 1) {
            cat(
                warn(
                    "  Warning: genlight object contains monomorphic loci which will be factored into heterozygosity estimates\n"
                )
            )
        }
    }
    
    # FLAG SCRIPT START
    
    if (verbose >= 1) {
        if (verbose == 5) {
            cat(report("\n\nStarting", funname, "[ Build =", build, "]\n\n"))
        } else {
            cat(report("\n\nStarting", funname, "\n\n"))
        }
    }
    
    # DO THE JOB
    
    # Bootstrap the unbiased expected heterozygosity (uHe) of each population
    # with the shared engine (pop.het + boot::boot, as used by
    # gl.report.heterozygosity), so this function resamples the same way. The
    # unit of resampling is set by boot.method: individuals ("ind", the default,
    # which captures the finite- and unequal-sample-size uncertainty) or loci
    # ("loc").
    if (verbose >= 2) {
        cat(report(
            paste0(
                "  Bootstrapping the heterozygosity of each population ",
                "(boot.method = '", boot.method,
                "') to build the pairwise difference distributions\n"
            )
        ))
        cat(important("  Please be patient .... go have a coffee\n"))
    }

    sgl <- seppop(x)
    npop <- nPop(x)

    # he_boot: nreps bootstrap replicates of uHe for each population (rows = pops)
    he_boot <- matrix(NA_real_, nrow = npop, ncol = nreps)
    out <- numeric(npop)                 # observed uHe per population
    for (i in seq_len(npop)) {
        pop_mat <- as.matrix(sgl[[i]])
        if (boot.method == "loc") {
            pop_mat <- t(pop_mat)
        }
        res_boots <- boot::boot(
            data = pop_mat,
            statistic = pop.het,
            n.invariant = 0,
            aHet = FALSE,
            boot_method = boot.method,
            R = nreps
        )
        # pop.het returns c(Ho, He, uHe, FIS); uHe is column 3
        out[i] <- res_boots$t0[3]
        he_boot[i, ] <- res_boots$t[, 3]
    }
    names(out) <- popNames(x)

    # Observed pairwise differences (pop1 - pop2) and their bootstrap
    # distributions. mat[y, z, ] holds the resampled difference uHe(y) - uHe(z);
    # individuals (or loci) in different populations are independent samples, so
    # the per-population bootstraps are paired replicate-wise to form it.
    D <- outer(out, out, "-")
    mat <- array(data = NA_real_, dim = c(npop, npop, nreps))
    for (y in seq_len(npop)) {
        for (z in seq_len(npop)) {
            mat[y, z, ] <- he_boot[y, ] - he_boot[z, ]
        }
    }
    
    # Calculate the p values, significance for pairwise differences Dataframe to store output
    df <-
        data.frame(matrix(ncol = 7, nrow = (nPop(x) * nPop(x) - nPop(x)) / 2))
    colnames(df) <-
        c("pop1", "pop2", "diff", "CIlower", "CIupper", "significance", "pval")
    
    count <- 0
    
    if (plot.out) {
        if (verbose >= 2) {
            cat(report(
                "  Cycling through the",
                (nPop(x) * nPop(x) - nPop(x)) / 2,
                "pairs of populations\n"
            ))
        }
        if (verbose >= 2) {
            cat(report("  Plotting sampling distributions\n"))
        }
        
        # set the limit values for the x axis in the plots
        x_axis_limits_lots <- c(min(mat, na.rm = TRUE), max(mat, na.rm = TRUE))
        # count how many plots are going to be created
        total_number_plots <- choose(nPop(x), 2)
        # create list to contain plots
        p_list <- list()
    }

        for (y in 1:(nPop(x) - 1)) {
            for (z in (y + 1):nPop(x)) {
                count <- count + 1
                
                # Prepare the parameters for the plot
                
                # Observed He
                obs <- D[y, z]
                
                # Upper and lower significance limits
                u1quantile <-
                    as.numeric(quantile(mat[y, z,], upper1))
                u2quantile <-
                    as.numeric(quantile(mat[y, z,], upper2))
                l1quantile <-
                    as.numeric(quantile(mat[y, z,], lower1))
                l2quantile <-
                    as.numeric(quantile(mat[y, z,], lower2))
                
                # Is zero within the range specified by the upper and lower limits
                signif_res <- NA_character_
                if (0 < u2quantile && 0 > l2quantile) {
                    signif_res <- paste0("non-sig @", lower2)
                }
                if (0 < u1quantile && 0 > l1quantile) {
                    signif_res <- paste0("non-sig @", lower1)
                }
                if (0 > u1quantile || 0 < l1quantile) {
                    signif_res <- paste0("sig @", lower1)
                }
                if (0 > u2quantile || 0 < l2quantile) {
                    signif_res <- paste0("sig @", lower2)
                }
                
                # Two-sided bootstrap p value that the difference is zero, using
                # the (1 + b) / (nreps + 1) convention so p is never exactly 0.
                values <- mat[y, z,]
                nb <- length(values)
                p <- 2 * min((1 + sum(values <= 0)) / (nb + 1),
                             (1 + sum(values >= 0)) / (nb + 1))
                p <- signif(min(p, 1), digits = 6)

                # Percentile bootstrap confidence interval for the difference
                ci <- as.numeric(quantile(values,
                                          c((1 - conf) / 2, 1 - (1 - conf) / 2),
                                          na.rm = TRUE))

                # Store results in a df
                df$pop1[count] <- popNames(x)[y]
                df$pop2[count] <- popNames(x)[z]
                df$diff[count] <- obs
                df$CIlower[count] <- signif(ci[1], 6)
                df$CIupper[count] <- signif(ci[2], 6)
                df$significance[count] <- signif_res
                df$pval[count] <- p

                # Plot the histogram of pairwise differences (only when plot.out = TRUE)
                if (plot.out) {
                
                # Construct the label
                title <-
                    paste0(popNames(x)[y], " vs ", popNames(x)[z])
                subtitle <- paste0(signif_res, " (p=", p, ")")
                
                plot_values <- as.data.frame(values)
                colnames(plot_values) <- "values"
                
                # these plots do not contain legends
                if (count < total_number_plots |
                    count != max_plots) {
                    # Add lines for the observed value of He, Zero, upper and lower levels of significance
                    suppressWarnings(
                        p_temp <-
                            ggplot(plot_values, aes(x = values)) + geom_histogram(
                                bins = 50,
                                color = plot.colors[1],
                                fill = plot.colors[2]
                            ) +
                            geom_vline(
                                xintercept = u1quantile,
                                color = "firebrick4",
                                linewidth = 1
                            ) + geom_vline(
                                xintercept = u2quantile,
                                color = "firebrick1",
                                linewidth = 1
                            ) + geom_vline(
                                xintercept = l1quantile,
                                color = "firebrick4",
                                linewidth = 1
                            ) + geom_vline(
                                xintercept = l2quantile,
                                color = "firebrick1",
                                linewidth = 1
                            ) + geom_vline(
                                xintercept = D[y, z],
                                color = "green",
                                linewidth = 2
                            ) + geom_vline(
                                xintercept = 0,
                                color = "blue",
                                linewidth = 1
                            ) +
                            coord_cartesian(xlim = x_axis_limits_lots) + xlab("Difference") + ylab("Count") + plot.theme + theme(plot.title = element_text(size = 12)) +
                            labs(title = title, subtitle = subtitle)
                    )
                    
                    p_list[[count]] <- p_temp
                    
                }
                
                # this plot contains the plot legends. it is just created every max_plots and in the last plot
                if (count == total_number_plots |
                    (count %% max_plots) == 0) {
                    suppressWarnings(
                        p_temp <-
                            ggplot(plot_values, aes(x = values)) + geom_histogram(
                                bins = 50,
                                color = plot.colors[1],
                                fill = plot.colors[2]
                            ) +
                            coord_cartesian(xlim = x_axis_limits_lots) + xlab("Difference") + ylab("Count") + plot.theme + theme(plot.title = element_text(size = 12)) +
                            labs(title = title, subtitle = subtitle) + geom_vline(
                                aes(xintercept = u1quantile, color = "alpha1"),
                                linewidth = 1
                            ) + geom_vline(
                                aes(xintercept = u2quantile,
                                    color = "alpha2"),
                                linewidth = 1
                            ) + geom_vline(
                                xintercept = l1quantile,
                                color = "firebrick4",
                                linewidth = 1
                            ) + geom_vline(
                                xintercept = l2quantile,
                                color = "firebrick1",
                                linewidth = 1
                            ) + geom_vline(
                                aes(
                                    xintercept = D[y, z],
                                    color = "Observed"
                                ),
                                linewidth = 2
                            ) + geom_vline(
                                aes(xintercept = 0,
                                    color = "Zero_value"),
                                linewidth = 1
                            ) + scale_color_manual(
                                name = "Values",
                                breaks = c("alpha1", "alpha2", "Observed", "Zero_value"),
                                values = c(
                                    Zero_value = "blue",
                                    Observed = "green",
                                    alpha1 = "firebrick4",
                                    alpha2 = "firebrick1"
                                ),
                                labels = c(
                                    paste("Sig. ", alpha1),
                                    paste("Sig. ", alpha2),
                                    "Observed",
                                    "Zero value"
                                )
                            ) + guides(color = guide_legend(
                                override.aes = list(size = 5),
                                ncol = 2
                            )) + theme(
                                legend.position = "bottom",
                                legend.title = element_text(face = "bold"),
                                legend.text = element_text(size = 10)
                            )
                    )
                    
                    p_list[[count]] <- suppressWarnings(p_temp)
                }
                }

            }
            
        }
        
    if (plot.out) {
        # PRINTING OUTPUTS
        match_call <-
            paste0(names(match.call()),
                   "_",
                   as.character(match.call()),
                   collapse = "_")
        # using package patchwork
        seq_1 <- seq(1, length(p_list), max_plots)
        seq_2 <- seq(1, length(p_list), max_plots) - 1
        seq_2 <- seq_2[-1]
        seq_2 <- c(seq_2, length(p_list))
        for (i in 1:ceiling((length(p_list) / max_plots))) {
            p_final <-
                suppressWarnings(wrap_plots(p_list[seq_1[i]:seq_2[i]], ncol = 2))
            
            suppressWarnings(print(p_final))
            # SAVE INTERMEDIATES TO TEMPDIR
            temp_plot <- paste0(plot.file,"_", seq_1[i], "_to_", seq_2[i])
            if(!is.null(plot.file)){
              tmp <- utils.plot.save(p_final,
                                     dir=plot.dir,
                                     file=temp_plot,
                                     verbose=verbose)
            }    
        }
            
                if (verbose >= 2) {
                    cat(report("  Saving the ggplot to session tempfile\n"))
                }
        }
        
    
    
    # Multiple-testing correction across the pairwise comparisons
    df$pval.adj <- signif(p.adjust(df$pval, method = "BH"), 6)

    print(df)

    # Optionally save the plot ---------------------
    
    if(!is.null(plot.file)){
      tmp <- utils.plot.save(df,
                             dir=plot.dir,
                             file=paste0("table_",plot.file),
                             verbose=verbose)
    }
    
    # FLAG SCRIPT END
    
    if (verbose >= 1) {
        cat(report("Completed:", funname, "\n"))
    }
    
    # RETURN
    
    invisible(df)
}
