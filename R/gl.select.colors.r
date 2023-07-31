#' @name gl.select.colors
# Preliminaries -- set parameter definitions --------------------
#' @title Selects colors from one of several palettes and outputs as a vector
#' @family graphics

#' @description
#' This function draws upon a number of specified color libraries to extract a
#' vector of colors for plotting. For use where the function that follows has a color
#' parameter expecting a vector of colors.

#' @param x Optionally, provide a gl object from which to determine the number
#' of populations [default NULL].
#' @param library Name of the color library to be used, one of 'brewer'
#' 'gr.palette', 'r.hcl' or 'baseR' [default scales::hue_pl].
#' @param palette Name of the color palette to be pulled from the specified
#' library, refer function help [default is library specific].
#' @param ncolors number of colors to be displayed and returned [default 9 or nPop(gl)].
#' @param select select bu number the colors to retain in the output vector; 
#' can repeat colors. [default NULL].
#' @param plot.display if TRUE, plot the colours in the plot window [default=TRUE]
#' @param verbose -- verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' progress log; 3, progress and results summary; 5, full report
#' [default 2 or as specified using gl.set.verbosity].
#' 
#' @details
#' Colors are chosen by specifying a library (one of 'brewer'
#' 'gr.palette', 'r.hcl' or 'baseR') and a palette within that library. Each
#' library has its own array of palettes, which can be listed as outlined below.
#' Alternatively, if you specify an incorrect palette, the list of available
#' palettes for the specified library will be listed.
#' 
#' The available color libraries and their palettes include:
#' \itemize{
#' \item library 'brewer' and the palettes available can be listed by
#' RColorBrewer::display.brewer.all() and RColorBrewer::brewer.pal.info.
#' \item library 'gr.palette' and the palettes available can be listed by
#' grDevices::palette.pals()
#' \item library 'r.hcl' and the palettes available can be listed by
#' grDevices::hcl.pals()
#' \item library 'baseR' and the palettes available are: 'rainbow','heat',
#' 'topo.colors','terrain.colors','cm.colors'.
#' }
#' 
#' If the library is not specified, then the default library 'scales' is set and
#' the default palette of 'hue_pal is set.
#' 
#' If the library is set but the palette is not specified, all palettes for that
#' library will be listed and a default palette will then be chosen.

#' The color palette will be displayed in the graphics window for the requested
#' number of colors (or 9 if not specified or nPop(gl) if a genlight object is 
#' specified),and the vector of colors returned by assignment for later use.

#' The select parameter can be used to select colors from the specified ncolors.
#' For example, select=c(1,1,3) will select color 1, 1 again and 3 to retain in
#' the final vector. This can be useful for fine-tuning color selection, and
#' matching colors and shapes.

#' @author Custodian: Arthur Georges -- Post to
#' \url{https://groups.google.com/d/forum/dartr}
#' 
#' @examples
#' # SET UP DATASET
#' gl <- testset.gl
#' levels(pop(gl))<-c(rep('Coast',5),rep('Cooper',3),rep('Coast',5),
#' rep('MDB',8),rep('Coast',7),'Em.subglobosa','Em.victoriae')
#' # EXAMPLES -- SIMPLE
#' colors <- gl.select.colors()
#' colors <- gl.select.colors(library='brewer',palette='Spectral',ncolors=6)
#' colors <- gl.select.colors(library='baseR',palette='terrain.colors',ncolors=6)
#' colors <- gl.select.colors(library='baseR',palette='rainbow',ncolors=12)
#' colors <- gl.select.colors(library='gr.hcl',palette='RdBu',ncolors=12)
#' colors <- gl.select.colors(library='gr.palette',palette='Pastel 1',ncolors=6)
#' # EXAMPLES -- SELECTING colorS
#' colors <- gl.select.colors(library='baseR',palette='rainbow',ncolors=12,select=c(1,1,1,5,8))
#' # EXAMPLES -- CROSS-CHECKING WITH A GENLIGHT OBJECT
#' colors <- gl.select.colors(x=gl,library='baseR',palette='rainbow',ncolors=12,select=c(1,1,1,5,8))
#' 
#' @seealso \code{\link{gl.select.shapes}}
#' 
#' @importFrom grDevices cm.colors hcl.pals palette.pals terrain.colors topo.colors rainbow
#' @export
#' 
#' @return A vector with the required number of colors

# Testing scripts
# cols <- gl.select.colors()
# cols <- gl.select.colors(x=testset.gl)
# cols <- gl.select.colors(x=testset.gl,ncolors=3)
# cols <- gl.select.colors(library="brewer")
# cols <- gl.select.colors(library="baseR")
# cols <- gl.select.colors(library="gr.hcl")
# cols <- gl.select.colors(library="gr.palette")
# cols <- gl.select.colors(library='brewer',palette='Spectral',ncolors=6)
# cols <- gl.select.colors(library='baseR',palette='Spectral',ncolors=6)
# cols <- gl.select.colors(library='gr.palette',palette='Spectral',ncolors=6)
# cols <- gl.select.colors(palette='Spectral',ncolors=6)

# Function -----------------
gl.select.colors <- function(x = NULL,
                             library = NULL,
                             palette = NULL,
                             ncolors = NULL,
                             select = NULL,
                             plot.display=TRUE,
                             verbose = NULL) {
# Preliminaries -----------------
    # SET VERBOSITY
    verbose <- gl.check.verbosity(verbose)
    
    # FLAG SCRIPT START
    funname <- match.call()[[1]]
    utils.flag.start(func = funname,
                     build = "v.2023.3",
                     verbose = verbose)
    
    # CHECK PACKAGES
    pkg <- "RColorBrewer"
    if (!(requireNamespace(pkg, quietly = TRUE))) {
      cat(error(
        "Package",
        pkg,
        " needed for this function to work. Please install it.\n"
      ))
      return(-1)
    }
    
    pkg <- "scales"
    if (!(requireNamespace(pkg, quietly = TRUE))) {
      cat(error(
        "Package",
        pkg,
        " needed for this function to work. Please install it.\n"
      ))
      return(-1)
    }
    
    # SCRIPT SPECIFIC ERROR CHECKING
    
    if (!is.null(x)) {
      datatype <- utils.check.datatype(x)
    }
    
    if (is.null(ncolors)) {
      if (!is.null(x)) {
        ncolors <- nPop(x)
        if (verbose >= 2) {
          cat(
            warn(
              "  Warning: Number of required colors not specified, set to number of pops",
              nPop(x),"in gl object\n"
            )
          )
        }
      } else {
        ncolors <- 9  
        if (verbose >= 2) {
          cat(
            warn(
              "  Warning: Number of required colors not specified, set to 9\n"
            )
          )
        }
        
      }
    }
    
    if (!is.null(select)) {
      if (!is.null(x)) {
        if (nPop(x) != length(select)) {
          stop(
            error(
              "Fatal Error: Number of specified colors",
              length(select),
              "does not correspond to number of populations",
              nPop(x),
              "in supplied genlight object\n"
            )
          )
        } else {
          if (verbose >= 2) {
            cat(
              report(
                "  Number of specified colors",
                length(select),
                "corresponds to number of populations in supplied genlight object\n"
              )
            )
          }
        }
      }
    }
    
    if (is.null(library)) {
        library <- "scales"
        palette <- "hue_pal"
        colors <- (scales::hue_pal())(ncolors)
        if (verbose >= 2) {
            cat(warn(
                "  Warning: No color library specified, set to",library,"and palette set to",palette,"\n"
            ))
            cat(warn("    Need to select one of baseR, brewer, scales, gr.palette or gr.hcl\n"))
        }
    } else {
        if (library == "brewer") {
            if (is.null(palette)) {
              palette <- "Spectral"
                if (verbose >= 2) {
                    cat(warn(
                        "  Warning: Palette not specified, set to",palette,"\n"
                    ))
                }

            }
            if (!(palette %in% row.names(RColorBrewer::brewer.pal.info))) {
                 palette <- "Spectral"
                 if (verbose >= 2) {
                    cat(
                        warn(
                            "  Warning: Nominated palette not available in RColorBrewer, should be one of\n"
                        )
                    )
                    cat(warn(paste(
                        row.names(RColorBrewer::brewer.pal.info),
                        collapse = ", "
                    ), "\n"))
                    cat(warn("  Set to Spectral\n"))
                }
            }
            colors <- RColorBrewer::brewer.pal(ncolors, palette)
            if (verbose >= 2) {
            cat(report(("  Library: RColorBrewer\n")))
            cat(report(("  Palette: brewer.pal\n")))
            }
            
        } else if (library == "gr.palette") {
          if (is.null(palette)) {
            if (verbose >= 2) {
              cat(warn(
                "  Warning: Palette not specified, set to Tableau 10\n"
              ))
            }
            palette <- "Tableau 10"
          }
          if (!(palette %in% grDevices::palette.pals())) {
            if (verbose >= 2) {
              cat(
                warn(
                  "  Warning: Nominated palette not available in grDevices::palette, should be one of\n"
                )
              )
              cat(warn(paste(
                palette.pals(), collapse = ", "
              ), "\n"))
              cat(warn("  Set to Tableau 10\n"))
            }
            palette <- "Tableau 10"
          }
          colors <-
            grDevices::palette.colors(n = ncolors, palette = palette)
          if (verbose >= 2) {
            cat(report(("  Library: grDevices\n")))
            cat(report(("  Palette: palette.pals\n")))
          }
        } else if (library == "gr.hcl") {
            if (is.null(palette)) {
                palette <- "Spectral"
                if (verbose >= 2) {
                    cat(warn(
                        "  Warning: Palette not specified, set to Spectral\n"
                    ))
                }
            }
            if (!(palette %in% grDevices::hcl.pals())) {
                if (verbose >= 2) {
                    cat(
                        warn(
                            "  Warning: Nominated palette not available in grDevices::hcl, should be one of\n"
                        )
                    )
                    cat(warn(paste(
                        hcl.pals(), collapse = ", "
                    ), "\n"))
                    cat(warn("  Set to Spectral\n"))
                }
                palette <- "Spectral"
            }
            colors <-
                grDevices::hcl.colors(n = ncolors, palette = palette)
            if (verbose >= 2) {
              cat(report(("  Library: grDevices\n")))
              cat(report(("  Palette: hcl.pals\n")))
            }
            
        } else if (library == "baseR") {
            if (is.null(palette)) {
                if (verbose >= 2) {
                    cat(warn(
                        "  Warning: Palette not specified, set to rainbow\n"
                    ))
                }
                palette <- "rainbow"
            }
            if (!(
                palette %in% c(
                    "rainbow",
                    "heat",
                    "topo.colors",
                    "terrain.colors",
                    "cm.colors"
                )
            )) {
                if (verbose >= 2) {
                    cat(
                        warn(
                            "  Warning: Nominated palette not available in base R, should be one of\n"
                        )
                    )
                    cat(
                        warn(
                            "  rainbow, heat.colors, topo.colors, terrain.colors, cm.colors\n"
                        )
                    )
                    cat(warn("  Set to rainbow\n"))
                }
                palette <- "rainbow"
            }
            if (palette == "rainbow") {
                colors <- rainbow(n = ncolors)
                if (verbose >= 2) {
                cat(report(("  Library: baseR\n")))
                cat(report(("  Palette: rainbow\n")))
                }
            } else if (palette == "topo.colors") {
                colors <- topo.colors(n = ncolors)
                if (verbose >= 2) {
                cat(report(("  Library: baseR\n")))
                cat(report((
                    "  Palette: topo.colors\n"
                )))
                }
            } else if (palette == "terrain.colors") {
                colors <- terrain.colors(n = ncolors)
                if (verbose >= 2) {
                cat(report(("  Library: baseR\n")))
                cat(report((
                    "  Palette: terrain.colors\n"
                )))
                }
            } else if (palette == "cm.colors") {
                colors <- cm.colors(n = ncolors)
                if (verbose >= 2) {
                cat(report(("  Library: baseR\n")))
                cat(report(("  Palette: cm.colors\n")))
                }
            } else {
                if (verbose >= 2) {
                    cat(
                        warn(
                            "  Warning: nominated palette not in Base R, selecting rainbow\n"
                        )
                    )
                }
                colors <- rainbow(n = ncolors)
            }
        }
    }

        if (library == "brewer") {
            library <- "RColorBrewer"
        }
        if (library == "gr.hcl") {
            library <- "grDevice-hcl"
        }
        if (library == "gr.palette") {
            library <- "grDevice-palette"
        }
        if (!is.null(select)) {
            colors <- colors[c(select)]
            if (verbose >= 2) {
              cat(
                "  Showing and returning",
                length(select),
                "of",
                ncolors,
                "colors for library",
                library,
                ": palette",
                palette,
                "\n"
            )
            }
        } else {
          if (verbose >= 2) {
            cat(
                "  Showing and returning",
                ncolors,
                "colors for library",
                library,
                ": palette",
                palette,
                "\n"
            )
          }
        }
        
        if(verbose >=1 && plot.display==TRUE){scales::show_col(colors)}

    # FLAG SCRIPT END ----------------
    
    if (verbose >= 1) {
        cat(report("Completed:", funname, "\n"))
    }
    # End block------------------
    
    return(colors)
}
