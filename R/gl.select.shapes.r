#' @name gl.select.shapes
# Preliminaries -- specify parameter definitions -------------------
#' @title Selects shapes from the base R shape palette and outputs as a vector
#' @family graphics
#'
#' @description
#' This script draws upon the standard R shape palette to extract a vector of
#' shapes for plotting, where the script that follows has a shape parameter
#' expecting a vector of shapes.

#' @param x Optionally, provide a gl object from which to determine the number
#' of populations [default NULL].
#' @param select Select the shapes to retain in the output vector
#' [default NULL, all shapes shown and returned].
#' @param plot.display If TRUE, the shape palette is displayed in the graphics
#' window [default TRUE; display is suppressed when verbose is 0].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' progress log; 3, progress and results summary; 5, full report
#' [default 2 or as specified using gl.set.verbosity].
#'
#' @details
#' By default the shape palette will be displayed in full in the graphics window
#' from which shapes can be selected in a subsequent run, and the vector of
#' shapes returned for later use.

#' The select parameter can be used to select shapes from the specified 26
#' shapes available (0-25). For example, select=c(1,1,3) will select shape 1, 1
#' again and 3 to retain in the final vector. This can be useful for fine-tuning
#' shape selection, and matching colors and shapes.

#' If a gl object is provided, the number of specified shapes must agree with
#' the number of populations; if select is not specified, one shape is assigned
#' per population.

#' The palette is displayed only when plot.display=TRUE (the default) and
#' verbose >= 1.

#' @author Custodian: Arthur Georges -- Post to
#' \url{https://groups.google.com/d/forum/dartr}

#' @examples
#' # SET UP DATASET
#' if (isTRUE(getOption("dartR_fbm"))) testset.gl <- gl.gen2fbm(testset.gl)
#' gl <- testset.gl
#' levels(pop(gl))<-c(rep('Coast',5),rep('Cooper',3),rep('Coast',5),
#' rep('MDB',8),rep('Coast',7),'Em.subglobosa','Em.victoriae')
#' # EXAMPLES
#' shapes <- gl.select.shapes() # Select and display available shapes
#' # Select and display a restricted set of shapes
#' shapes <- gl.select.shapes(select=c(1,1,1,5,8))
#'  # Select set of shapes and check with no. of pops.
#' shapes <- gl.select.shapes(x=gl,select=c(1,1,1,5,8))

#' @seealso \code{\link{gl.select.colors}}

#' @return A vector with the required number of shapes
#' @export
#'
# Function --------------------
gl.select.shapes <- function(x = NULL,
                             select = NULL,
                             plot.display = TRUE,
                             verbose = NULL) {
  # Preliminaries -------------------------------
    # SET VERBOSITY
    verbose <- gl.check.verbosity(verbose)

    if (verbose == 0) {
        plot.display <- FALSE
    }

    # FLAG SCRIPT START
    funname <- match.call()[[1]]
    utils.flag.start(func = funname,
                     build = "Jody",
                     verbose = verbose)

    # SCRIPT SPECIFIC ERROR CHECKING

    if (!is.null(select)) {
        if (min(select) < 0 | max(select) > 25) {
            stop(error(
                "Fatal Error: specified shapes must be in the range 0-25\n"
            ))
        }
    }

    if (!is.null(x)) {
        datatype <- utils.check.datatype(x, verbose = verbose)
        if (!is.null(select)) {
            if (length(select) != nPop(x)) {
                stop(error(
                    "Fatal Error: the number of specified shapes (",
                    length(select),
                    ") must agree with the number of populations (",
                    nPop(x),
                    ") in the gl object\n"
                ))
            }
        } else {
            if (nPop(x) > 26) {
                stop(error(
                    "Fatal Error: only 26 shapes (0-25) are available but the",
                    " gl object has ",
                    nPop(x),
                    " populations; specify select manually (repeats allowed)\n"
                ))
            }
            if (verbose >= 2) {
                cat(report(
                    "  Setting the number of shapes to the number of populations\n"
                ))
            }
            select <- 0:(nPop(x) - 1)
        }
    } else if (is.null(select)) {
        if (verbose >= 2) {
            cat(
                warn(
                    "  Warning: Required shapes not specified, displaying and returning all available 26 shapes\n"
                )
            )
        }
        select <- 0:25
    }

    # DO THE JOB -----------------

    if (plot.display) {
        y.coord <-
            rev(c(rep(1, 6), rep(2, 5), rep(3, 5), rep(4, 5), rep(5, 5)))
        y.coord <- y.coord[1:length(select)]
        x.coord <- c(rep(1:5, 5), 6)
        x.coord <- x.coord[1:length(select)]
        plot(
            x.coord,
            y.coord,
            pch = select,
            cex = 1.5,
            ylim = c(1, 5.5),
            xlim = c(1, 6.5),
            axes = FALSE,
            xlab = "",
            ylab = "",
            bg = "blue"
        )
        text(x.coord, y.coord, labels = select, pos = 3)
    }

    if (verbose >= 2) {
        cat(report(
            "  Displaying and returning shapes",
            paste(select, collapse = ", "),
            "\n"
        ))
    }

    # FLAG SCRIPT END ----------------------

    if (verbose >= 1) {
        cat(report("Completed:", funname, "\n"))
    }
    # End block -------------------

    return(select)
}
