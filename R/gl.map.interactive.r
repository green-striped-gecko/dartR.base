#' @name gl.map.interactive
#' @title Creates an interactive map (based on latlon) from a genlight object
#' @family graphics
#' @description
#' Plots the individuals of a genlight object on an interactive leaflet map,
#' using the coordinates stored in \code{x@other$latlon}, with optional
#' population labels and optional links between populations or individuals
#' drawn from a distance matrix.
#' @param x A genlight object (including coordinates within the latlon slot) 
#' [required].
#' @param matrix A distance matrix between populations or individuals, either
#' a square matrix or a \code{dist} object such as returned by
#' \code{\link{gl.dist.pop}} or \code{\link{gl.dist.ind}}. The matrix is
#' visualised as links between individuals/populations. If the matrix has row
#' and column names matching the population or individual names they are used
#' to align it to the object; otherwise rows and columns are assumed to follow
#' the order of the populations/individuals in the object. If matrix is
#' asymmetric two lines with arrows are plotted [default NULL].
#' @param standard If a matrix is provided, values are standardised to be
#' between 1 and 10 if set to TRUE, otherwise taken as given. Symmetric links
#' use the value for line colour and width; asymmetric links use it for arrow
#' thickness [default TRUE].
#' @param symmetric If a symmetric matrix is provided only one line is drawn
#' based on the lower triangle of the matrix. If set to false arrows indicating
#' the direction are used instead [default TRUE].
#' @param pop.labels Population labels at the center of the individuals of
#'  populations [default TRUE].
#' @param pop.labels.cex Size of population labels [default 12].
#' @param ind.circles Should individuals be plotted as circles [default TRUE].
#' @param ind.circle.cols Colors of circles. A color palette or a vector with
#' as many colors as there are populations in the dataset [default rainbow].
#' @param ind.circle.cex Size of circles in pixels [default 10].
#' @param ind.circle.transparency Transparency of circles between 0=invisible 
#' and 1=no transparency [default 0.8].
#' @param palette.links Color palette for the symmetric links in case a matrix
#'  is provided. Ignored when \code{symmetric = FALSE}, where arrow colours
#'  encode the direction of the larger value [default NULL].
#' @param legend.title Legend's title for the symmetric links in case a matrix
#'  is provided [default NULL].
#' @param provider Passed to leaflet [default "Esri.NatGeoWorldMap"].
#' @param scale.bar Whether to add a scale bar [default TRUE].
#' @param raster.image Path to a georeferenced raster image to plot 
#' [default NULL].
#' @param raster.opacity The opacity of the raster, expressed from 0 to 1 
#' [default 0.5].
#' @param raster.colors The color palette to use to color the raster values
#'  [default scales::viridis_pal(option = "D")(255)].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' brief progress messages; 3, progress and results summary; 5, full report
#' [default 2, unless specified using gl.set.verbosity].
#' 
#' @details 
#' A wrapper around the \pkg{leaflet} package. For possible background 
#' maps check as specified via the provider:
#' \url{http://leaflet-extras.github.io/leaflet-providers/preview/index.html}
#' 
#' The palette.links argument can be any of the following:
#' A character vector of RGB or named colors. Examples: palette(), 
#' c("#000000", "#0000FF", "#FFFFFF"), topo.colors(10)
#' 
#' The name of an RColorBrewer palette, e.g. "BuPu" or "Greens".
#' 
#' The full name of a viridis palette: "viridis", "magma", "inferno", 
#' or "plasma".
#' 
#' A function that receives a single value between 0 and 1 and returns a color.
#'  Examples: colorRamp(c("#000000", "#FFFFFF"), interpolate = "spline").
#'
#' Symmetric links are drawn only for pairs whose value is above 0 and not
#' missing; a pair with a missing value in an asymmetric matrix is skipped
#' in that direction and drawn in grey in the other.
#' 
#' @author Author(s): Bernd Gruber. Custodian: Bernd Gruber -- Post to 
#' \url{https://groups.google.com/d/forum/dartr}
#' 
#' @examples
#' require("dartR.data")
#' if (isTRUE(getOption("dartR_fbm"))) bandicoot.gl <- gl.gen2fbm(bandicoot.gl)
#' gl.map.interactive(bandicoot.gl)
#' cols <- c("red","blue","yellow")
#' if (isTRUE(getOption("dartR_fbm"))) platypus.gl <- gl.gen2fbm(platypus.gl)
#' gl.map.interactive(platypus.gl, ind.circle.cols=cols, ind.circle.cex=10, 
#' ind.circle.transparency=0.5)
#' 
#' @importFrom methods is
#' @export
#' @return A leaflet map (htmlwidget), displayed when printed.

gl.map.interactive <- function(x,
                               matrix = NULL,
                               standard = TRUE,
                               symmetric = TRUE,
                               pop.labels = TRUE,
                               pop.labels.cex = 12,
                               ind.circles = TRUE,
                               ind.circle.cols = rainbow,
                               ind.circle.cex = 10,
                               ind.circle.transparency = 0.8,        
                               palette.links = NULL,
                               legend.title = NULL,
                               provider = "Esri.NatGeoWorldMap",
                               scale.bar = TRUE, 
                               raster.image = NULL,
                               raster.opacity = 0.5,
                               raster.colors = scales::viridis_pal(option = "D")(255),
                               verbose = NULL) {
    
    # SET VERBOSITY
    verbose <- gl.check.verbosity(verbose)
    
    # FLAG SCRIPT START
    funname <- match.call()[[1]]
    utils.flag.start(func = funname,
                     verbose = verbose)
    
    # CHECK DATATYPE
    datatype <- utils.check.datatype(x, verbose = verbose)
    
    # FUNCTION SPECIFIC ERROR CHECKING
    
    # CHECK IF PACKAGES ARE INSTALLED
    pkgs <- c("leaflet", "leaflet.minicharts")
    if (!is.null(raster.image)) {
        pkgs <- c(pkgs, "terra", "scales")
    }
    for (pkg in pkgs) {
        if (!(requireNamespace(pkg, quietly = TRUE))) {
            stop(error(
                "Package",
                pkg,
                " needed for this function to work. Please install it.\n"
            ))
        }
    }
    
    if (is.null(x@other$latlon)) {
        stop(error(
            "No valid coordinates are supplied at gl@other$latlon"
        ))
    }
    
    if (sum(colnames(x@other$latlon) %in% c("lat", "lon")) != 2) {
        stop(error(
            "Coordinates under gl@other$latlon are not named 'lat' and 'lon'."
        ))
    }
    
    if (!is.null(matrix)) {
        if (inherits(matrix, "dist")) {
            matrix <- as.matrix(matrix)
        }
        if (!is.matrix(matrix) || nrow(matrix) != ncol(matrix)) {
            stop(error(
                "The matrix argument must be a square matrix or a dist object."
            ))
        }
        if (nrow(matrix) != nInd(x) & nrow(matrix) != nPop(x)) {
            stop(
                error(
"The dimension of the provided matrix does neither match the number of 
individuals nor the number of populations."
                )
            )
        }
    }
    
    # if pop colors is a palette
    if (is(ind.circle.cols, "function")) {
        cols <- ind.circle.cols(length(levels(pop(x))))
    }
    # if pop colors is a vector
    if (!is(ind.circle.cols, "function")) {
        cols <- ind.circle.cols
    }
    if (length(cols) < nPop(x)) {
        stop(error(
            "ind.circle.cols has", length(cols), "colours but the dataset has",
            nPop(x), "populations."
        ))
    }
    ic <- cols[as.numeric(pop(x))]
    
    # DO THE JOB
    
    # coordinates may be stored as a matrix; work on a data frame
    df <- as.data.frame(x@other$latlon)
    # population centres, one row per population in popNames order
    centers <- data.frame(
        lon = as.vector(tapply(df$lon, pop(x), mean, na.rm = TRUE)),
        lat = as.vector(tapply(df$lat, pop(x), mean, na.rm = TRUE)),
        row.names = popNames(x)
    )
    
    m <- leaflet::leaflet() %>%
        leaflet::addProviderTiles(provider)
    
    if (ind.circles) {
        m <- m %>%
            leaflet::addCircles(
                lng = df$lon,
                lat = df$lat,
                popup = indNames(x),
                color = ic,
                opacity = ind.circle.transparency,
                weight = ind.circle.cex
            )
    }
    
    if (pop.labels) {
        m <- m %>%
            leaflet::addLabelOnlyMarkers(
                lng = centers[, "lon"],
                lat = centers[, "lat"],
                label = popNames(x),
                labelOptions = leaflet::labelOptions(
                    noHide = TRUE,
                    direction = "top",
                    textOnly = TRUE,
                    textsize = paste0(pop.labels.cex, "px")
                )
            )
    }
    
    if (scale.bar) {
        m <- m %>%
            leaflet::addScaleBar(
                position = "bottomright",
                options = leaflet::scaleBarOptions(
                    metric        = TRUE,
                    imperial      = FALSE,
                    maxWidth      = 300,
                    updateWhenIdle = TRUE 
                )
            )
    }
    
    if (!is.null(matrix)) {
        
        # align the matrix to the object: by name when dimnames match the
        # population/individual names, otherwise in object order
        if (nrow(matrix) == nPop(x)) {
            nms <- popNames(x)
            xys <- centers
        } else {
            nms <- indNames(x)
            xys <- df
        }
        named <- !is.null(rownames(matrix)) && !is.null(colnames(matrix))
        if (named && !anyDuplicated(nms) &&
            all(nms %in% rownames(matrix)) && all(nms %in% colnames(matrix))) {
            matrix <- matrix[nms, nms, drop = FALSE]
        } else if (named && verbose >= 2) {
            cat(warn(
                "  Warning: matrix dimnames do not match the population or",
                "individual names; rows and columns are taken in object order.\n"
            ))
        }
        
        # standardize
        if (standard) {
            matrix[, ] <-
                ((matrix[, ] - min(matrix, na.rm = TRUE)) / 
                   (max(matrix, na.rm = TRUE) - 
                      min(matrix, na.rm = TRUE))) * 9 + 1
        }
        
        vals <- unique(as.vector(matrix))
        vals <- vals[!is.na(vals)]
        
        if (is.null(palette.links)) {
            palette.links <- gl.colors("div", verbose = 0)(length(vals))
        }
        
        qpal <- leaflet::colorNumeric(palette = palette.links, domain = vals)
        
        if (symmetric) {
            for (ii in 1:nrow(matrix)) {
                for (i in ii:nrow(matrix)) {
                    v <- matrix[i, ii]
                    # lower triangle only, no self links, no missing or zero
                    if (i != ii && !is.na(v) && v > 0) {
                        m <- m %>%
                            leaflet::addPolylines(
                                lng = c(xys[i, "lon"], xys[ii, "lon"]),
                                lat = c(xys[i, "lat"], xys[ii, "lat"]),
                                color = qpal(v),
                                weight = v,
                                opacity = 1
                            )
                    }
                }
            }
            m <- m %>% leaflet::addLegend(
                pal = qpal, 
                values = vals, 
                group = "addPolylines", 
                position = "bottomleft",
                title = legend.title) 
        }
        
        if (!symmetric) {
            for (i in 1:nrow(matrix)) {
                for (ii in 1:nrow(matrix)) {
                    if (i != ii) {
                        from <- xys[i, ]
                        to <- xys[ii, ]
                        v1 <- matrix[i, ii]
                        v2 <- matrix[ii, i]
                        # nothing to draw in this direction
                        if (is.na(v1)) {
                            next
                        }
                        if (is.na(v2)) {
                            lcols <- "#333333"
                        } else if (v1 > v2) {
                            lcols <- "#FFAA00"
                        } else if (v1 < v2) {
                            lcols <- "#00AAFF"
                        } else {
                            lcols <- "#00AA00"
                        }
                        m <- m %>%
                            leaflet.minicharts::addFlows(
                                lng0 = as.numeric(from["lon"]),
                                lng1 = as.numeric(to["lon"]),
                                lat0 = as.numeric(from["lat"]),
                                lat1 = as.numeric(to["lat"]),
                                flow = v1,
                                color = lcols,
                                maxThickness = 10,
                                minThickness = 0,
                                maxFlow = max(matrix, na.rm = TRUE),
                                opacity = 0.8
                            )
                    }
                }
            }
        }
    }
    
    if (!is.null(raster.image)) {
        r <- terra::rast(raster.image)
        m <- m %>% 
            leaflet::addRasterImage(r, 
                                    opacity = raster.opacity,
                                    colors = raster.colors)
    } 
    
    # FLAG SCRIPT END
    
    if (verbose >= 1) {
        cat(report("Completed:", funname, "\n"))
    }
    
    # RETURN
    
    return(m)
    
}
