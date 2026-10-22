# Shared by the population-colour regression tests (one per plotting
# function). platypus.gl with levels that differ from both the order of
# first appearance (TENTERFIELD first) and alphabetical order, so a colour
# keyed by position to either order lands on the wrong population.
pop_colours_gl <- function() {
  x <- platypus.gl
  pop(x) <- factor(as.character(pop(x)),
                   levels = c("SEVERN_BELOW", "TENTERFIELD", "SEVERN_ABOVE"))
  x
}
pop_colours_named <- c(TENTERFIELD = "green", SEVERN_ABOVE = "red",
                       SEVERN_BELOW = "blue")
# an unnamed vector is assigned in the order of levels(pop(x))
pop_colours_unnamed <- c("red", "blue", "green")
pop_colours_expected <- c(SEVERN_BELOW = "blue", TENTERFIELD = "green",
                          SEVERN_ABOVE = "red")
pop_colours_expected_unnamed <- c(SEVERN_BELOW = "red", TENTERFIELD = "blue",
                                  SEVERN_ABOVE = "green")

# colour the discrete `aes` scale of a ggplot gives each population in `pops`;
# for a plot with no such scale (colours set as a constant per row), the
# colours drawn in each facet panel named by population
pop_colours_of <- function(p, pops, aes = "fill") {
  b <- ggplot2::ggplot_build(p)
  sc <- b$plot$scales$get_scales(aes)
  if (!is.null(sc) && inherits(sc, "ScaleDiscrete")) {
    return(setNames(as.character(sc$map(pops)), pops))
  }
  lay <- b$layout$layout
  d <- b$data[[1]]
  pop <- as.character(lay$pop[match(d$PANEL, lay$PANEL)])
  unlist(tapply(d[[aes]], pop, function(v) paste(unique(v), collapse = "/"),
                simplify = FALSE))[pops]
}

# the plots a function saved to plot.dir as plot.file (.RDS), split into
# patchwork panels
pop_colours_plots <- function(dir) {
  f <- list.files(dir, "\\.RDS$", full.names = TRUE, ignore.case = TRUE)
  unlist(lapply(f, function(ff) {
    p <- readRDS(ff)
    if (inherits(p, "patchwork")) {
      lapply(seq_len(length(p$patches$plots) + 1), function(i) p[[i]])
    } else {
      list(p)
    }
  }), recursive = FALSE)
}
