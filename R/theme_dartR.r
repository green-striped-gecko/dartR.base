# SET THEME FOR PLOTS
#' @name theme_dartR
#' @title Default theme for dartR plots
#' @family environment
#' 
#' @description
#' This is the theme used as default for dartR plots.
#' This function controls all non-data display elements in the plots.

#' @param base_size base font size, given in pts.
#' @param base_family base font family
#' @param base_line_size base size for line elements
#' @param base_rect_size base size for rect elements
#' @return a the standard dartR theme to be used in ggplots
#' @examples
#' ggplot(data.frame(dummy=rnorm(1000)),aes(dummy)) +
#' geom_histogram(binwidth=0.1) + theme_dartR()
#' 
#' @export

# Version v.2023.2

# The half-line (base-fontsize / 2) sets up the basic vertical rhythm of the
# theme. Most margins will be set to this value.  However, when we work with
# relative sizes, we may want to multiply `half_line` with the appropriate
# relative size. This applies in particular for axis tick sizes. And also, for
# axis ticks and axis titles, `half_size` is too large a distance, and we use
# `half_size/2` instead.
theme_dartR <- function(base_size = 11,
                        base_family = "",
                        base_line_size = base_size / 22,
                        base_rect_size = base_size / 22) {
  half_line <- base_size / 2
  
  # Throughout the theme, we use three font sizes, `base_size` (`rel(1)`)
  #for
  # normal, `rel(0.8)` for small, and `rel(1.2)` for large.
  
  # Elements in this first block aren't used directly, but are inherited by
  # others
  t <-
    theme(
      line = element_line(
        color = "black",
        size = base_line_size,
        linetype = 1,
        lineend = "butt"
      ),
      rect = element_rect(
        fill = "white",
        color = "black",
        size = base_rect_size,
        linetype = 1
      ),
      text = element_text(
        family = base_family,
        face = "plain",
        color = "black",
        size = base_size,
        lineheight = 0.9,
        hjust = 0.5,
        vjust = 0.5,
        angle = 0,
        margin = margin(),
        debug = FALSE
      ),
      axis.line = element_blank(),
      axis.line.x = NULL,
      axis.line.y = NULL,
      axis.text = element_text(size = rel(1), color = "black"),
      axis.text.x = element_text(margin = margin(t = 0.8 * half_line / 2),
                                 vjust = 1),
      axis.text.x.top = element_text(margin = margin(b = 0.8 * half_line / 2), vjust = 0),
      axis.text.y = element_text(margin = margin(r = 0.8 *
                                                   half_line / 2), hjust = 1),
      axis.text.y.right = element_text(margin = margin(l = 0.8 * half_line /
                                                         2), hjust = 0),
      axis.ticks = element_line(color = "gray80"),
      axis.ticks.length = unit(half_line / 2, "pt"),
      axis.ticks.length.x = NULL,
      axis.ticks.length.x.top = NULL,
      axis.ticks.length.x.bottom = NULL,
      axis.ticks.length.y = NULL,
      axis.ticks.length.y.left = NULL,
      axis.ticks.length.y.right = NULL,
      axis.title.x = element_text(
        margin = margin(t = half_line / 2),
        vjust = 1,
        face = "bold"
      ),
      axis.title.x.top = element_text(margin = margin(b = half_line / 2), vjust = 0),
      axis.title.y = element_text(
        angle = 90,
        margin = margin(r = half_line / 2),
        vjust = 1,
        face = "bold"
      ),
      axis.title.y.right = element_text(
        angle = -90,
        margin = margin(l = half_line / 2),
        vjust = 0
      ),
      legend.background = element_rect(color = "transparent"),
      legend.spacing = unit(2 * half_line, "pt"),
      legend.spacing.x = NULL,
      legend.spacing.y = NULL,
      legend.margin = margin(half_line, half_line, half_line, half_line),
      legend.key = element_rect(fill = "white",
                                color = NA),
      legend.key.size = unit(1.2, "lines"),
      legend.key.height = NULL,
      legend.key.width = NULL,
      legend.text = element_text(size = rel(1.2),
                                 face = "bold"),
      legend.text.align = NULL,
      legend.title = element_text(
        size = rel(1),
        face = "bold",
        hjust = 0
      ),
      legend.title.align = NULL,
      legend.position = "right",
      legend.direction = NULL,
      legend.justification = "center",
      legend.box = NULL,
      legend.box.margin = margin(0,
                                 0, 0, 0, "cm"),
      legend.box.background = element_blank(),
      legend.box.spacing = unit(2 * half_line, "pt"),
      panel.background = element_rect(fill = "white",
                                      color = NA),
      panel.border = element_blank(),
      panel.grid = element_line(color = "gray80"),
      panel.grid.minor = element_line(size = rel(0.5)),
      panel.spacing = unit(half_line, "pt"),
      panel.spacing.x = NULL,
      panel.spacing.y = NULL,
      panel.ontop = FALSE,
      strip.background = element_rect(fill = "white",
                                      color = "black"),
      strip.text = element_text(
        color = "black",
        size = rel(1),
        face = "bold",
        margin = margin(
          0.8 * half_line,
          0.8 *
            half_line,
          0.8 * half_line,
          0.8 * half_line
        )
      ),
      strip.text.x = element_text(size = 14, face = "bold"),
      strip.text.y = element_text(angle = -90),
      strip.text.y.left = element_text(angle = 90),
      strip.placement = "inside",
      strip.placement.x = NULL,
      strip.placement.y = NULL,
      strip.switch.pad.grid = unit(half_line / 2,
                                   "pt"),
      strip.switch.pad.wrap = unit(half_line / 2, "pt"),
      plot.background = element_rect(color = "white"),
      plot.title = element_text(
        face = "bold",
        size = rel(1.2),
        hjust = 0.5,
        vjust = 1,
        margin = margin(b = half_line)
      ),
      plot.title.position = "panel",
      plot.subtitle = element_text(
        hjust = 0.5,
        vjust = 1,
        margin = margin(b = half_line)
      ),
      plot.caption = element_text(
        size = rel(0.8),
        hjust = 1,
        vjust = 1,
        margin = margin(t = half_line)
      ),
      plot.caption.position = "panel",
      plot.tag = element_text(
        size = rel(1.2),
        hjust = 0.5,
        vjust = 0.5
      ),
      plot.tag.position = "topleft",
      plot.margin = margin(half_line, half_line, half_line, half_line),
      complete = TRUE
    )
}
