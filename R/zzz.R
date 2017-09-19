.onLoad <- function(libname="BBSclim", pkgname="BBSclim"){
  options(digits=4)
  library(ggplot2)
  theme_crushing <- function(base_size = 12, base_family = "") {
    half_line <- base_size/2
    theme(
      # Elements in this first block aren't used directly, but are inherited
      # by others
      line =               element_line(size = 0.5, linetype = 1, colour = "black",
                                        lineend = "butt"),
      rect =               element_rect(fill = "white", colour = "black", size = 0.5, linetype = 1),
      text =               element_text(family = base_family, face = "plain", colour = "black",
                                        size = base_size,
                                        hjust = 0.5, vjust = 0.5, angle = 0, lineheight = 0.9, margin = margin(),
                                        debug = FALSE),
      axis.text =          element_text(colour = "grey40"),
      axis.title =         element_text(colour = "grey20", vjust = 0.35),
      strip.text =         element_text(size = rel(0.8)),

      axis.line =          element_line(),
      axis.line.x =        element_line(size=.7, color = "grey60"),
      axis.line.y =        element_line(size=.7, color = "grey60"),
      axis.text.x =        element_text(size = base_size*1.1, lineheight = 0.9,
                                        margin = margin(t = 0.8 * half_line/2), vjust = 1),
      axis.text.y =        element_text(size = base_size*1.1, lineheight = 0.9,
                                        margin = margin(r = 0.8 * half_line/2), vjust = 0.5),
      axis.ticks =         element_line(colour = "grey60", size = 0.2),
      axis.title.x =       element_text(size = base_size*1.4, vjust = 0.3,
                                        margin = margin(t = 10, b = 0.8 * half_line/2)),
      axis.title.y =       element_text(size = base_size*1.4, angle = 90, vjust = 1,
                                        margin = margin(r = 10, l = 0.8 * half_line/2)),
      axis.ticks.length =  grid::unit(0.3, "lines"),

      legend.background =  element_rect(colour = NA),
      legend.margin =      grid::unit(0.2, "cm"),
      legend.key =         element_rect(colour = "grey80"),
      legend.key.size =    grid::unit(1.2, "lines"),
      legend.key.height =  NULL,
      legend.key.width =   NULL,
      legend.text =        element_text(size = base_size * 0.8),
      legend.text.align =  NULL,
      legend.title =      element_blank(),
      legend.title.align = NULL,
      legend.position =    "right",
      legend.direction =   NULL,
      legend.justification = "center",
      legend.box =         NULL,

      panel.background =   element_rect(fill = "white", colour = NA),
      panel.border =       element_rect(fill = NA, color = NA,size=.5),
      panel.grid.major =   element_blank(),
      panel.grid.minor =   element_blank(),
      panel.margin =       grid::unit(half_line, "pt"),
      panel.margin.x =     NULL,
      panel.margin.y =     NULL,
      panel.ontop =        FALSE,

      strip.background =   element_rect(fill = NA, colour = NA),
      strip.text.x =       element_text(size = base_size, margin = margin(t = half_line, b = half_line)),
      strip.text.y =       element_text(size = base_size, angle = -90, margin = margin(l = half_line, r = half_line)),
      strip.switch.pad.grid = unit(0.1, "cm"),
      strip.switch.pad.wrap = unit(0.1, "cm"),

      plot.background =   element_rect(colour = NA),
      plot.title =        element_text(size = base_size * 1.7, face="bold",vjust=2, margin = margin(b = half_line * 1.2)),
      plot.margin =       grid::unit(c(1, 1.5, 0.8, 0.8), "lines"),

      complete = TRUE
    )
  }

  theme_set(theme_crushing())
  scale_colour_discrete <- function(...) ggthemes::scale_color_solarized()
  update_geom_defaults("point", list(size = 3))
  update_geom_defaults("line", list(size = 0.8))
}
