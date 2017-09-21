#' PlotRoutes
#'
#' Plot map with BBS route locations
#' @param alpha 4-letter alpha code for species of interest
#' @param psi Should occupancy probability (in year 1) be plotted? (*doesn't work yet*)
#' @export

PlotRoutes <- function(alpha){
  dat <-  readRDS(here::here(paste0("inst/output/", alpha, "/bbs_data.rds")))

  rts <- data.frame(Longitude = dat$lon, Latitude = dat$lat, z = apply(dat$h, 1, function(x) max(x, na.rm = TRUE)))
  obs <- dplyr::filter(rts, z == 1)
  buff <- dplyr::filter(rts, z == 0)
  usa <- ggplot2::map_data("state")
  canada <- ggplot2::map_data("world", "Canada")
  mexico <- ggplot2::map_data("world", "Mexico")


  xmin <- min(dat$lon) - 2
  xmax <- max(dat$lon) + 2

  ymin <- min(dat$lat) - 2
  ymax <- max(dat$lat) + 2


  p <- ggplot() + theme(panel.background = element_rect(fill = "#CDD2D4", color = NA), axis.title = element_blank())
  p <- p + coord_map(projection = "lambert", lat0 = 25, lat1 = 50, xlim = c(xmin, xmax), ylim = c(ymin, ymax))
  p <- p + geom_polygon(data = usa, aes(x=long, y = lat, group = group), fill = "#F0F0F1",  color="#909090") +
    geom_polygon(data = canada, aes(x=long, y = lat, group = group), fill = "#F0F0F1", color="#909090") +
    geom_polygon(data = mexico, aes(x=long, y = lat, group = group),  fill = "#F0F0F1", color="#909090")
  p <- p + geom_point(data = buff, aes(x = Longitude, y = Latitude), color = "grey70", size = 2)
  p <- p + geom_point(data = obs, aes(x = Longitude, y = Latitude), color = "#c45b4d", size = 3)
  p <- p + scale_x_continuous(breaks = seq(from = 10*(xmin%/%10 + as.logical(xmin%%10)),
                                           to = 10*(xmax%/%10 + as.logical(xmax%%10)), by = 10))
  p <- p + scale_y_continuous(breaks = seq(from = 5*(ymin%/%5 + as.logical(ymin%%5)),
                                           to = 5*(ymax%/%5 + as.logical(ymax%%5)), by = 10))
  p
}


#' PlotLat
#'
#' Plot graph of changes in latitude indices
#' @param alpha 4-letter alpha code for species of interest
#' @export

PlotLat <- function(alpha, ci = FALSE){
  indices <- read.csv(here::here(paste0('inst/output/', alpha, '/indices.csv')))
  lat.indices <- dplyr::filter(indices, ind == "n.lat" | ind == "s.lat" |
                               ind == "n.core" | ind == "s.core" | ind == "avg.lat")


  p <- ggplot(lat.indices, aes(x = Year, y = Est, group = ind, color = ind))
  p <- p + geom_line(aes(linetype = ind))
  p <- p + scale_y_continuous("Latitude")
  p <- p + scale_linetype_manual(values = c("solid", "dashed", "longdash", "dashed", "longdash"))
  p <- p + scale_color_manual(values = c("black", "grey35", "grey50", "grey35", "grey50"))
  p <- p + scale_x_continuous(breaks = seq(from = min(lat.indices$Year),
                                           to = max(lat.indices$Year),
                                           by = 2))
  p <- p + theme(legend.position = "none")

  if(ci){
    p <- p + geom_ribbon(aes(ymin = LCI, ymax = UCI), alpha = 0.2)
    p
  }else{
    p
  }
}

#' PlotLon
#'
#' Plot graph of changes in longitude indices
#' @param alpha 4-letter alpha code for species of interest
#' @export

PlotLon <- function(alpha, ci = FALSE){
  indices <- read.csv(here::here(paste0('inst/output/', alpha, '/indices.csv')))
  lon.indices <- dplyr::filter(indices, ind == "w.lon" | ind == "e.lon" |
                                 ind == "w.core" | ind == "e.core" | ind == "avg.lon")


  p <- ggplot(lon.indices, aes(x = Est, y = Year, group = ind, color = ind))
  p <- p + geom_path(aes(linetype = ind))
  p <- p + scale_x_continuous("(West)               Longitude                 (East)")
  p <- p + scale_linetype_manual(values = c("solid", "dashed", "longdash", "dashed", "longdash"))
  p <- p + scale_color_manual(values = c("black", "grey35", "grey50", "grey35", "grey50"))
  p <- p + scale_y_continuous(breaks = seq(from = min(lon.indices$Year),
                                           to = max(lon.indices$Year),
                                           by = 2))
  p <- p + theme(legend.position = "none")

  if(ci){
    p <- p + geom_ribbon(aes(xmin = LCI, xmax = UCI), alpha = 0.2)
    p
  }else{
    p
  }
}

#' PlotPsi
#'
#' Plot graph of changes in range-wide occupancy
#' @param alpha 4-letter alpha code for species of interest
#' @export

PlotPsi <- function(alpha){
  indices <- read.csv(here::here(paste0('inst/output/', alpha, '/indices.csv')))
  psi.indices <- dplyr::filter(indices, ind == "avg.psi")
  y.min <- 0
  y.max <- plyr::round_any(max(psi.indices$UCI), 0.1, f = ceiling)

  p <- ggplot(psi.indices, aes(x = Year, y = Est, group = ind))
  p <- p + geom_point()
  p <- p + geom_line()
  p <- p + geom_point(color = "white", size = 6)
  p <- p + geom_point(size = 4)
  p <- p + geom_ribbon(aes(ymin = LCI, ymax = UCI), alpha = 0.2)
  p <- p + scale_y_continuous("Occupancy", limits = c(y.min, y.max))
  p <- p + scale_x_continuous(breaks = seq(from = min(psi.indices$Year),
                                           to = max(psi.indices$Year),
                                           by = 2))
  p
}

#' MapPsi
#'
#' Annual occupancy probability maps
#' @param alpha 4-letter alpha code for species of interest
#' @param proj Should maps be projected? (Probably but much slower)
#' @export


MapPsi <- function(alpha, proj = TRUE){
  dat <-  readRDS(here::here(paste0("inst/output/", alpha, "/bbs_data.rds")))
  psi <- readRDS(here::here(paste0('inst/output/', alpha, '/psi.rds')))
  indices <- read.csv(here::here(paste0('inst/output/', alpha, '/indices.csv')))
  limits <- dplyr::filter(indices, ind %in% c("s.lat", "n.lat"))
  core <- dplyr::filter(indices, ind %in% c("s.core", "n.core"))
  center <- dplyr::filter(indices, ind %in% c("avg.lat"))

  usa <- ggplot2::map_data("state")
  canada <- ggplot2::map_data("world", "Canada")
  mexico <- ggplot2::map_data("world", "Mexico")

  xmin <- min(psi$Longitude) - 1
  xmax <- max(psi$Longitude) + 1

  ymin <- min(psi$Latitude) - 1
  ymax <- max(psi$Latitude) + 1


  if(proj){
    p <- ggplot() + facet_wrap(~Year)
    p <- p + geom_polygon(data = usa, aes(x=long, y = lat, group = group), fill = "#F0F0F1",  color="#909090") +
      geom_polygon(data = canada, aes(x=long, y = lat, group = group), fill = "#F0F0F1", color="#909090") +
      geom_polygon(data = mexico, aes(x=long, y = lat, group = group),  fill = "#F0F0F1", color="#909090")
    p <- p + geom_tile(data = psi, aes(x = Longitude, y = Latitude, fill = Psi))
    p <- p + geom_polygon(data = usa, aes(x=long, y = lat, group = group), fill = NA,  color="#909090") +
      geom_polygon(data = canada, aes(x=long, y = lat, group = group), fill = "#F0F0F1", color="#909090")
    p <- p + geom_polygon(data = mexico, aes(x=long, y = lat, group = group),  fill = NA, color="#909090")
    p <- p + coord_map(projection = "lambert", lat0 = 25, lat1 = 50, xlim = c(xmin, xmax), ylim = c(ymin, ymax))
    p <- p + geom_hline(data = limits, aes(yintercept = Est, group = ind), linetype = "dashed", color = "grey70")
    p <- p + geom_hline(data = core, aes(yintercept = Est, group = ind), linetype = "longdash", color = "grey70")
    p <- p + geom_hline(data = center, aes(yintercept = Est), color = "grey70")
    p <- p + scale_fill_continuous(low = "#F0F0F1", high = "#ef6e0b")
    p <- p + scale_x_continuous("Longitude",breaks = seq(from = 10*(xmin%/%10 + as.logical(xmin%%10)),
                                                         to = 10*(xmax%/%10 + as.logical(xmax%%10)), by = 10))
    p <- p + scale_y_continuous("Latitude", breaks = seq(from = 5*(ymin%/%5 + as.logical(ymin%%5)),
                                                         to = 5*(ymax%/%5 + as.logical(ymax%%5)), by = 10))
    p <- p + theme(panel.background = element_rect(fill = "#CDD2D4", color = NA), axis.title = element_blank())
    p
  }else{
    p <- ggplot() + facet_wrap(~Year) + coord_fixed(ratio = 1.3, xlim = c(xmin, xmax), ylim = c(ymin, ymax))
    p <- p + geom_polygon(data = usa, aes(x=long, y = lat, group = group), fill = "#F0F0F1",  color="#909090") +
      geom_polygon(data = canada, aes(x=long, y = lat, group = group), fill = "#F0F0F1", color="#909090") +
      geom_polygon(data = mexico, aes(x=long, y = lat, group = group),  fill = "#F0F0F1", color="#909090")
    p <- p + geom_raster(data = psi, aes(x = Longitude, y = Latitude, fill = Psi)) + facet_wrap(~Year)
    p <- p + geom_polygon(data = usa, aes(x=long, y = lat, group = group), fill = NA,  color="#909090") +
      geom_polygon(data = canada, aes(x=long, y = lat, group = group), fill = NA, color="#909090")
    p <- p + geom_polygon(data = mexico, aes(x=long, y = lat, group = group),  fill = "#F0F0F1", color="#909090")
    p <- p + geom_hline(data = limits, aes(yintercept = Est, group = ind), linetype = "dashed", color = "grey70")
    p <- p + geom_hline(data = core, aes(yintercept = Est, group = ind), linetype = "longdash", color = "grey70")
    p <- p + geom_hline(data = center, aes(yintercept = Est), color = "grey70")
    p <- p + scale_fill_continuous(low = "#F0F0F1", high = "#ef6e0b")
    p <- p + theme(panel.background = element_rect(fill = "#CDD2D4", color = NA), axis.title = element_blank())
    p <- p + scale_x_continuous("Longitude", breaks = seq(from = 10*(xmin%/%10 + as.logical(xmin%%10)),
                                                          to = 10*(xmax%/%10 + as.logical(xmax%%10)), by = 10))
    p <- p + scale_y_continuous("Latitude", breaks = seq(from = 5*(ymin%/%5 + as.logical(ymin%%5)),
                                                         to = 5*(ymax%/%5 + as.logical(ymax%%5)), by = 10))
    p
  }


}


#' MapDiff
#'
#' Initial/final & difference occupancy maps
#' @param alpha 4-letter alpha code for species of interest
#' @export


MapDiff <- function(alpha){
  psi <- readRDS(here::here(paste0('inst/output/', alpha, '/psi.rds')))

  psi1 <- dplyr::filter(psi, Year == min(Year))
  psil <- dplyr::filter(psi, Year == max(Year))
  psid <- dplyr::mutate(psil, Diff = Psi - psi1$Psi)
  psid <- dplyr::select(psid, -Psi)
  psid <- dplyr::rename(psid, Psi = Diff)

  psi2 <- dplyr::bind_rows(psi1, psid)
  psi2 <- dplyr::mutate(psi2, Period = rep(c("Initial", "Difference"), each = nrow(psi1)))

  usa <- ggplot2::map_data("state")
  canada <- ggplot2::map_data("world", "Canada")
  mexico <- ggplot2::map_data("world", "Mexico")

  xmin <- min(psi$Longitude) - 1
  xmax <- max(psi$Longitude) + 1

  ymin <- min(psi$Latitude) - 1
  ymax <- max(psi$Latitude) + 1

  p <- ggplot()
  p <- p + geom_polygon(data = usa, aes(x=long, y = lat, group = group), fill = "#F0F0F1",  color="#909090") +
    geom_polygon(data = canada, aes(x=long, y = lat, group = group), fill = "#F0F0F1", color="#909090") +
    geom_polygon(data = mexico, aes(x=long, y = lat, group = group),  fill = "#F0F0F1", color="#909090")
  p <- p + geom_tile(data = psi1, aes(x = Longitude, y = Latitude, fill = Psi))
  p <- p + geom_polygon(data = usa, aes(x=long, y = lat, group = group), fill = NA,  color="#909090") +
    geom_polygon(data = canada, aes(x=long, y = lat, group = group), fill = NA, color="#909090")
  p <- p + geom_polygon(data = mexico, aes(x=long, y = lat, group = group),  fill = "#F0F0F1", color="#909090")
  p <- p + coord_map(projection = "lambert", lat0 = 25, lat1 = 50, xlim = c(xmin, xmax), ylim = c(ymin, ymax))
  p <- p + scale_fill_continuous(low = "#F0F0F1", high = "#ef6e0b",  name =  "Occupancy \nProb.")
  p <- p + scale_x_continuous("Longitude",breaks = seq(from = 10*(xmin%/%10 + as.logical(xmin%%10)),
                                                       to = 10*(xmax%/%10 + as.logical(xmax%%10)), by = 10))
  p <- p + scale_y_continuous("Latitude", breaks = seq(from = 5*(ymin%/%5 + as.logical(ymin%%5)),
                                                       to = 5*(ymax%/%5 + as.logical(ymax%%5)), by = 10))
  p <- p + theme(panel.background = element_rect(fill = "#CDD2D4", color = NA),
                 axis.title = element_blank(), legend.title = element_text(size = 10))
  p <- p + labs(title = paste0(unique(psi1$Year), " (Initial)"))

  q <- ggplot()
  q <- q + geom_polygon(data = usa, aes(x=long, y = lat, group = group), fill = "#F0F0F1",  color="#909090") +
    geom_polygon(data = canada, aes(x=long, y = lat, group = group), fill = "#F0F0F1", color="#909090") +
    geom_polygon(data = mexico, aes(x=long, y = lat, group = group),  fill = "#F0F0F1", color="#909090")
  q <- q + geom_tile(data = psil, aes(x = Longitude, y = Latitude, fill = Psi))
  q <- q + geom_polygon(data = usa, aes(x=long, y = lat, group = group), fill = NA,  color="#909090") +
    geom_polygon(data = canada, aes(x=long, y = lat, group = group), fill = NA, color="#909090")
  q <- q + geom_polygon(data = mexico, aes(x=long, y = lat, group = group),  fill = "#F0F0F1", color="#909090")
  q <- q + coord_map(projection = "lambert", lat0 = 25, lat1 = 50, xlim = c(xmin, xmax), ylim = c(ymin, ymax))
  q <- q + scale_fill_continuous(low = "#F0F0F1", high = "#ef6e0b", name =  "Occupancy \nProb.")
  q <- q + scale_x_continuous("Longitude",breaks = seq(from = 10*(xmin%/%10 + as.logical(xmin%%10)),
                                                       to = 10*(xmax%/%10 + as.logical(xmax%%10)), by = 10))
  q <- q + scale_y_continuous("Latitude", breaks = seq(from = 5*(ymin%/%5 + as.logical(ymin%%5)),
                                                       to = 5*(ymax%/%5 + as.logical(ymax%%5)), by = 10))
  q <- q + theme(panel.background = element_rect(fill = "#CDD2D4", color = NA),
                 axis.title = element_blank(), legend.title = element_text(size = 10))
  q <- q + labs(title = unique(psil$Year))


  r <- ggplot()
  r <- r + geom_polygon(data = usa, aes(x=long, y = lat, group = group), fill = "#F0F0F1",  color="#909090") +
    geom_polygon(data = canada, aes(x=long, y = lat, group = group), fill = "#F0F0F1", color="#909090") +
    geom_polygon(data = mexico, aes(x=long, y = lat, group = group),  fill = "#F0F0F1", color="#909090")
  r <- r + geom_tile(data = psid, aes(x = Longitude, y = Latitude, fill = Psi))
  r <- r + geom_polygon(data = usa, aes(x=long, y = lat, group = group), fill = NA,  color="#909090") +
    geom_polygon(data = canada, aes(x=long, y = lat, group = group), fill = NA, color="#909090")
  r <- r + geom_polygon(data = mexico, aes(x=long, y = lat, group = group),  fill = "#F0F0F1", color="#909090")
  r <- r + coord_map(projection = "lambert", lat0 = 25, lat1 = 50, xlim = c(xmin, xmax), ylim = c(ymin, ymax))
  r <- r + scale_fill_gradient2(low = "#c45b4d", mid ="#F0F0F1", high = "#268bd2", name = "Difference")
  r <- r + scale_x_continuous("Longitude",breaks = seq(from = 10*(xmin%/%10 + as.logical(xmin%%10)),
                                                       to = 10*(xmax%/%10 + as.logical(xmax%%10)), by = 10))
  r <- r + scale_y_continuous("Latitude", breaks = seq(from = 5*(ymin%/%5 + as.logical(ymin%%5)),
                                                       to = 5*(ymax%/%5 + as.logical(ymax%%5)), by = 10))
  r <- r + theme(panel.background = element_rect(fill = "#CDD2D4", color = NA),
                 axis.title = element_blank(), legend.title = element_text(size = 10))

  s <- gridExtra::grid.arrange(gridExtra::arrangeGrob(p, q, ncol = 2), r, ncol = 1)
  grid::grid.draw(s)
}

#' md2html
#'
#' Function to covert .md file to html (ht https://github.com/richfitz/modeladequacy)
#' @export

md2html <- function(filename, dest = NULL) {
  if(is.null(dest)){ dest <- paste0(tools::file_path_sans_ext(filename), ".html")}
  opts <- setdiff(markdownHTMLOptions(TRUE), 'base64_images')
  markdownToHTML(filename, dest, options = opts)
}


#' PlotIndex
#'
#' Function to plot multi-species indices of range change
#' @export

PlotIndex <- function(group_name, lat = TRUE, start.year = 1972){
  indices <- read.csv(here::here(paste0("inst/output/indices/", group_name, ".csv")))

  if(lat){
    indices <- dplyr::filter(indices, grepl("lat", ind) & Year >= start.year)
    indices$ind <- droplevels(indices$ind)

    indices$ind <- dplyr::recode(indices$ind, s.lat = "Southern range limit",
                                 avg.lat = "Mean breeding latitude",
                                 n.lat = "Northern range limit")
    indices$ind <- factor(indices$ind, levels = c("Northern range limit", "Mean breeding latitude",  "Southern range limit"))
    tot_index <- dplyr::filter(indices, Species == "Composite")
    spp_index <- dplyr::filter(indices, Species != "Composite")

    p <- ggplot()
    p <- p + geom_line(data = spp_index, aes(x = Year, y = value,
                   group = interaction(ind, Species)), alpha = 0.15, size = 0.7)
    p <- p + scale_color_manual(values = rep("black", 5), guide = FALSE)
    p <- p + geom_line(data = tot_index, aes(x = Year, y = value, group = ind),
                       color = "#cb4b16", size = 1.25)
    p <- p + scale_y_continuous("Latitude")
    p <- p + facet_wrap(~ind, ncol = 1, scales = "free")

  }else{
    indices <- dplyr::filter(indices, grepl("lon", ind) & Year >= start.year)
    indices$ind <- droplevels(indices$ind)

    indices$ind <- dplyr::recode(indices$ind, w.lon = "Western range limit", avg.lon = "Mean breeding longitude", e.lon = "Eastern range limit")
    indices$ind <- factor(indices$ind, levels = c("Western range limit", "Mean breeding longitude",  "Eastern range limit"))

    tot_index <- dplyr::filter(indices, Species == "Composite")
    spp_index <- dplyr::filter(indices, Species != "Composite")

    p <- ggplot()
    p <- p + geom_path(data = spp_index, aes(x = value, y = Year,
                                             group = interaction(ind, Species)), alpha = 0.15, size = 0.7)
    p <- p + scale_color_manual(values = rep("black", 5), guide = FALSE)
    p <- p + geom_path(data = tot_index, aes(x = value, y = Year, group = ind),
                       color = "#cb4b16", size = 1.25)
    p <- p + scale_x_continuous("Longitude")
    p <- p + facet_wrap(~ind, scales = "free")
  }

  p

}
