#' GetIndices
#'
#' Estimate indices of range dynamics from annual occupancy estimates
#' @param alpha Vector of alpha codes for species of interest
#' @return Data frame containing annual estimates of the following indices:
#' @return    avg.psi = proportion of area occupied (i.e., range size)
#' @return    s.lat = southern range limit
#' @return    s.core = southern core range limit
#' @return    n.lat = northern range limit
#' @return    n.core = northern core range limit
#' @return    w.lat = western range limit
#' @return    w.core = western core range limit
#' @return    e.lat = eastern range limit
#' @return    e.core = eastern core range limit
#' @return    avg.lat = occupancy-weighted mean breeding latitude
#' @return    avg.lon = occupancy-weighted mean breeding longitude
#' @export

GetIndices <- function(spp = NULL, alpha = NULL){
  if(!is.null(spp)){
    cores <- parallel::detectCores()
    if(length(spp) < cores) cores <- length(spp)
    doParallel::registerDoParallel(cores = cores)

    ### Run posterior predictive checks in parallel
    indices <- foreach::foreach(i = 1:length(spp), .combine = c,
                                    .packages = c("dplyr", "BayesCorrOcc")) %dopar%{
                                      occ <- readRDS(paste0('inst/output/', spp[i], '/occ.rds'))
                                      dat <- readRDS(paste0("inst/output/", spp[i], "/bbs_data.rds"))

                                      years <- seq(from = dat$start_year, to = dat$end_year)

                                      avg.psi <- s.lat <- s.core <- n.lat <- n.core <- w.lon <- w.core <- e.lon <- e.core <- avg.lat <- avg.lon <- data.frame(Est = rep(NA, dat$nYears) , LCI = rep(NA, dat$nYears), UCI = rep(NA, dat$nYears), ind = rep(NA, dat$nYears), Year = years)

                                      for(tt in 1:dat$nYears){
                                        avg.psi$Est[tt] <- mean(apply(occ$occ[,,tt], 1, mean))
                                        avg.psi$LCI[tt] <- quantile(apply(occ$occ[,,tt], 1, mean), probs = 0.025)
                                        avg.psi$UCI[tt] <- quantile(apply(occ$occ[,,tt], 1, mean), probs = 0.975)
                                        avg.psi$ind[tt] <- "avg.psi"

                                        s.lat$Est[tt] <- mean(apply(occ$occ[,,tt], 1, function(x) range.limit(cell.probs = x, prob = 0.975, coord = occ$xy$lat, limit = "south")))
                                        s.lat$LCI[tt] <- quantile(apply(occ$occ[,,tt], 1, function(x) range.limit(cell.probs = x, prob = 0.975, coord = occ$xy$lat, limit = "south")), probs = 0.025)
                                        s.lat$UCI[tt] <- quantile(apply(occ$occ[,,tt], 1, function(x) range.limit(cell.probs = x, prob = 0.975, coord = occ$xy$lat, limit = "south")), probs = 0.975)
                                        s.lat$ind[tt] <- "s.lat"

                                        s.core$Est[tt] <- mean(apply(occ$occ[,,tt], 1, function(x) range.limit(cell.probs = x, prob = 0.75, coord = occ$xy$lat, limit = "south")))
                                        s.core$LCI[tt] <- quantile(apply(occ$occ[,,tt], 1, function(x) range.limit(cell.probs = x, prob = 0.75, coord = occ$xy$lat, limit = "south")), probs = 0.025)
                                        s.core$UCI[tt] <- quantile(apply(occ$occ[,,tt], 1, function(x) range.limit(cell.probs = x, prob = 0.75, coord = occ$xy$lat, limit = "south")), probs = 0.975)
                                        s.core$ind[tt] <- "s.core"

                                        n.lat$Est[tt] <- mean(apply(occ$occ[,,tt], 1, function(x) range.limit(cell.probs = x, prob = 0.975, coord = occ$xy$lat, limit = "north")))
                                        n.lat$LCI[tt] <- quantile(apply(occ$occ[,,tt], 1, function(x) range.limit(cell.probs = x, prob = 0.975, coord = occ$xy$lat, limit = "north")), probs = 0.025)
                                        n.lat$UCI[tt] <- quantile(apply(occ$occ[,,tt], 1, function(x) range.limit(cell.probs = x, prob = 0.975, coord = occ$xy$lat, limit = "north")), probs = 0.975)
                                        n.lat$ind[tt] <- "n.lat"

                                        n.core$Est[tt] <- mean(apply(occ$occ[,,tt], 1, function(x) range.limit(cell.probs = x, prob = 0.75, coord = occ$xy$lat, limit = "north")))
                                        n.core$LCI[tt] <- quantile(apply(occ$occ[,,tt], 1, function(x) range.limit(cell.probs = x, prob = 0.75, coord = occ$xy$lat, limit = "north")), probs = 0.025)
                                        n.core$UCI[tt] <- quantile(apply(occ$occ[,,tt], 1, function(x) range.limit(cell.probs = x, prob = 0.75, coord = occ$xy$lat, limit = "north")), probs = 0.975)
                                        n.core$ind[tt] <- "n.core"

                                        w.lon$Est[tt] <- mean(apply(occ$occ[,,tt], 1, function(x) range.limit(cell.probs = x, prob = 0.975, coord = occ$xy$lon, limit = "west")))
                                        w.lon$LCI[tt] <- quantile(apply(occ$occ[,,tt], 1, function(x) range.limit(cell.probs = x, prob = 0.975, coord = occ$xy$lon, limit = "west")), probs = 0.025)
                                        w.lon$UCI[tt] <- quantile(apply(occ$occ[,,tt], 1, function(x) range.limit(cell.probs = x, prob = 0.975, coord = occ$xy$lon, limit = "west")), probs = 0.975)
                                        w.lon$ind[tt] <- "w.lon"

                                        w.core$Est[tt] <- mean(apply(occ$occ[,,tt], 1, function(x) range.limit(cell.probs = x, prob = 0.75, coord = occ$xy$lon, limit = "west")))
                                        w.core$LCI[tt] <- quantile(apply(occ$occ[,,tt], 1, function(x) range.limit(cell.probs = x, prob = 0.75, coord = occ$xy$lon, limit = "west")), probs = 0.025)
                                        w.core$UCI[tt] <- quantile(apply(occ$occ[,,tt], 1, function(x) range.limit(cell.probs = x, prob = 0.75, coord = occ$xy$lon, limit = "west")), probs = 0.975)
                                        w.core$ind[tt] <- "w.core"

                                        e.lon$Est[tt] <- mean(apply(occ$occ[,,tt], 1, function(x) range.limit(cell.probs = x, prob = 0.975, coord = occ$xy$lon, limit = "east")))
                                        e.lon$LCI[tt] <- quantile(apply(occ$occ[,,tt], 1, function(x) range.limit(cell.probs = x, prob = 0.975, coord = occ$xy$lon, limit = "east")), probs = 0.025)
                                        e.lon$UCI[tt] <- quantile(apply(occ$occ[,,tt], 1, function(x) range.limit(cell.probs = x, prob = 0.975, coord = occ$xy$lon, limit = "east")), probs = 0.975)
                                        e.lon$ind[tt] <- "e.lon"

                                        e.core$Est[tt] <- mean(apply(occ$occ[,,tt], 1, function(x) range.limit(cell.probs = x, prob = 0.75, coord = occ$xy$lon, limit = "east")))
                                        e.core$LCI[tt] <- quantile(apply(occ$occ[,,tt], 1, function(x) range.limit(cell.probs = x, prob = 0.75, coord = occ$xy$lon, limit = "east")), probs = 0.025)
                                        e.core$UCI[tt] <- quantile(apply(occ$occ[,,tt], 1, function(x) range.limit(cell.probs = x, prob = 0.75, coord = occ$xy$lon, limit = "east")), probs = 0.975)
                                        e.core$ind[tt] <- "e.core"

                                        avg.lat$Est[tt] <- mean(apply(occ$occ[,,tt], 1, function(x) sum(occ$xy$lat * x, na.rm = TRUE)/sum(x, na.rm = TRUE)))
                                        avg.lat$LCI[tt] <- quantile(apply(occ$occ[,,tt], 1, function(x) sum(occ$xy$lat * x, na.rm = TRUE)/sum(x, na.rm = TRUE)), probs = 0.025)
                                        avg.lat$UCI[tt] <- quantile(apply(occ$occ[,,tt], 1, function(x) sum(occ$xy$lat * x, na.rm = TRUE)/sum(x, na.rm = TRUE)), probs = 0.975)
                                        avg.lat$ind[tt] <- "avg.lat"

                                        avg.lon$Est[tt] <- mean(apply(occ$occ[,,tt], 1, function(x) sum(occ$xy$lon * x, na.rm = TRUE)/sum(x, na.rm = TRUE)))
                                        avg.lon$LCI[tt] <- quantile(apply(occ$occ[,,tt], 1, function(x) sum(occ$xy$lon * x, na.rm = TRUE)/sum(x, na.rm = TRUE)), probs = 0.025)
                                        avg.lon$UCI[tt] <- quantile(apply(occ$occ[,,tt], 1, function(x) sum(occ$xy$lon * x, na.rm = TRUE)/sum(x, na.rm = TRUE)), probs = 0.975)
                                        avg.lon$ind[tt] <- "avg.lon"
                                      }


                                      indices <- dplyr::bind_rows(avg.psi, s.lat)
                                      indices <- dplyr::bind_rows(indices, s.core)
                                      indices <- dplyr::bind_rows(indices, n.lat)
                                      indices <- dplyr::bind_rows(indices, n.core)
                                      indices <- dplyr::bind_rows(indices, w.lon)
                                      indices <- dplyr::bind_rows(indices, w.core)
                                      indices <- dplyr::bind_rows(indices, e.lon)
                                      indices <- dplyr::bind_rows(indices, e.core)
                                      indices <- dplyr::bind_rows(indices, avg.lat)
                                      indices <- dplyr::bind_rows(indices, avg.lon)

                                      write.csv(indices, file = paste0("inst/output/", spp[i], "/indices.csv"), row.names = FALSE)
                                      return(spp[i])
                                    }
    return(indices)
  }

  if(!is.null(alpha)){
    occ <- readRDS(paste0('inst/output/', alpha, '/occ.rds'))
    dat <- readRDS(paste0("inst/output/", alpha, "/bbs_data.rds"))

    years <- seq(from = dat$start_year, to = dat$end_year)

    avg.psi <- s.lat <- s.core <- n.lat <- n.core <- w.lon <- w.core <- e.lon <- e.core <- avg.lat <- avg.lon <- data.frame(Est = rep(NA, dat$nYears) , LCI = rep(NA, dat$nYears), UCI = rep(NA, dat$nYears), ind = rep(NA, dat$nYears), Year = years)

    for(tt in 1:dat$nYears){
      avg.psi$Est[tt] <- mean(apply(occ$occ[,,tt], 1, mean))
      avg.psi$LCI[tt] <- quantile(apply(occ$occ[,,tt], 1, mean), probs = 0.025)
      avg.psi$UCI[tt] <- quantile(apply(occ$occ[,,tt], 1, mean), probs = 0.975)
      avg.psi$ind[tt] <- "avg.psi"

      s.lat$Est[tt] <- mean(apply(occ$occ[,,tt], 1, function(x) range.limit(cell.probs = x, prob = 0.975, coord = occ$xy$lat, limit = "south")))
      s.lat$LCI[tt] <- quantile(apply(occ$occ[,,tt], 1, function(x) range.limit(cell.probs = x, prob = 0.975, coord = occ$xy$lat, limit = "south")), probs = 0.025)
      s.lat$UCI[tt] <- quantile(apply(occ$occ[,,tt], 1, function(x) range.limit(cell.probs = x, prob = 0.975, coord = occ$xy$lat, limit = "south")), probs = 0.975)
      s.lat$ind[tt] <- "s.lat"

      s.core$Est[tt] <- mean(apply(occ$occ[,,tt], 1, function(x) range.limit(cell.probs = x, prob = 0.75, coord = occ$xy$lat, limit = "south")))
      s.core$LCI[tt] <- quantile(apply(occ$occ[,,tt], 1, function(x) range.limit(cell.probs = x, prob = 0.75, coord = occ$xy$lat, limit = "south")), probs = 0.025)
      s.core$UCI[tt] <- quantile(apply(occ$occ[,,tt], 1, function(x) range.limit(cell.probs = x, prob = 0.75, coord = occ$xy$lat, limit = "south")), probs = 0.975)
      s.core$ind[tt] <- "s.core"

      n.lat$Est[tt] <- mean(apply(occ$occ[,,tt], 1, function(x) range.limit(cell.probs = x, prob = 0.975, coord = occ$xy$lat, limit = "north")))
      n.lat$LCI[tt] <- quantile(apply(occ$occ[,,tt], 1, function(x) range.limit(cell.probs = x, prob = 0.975, coord = occ$xy$lat, limit = "north")), probs = 0.025)
      n.lat$UCI[tt] <- quantile(apply(occ$occ[,,tt], 1, function(x) range.limit(cell.probs = x, prob = 0.975, coord = occ$xy$lat, limit = "north")), probs = 0.975)
      n.lat$ind[tt] <- "n.lat"

      n.core$Est[tt] <- mean(apply(occ$occ[,,tt], 1, function(x) range.limit(cell.probs = x, prob = 0.75, coord = occ$xy$lat, limit = "north")))
      n.core$LCI[tt] <- quantile(apply(occ$occ[,,tt], 1, function(x) range.limit(cell.probs = x, prob = 0.75, coord = occ$xy$lat, limit = "north")), probs = 0.025)
      n.core$UCI[tt] <- quantile(apply(occ$occ[,,tt], 1, function(x) range.limit(cell.probs = x, prob = 0.75, coord = occ$xy$lat, limit = "north")), probs = 0.975)
      n.core$ind[tt] <- "n.core"

      w.lon$Est[tt] <- mean(apply(occ$occ[,,tt], 1, function(x) range.limit(cell.probs = x, prob = 0.975, coord = occ$xy$lon, limit = "west")))
      w.lon$LCI[tt] <- quantile(apply(occ$occ[,,tt], 1, function(x) range.limit(cell.probs = x, prob = 0.975, coord = occ$xy$lon, limit = "west")), probs = 0.025)
      w.lon$UCI[tt] <- quantile(apply(occ$occ[,,tt], 1, function(x) range.limit(cell.probs = x, prob = 0.975, coord = occ$xy$lon, limit = "west")), probs = 0.975)
      w.lon$ind[tt] <- "w.lon"

      w.core$Est[tt] <- mean(apply(occ$occ[,,tt], 1, function(x) range.limit(cell.probs = x, prob = 0.75, coord = occ$xy$lon, limit = "west")))
      w.core$LCI[tt] <- quantile(apply(occ$occ[,,tt], 1, function(x) range.limit(cell.probs = x, prob = 0.75, coord = occ$xy$lon, limit = "west")), probs = 0.025)
      w.core$UCI[tt] <- quantile(apply(occ$occ[,,tt], 1, function(x) range.limit(cell.probs = x, prob = 0.75, coord = occ$xy$lon, limit = "west")), probs = 0.975)
      w.core$ind[tt] <- "w.core"

      e.lon$Est[tt] <- mean(apply(occ$occ[,,tt], 1, function(x) range.limit(cell.probs = x, prob = 0.975, coord = occ$xy$lon, limit = "east")))
      e.lon$LCI[tt] <- quantile(apply(occ$occ[,,tt], 1, function(x) range.limit(cell.probs = x, prob = 0.975, coord = occ$xy$lon, limit = "east")), probs = 0.025)
      e.lon$UCI[tt] <- quantile(apply(occ$occ[,,tt], 1, function(x) range.limit(cell.probs = x, prob = 0.975, coord = occ$xy$lon, limit = "east")), probs = 0.975)
      e.lon$ind[tt] <- "e.lon"

      e.core$Est[tt] <- mean(apply(occ$occ[,,tt], 1, function(x) range.limit(cell.probs = x, prob = 0.75, coord = occ$xy$lon, limit = "east")))
      e.core$LCI[tt] <- quantile(apply(occ$occ[,,tt], 1, function(x) range.limit(cell.probs = x, prob = 0.75, coord = occ$xy$lon, limit = "east")), probs = 0.025)
      e.core$UCI[tt] <- quantile(apply(occ$occ[,,tt], 1, function(x) range.limit(cell.probs = x, prob = 0.75, coord = occ$xy$lon, limit = "east")), probs = 0.975)
      e.core$ind[tt] <- "e.core"

      avg.lat$Est[tt] <- mean(apply(occ$occ[,,tt], 1, function(x) sum(occ$xy$lat * x, na.rm = TRUE)/sum(x, na.rm = TRUE)))
      avg.lat$LCI[tt] <- quantile(apply(occ$occ[,,tt], 1, function(x) sum(occ$xy$lat * x, na.rm = TRUE)/sum(x, na.rm = TRUE)), probs = 0.025)
      avg.lat$UCI[tt] <- quantile(apply(occ$occ[,,tt], 1, function(x) sum(occ$xy$lat * x, na.rm = TRUE)/sum(x, na.rm = TRUE)), probs = 0.975)
      avg.lat$ind[tt] <- "avg.lat"

      avg.lon$Est[tt] <- mean(apply(occ$occ[,,tt], 1, function(x) sum(occ$xy$lon * x, na.rm = TRUE)/sum(x, na.rm = TRUE)))
      avg.lon$LCI[tt] <- quantile(apply(occ$occ[,,tt], 1, function(x) sum(occ$xy$lon * x, na.rm = TRUE)/sum(x, na.rm = TRUE)), probs = 0.025)
      avg.lon$UCI[tt] <- quantile(apply(occ$occ[,,tt], 1, function(x) sum(occ$xy$lon * x, na.rm = TRUE)/sum(x, na.rm = TRUE)), probs = 0.975)
      avg.lon$ind[tt] <- "avg.lon"
    }


    indices <- dplyr::bind_rows(avg.psi, s.lat)
    indices <- dplyr::bind_rows(indices, s.core)
    indices <- dplyr::bind_rows(indices, n.lat)
    indices <- dplyr::bind_rows(indices, n.core)
    indices <- dplyr::bind_rows(indices, w.lon)
    indices <- dplyr::bind_rows(indices, w.core)
    indices <- dplyr::bind_rows(indices, e.lon)
    indices <- dplyr::bind_rows(indices, e.core)
    indices <- dplyr::bind_rows(indices, avg.lat)
    indices <- dplyr::bind_rows(indices, avg.lon)

    write.csv(indices, file = paste0("inst/output/", alpha, "/indices.csv"), row.names = FALSE)
  }
}


#' range.limit
#'
#' Estimate range limit using cumulative probability method


range.limit <- function(cell.probs, prob, coord, limit){
  xy <- data.frame(x = cell.probs, y = coord)
  x <- dplyr::arrange(xy, y)$x/sum(xy$x, na.rm = TRUE)
  y <- dplyr::arrange(xy, y)$y
  y <- y[!is.na(x)]
  x <- x[!is.na(x)]


    if(limit == "south"){
      lim <- predict(smooth.spline(cumsum(x), y, spar=0.1), (1 - prob))$y
    }

    if(limit == "north"){
      lim <- predict(smooth.spline(cumsum(x), y, spar=0.1), prob)$y
    }

  if(limit == "west"){
    lim <- predict(smooth.spline(cumsum(x), y, spar=0.1), (1 - prob))$y
  }

  if(limit == "east"){
    lim <- predict(smooth.spline(cumsum(x), y, spar=0.1), prob)$y
  }

  lim
}
