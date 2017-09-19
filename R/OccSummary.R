#' OccSummary
#'
#' Estimate mean & CI psi for all raster cells
#' @param alpha Vector of alpha codes for species of interest
#' @return Data frame containing annual estimates of the following indices:
#' @export

OccSummary <- function(alpha){
    occ <- readRDS(paste0('inst/output/', alpha, '/occ.rds'))
    dat <- readRDS(paste0("inst/output/", alpha, "/bbs_data.rds"))

    years <- seq(from = dat$start_year, to = dat$end_year)


    Psi <- LCI <- UCI <- NULL
    for(tt in 1:dat$nYears){
      tPsi <- apply(occ$occ[,,tt], 2, mean)
      tLCI <- apply(occ$occ[,,tt], 2, function(x) quantile(x, probs = 0.025))
      tUCI <- apply(occ$occ[,,tt], 2, function(x) quantile(x, probs = 0.975))

      Psi <- c(Psi, tPsi)
      LCI <- c(LCI, tLCI)
      UCI <- c(UCI, tUCI)
    }

    psi <- data.frame(Year = rep(years, each = dim(occ$occ)[2]),
                      Latitude = rep(occ$xy$lat, dat$nYears),
                      Longitude = rep(occ$xy$lon, dat$nYears),
                      Psi = Psi, LCI = LCI, UCI = UCI)

    saveRDS(psi, file = paste0("inst/output/", alpha, "/psi.rds"))
}
