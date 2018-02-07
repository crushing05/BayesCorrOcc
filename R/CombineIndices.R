#' CombineIndices
#'
#' Estimate summary indices for single species or group of species
#' @param alpha 4-letter alpha code for single species of interest
#' @param spp Vector containing the 4-letter alpha codes for species group of interest
#' @return Data frame containing the annual indices
#' @export

CombineIndices <- function(alpha = NULL, spp = NULL, group_name = NULL){
  indices <- c("avg.psi", "s.lat", "s.core", "n.lat", "n.core",
          "w.lon", "w.core", "e.lon", "e.core", "avg.lat", "avg.lon")

  if(!is.null(alpha)){
    ind_post <- readRDS(paste0('inst/output/', alpha, '/indices_post.rds'))
    dat <- readRDS(paste0("inst/output/", alpha, "/bbs_data.rds"))
    years <- seq(from = dat$start_year, to = dat$end_year)

    Est <- LCI <- UCI <- NULL
      for(tt in 1:dat$nYears){
        Est.temp <- apply(ind_post[,,tt], 1, mean)
        Est <- c(Est, Est.temp)

        LCI.temp <- apply(ind_post[,,tt], 1, function(x) quantile(x, probs = 0.025))
        LCI <- c(LCI, LCI.temp)

        UCI.temp <- apply(ind_post[,,tt], 1, function(x) quantile(x, probs = 0.975))
        UCI <- c(UCI, UCI.temp)
      }

      indices <- data.frame(ind = c("avg.psi", "s.lat", "s.core", "n.lat", "n.core",
                                    "w.lon", "w.core", "e.lon", "e.core", "avg.lat", "avg.lon"),
                            Year = rep(years, each = 11),
                            Est = Est, LCI = LCI, UCI = UCI)

    write.csv(indices, file = paste0("inst/output/", alpha, "/indices.csv"), row.names = FALSE)
  }

  if(!is.null(spp)){
    spp_test <- dir.exists("inst/output/indices")

    if(!spp_test){
      dir.create("inst/output/indices")
    }

    dat <- readRDS(paste0("inst/output/", spp[1], "/bbs_data.rds"))
    years <- seq(from = dat$start_year, to = dat$end_year)

    for(jj in 1:11){
      for(ii in 1:length(spp)){
        temp <- readRDS(paste0('inst/output/', spp[ii], '/indices_post.rds'))
        temp_df <- data.frame(Species = spp[ii], est = c(temp[jj,,]),
                              Year = rep(years, each = dim(temp)[2]), it = rep(seq(from = 1, to = dim(temp)[2]), dim(temp)[3]))
        temp_df <- dplyr::group_by(temp_df, it)
        temp_df <- dplyr::mutate(temp_df, est2 =  est - est[1] + 1)
        temp_df <- dplyr::ungroup(temp_df)
        if(ii == 1){ind_df <- temp_df}else{ind_df <- suppressWarnings(dplyr::bind_rows(ind_df, temp_df))}
      }
      ind_summ <- dplyr::group_by(ind_df, Year)
      ind_summ <- dplyr::summarise(ind_summ,
                                    Est = mean(est2),
                                    LCI = quantile(est2, probs = 0.025),
                                    UCI = quantile(est2, probs = 0.975))
      ind_summ <- dplyr::mutate(ind_summ, ind = indices[jj])
      if(jj == 1){indices2 <- ind_summ}else{indices2  <- suppressWarnings(dplyr::bind_rows(indices2 , ind_summ))}
    }

    write.csv(indices2, paste0("inst/output/indices/", group_name, ".csv"), row.names = FALSE)
  }
}



