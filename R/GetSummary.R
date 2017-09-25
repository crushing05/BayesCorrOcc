#' GetSummary
#'
#' Generate summary information for report
#' @export

GetSummary <- function(alpha){
  common <- BayesCorrOcc::code_lookup$common[BayesCorrOcc::code_lookup$alpha == toupper(alpha)]

  ## Data, JAGS object, and occupancy estimates
  occ <- readRDS(paste0('inst/output/', alpha, '/occ.rds'))
  dat <- readRDS(paste0("inst/output/", alpha, "/bbs_data.rds"))
  ppc <- readRDS(paste0("inst/output/", alpha, "/ppc.rds"))

  years <- seq(from = dat$start_year, to = dat$end_year)
  nYears <- dat$nYears

  ## Extract beta coef estimates and CI's
    ## Beta coefficients for psi, gamma, & epsilon
    psi_beta_tab <- MakeBetatab(alpha)

    ## Beta coefficients for th, th0, p, & omega
    p_beta_tab <- MakeBetatab(alpha, nuisance = TRUE)


    ## Indicator variable posterior means
    psi_g_tab <- MakeGtab(alpha)

    summ <- list(spp_name = common,
                 spp_alpha = alpha,
                 n.routes = sum(apply(dat$h, 1, max, na.rm = TRUE)),
                 n.buffer = dim(dat$h)[1] - sum(apply(dat$h, 1, max, na.rm = TRUE)),
                 psi.betas = psi_beta_tab,
                 p.betas = p_beta_tab,
                 g.tab = psi_g_tab,
                 p = ppc$p)

  summ
}

#' MakeBetatab
#'
#' Make Beta table for including in report
#' @param nuisance If true, returns AIC table of theta, theta', p, & omega; if false, returns AIC table for psi, gamma & epsilon


MakeBetatab <- function(alpha, nuisance = FALSE){
  sim_list <- readRDS(paste0("inst/output/", alpha, "/jags_fit.rds"))

  if(!nuisance){
    ## Covert covs included in psi, gam, & eps models to factor, set levels
    covs_use <- factor(c("tmp", "Twet", "Prec", "Pwarm", "dtr",
                          "sq_tmp", "sq_Twet", "sq_Prec", "sq_Pwarm", "sq_dtr"),
                        levels = c("tmp", "Twet", "Prec", "Pwarm", "dtr",
                                   "sq_tmp", "sq_Twet", "sq_Prec", "sq_Pwarm", "sq_dtr"))


    ## Data frame containing beta coeffecients and se
    beta_est <- matrix(NA, nrow = 11, ncol = 3)
    colnames(beta_est) <- c("$\\psi$", "$\\gamma$", "$\\epsilon$")


    ## Fill in intercept values
    beta_est[1, 1] <- paste(trunc(sim_list$mean$b[1]*100)/100,
                            " (", trunc(sim_list$q2.5$b[1]*100)/100,
                            " -- ", trunc(sim_list$q97.5$b[1]*100)/100, ")", sep = "")
    beta_est[1, 2] <- paste(trunc(sim_list$mean$beta.gam0*100)/100,
                            " (", trunc(sim_list$q2.5$beta.gam0*100)/100,
                            " -- ", trunc(sim_list$q97.5$beta.gam0*100)/100, ")", sep = "")
    beta_est[1, 3] <- paste(trunc(sim_list$mean$beta.eps0*100)/100,
                            " (", trunc(sim_list$q2.5$beta.eps0*100)/100,
                            " -- ", trunc(sim_list$q97.5$beta.eps0*100)/100, ")", sep = "")

    ## Fill in coefficients for climate covariates
    for(i in 1:10){
      # Psi
        beta_est[i + 1, 1] <- paste(trunc(sim_list$mean$betaT.psi[i]*100)/100,
                                    " (", trunc(sim_list$q2.5$betaT.psi[i]*100)/100,
                                    " -- ", trunc(sim_list$q97.5$betaT.psi[i]*100)/100, ")", sep = "")


      # Gamma
        beta_est[i + 1, 2] <- paste(trunc(sim_list$mean$betaT.gam[i]*100)/100,
                                    " (", trunc(sim_list$q2.5$betaT.gam[i]*100)/100,
                                    " -- ", trunc(sim_list$q97.5$betaT.gam[i]*100)/100, ")", sep = "")



      # Epsilon
        beta_est[i + 1, 3] <- paste(trunc(sim_list$mean$betaT.eps[i]*100)/100,
                                    " (", trunc(sim_list$q2.5$betaT.eps[i]*100)/100,
                                    " -- ", trunc(sim_list$q97.5$betaT.eps[i]*100)/100, ")", sep = "")


    }


    ## Rename climate covariates
    covs_use <- dplyr::recode(covs_use, tmp = "Temp", sq_tmp = "Temp^2",
                               Twet = "Temp, Wettest Qrt", sq_Twet = "Temp, Wettest Qrt^2",
                               Prec = "Precip", sq_Prec = "Precip^2",
                               Pwarm = "Precip, Warmest Qrt", sq_Pwarm = "Precip, Warmest Qrt^2",
                               dtr = "Diurnal temp range", sq_dtr = "Diurnal temp range^2")

    ## Covert to data frame, add intercept, covert to character, replace NA with "-"
    beta_df <- as.data.frame(beta_est)
    covs <- data.frame(cov = c("Intercept", as.character(levels(covs_use[order(covs_use)]))[1:10]))
    beta_df <- dplyr::bind_cols(covs, beta_df)
    beta_df[, 2:4] <- as.character(unlist(beta_df[, 2:4]))
    names(beta_df)[1] <- ""

  }else{
    ## Data frame containing beta coeffecients and se for theta, theta', p, omega ----
    beta_est <- matrix(NA, nrow = 4, ncol = 3)
    colnames(beta_est) <- c( "$p$", "$\\theta$", "$\\theta'$")


    # Fill in intercept values (if annual == TRUE, different intercept for each year so don't include)
    beta_est[1, 1] <- paste(trunc(sim_list$mean$alpha0*100)/100,
                            " (", trunc(sim_list$q2.5$alpha0*100)/100, " -- ",
                            trunc(sim_list$q97.5$alpha0*100)/100,")", sep = "")

    beta_est[1, 2] <- paste(trunc(sim_list$mean$xpsi[1]*100)/100,
                             " (", trunc(sim_list$q2.5$xpsi[1]*100)/100, " -- ",
                            trunc(sim_list$q97.5$xpsi[1]*100)/100,")", sep = "")
    beta_est[1, 3] <- paste(trunc(sim_list$mean$xpsi[2]*100)/100,
                            " (", trunc(sim_list$q2.5$xpsi[2]*100)/100, " -- ",
                            trunc(sim_list$q97.5$xpsi[2]*100)/100,")", sep = "")



    # Fill in novice observer & wind effects, random observer variance
    beta_est[2, 1] <- paste(trunc(sim_list$mean$alpha1*100)/100,
                            " (", trunc(sim_list$q2.5$alpha1*100)/100, " -- ",
                            trunc(sim_list$q97.5$alpha1*100)/100,")", sep = "")

    beta_est[3, 1] <- paste(trunc(sim_list$mean$alpha2*100)/100,
                            " (", trunc(sim_list$q2.5$alpha2*100)/100, " -- ",
                            trunc(sim_list$q97.5$alpha2*100)/100,")", sep = "")

    beta_est[4, 1] <- paste(trunc(sim_list$mean$sigma.obs*100)/100,
                            " (", trunc(sim_list$q2.5$sigma.obs*100)/100, " -- ",
                            trunc(sim_list$q97.5$sigma.obs*100)/100,")", sep = "")


    ## Covert to data frame, add intercept, covert to character, replace NA with "-"
    beta_df <- as.data.frame(beta_est)
    covs <- data.frame(cov = c("Intercept", "Novice Observer", "Wind", "/sigma^2_{Obs}"))
    beta_df <- dplyr::bind_cols(covs, beta_df)
    beta_df[, 2:4] <- as.character(unlist(beta_df[, 2:4]))
    beta_df[is.na(beta_df)] <- "-"
    names(beta_df)[1] <- ""
  }

  beta_df
}


#' MakeGtab
#'
#' Make indicator variable table for including in report
#' @param alpha 4-letter species code


MakeGtab <- function(alpha){
  sim_list <- readRDS(paste0("inst/output/", alpha, "/jags_fit.rds"))

    ## Covert covs included in psi, gam, & eps models to factor, set levels
    covs_use <- factor(c("tmp", "Twet", "Prec", "Pwarm", "dtr",
                         "sq_tmp", "sq_Twet", "sq_Prec", "sq_Pwarm", "sq_dtr"),
                       levels = c("tmp", "Twet", "Prec", "Pwarm", "dtr",
                                  "sq_tmp", "sq_Twet", "sq_Prec", "sq_Pwarm", "sq_dtr"))


    ## Data frame containing beta coeffecients and se
    beta_est <- matrix(NA, nrow = 10, ncol = 3)
    colnames(beta_est) <- c("$\\psi$", "$\\gamma$", "$\\epsilon$")


    ## Fill in coefficients for climate covariates
    for(i in 1:10){
      # Psi
      beta_est[i, 1] <- trunc(sim_list$mean$g.psi[i]*100)/100


      # Gamma
      beta_est[i, 2] <- trunc(sim_list$mean$g.gam[i]*100)/100


      # Epsilon
      beta_est[i, 3] <- trunc(sim_list$mean$g.eps[i]*100)/100


    }


    ## Rename climate covariates
    covs_use <- dplyr::recode(covs_use, tmp = "Temp", sq_tmp = "Temp^2",
                              Twet = "Temp, Wettest Qrt", sq_Twet = "Temp, Wettest Qrt^2",
                              Prec = "Precip", sq_Prec = "Precip^2",
                              Pwarm = "Precip, Warmest Qrt", sq_Pwarm = "Precip, Warmest Qrt^2",
                              dtr = "Diurnal temp range", sq_dtr = "Diurnal temp range^2")

    ## Covert to data frame, add intercept, covert to character, replace NA with "-"
    beta_df <- as.data.frame(beta_est)
    covs <- data.frame(cov =  as.character(levels(covs_use[order(covs_use)]))[1:10])
    beta_df <- dplyr::bind_cols(covs, beta_df)
    beta_df[, 2:4] <- as.character(unlist(beta_df[, 2:4]))
    names(beta_df)[1] <- ""

beta_df
}
