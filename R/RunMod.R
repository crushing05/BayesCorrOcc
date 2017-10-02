#' RunMod
#'
#' Run model
#' @param alpha alpha code for species of interest; if NULL, runs model in parallel for all species
#' @export

RunMod <- function(spp = NULL, alpha = NULL, nI = 250, nA = 100, nC = 2, nT = 5, Parallel = FALSE){
  if(!is.null(spp)){
    ## Check if model has already been run for species
    spp_run <- NULL
    for(i in 1:length(spp)){
      spp_test <- file.exists(paste0("inst/output/", spp[i], "/jags_fit.rds"))
      if(!spp_test) spp_run <- c(spp_run, as.character(spp[i]))
    }

    if(is.null(spp_run)){
      return(spp)
    }else{
      ### Register core
      cores <- parallel::detectCores()
      if(length(spp_run) < cores) cores <- length(spp_run)
      doParallel::registerDoParallel(cores = cores)

      ### Run models in parallel
      mods_run <- foreach::foreach(i = 1:length(spp_run), .combine = c,
                                   .packages = c("dplyr", "BayesCorrOcc", "jagsUI")) %dopar%{

                                     ### Read data
                                     dat <-  readRDS(paste0("inst/output/", spp_run[i], "/bbs_data.rds"))
                                     covs <- readRDS(paste0("inst/output/", spp_run[i], "/biovars.rds"))
                                     inits <- readRDS(paste0("inst/output/", spp_run[i], "/inits.rds"))

                                     ### For inits w/ very large (i.e., likely poorly estimated), change to 0 to ensure model doesn't choke
                                     inits$psi.se[abs(inits$psi.betas) > 8] <- 1
                                     inits$psi.betas[abs(inits$psi.betas) > 8] <- 0
                                     inits$gam.se[abs(inits$gam.betas) > 8] <- 1
                                     inits$gam.betas[abs(inits$gam.betas) > 8] <- 0
                                     inits$eps.se[abs(inits$eps.betas) > 8] <- 1
                                     inits$eps.betas[abs(inits$eps.betas) > 8] <- 0
                                     inits$p.betas[abs(inits$p.betas) > 8] <- 0
                                     inits$th.betas[abs(inits$th.betas) > 8] <- 0

                                     ### Get data for GAM JAGS model
                                     jagam.data <- data.frame(z = rep(1, length(dat$lat)), x = dat$lon, y = dat$lat)
                                     jagam.mod <- mgcv::jagam(z ~ s(x, y), data = jagam.data, family = "binomial", file = "inst/jags/jagam.jags")


                                     ### Data for JAGS
                                     jags.data <- list(h = dat$h, nStops = dat$nStops, nRoutes = dat$nRoutes, nYears = dat$nYears,
                                                       Xp = dat$wind, nov = dat$nov, obs = dat$obs, nObs = max(dat$obs) - 1,
                                                       Xclim = covs, nPred = dim(covs)[2]/2,
                                                       X = jagam.mod$jags.data$X, S1 = jagam.mod$jags.data$S1, zero = jagam.mod$jags.data$zero,
                                                       mu.psi = inits$psi.betas, se.psi = inits$psi.se,
                                                       mu.gam = inits$gam.betas, se.gam = inits$gam.se,
                                                       mu.eps = inits$eps.betas, se.eps = inits$eps.se)


                                     ### Parameters to monitor
                                     jags.params <- c("xpsi", "lambda", "betaT.psi", "g.psi", "beta.gam0", "betaT.gam", "g.gam", "beta.eps0", "betaT.eps",
                                                      "g.eps", "b", "alpha0", "alpha1", "alpha2", "sigma.obs", "sigma.beta", "rho", "omega")


                                     ### Initial values
                                     y <- dat$h
                                     y[is.na(y)] <- rbinom(n = length(y[is.na(y)]), size = 1, prob = 0.5)
                                     jags.inits <- function(){list(y = y, z = apply(dat$h, c(1, 3), max),
                                                                   alpha0 = inits$p.betas[1], alpha1 = rnorm(1), alpha2 = rnorm(1),
                                                                   omega = c(rnorm(max(dat$obs) - 1), NA), sigma.obs = runif(1, 0, 5),
                                                                   betaT.psi = inits$psi.betas, sigma.beta = runif(1, 0, 5),
                                                                   betaT.gam = inits$gam.betas,
                                                                   betaT.eps = inits$eps.betas,
                                                                   pind = runif(1, 0, 1),
                                                                   xpsi = plogis(inits$th.betas), b = jagam.mod$jags.ini$b,
                                                                   lambda = jagam.mod$jags.ini$lambda,
                                                                   g.psi = rbinom(dim(covs)[2], size = 1, prob = 0.5),
                                                                   g.gam = rbinom(dim(covs)[2], size = 1, prob = 0.5),
                                                                   g.eps = rbinom(dim(covs)[2], size = 1, prob = 0.5))}


                                     ### Fit model
                                     mod <- system.file("jags", "cor_Occ_dyn.jags", package = "BayesCorrOcc")
                                     jags.fit <- jagsUI::jags(data = jags.data, parameters.to.save = jags.params,
                                                              inits = jags.inits, model.file = mod,
                                                              n.chains = nC, n.iter = nI, n.adapt = nA, n.burnin = nA/2, n.thin = nT,
                                                              parallel = Parallel, verbose = FALSE)

                                     ### Save output
                                     saveRDS(jags.fit, paste0("inst/output/", spp_run[i], "/jags_fit.rds"))
                                     return(spp_run[i])
                                   }
      return(mods_run)
    }
  }

  if(!is.null(alpha)){
    ### Read data
    dat <-  readRDS(paste0("inst/output/", alpha, "/bbs_data.rds"))
    covs <- readRDS(paste0("inst/output/", alpha, "/biovars.rds"))
    inits <- readRDS(paste0("inst/output/", alpha, "/inits.rds"))

    ### For inits w/ very large (i.e., likely poorly estimated), change to 0 to ensure model doesn't choke
    inits$psi.se[abs(inits$psi.betas) > 8] <- 1
    inits$psi.betas[abs(inits$psi.betas) > 8] <- 0
    inits$gam.se[abs(inits$gam.betas) > 8] <- 1
    inits$gam.betas[abs(inits$gam.betas) > 8] <- 0
    inits$eps.se[abs(inits$eps.betas) > 8] <- 1
    inits$eps.betas[abs(inits$eps.betas) > 8] <- 0
    inits$p.betas[abs(inits$p.betas) > 8] <- 0
    inits$th.betas[abs(inits$th.betas) > 8] <- 0

    ### Get data for GAM JAGS model
    jagam.data <- data.frame(z = rep(1, length(dat$lat)), x = dat$lon, y = dat$lat)
    jagam.mod <- mgcv::jagam(z ~ s(x, y), data = jagam.data, family = "binomial", file = "inst/jags/jagam.jags")


    ### Data for JAGS
    jags.data <- list(h = dat$h, nStops = dat$nStops, nRoutes = dat$nRoutes, nYears = dat$nYears,
                      Xp = dat$wind, nov = dat$nov, obs = dat$obs, nObs = max(dat$obs) - 1,
                      Xclim = covs, nPred = dim(covs)[2]/2,
                      X = jagam.mod$jags.data$X, S1 = jagam.mod$jags.data$S1, zero = jagam.mod$jags.data$zero,
                      mu.psi = inits$psi.betas, se.psi = inits$psi.se,
                      mu.gam = inits$gam.betas, se.gam = inits$gam.se,
                      mu.eps = inits$eps.betas, se.eps = inits$eps.se)


    ### Parameters to monitor
    jags.params <- c("xpsi", "lambda", "betaT.psi", "g.psi", "beta.gam0", "betaT.gam", "g.gam", "beta.eps0", "betaT.eps",
                     "g.eps", "b", "alpha0", "alpha1", "alpha2", "sigma.obs", "sigma.beta", "rho", "omega")


    ### Initial values
    y <- dat$h
    y[is.na(y)] <- rbinom(n = length(y[is.na(y)]), size = 1, prob = 0.5)
    jags.inits <- function(){list(y = y, z = apply(dat$h, c(1, 3), max),
                                  alpha0 = inits$p.betas[1], alpha1 = rnorm(1), alpha2 = rnorm(1),
                                  omega = c(rnorm(max(dat$obs) - 1), NA), sigma.obs = runif(1, 0, 5),
                                  betaT.psi = inits$psi.betas, sigma.beta = runif(1, 0, 5),
                                  betaT.gam = inits$gam.betas,
                                  betaT.eps = inits$eps.betas,
                                  pind = runif(1, 0, 1),
                                  xpsi = plogis(inits$th.betas), b = jagam.mod$jags.ini$b,
                                  lambda = jagam.mod$jags.ini$lambda,
                                  g.psi = rbinom(dim(covs)[2], size = 1, prob = 0.5),
                                  g.gam = rbinom(dim(covs)[2], size = 1, prob = 0.5),
                                  g.eps = rbinom(dim(covs)[2], size = 1, prob = 0.5))}


    ### Fit model
    mod <- system.file("jags", "cor_Occ_dyn.jags", package = "BayesCorrOcc")
    jags.fit <- jagsUI::jags(data = jags.data, parameters.to.save = jags.params,
                             inits = jags.inits, model.file = mod,
                             n.chains = nC, n.iter = nI, n.adapt = nA, n.burnin = nA/2, n.thin = nT,
                             parallel = Parallel, verbose = FALSE)

    ### Save output
    saveRDS(jags.fit, paste0("inst/output/", alpha, "/jags_fit.rds"))
  }
}

