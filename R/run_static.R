#' run_static
#'
#' Run single-season correlated detection model
#' @param sim_data Data object generated from sim_data_static()
#' @param jagam.mod Model object from mgcv::jagam()
#' @param nC Number of chains (default = 3)
#' @param nI Number of interations (default = 25000)
#' @param nB Number of adaption and burn-in iterations (default = 5000)
#' @param nT Thinning rate (default = 10)
#' @param Parallel Should chains be run in parallel (default = TRUE)
#' @export

run_static <- function(sim_data, jagam.mod, nC = 3, nI = 25000, nB = 5000, nT = 10, Parallel = TRUE){

    ## Mu/se of priors when g == 0
    mu.est <- sim_data$beta
    se.est <- rep(0.25, dim(sim_data$X)[2])

    jags.data <- list(h = sim_data$h, nStops = sim_data$nStops, nRoutes = sim_data$nRoutes,
                      stop = sim_data$stop, stop2 = sim_data$stop2, Xp = sim_data$Xp,
                      X1 = sim_data$X, nPred = dim(sim_data$X)[2]/2,
                      X = jagam.mod$jags.data$X, S1 = jagam.mod$jags.data$S1, zero = jagam.mod$jags.data$zero,
                      mu = mu.est, se = se.est)

    jags.params <- c("xpsi", "lambda", "beta", "g", "alpha0", "alpha1", "pind", "psi", "p", "z")

    jags.inits <- function(){list(y = sim_data$y, z = apply(sim_data$y, 1, max),
                                  alpha0 = rnorm(1), alpha1 = rnorm(1), alpha2 = rnorm(1),
                                  betaT = rnorm(dim(sim_data$X)[2]), sigma.beta = runif(1, 0, 5),
                                  pind = runif(1, 0, 1),
                                  xpsi = runif(2, 0, 1), b = jagam.mod$jags.ini$b,
                                  lambda = jagam.mod$jags.ini$lambda)}

    jags.fit <- jagsUI::jags(data = jags.data, parameters.to.save = jags.params,
                             inits = jags.inits, model.file = "inst/jags/cor_Occ_static.jags",
                             n.chains = nC, n.iter = nI, n.adapt = nB, n.burnin = nB, n.thin = nT,
                             parallel = Parallel)
    return(jags.fit)
}
