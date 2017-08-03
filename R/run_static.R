#' run_static
#'
#' Run single-season correlated detection model
#' @param sim_data Data object generated from sim_data_static()
#' @param nC Number of chains (default = 3)
#' @param nI Number of interations (default = 25000)
#' @param nB Number of adaption and burn-in iterations (default = 5000)
#' @param nT Thinning rate (default = 10)
#' @param Parallel Should chains be run in parallel (default = TRUE)
#' @export

run_static <- function(sim_data, nC = 3, nI = 25000, nB = 5000, nT = 10, Parallel = TRUE){

    jags.data <- list(h = sim_data$h, nStops = sim_data$nStops, nRoutes = sim_data$nRoutes,
                      stop = sim_data$stop, stop2 = sim_data$stop2)

    jags.params <- c("xpsi", "psi", "alpha0", "alpha1", "alpha2", "p", "z")

    jags.inits <- function(){list(y = sim_data$y, z = apply(sim_data$y, 1, max),
                                  alpha0 = rnorm(1), alpha1 = rnorm(1), alpha2 = rnorm(1),
                                  psi = runif(1, 0, 1), xpsi = runif(2, 0, 1))}

    jags.fit <- jagsUI::jags(data = jags.data, parameters.to.save = jags.params,
                             inits = jags.inits, model.file = "inst/jags/cor_Occ_static.jags",
                             n.chains = nC, n.iter = nI, n.adapt = nB, n.burnin = nB, n.thin = nT,
                             parallel = Parallel)
    return(jags.fit)
}
