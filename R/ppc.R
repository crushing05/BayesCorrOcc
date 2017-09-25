#' ppc
#'
#' Perform posterior predictive check on fitted model
#' @param alpha alpha code for species of interest
#' @return Bayesian p-value
#' @export

ppc <- function(alpha){
  ### Read data & fitted model object
  dat <-  readRDS(paste0("inst/output/", alpha, "/bbs_data.rds"))
  mod <- readRDS(paste0("inst/output/", alpha, "/jags_fit.rds"))
  bio <- readRDS(paste0("inst/output/", alpha, "/biovars.rds"))

  ### Get lat/lon splines
  jagam.data <- data.frame(z = rep(1, length(dat$lat)), x = dat$lon, y = dat$lat)
  jagam.mod <- mgcv::jagam(z ~ s(x, y), data = jagam.data, family = "binomial", file = "inst/jags/jagam.jags")

  ### Empty matrices to store observed and simulated X2 estimates
  fit <- matrix(numeric(length = mod$mcmc.info$n.samples*dat$nYears), nrow = dat$nYears)
  fit.new <- matrix(numeric(length = mod$mcmc.info$n.samples*dat$nYears), nrow = dat$nYears)

  ### Store annual psi estimats for each route
  PSI <- array(dim = c(dat$nRoutes, dat$nYears, mod$mcmc.info$n.samples))

  ### Fit in year 1
  ## Observed # of routes w/ each possible detection history
  obs_hist <- BayesCorrOcc::obs_h(dat$h, 1)

  ### For each posterior sample...
  for(i in 1:mod$mcmc.info$n.samples){
    ### Matrix to store simulated detection histories
    sim_h <- matrix(NA, nrow = dat$nRoutes, ncol = 5)

      ### For each route, estimate psi and p
      PSI[, 1, i] <- plogis(jagam.mod$jags.data$X %*% mod$sims.list$b[i,] +
                       bio[,,1] %*% (mod$sims.list$g.psi[i,] * mod$sims.list$betaT.psi[i,]))
      p1 <- plogis(mod$sims.list$alpha0[i] + mod$sims.list$alpha1[i] * dat$wind[, 1] +
                     mod$sims.list$alpha2[i] * dat$nov[, 1] + mod$sims.list$omega[dat$obs[, 1]])

      ### Estimate expected prob for each possible detection history|psi, p, xpsi
      exp_pr <- BayesCorrOcc::y_probs(psi = PSI[obs_hist$no_hist, 1, i], xpsi = mod$sims.list$xpsi[i,], p = p1[obs_hist$no_hist])

      ### Equilibrium stop-level availability
      pi <- mod$sims.list$xpsi[i, 1] / (mod$sims.list$xpsi[i, 1] + 1 - mod$sims.list$xpsi[i, 2])

      ### Simulated availability and detection histories
      sim_y <- matrix(nrow = dat$nRoutes, ncol = 5)
      sim_h <- matrix(nrow = dat$nRoutes, ncol = 5)

      ### Simulate route-level occupancy|psi
      z.new <- rbinom(n = dat$nRoutes, size = 1, prob = PSI[, 1, i])

      ### Simulate stop-level availability|z.new, xpsi
        sim_y[, 1] <- rbinom(n = dat$nRoutes, size = 1, prob = z.new * pi)
        sim_h[, 1] <- rbinom(n = dat$nRoutes, size = 1, prob = sim_y[, 1] * p1)
        for(k in 2:5){
          sim_y[, k] <- rbinom(n = dat$nRoutes, size = 1, prob = z.new * mod$sims.list$xpsi[i, (sim_y[, k - 1] + 1)])
          sim_h[, k] <- rbinom(n = dat$nRoutes, size = 1, prob = sim_y[, k] * p1)
        }

    ### Counts of each possible detection history of simulated histories
    sim_hist <- BayesCorrOcc::obs_h(sim_h[obs_hist$no_hist,], sim = TRUE)

    ### Expected counts of each possible detection history
    exp_h <- apply(exp_pr, 2, sum)

    ### Observed & simulated X2 statistics
    fit[1, i] <- sum((obs_hist$obs_h - exp_h)^2/exp_h)
    fit.new[1, i] <- sum((sim_h2 - exp_h)^2/exp_h)
  }


  ### Fit for years 2 - nYears
  for(t in 2:dat$nYears){
    ## Observed # of routes w/ each possible detection history
    obs_hist <- BayesCorrOcc::obs_h(dat$h, t)


    ### For each posterior sample...
    for(i in 1:mod$mcmc.info$n.samples){
      ### Matrices to store expected and simulated detection histories
      exp_pr <- matrix(NA, nrow = dat$nRoutes, ncol = 32)
      sim_h <- matrix(NA, nrow = dat$nRoutes, ncol = 5)

        ### For each route, estimate gamma, epsilon, psi and p
        gamma <- plogis(mod$sims.list$beta.gam0[i] + bio[,,t] %*% (mod$sims.list$g.gam[i,] * mod$sims.list$betaT.gam[i,]))
        epsilon <- plogis(mod$sims.list$beta.eps[i] + bio[,,t] %*% (mod$sims.list$g.eps[i,] * mod$sims.list$betaT.eps[i,]))

        PSI[, t, i] <- PSI[, t - 1, i] * epsilon + (1 - PSI[, t - 1, i]) * gamma
        p1 <- plogis(mod$sims.list$alpha0[i] + mod$sims.list$alpha1[i] * dat$wind[, t] +
                       mod$sims.list$alpha2[i] * dat$nov[, t] + mod$sims.list$omega[dat$obs[, t]])

        ### Estimate expected prob for each possible detection history|psi, p, xpsi
        exp_pr <- BayesCorrOcc::y_probs(psi = PSI[obs_hist$no_hist, t, i], xpsi = mod$sims.list$xpsi[i,], p = p1[obs_hist$no_hist])

        ### Equilibrium stop-level availability
        pi <- mod$sims.list$xpsi[i, 1] / (mod$sims.list$xpsi[i, 1] + 1 - mod$sims.list$xpsi[i, 2])

        ### Simulated availability and detection histories
        sim_y <- matrix(nrow = dat$nRoutes, ncol = 5)
        sim_h <- matrix(nrow = dat$nRoutes, ncol = 5)

        ### Simulate route-level occupancy|psi
        z.new <- rbinom(n = dat$nRoutes, size = 1, prob = PSI[, t, i])

        ### Simulate stop-level availability|z.new, xpsi
        sim_y[, 1] <- rbinom(n = dat$nRoutes, size = 1, prob = z.new * pi)
        sim_h[, 1] <- rbinom(n = dat$nRoutes, size = 1, prob = sim_y[, 1] * p1)
        for(k in 2:5){
          sim_y[, k] <- rbinom(n = dat$nRoutes, size = 1, prob = z.new * mod$sims.list$xpsi[i, (sim_y[, k - 1] + 1)])
          sim_h[, k] <- rbinom(n = dat$nRoutes, size = 1, prob = sim_y[, k] * p1)
        }


      ### Counts of each possible detection history of simulated histories
      sim_hist <- BayesCorrOcc::obs_h(sim_h[obs_hist$no_hist,], sim = TRUE)

      ### Expected counts of each possible detection history
      exp_h <- apply(exp_pr, 2, sum)

      ### Observed & simulated X2 statistics
      fit[t, i] <- sum((obs_hist$obs_h - exp_h)^2/exp_h)
      fit.new[t, i] <- sum((sim_h2 - exp_h)^2/exp_h)
    }

  }

  ppc <- list(fit = fit, fit.new = fit.new, p = sum(fit > fit.new)/length(fit))
  saveRDS(ppc, file = paste0("inst/output/", dat$alpha, "/ppc.rds"))
}


#' obs_h
#'
#' Count observed number of routes with each detection history
#' @param h detection history matrix
#' @param year Year for with observed counts needed
#' @return List containing observed counts for each detection history & index of NA counts
#' @export


obs_h <- function(h, year, sim = FALSE){
  hists <- c("00000", "00001", "00010", "00011", "00100", "00101", "00110", "00111",
             "01000", "01001", "01010", "01011", "01100", "01101", "01110", "01111",
             "10000", "10001", "10010", "10011", "10100", "10101", "10110", "10111",
             "11000", "11001", "11010", "11011", "11100", "11101", "11110", "11111")

  if(sim){
    x <- as.data.frame(h)
    x2 <- tidyr::unite(x, hist, 1:5, sep = "")
    y <- table(x2)

    non_zero <- hists[which(hists %in% names(y))]

    obs <- integer(length = 32)
    names(obs) <- hists

    for(i in 1:length(non_zero)){
      obs[names(obs) == non_zero[i]] <- y[names(y) == non_zero[i]]
    }
    return(obs)
  }else{
    x <- as.data.frame(h[,,year])
    x2 <- tidyr::unite(x, hist, 1:5, sep = "")
    y <- table(x2)

    non_zero <- hists[which(hists %in% names(y))]

    obs <- integer(length = 32)
    names(obs) <- hists

    for(i in 1:length(non_zero)){
      obs[names(obs) == non_zero[i]] <- y[names(y) == non_zero[i]]
    }

    no_hist <- which(x2 != "NANANANANA")

    obs_H <- list(obs_h = obs, no_hist = no_hist)
    return(obs_H)
  }

}
