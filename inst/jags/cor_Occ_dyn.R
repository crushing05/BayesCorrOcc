sink(file="inst/jags/cor_Occ_dyn.jags")
cat("
    model {

    #### Prior distributions

    ## Priors for linear indicator variables -- prob = 1 if quadratic term in model, 0.5 otherwise
    # for(ii in 1:(nPred * 2)){
    #   g.psi[ii] ~ dbern(0.5)
    #   g.gam[ii] ~ dbern(0.5)
    #   g.eps[ii] ~ dbern(0.5)
    # }


    ## Priors for quadratic indicator variables - prob = 0.5
    # for(ii in (nPred + 1):(nPred * 2)){
    #   g.psi[ii] ~ dbern(g.psi[ii - nPred] * 0.5)
    #   g.gam[ii] ~ dbern(g.psi[ii - nPred] * 0.5)
    #   g.eps[ii] ~ dbern(g.psi[ii - nPred] * 0.5)
    # }

    for(ii in 1:(nPred*2)){
      beta.psi[ii] ~ dnorm(0, 0.01)
      beta.gam[ii] ~ dnorm(0, 0.01)
      beta.eps[ii] ~ dnorm(0, 0.01)
    }

    ## Priors for betas - Normal(0, tau.beta) if g == 1; Normal(mu, 1/se^2) if g == 0
    # for(ii in 1:nPred){
    #   ## Initial occupancy
    #   betaT.psi[ii] ~ dnorm(bpriorm.psi[ii], tprior.psi[ii])T(-10, 10)
    #   beta.psi[ii] <- g.psi[ii] * betaT.psi[ii]
    #
    #   bpriorm.psi[ii] <- (1 - g.psi[ii]) * mu.psi[ii]
    #   tprior.psi[ii] <- g.psi[ii] * 0.01 + (1 - g.psi[ii]) * pow(se.psi[ii], -2)
    #
    #   ## Colonization
    #   betaT.gam[ii] ~ dnorm(bpriorm.gam[ii], tprior.gam[ii])T(-10, 10)
    #   beta.gam[ii] <- g.gam[ii] * betaT.gam[ii]
    #
    #   bpriorm.gam[ii] <- (1 - g.gam[ii]) * mu.gam[ii]
    #   tprior.gam[ii] <- g.gam[ii] * 0.01 + (1 - g.gam[ii]) * pow(se.gam[ii], -2)
    #
    #   ## Extinction
    #   betaT.eps[ii] ~ dnorm(bpriorm.eps[ii], tprior.eps[ii])T(-10, 10)
    #   beta.eps[ii] <- g.eps[ii] * betaT.eps[ii]
    #
    #   bpriorm.eps[ii] <- (1 - g.eps[ii]) * mu.eps[ii]
    #   tprior.eps[ii] <- g.eps[ii] * 0.01 + (1 - g.eps[ii]) * pow(se.eps[ii], -2)
    # }
    #
    # for(ii in (nPred + 1):(nPred * 2)){
    #   ## Initial occupancy
    #   betaT.psi[ii] ~ dnorm(bpriorm.psi[ii], tprior.psi[ii])T(-10, 10)
    #   beta.psi[ii] <- g.psi[ii] * g.psi[ii - nPred] * betaT.psi[ii]
    #
    #   bpriorm.psi[ii] <- (1 - g.psi[ii]) * mu.psi[ii]
    #   tprior.psi[ii] <- g.psi[ii] * 0.01 + (1 - g.psi[ii]) * pow(se.psi[ii], -2)
    #
    #   ## Colonization
    #   betaT.gam[ii] ~ dnorm(bpriorm.gam[ii], tprior.gam[ii])T(-10, 10)
    #   beta.gam[ii] <- g.gam[ii] * g.gam[ii - nPred] * betaT.gam[ii]
    #
    #   bpriorm.gam[ii] <- (1 - g.gam[ii]) * mu.gam[ii]
    #   tprior.gam[ii] <- g.gam[ii] * 0.01 + (1 - g.gam[ii]) * pow(se.gam[ii], -2)
    #
    #   ## Extinction
    #   betaT.eps[ii] ~ dnorm(bpriorm.eps[ii], tprior.eps[ii])T(-10, 10)
    #   beta.eps[ii] <- g.eps[ii] * g.eps[ii - nPred] * betaT.eps[ii]
    #
    #   bpriorm.eps[ii] <- (1 - g.eps[ii]) * mu.eps[ii]
    #   tprior.eps[ii] <- g.eps[ii] * 0.01 + (1 - g.eps[ii]) * pow(se.eps[ii], -2)
    # }
    #
    #
    # beta.gam0 ~ dnorm(0, 0.1)T(-10, 10)
    # beta.eps0 ~ dnorm(0, 0.1)T(-10, 10)


    ## Detection priors
    alpha0 ~ dnorm(0, 0.1)T(-10, 10)
    alpha1 ~ dnorm(0, 0.1)T(-10, 10)
    alpha2 ~ dnorm(0, 0.1)T(-10, 10)

    for(ii in 1:nObs){
      omega[ii] ~ dnorm(0, tau.obs)
    }
    omega[(nObs + 1)] <- 0

    tau.obs <- pow(sigma.obs, -2)
    sigma.obs ~ dunif(0, 10)


    ## Spatial correlation priors
    xpsi[1] ~ dunif(0, 1)
    xpsi[2] ~ dunif(0, 1)


    #### GAM priors from mgcv::jagam()

    ## Parametric effect priors
    for (ii in 1:1) {
      b.psi[ii] ~ dnorm(0,0.033)
      b.gam[ii] ~ dnorm(0,0.033)
      b.eps[ii] ~ dnorm(0,0.033)
    }


    ## prior for s(sim_dat$xy$x,sim_dat$xy$y)...
    K1 <- S1[1:59,1:59] * lambda.psi[1]
    b.psi[2:60] ~ dmnorm(zero[2:60], K1)

    K2 <- S1[1:59,1:59] * lambda.gam[1]
    b.gam[2:60] ~ dmnorm(zero[2:60], K2)

    K3 <- S1[1:59,1:59] * lambda.eps[1]
    b.eps[2:60] ~ dmnorm(zero[2:60], K3)


    ## smoothing parameter priors
    for (ii in 1:1) {
      lambda.psi[ii] ~ dgamma(.05,.005)
      rho.psi[ii] <- log(lambda.psi[ii])

      lambda.gam[ii] ~ dgamma(.05,.005)
      rho.gam[ii] <- log(lambda.gam[ii])

      lambda.eps[ii] ~ dgamma(.05,.005)
      rho.eps[ii] <- log(lambda.eps[ii])
    }


    ## Detection probability
    p[1:nRoutes, 1, 1:nYears] <- 0
    
    for(tt in 1:nYears){
     lp[1:nRoutes, 2, tt] <-  alpha0 + alpha1 %*% Xp[, tt] + alpha2 %*% nov[, tt] + omega[obs[, tt]]
    }


    ## Initial occupancy probability
    lpsi <- X %*% b.psi + Xclim[,,1] %*% beta.psi


    ## Colonization/extinction probability
    for(tt in 2:nYears){
      lgam[1:nRoutes, tt - 1] <- X %*% b.gam + Xclim[,,tt] %*% beta.gam
      leps[1:nRoutes, tt - 1] <- X %*% b.eps + Xclim[,,tt] %*% beta.eps
    }


    ### Likelihood
    for (ii in 1:nRoutes) {
    ## Detection probability
      logit(p[ii, 2, 1]) <- lp[ii, 2, 1]

      ## Initial occupancy
      logit(psi[ii, 1]) <- lpsi[ii]
      z[ii, 1] ~ dbern(psi[ii, 1])


      ##  y: local availability at stop 1 -- 0 = locally unavailable,  1 = locally available
      y[ii, 1, 1] ~ dbern(xpsi[1]/(xpsi[1] + (1 - xpsi[2])))
      h[ii, 1, 1] ~ dbern(p[ii, (z[ii, 1] * y[ii, 1, 1] + 1), 1])


      ## Availability at stops 2-nStops
      for (jj in 2:nStops) {
        y[ii, jj, 1] ~ dbern(xpsi[(z[ii, 1] * y[ii, jj - 1, 1] + 1)])
        h[ii, jj, 1] ~ dbern(p[ii, (z[ii, 1] * y[ii, jj, 1] + 1), 1])
      } # jj


     for(tt in 2:nYears){
      ## Detection probability
      logit(p[ii, 2, tt]) <- lp[ii, 2, tt]

      ## Colonization/extinction probability
      logit(gam[ii, tt - 1]) <- lgam[ii, tt - 1]
      logit(eps[ii, tt - 1]) <- leps[ii, tt - 1]

      ## Occupancy
      psi[ii, tt] <- z[ii, tt - 1] * (1 - eps[ii, tt - 1]) + (1 - z[ii, tt - 1]) * gam[ii, tt - 1]
      z[ii, tt] ~ dbern(psi[ii, tt])


      ##  y: local availability at stop 1 -- 0 = locally unavailable,  1 = locally available
      y[ii, 1, tt] ~ dbern(xpsi[1]/(xpsi[1] + (1 - xpsi[2])))
      h[ii, 1, tt] ~ dbern(p[ii, (z[ii, tt] * y[ii, 1, tt] + 1), tt])


      ## Availability at stops 2-nStops
      for (jj in 2:nStops) {
        y[ii, jj, tt] ~ dbern(xpsi[(z[ii, tt] * y[ii, jj - 1, tt] + 1)])
        h[ii, jj, tt] ~ dbern(p[ii,( z[ii, tt] * y[ii, jj, tt] + 1), tt])
      } # jj
     }
    } # ii

    }
    ", fill=TRUE)
sink()
