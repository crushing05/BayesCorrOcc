sink(file="inst/jags/cor_Occ_static.jags")
cat("
    model {

    #### Prior distributions
    for(ii in 1:nPred){
      g[ii] ~ dbern(pind)
      betaT[ii] ~ dnorm(0, 0.01)T(-10, 10)
      beta[ii] <- g[ii] * betaT[ii]
    }

    tau.beta <- pow(sigma.beta, -2)
    sigma.beta ~ dunif(0, 10)
    pind ~ dbeta(2, 8)

    alpha0 ~ dnorm(0, 0.1)T(-10, 10)
    alpha1 ~ dnorm(0, 0.1)T(-10, 10)
    alpha2 ~ dnorm(0, 0.1)T(-10, 10)

    xpsi[1] ~ dunif(0, 1)
    xpsi[2] ~ dunif(0, 1)


    ## Parametric effect priors
    for (ii in 1:1) { b[ii] ~ dnorm(0,0.033) }


    ## prior for s(sim_dat$xy$x,sim_dat$xy$y)...
    K1 <- S1[1:29,1:29] * lambda[1]  + S1[1:29,30:58] * lambda[2]
    b[2:30] ~ dmnorm(zero[2:30],K1)


    ## smoothing parameter priors
    for (ii in 1:2) {
      lambda[ii] ~ dgamma(.05,.005)
      rho[ii] <- log(lambda[ii])
    }

    eta <- X %*% b + X1 %*% beta

    for (ii in 1:nRoutes) {
      ## Detection probability
      p[ii, 1] <- 0
      logit(p[ii, 2]) <- alpha0 + alpha1 * Xp[ii]


      ## Route-level occupancy probability
      logit(psi[ii]) <- eta[ii]


      ## Route-level occupancy
      z[ii] ~ dbern(psi[ii])


      ##  y: local availability at stop 1 -- 0 = locally unavailable,  1 = locally available
      y[ii, 1] ~ dbern(xpsi[1]/(xpsi[1] + (1 - xpsi[2])))
      h[ii, 1] ~ dbern(p[ii, z[ii] * y[ii, 1] + 1])


      ## Availability at stops 2-nStops
      for (jj in 2:nStops) {
        y[ii, jj] ~ dbern(xpsi[(z[ii] * y[ii, jj - 1] + 1)])
        h[ii, jj] ~ dbern(p[ii, z[ii] * y[ii, jj] + 1])
      } # jj
    } # ii

    }
    ", fill=TRUE)
sink()
