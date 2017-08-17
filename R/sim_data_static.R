#' sim_data_static
#'
#' Simulation of single-season count data for the static correlated-detection N-mixture model
#' @param nRoutes Number of sample units
#' @param nStops Number of stops within each sample unit
#' @param alpha0 Intercept for detection model
#' @param alpha1 Linear effect of Xp on detection probability
#' @param sigma1 SD of Xp
#' @param sigma2 SD of X
#' @param sigma3 Noise
#' @param theta Vector containing the correlation terms
#' @param nSum Optional number of stops to aggregate
#' @export

sim_data_static <- function(nRoutes = 150, nStops = 50, nPred = 5,
                             alpha0 = -0.5, alpha1 = 0.75, sigma1 = 0.5,
                             sigma2 = 0.5, sigma3 = 0.05,
                             theta = c(0.3, 0.75), sx = 0.2, sy = 0.2,
                             nSum = NULL){

  ## Lat/lon
  xy <- data.frame(x = runif(nRoutes, -100, -70),
                   y = runif(nRoutes, 30, 50))

  sxy <- data.frame(x = (xy$x - min(xy$x))/(max(xy$x) - min(xy$x)),
                    y = (xy$y - min(xy$y))/(max(xy$y) - min(xy$y)))
  sigma <- cov(sxy)


  ## Standardize stop number
  stop <- scale(seq(1:nStops))[,1]
  stop2 <- stop^2


  ## Stop-level detection probability w/ 1 covariate
  Xp <- rnorm(nRoutes, sd = sigma1)
  lp <- alpha0 + alpha1 * Xp
  p <- exp(lp)/(1 + exp(lp))


  ## Route-level occupancy probability w/ lat/lon spline & nPred covariates (linear + quadratic effects)
  beta1 <- rgamma(nPred, 1, 1)    # Linear slopes
  beta2 <- -rgamma(nPred, 1, 1)   # Quadratic slopes
  beta <- c(beta1, beta2)

  X1 <- matrix(rnorm(nRoutes * nPred, sd = sigma2), nrow = nRoutes, ncol = nPred) # Covariate values
  X2 <- X1^2
  X <- cbind(X1, X2)

  eta <- rnorm(nRoutes, sd = sigma3)   # Noise

  f <- (1/(2*pi*sx*sy*sqrt(1 - sigma[1,2]))) *
        exp(-(1/(2*(1-sigma[1,2]^2))*
        (sxy$x-0.5)^2/sx^2 + (sxy$y - 0.5)^2/sy^2 - (2*sigma[1,2]*(sxy$x - 0.5)*(sxy$y - 0.5))/sx*sy))

  l.psi <- f + X %*% beta + eta
  psi <- exp(l.psi)/(1 + exp(l.psi))


  # Route-level occupancy
  z <- rbinom(n = nRoutes, size = 1, prob = psi)


  ## Equilibrium proportion of available sites
  theta1 <- theta[1] / (theta[1] + (1 - theta[2]))

  ## Stop level data
  y <- h <- matrix(NA, nrow = nRoutes, ncol = nStops)  # Stop availability & observed occupancy

  for(i in 1:nRoutes){
    y[i, 1] <- rbinom(n = 1, size = 1, prob = z[i] * theta1)
    h[i, 1] <- rbinom(n = 1, size = 1, prob = y[i, 1] * p[i])
    for(j in 2:nStops){
      y[i, j] <- rbinom(n = 1, size = 1, prob = z[i] * theta[y[i, j - 1] + 1])
      h[i, j] <- rbinom(n = 1, size = 1, prob = y[i, j] * p[i])
    }
  }



  if(!is.null(nSum)){
    ## Sum stop-level counts
    for(i in 1:nRoutes){
      yt <- unname(tapply(y[i,], (seq_along(y[i,])-1) %/% nSum, max))
      ht <- unname(tapply(h[i,], (seq_along(h[i,])-1) %/% nSum, max))
      if(i == 1){
        y2 <- yt
        h2 <- ht
      }else{
        y2 <- rbind(y2, yt)
        h2 <- rbind(h2, ht)
      }
    }

    ## New number of stops
    nStops <- dim(h2)[2]

    ## New scaled stop number
    stop <- scale(seq(1:nStops))[,1]
    stop2 <- stop^2

    ## New stop-level availability
    y <- y2
    h <- h2
  }


  sim_data <- list(psi = psi, z = z, h = h, y = y, p = p, Xp = Xp, X = X, xy = xy, nRoutes = nRoutes, nStops = nStops,
                   beta0 = beta0, beta = beta, sigma1 = sigma1, sigma2 = sigma2, sigma3 = sigma3, alpha0 = alpha0, alpha1 = alpha1, #alpha2 = alpha2,
                   stop = stop, stop2 = stop2, theta = theta)

  return(sim_data)
}
