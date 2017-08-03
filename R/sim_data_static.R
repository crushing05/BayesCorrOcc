#' sim_data_static
#'
#' Simulation of single-season count data for the static correlated-detection N-mixture model
#' @param nRoutes Number of sample units
#' @param nStops Number of stops within each sample unit
#' @param alpha0 Intercept for detection model
#' @param alpha1 Linear effect of stop number on detection probability
#' @param alpha2 Quadratic effect of stop number on detection probability
#' @param beta0 Log expected number of individuals at each sampling unit
#' @param sigma Standard deviation around mean count
#' @param theta Vector containing the correlation terms
#' @param nSum Optional number of stops to aggregate
#' @export

sim_data_static <- function(nRoutes = 150, nStops = 50,
                             alpha0 = -0.5, alpha1 = -0.3, alpha2 = 0.2,
                             beta0 = 2, beta1 = 1.5, sigma = 0.5,
                             theta = c(0.3, 0.75),
                             nSum = NULL){


  ## Standardize stop number
  stop <- scale(seq(1:nStops))[,1]
  stop2 <- stop^2

  ## Stop detection probability
  lp <- alpha0 + alpha1 * stop + alpha2 * stop2
  p <- exp(lp)/(1 + exp(lp))

  ## Route-level abundance
  eta <- rnorm(nRoutes, 0, sigma)
  X <- rnorm(nRoutes)
  l.psi <- beta0 + beta1 * X + eta
  psi <- exp(l.psi)/(1 + exp(l.psi))

  # Route-level occupancy
  z <- rbinom(n = nRoutes, size = 1, prob = psi)

  ## Equilibrium proportion of available sites
  theta1 <- theta[1] / (theta[1] + (1 - theta[2]))

  ## Stop level data
  y <- matrix(NA, nrow = nRoutes, ncol = nStops)  # Stop availability
  h <- matrix(NA, nrow = nRoutes, ncol = nStops)  # Observed presence/absence

  for(i in 1:nRoutes){
    y[i, 1] <- rbinom(n = 1, size = 1, prob = z[i] * theta1)
    h[i, 1] <- rbinom(n = 1, size = 1, prob = y[i, 1] * p[1])
    for(j in 2:nStops){
      y[i, j] <- rbinom(n = 1, size = 1, prob = z[i] * theta[y[i, j - 1] + 1])
      h[i, j] <- rbinom(n = 1, size = 1, prob = y[i, j] * p[j])
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

  sim_data <- list(h = h, y = y, p = p, X = X, nRoutes = nRoutes, nStops = nStops,
                   beta0 = beta0, sigma = sigma, alpha0 = alpha0, alpha1 = alpha1, alpha2 = alpha2,
                   stop = stop, stop2 = stop2, theta = theta)

  return(sim_data)
}
