#' y_probs
#'
#' Detection history probabilities, given model parameters
#' @param psi Predicted route-level occupancy
#' @param xpsi spatial correlation parameters
#' @param p detection probability
#' @export

y_probs <- function(psi, xpsi, p){
  pi <- xpsi[1] / (xpsi[1] + 1 - xpsi[2])
  pr_y1 <- pi
  pr_y0 <- 1 - pi

  pr_h1 <- pi * p
  pr_h0 <- pi * (1 - p)

  pr <- matrix(NA, nrow = length(psi), ncol = 32)

  #h00000
  pr[,1] <- psi * (pr_y0 * (1 - xpsi[1]) ^ 4 +#00000
                     pr_h0 * (1 - xpsi[2]) * (1 - xpsi[1]) ^ 3 +#10000
                     pr_y0 * xpsi[1] * (1 - p) * (1 - xpsi[2]) * (1 - xpsi[1]) ^ 2 + #01000
                     pr_y0 * (1 - xpsi[1]) * xpsi[1] * (1 - p) * (1 - xpsi[2]) * (1 - xpsi[1]) + #00100
                     pr_y0 * (1 - xpsi[1]) * (1 - xpsi[1]) * xpsi[1] * (1 - p) * (1 - xpsi[2]) + #00010
                     pr_y0 * (1 - xpsi[1]) * (1 - xpsi[1]) * (1 - xpsi[1]) * xpsi[1] * (1 - p) + #00001
                     pr_h0 * xpsi[2] * (1 - p) * (1 - xpsi[2]) * (1 - xpsi[1]) * (1 - xpsi[1]) + #11000
                     pr_h0 * (1 - xpsi[2]) * xpsi[1] * (1 - p) * (1 - xpsi[2]) * (1 - xpsi[1]) + #10100
                     pr_h0 * (1 - xpsi[2]) * (1 - xpsi[1]) * xpsi[1] * (1 - p) * (1 - xpsi[2]) + #10010
                     pr_h0 * (1 - xpsi[2]) * (1 - xpsi[1]) * (1 - xpsi[1]) * xpsi[1] * (1 - p) + #10001
                     pr_y0 * xpsi[1] * (1 - p) * xpsi[2] * (1 - p) * (1 - xpsi[2]) * (1 - xpsi[1]) + #01100
                     pr_y0 * xpsi[1] * (1 - p) * (1 - xpsi[2]) * xpsi[1] * (1 - p) * (1 - xpsi[2]) + #01010
                     pr_y0 * xpsi[1] * (1 - p) * (1 - xpsi[2]) * (1 - xpsi[1]) * xpsi[1] * (1 - p) + #01001
                     pr_y0 * (1 - xpsi[1]) * xpsi[1] * (1 - p) * xpsi[2] * (1 - p) * (1 - xpsi[2]) + #00110
                     pr_y0 * (1 - xpsi[1]) * xpsi[1] * (1 - p) * (1 - xpsi[2]) * xpsi[1] * (1 - p) + #00101
                     pr_y0 * (1 - xpsi[1]) * (1 - xpsi[1]) * xpsi[1] * (1 - p) * xpsi[2] * (1 - p) + #00011
                     pr_y0 * (1 - xpsi[1]) * xpsi[1] * (1 - p) * xpsi[2] * (1 - p) * xpsi[2] * (1 - p) + #00111
                     pr_y0 * xpsi[1] * (1 - p) * (1 - xpsi[2]) * xpsi[1] * (1 - p) * xpsi[2] * (1 - p) + #01011
                     pr_y0 * xpsi[1] * (1 - p) * xpsi[2] * (1 - p) * (1 - xpsi[2]) * xpsi[1] * (1 - p) + #01101
                     pr_y0 * xpsi[1] * (1 - p) * xpsi[2] * (1 - p) * xpsi[2] * (1 - p) * (1 - xpsi[2]) + #01110
                     pr_h0 * (1 - xpsi[2]) * (1 - xpsi[1]) * xpsi[1] * (1 - p) * xpsi[2] * (1 - p) + #10011
                     pr_h0 * (1 - xpsi[2]) * xpsi[1] * (1 - p) * (1 - xpsi[2]) * xpsi[1] * (1 - p) + #10101
                     pr_h0 * (1 - xpsi[2]) * xpsi[1] * (1 - p) * xpsi[2] * (1 - p) * (1 - xpsi[2]) + #10110
                     pr_h0 * xpsi[2] * (1 - p) * (1 - xpsi[2]) * (1 - xpsi[1]) * xpsi[1] * (1 - p) + #11001
                     pr_h0 * xpsi[2] * (1 - p) * (1 - xpsi[2]) * xpsi[1] * (1 - p) * (1 - xpsi[2]) + #11010
                     pr_h0 * xpsi[2] * (1 - p) * xpsi[2] * (1 - p) * (1 - xpsi[2]) * (1 - xpsi[1]) + #11100
                     pr_y0 * xpsi[1] * (1 - p) * (xpsi[2] * (1 - p)) ^ 3 + #01111
                     pr_h0 * (1 - xpsi[2]) * xpsi[1] * (1 - p) * (xpsi[2] * (1 - p)) ^ 2 + #10111
                     pr_h0 * xpsi[2] * (1 - p) * (1 - xpsi[2]) * xpsi[1] * (1 - p) * xpsi[2] * (1 - p) + #11011
                     pr_h0 * xpsi[2] * (1 - p) * xpsi[2] * (1 - p) * (1 - xpsi[2]) * xpsi[1] * (1 - p) + #11101
                     pr_h0 * xpsi[2] * (1 - p) * xpsi[2] * (1 - p) * xpsi[2] * (1 - p) * (1 - xpsi[2]) + #11110
                     pr_h0 * (xpsi[2] * (1 - p)) ^ 4) + (1 - psi)

  #h00001
  pr[,2] <- psi * (pr_y0 * (1 - xpsi[1]) ^ 3 * xpsi[1] + #0000
                    pr_h0 * (1 - xpsi[2]) * (1 - xpsi[1])^2 * xpsi[1] + #1000
                    pr_y0 * xpsi[1] * (1 - p) * (1 - xpsi[2]) * (1 - xpsi[1]) * xpsi[1] + #0100
                    pr_y0 * (1 - xpsi[1]) * xpsi[1] * (1 - p) * (1 - xpsi[2]) * xpsi[1] + #0010
                    pr_y0 * (1 - xpsi[1]) * (1 - xpsi[1]) * xpsi[1] * (1- p) * xpsi[2] + #0001
                    pr_h0 * xpsi[2] * (1 - p) * (1 - xpsi[2]) * (1 - xpsi[1]) * xpsi[1] + #1100
                    pr_h0 * (1 - xpsi[2]) * xpsi[1] * (1 - p) * (1 - xpsi[2]) * xpsi[1] + #1010
                    pr_h0 * (1 - xpsi[2]) * (1 - xpsi[1]) * xpsi[1] * (1 - p) * xpsi[2] + #1001
                    pr_y0 * xpsi[1] * (1 - p) * xpsi[2] * (1 - p) * (1 - xpsi[2]) * xpsi[1] + #0110
                    pr_y0 * xpsi[1] * (1 - p) * (1 - xpsi[2]) * xpsi[1] * (1 - p) * xpsi[2] + #0101
                    pr_y0 * (1 - xpsi[1]) * xpsi[1] * (1 - p) * xpsi[2] * (1 - p) * xpsi[2] + #0011
                    pr_y0 * xpsi[1] * (1 - p) * xpsi[2] * (1 - p) * xpsi[2] * (1 - p) * xpsi[2] + #0111
                    pr_h0 * (1 - xpsi[2]) * xpsi[1] * (1 - p) * xpsi[2] * (1 - p) * xpsi[2] + #1011
                    pr_h0 * xpsi[2] * (1 - p) * (1 - xpsi[2]) * xpsi[1] * (1 - p) * xpsi[2] + #1101
                    pr_h0 * xpsi[2] * (1 - p) * xpsi[2] * (1 - p) * (1 - xpsi[2]) * xpsi[1] + #1110
                    pr_h0 * (xpsi[2] * (1 - p)) ^ 3 * xpsi[2]) * p

  #h00010
  pr[,3] <- psi * (pr_y0 * (1 - xpsi[1]) * (1 - xpsi[1]) * xpsi[1] + #000
                    pr_y0 * (1 - xpsi[1]) * xpsi[1] * (1 - p) * xpsi[2] + #001
                    pr_h0 * (1 - xpsi[2]) * (1 - xpsi[1]) * xpsi[1] + #100
                    pr_y0 * xpsi[1] * (1 - p) * (1 - xpsi[2]) * xpsi[1] + #010
                    pr_h0 * xpsi[2] * (1 - p) * (1 - xpsi[2]) * xpsi[1] + #110
                    pr_h0 * (1 - xpsi[2]) * xpsi[1] * (1 - p) * xpsi[2] + #101
                    pr_y0 * xpsi[1] * (1 - p) * xpsi[2] * (1 - p) * xpsi[2] + #011
                    pr_h0 * xpsi[2] * (1 - p) * xpsi[2] * (1 - p) * xpsi[2]) * p *  ((1 - xpsi[2]) + xpsi[2] * (1 - p))

  #h00011
  pr[,4] <- psi * (pr_y0 * (1 - xpsi[1]) * (1 - xpsi[1]) * xpsi[1] + #000
                     pr_y0 * (1 - xpsi[1]) * xpsi[1] * (1 - p) * xpsi[2] + #001
                     pr_h0 * (1 - xpsi[2]) * (1 - xpsi[1]) * xpsi[1] + #100
                     pr_y0 * xpsi[1] * (1 - p) * (1 - xpsi[2]) * xpsi[1] + #010
                     pr_h0 * xpsi[2] * (1 - p) * (1 - xpsi[2]) * xpsi[1] + #110
                     pr_h0 * (1 - xpsi[2]) * xpsi[1] * (1 - p) * xpsi[2] + #101
                     pr_y0 * xpsi[1] * (1 - p) * xpsi[2] * (1 - p) * xpsi[2] + #011
                     pr_h0 * xpsi[2] * (1 - p) * xpsi[2] * (1 - p) * xpsi[2]) * p * xpsi[2] * p

  #h00100
  pr[,5] <- psi * (pr_y0 * (1 - xpsi[1]) * xpsi[1] + #00
                    pr_y0 * xpsi[1] * (1 - p) * xpsi[2] + #01
                    pr_h0 * (1 - xpsi[2]) * xpsi[1] +#10
                    pr_h0 * xpsi[2] * (1 - p) * xpsi[2]) * p * ((1 - xpsi[2]) * (1 - xpsi[1]) +
                                                                (1 - xpsi[2]) * xpsi[1] * (1 - p) +
                                                                xpsi[2] * (1 - p) * (1 - xpsi[2]) +
                                                                (xpsi[2] * (1 - p))^2)
  #h00101
  pr[,6] <- psi * (pr_y0 * (1 - xpsi[1]) * xpsi[1] + #00
                     pr_y0 * xpsi[1] * (1 - p) * xpsi[2] + #01
                     pr_h0 * (1 - xpsi[2]) * xpsi[1] +#10
                     pr_h0 * xpsi[2] * (1 - p) * xpsi[2]) * p * ((1 - xpsi[2]) * xpsi[1] + xpsi[2] * (1 - p) * xpsi[2]) * p

  #h00110
  pr[,7] <- psi * (pr_y0 * (1 - xpsi[1]) * xpsi[1] + #00
                     pr_y0 * xpsi[1] * (1 - p) * xpsi[2] + #01
                     pr_h0 * (1 - xpsi[2]) * xpsi[1] +#10
                     pr_h0 * xpsi[2] * (1 - p) * xpsi[2]) * p * xpsi[2] * p * ((1 - xpsi[2]) + xpsi[2] * (1 - p))

  #h00111
  pr[,8] <- psi * (pr_y0 * (1 - xpsi[1]) * xpsi[1] + #00
                     pr_y0 * xpsi[1] * (1 - p) * xpsi[2] + #01
                     pr_h0 * (1 - xpsi[2]) * xpsi[1] +#10
                     pr_h0 * xpsi[2] * (1 - p) * xpsi[2]) * p * xpsi[2] * p * xpsi[2] * p

  #h01000
  pr[,9] <- psi * (pr_y0 * xpsi[1] +
                    pr_h0 * xpsi[2]) * p * ((1 - xpsi[2]) * (1 - xpsi[1]) ^ 2 + #000
                                             xpsi[2] * (1 - p) * (1 - xpsi[2]) * (1 - xpsi[1]) + #100
                                             (1 - xpsi[2]) * xpsi[1] * (1 - p) * (1 - xpsi[2]) + #010
                                             (1 - xpsi[2]) * (1 - xpsi[1]) * xpsi[1] * (1 - p) + #001
                                             xpsi[2] * (1 - p) * xpsi[2] * (1 - p) * (1 - xpsi[2]) + #110
                                             xpsi[2] * (1 - p) * (1 - xpsi[2]) * xpsi[1] * (1 - p) + #101
                                             (1 - xpsi[2]) * xpsi[1] * (1 - p) * xpsi[2] * (1 - p) + #011
                                             (xpsi[2] * (1 - p))^3)

  #h01001
  pr[,10] <- psi * (pr_y0 * xpsi[1] +
                     pr_h0 * xpsi[2]) * p * ((1 - xpsi[2]) * (1 - xpsi[1]) * xpsi[1] +
                                              (1 - xpsi[2]) * xpsi[1] * (1 - p) * xpsi[2] +
                                              xpsi[2] * (1 - p) * (1 - xpsi[2]) * xpsi[1] +
                                              xpsi[2] * (1 - p) * xpsi[2] * (1 - p) * xpsi[2]) * p

  #h01010
  pr[,11] <- psi * (pr_y0 * xpsi[1] +
                     pr_h0 * xpsi[2]) * p * ((1 - xpsi[2]) * xpsi[1] +
                                              xpsi[2] * (1 - p) * xpsi[2]) * p * ((1 - xpsi[2]) + xpsi[2] * (1 - p))


  #h01011
  pr[,12] <- psi * (pr_y0 * xpsi[1] +
                     pr_h0 * xpsi[2]) * p * ((1 - xpsi[2]) * xpsi[1] +
                                              xpsi[2] * (1 - p) * xpsi[2]) * p * xpsi[2] * p

  #h01100
  pr[,13] <- psi * (pr_y0 * xpsi[1] +
                     pr_h0 * xpsi[2]) * p * xpsi[2] * p * ((1 - xpsi[2]) * (1 - xpsi[1]) + #00
                                                           (1 - xpsi[2]) * xpsi[1] * (1 - p) + #01
                                                           xpsi[2] * (1 - p) * (1 - xpsi[2]) + #10
                                                           (xpsi[2] * (1 - p)) ^ 2) #11
  #h01101
  pr[,14] <- psi * (pr_y0 * xpsi[1] +
                     pr_h0 * xpsi[2]) * p * xpsi[2] * p  * ((1 - xpsi[2]) * xpsi[1] + xpsi[2] * (1 - p) * xpsi[2]) * p

  #h01110
  pr[,15] <- psi * (pr_y0 * xpsi[1] +
                     pr_h0 * xpsi[2]) * p * xpsi[2] * p * xpsi[2] * p * ((1 - xpsi[2]) + xpsi[2] * (1 - p))

  #h01111
  pr[,16] <- psi * (pr_y0 * xpsi[1] +
                     pr_h0 * xpsi[2]) * p * (xpsi[2] * p) ^ 3

  #h10000
  pr[,17] <- psi * pr_h1 * ((1 - xpsi[2]) * (1 - xpsi[1]) ^ 3 + #0000
                            xpsi[2] * (1 - p) * (1 - xpsi[2]) * (1 - xpsi[1])^2 + #1000
                            (1 - xpsi[2]) * xpsi[1] * (1 - p) * (1 - xpsi[2]) * (1 - xpsi[1]) + #0100
                            (1 - xpsi[2]) * (1 - xpsi[1]) * xpsi[1] * (1 - p) * (1 - xpsi[2]) + #0010
                            (1 - xpsi[2]) * (1 - xpsi[1]) * (1 - xpsi[1]) * xpsi[1] * (1- p) + #0001
                            xpsi[2] * (1 - p) * xpsi[2] * (1 - p) * (1 - xpsi[2]) * (1 - xpsi[1]) + #1100
                            xpsi[2] * (1 - p) * (1 - xpsi[2]) * xpsi[1] * (1 - p) * (1 - xpsi[2]) + #1010
                            xpsi[2] * (1 - p) * (1 - xpsi[2]) * (1 - xpsi[1]) * xpsi[1] * (1 - p) + #1001
                            (1 - xpsi[2]) * xpsi[1] * (1 - p) * xpsi[2] * (1 - p) * (1 - xpsi[2]) + #0110
                            (1 - xpsi[2]) * xpsi[1] * (1 - p) * (1 - xpsi[2]) * xpsi[1] * (1 - p) + #0101
                            (1 - xpsi[2]) * (1 - xpsi[1]) * xpsi[1] * (1 - p) * xpsi[2] * (1 - p) + #0011
                            (1 - xpsi[2]) * xpsi[1] * (1 - p) * xpsi[2] * (1 - p) * xpsi[2] * (1 - p) + #0111
                            xpsi[2] * (1 - p) * (1 - xpsi[2]) * xpsi[1] * (1 - p) * xpsi[2] * (1 - p) + #1011
                            xpsi[2] * (1 - p) * xpsi[2] * (1 - p) * (1 - xpsi[2]) * xpsi[1] * (1 - p) + #1101
                            xpsi[2] * (1 - p) * xpsi[2] * (1 - p) * xpsi[2] * (1 - p) * (1 - xpsi[2]) + #1110
                            (xpsi[2] * (1 - p)) ^ 4)

  #h10001
  pr[,18] <- psi * pr_h1 * ((1 - xpsi[2]) * (1 - xpsi[1]) ^ 2 * xpsi[1] + #0001
                             (1 - xpsi[2])* (1 - xpsi[1]) * xpsi[1] * (1 - p) * xpsi[2] + #0011
                             (1 - xpsi[2]) * xpsi[1] * (1 - p) * (1 - xpsi[2]) * xpsi[1] + #0101
                             xpsi[2] * (1 - p) * (1 - xpsi[2]) * (1 - xpsi[1]) * xpsi[1] + #1001
                             (1 - xpsi[2]) * xpsi[1] * (1 - p) * xpsi[2] * (1 - p) * xpsi[2] + #0111
                             xpsi[2] * (1 - p) * (1 - xpsi[2]) * xpsi[1] * (1 - p) * xpsi[2] + #1011
                             xpsi[2] * (1 - p) * xpsi[2] * (1 - p) * (1 - xpsi[2]) * xpsi[1] + #1101
                             xpsi[2] * (1 - p) * xpsi[2] * (1 - p) * xpsi[2] * (1 - p) * xpsi[2]) * p

  #h10010
  pr[,19] <- psi * pr_h1 * ((1 - xpsi[2]) * (1 - xpsi[1]) * xpsi[1] +
                            (1 - xpsi[2]) * xpsi[1] * (1 - p) * xpsi[2] +
                            xpsi[2] * (1 - p) * (1 - xpsi[2]) * xpsi[1] +
                            xpsi[2] * (1 - p) * xpsi[2] * (1 - p) * xpsi[2]) * p * ((1 - xpsi[2]) + xpsi[2] * (1 - p))

  #h10011
  pr[,20] <- psi * pr_h1 * ((1 - xpsi[2]) * (1 - xpsi[1]) * xpsi[1] +
                             (1 - xpsi[2]) * xpsi[1] * (1 - p) * xpsi[2] +
                             xpsi[2] * (1 - p) * (1 - xpsi[2]) * xpsi[1] +
                             xpsi[2] * (1 - p) * xpsi[2] * (1 - p) * xpsi[2]) * p * xpsi[2] * p

  #h10100
  pr[,21] <- psi * pr_h1 * ((1 - xpsi[2]) * xpsi[1] +
                            xpsi[2] * (1 - p) * xpsi[2]) * p * ((1 - xpsi[2]) * (1 - xpsi[1]) + #00
                                                                (1 - xpsi[2]) * xpsi[1] * (1 - p) + #01
                                                                xpsi[2] * (1 - p) * (1 - xpsi[2]) + #10
                                                                (xpsi[2] * (1 - p)) ^ 2) #11
  #h10101
  pr[,22] <- psi * pr_h1 * ((1 - xpsi[2]) * xpsi[1] +
                             xpsi[2] * (1 - p) * xpsi[2]) * p * ((1 - xpsi[2]) * xpsi[1] +
                                                                 xpsi[2] * (1 - p) * xpsi[2]) * p

  #h10110
  pr[,23] <- psi * pr_h1 * ((1 - xpsi[2]) * xpsi[1] +
                             xpsi[2] * (1 - p) * xpsi[2]) * p * xpsi[2] * p * ((1 - xpsi[2]) + xpsi[2] * (1 - p))

  #h10111
  pr[,24] <- psi * pr_h1 * ((1 - xpsi[2]) * xpsi[1] +
                             xpsi[2] * (1 - p) * xpsi[2]) * p * xpsi[2] * p * xpsi[2] * p

  #h11000
  pr[,25] <- psi * pr_h1 * xpsi[2] * p * ((1 - xpsi[2]) * (1 - xpsi[1]) ^ 2 + #000
                                         xpsi[2] * (1 - p) * (1 - xpsi[2]) * (1 - xpsi[1]) + #100
                                         (1 - xpsi[2]) * xpsi[1] * (1 - p) * (1 - xpsi[2]) + #010
                                         (1 - xpsi[2]) * (1 - xpsi[1]) * xpsi[1] * (1 - p) + #001
                                         xpsi[2] * (1 - p) * xpsi[2] * (1 - p) * (1 - xpsi[2]) + #110
                                         xpsi[2] * (1 - p) * (1 - xpsi[2]) * xpsi[1] * (1 - p) + #101
                                         (1 - xpsi[2]) * xpsi[1] * (1 - p) * xpsi[2] * (1 - p) + #011
                                         (xpsi[2] * (1 - p))^3)
  #h11001
  pr[,26] <- psi * pr_h1 * xpsi[2] * p * ((1 - xpsi[2]) * (1 - xpsi[1]) * xpsi[1] +
                                          (1 - xpsi[2]) * xpsi[1] * (1 - p) * xpsi[2] +
                                          xpsi[2] * (1 - p) * (1 - xpsi[2]) * xpsi[1] +
                                          xpsi[2] * (1 - p) * xpsi[2] * (1 - p) * xpsi[2]) * p
  #h11010
  pr[,27] <- psi * pr_h1 * xpsi[2] * p * ((1 - xpsi[2]) * xpsi[1] +
                                          xpsi[2] * (1 - p) * xpsi[2]) * p * ((1 - xpsi[2]) + xpsi[2] * (1 - p))

  #h11011
  pr[,28] <- psi * pr_h1 * xpsi[2] * p * ((1 - xpsi[2]) * xpsi[1] +
                                          xpsi[2] * (1 - p) * xpsi[2]) * p * xpsi[2] * p

  #h11100
  pr[,29] <- psi * pr_h1 * xpsi[2] * p * xpsi[2] * p * ((1 - xpsi[2]) * (1 - xpsi[1]) + #00
                                                       (1 - xpsi[2]) * xpsi[1] * (1 - p) + #01
                                                       xpsi[2] * (1 - p) * (1 - xpsi[2]) + #10
                                                       (xpsi[2] * (1 - p)) ^ 2) #11

  #h11101
  pr[,30] <- psi * pr_h1 * xpsi[2] * p * xpsi[2] * p * ((1 - xpsi[2]) * xpsi[1] +
                                                       xpsi[2] * (1 - p) * xpsi[2]) * p

  #h11110
  pr[,31] <- psi * pr_h1 * xpsi[2] * p * xpsi[2] * p * xpsi[2] * p * ((1 - xpsi[2]) + xpsi[2] * (1 - p))

  #h11111
  pr[,32] <- psi * pr_h1 * (xpsi[2] * p) ^ 4

  return(pr)
}

