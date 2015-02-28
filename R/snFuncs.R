



dSN <- function(y, mu = 0, sigma2 = 1, shape=1){
  dens <- 2*dnorm(y, mu, sqrt(sigma2))*pnorm(shape*((y - mu)/sqrt(sigma2)))
  return(dens)
}

d.mixedSN <- function(x, pi1, mu, sigma2, shape){
    g <- length(pi1)
    dens <- 0
    for (j in 1:g) dens <- dens + pi1[j]*dSN(x, mu[j], sigma2[j], shape[j])
    return(dens)
  }

SNmean <- function(mu, sigma2, shape){

   xi <- mu
   omega <- sqrt(sigma2)
   alpha <- shape

   C <- sqrt(2 / pi)
   delta <- alpha / sqrt(1 + alpha^2)

   mean <- xi + omega * delta * C
   mean
}

SNsd <- function(mu, sigma2, shape){

	xi <- mu
	omega2 <- sigma2
	alpha <- shape

	delta <- alpha / sqrt(1 + alpha^2)

	c <- (2 * delta^2) / pi

	variance <- omega2 * (1 - c)

	sd <- sqrt(variance)

}
