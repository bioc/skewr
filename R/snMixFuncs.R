############################################################
##
## All-in-one function that allows for subsetting
## of probes and skew normal modeling.
## Returns a Skew.normal object as produced by
## the smsn.mix function from the mixsmsn package, but which
## includes the means, modes, and densities of the
## individual components
##
############################################################



getSNparams <- function(MethyLumiSet, allele = c('M', 'U'),
                        type = c('I-red', 'I-green', 'II'), snps = TRUE, idmr = TRUE, ch = FALSE){

	object <- subsetProbes(MethyLumiSet, allele = allele, type = type, snps = snps, idmr = idmr, ch = ch)
	object <- prepSamp(object)

	if (type == 'II') g <- 2
	else g <- 3

	x <- seq(0, 18, length = 1000000)
	dens <- matrix(NA_real_, ncol = g, nrow = length(x))
	dens.list <- vector("list", g)
	means <- vector("numeric", g)
	modes <- vector("numeric", g)
	dens.max <- vector("numeric", g)

	if(is.matrix(object)){
		n.samp <- ncol(object)
		mix <- apply(object, 2, function(xx){
  											xx <- na.omit(xx)
												smsn.mix(xx, nu = 3, g = g, obs.prob = TRUE, calc.im = FALSE)
												})

		for(i in 1:ncol(object)){
			pii <- mix[[i]]$pii
			mu <- mix[[i]]$mu
			sigma2 <- mix[[i]]$sigma2
			shape <- mix[[i]]$shape
			for(j in 1:g){
				dens[,j] <- dSN(x, mu[j], sigma2[j], shape[j])
				dens.list[[j]] <- data.frame(x = x, y = dens[,j])
				means[j] <- SNmean(mu[j], sigma2[j], shape[j])
				dens.max[j] <- which.max(dens.list[[j]]$y)
				modes[j] <- dens.list[[j]]$x[dens.max[j]]
			}
			mix[[i]]$means <- means
			mix[[i]]$modes <- modes
			mix[[i]]$dens.list <- dens.list
		}
	}
	else{
		mix <- smsn.mix(na.omit(object), nu = 3, g = g, obs.prob = TRUE, calc.im  = FALSE)
		pii <- mix$pii
		mu <- mix$mu
		sigma2 <- mix$sigma2
		shape <- mix$shape
		for(j in 1:g){
			dens[,j] <- dSN(x, mu[j], sigma2[j], shape[j])
			dens.list[[j]] <- data.frame(x = x, y = dens[,j])
			means[j] <- SNmean(mu[j], sigma2[j], shape[j])
			dens.max[j] <- which.max(dens.list[[j]]$y)
			modes[j] <- dens.list[[j]]$x[dens.max[j]]
		}
		mix$means <- means
		mix$modes <- modes
		mix$dens.list <- dens.list
	}
	mix

}


