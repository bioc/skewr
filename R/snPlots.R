


plotMix <- function(sample, model){

	sample <- prepSamp(sample)

	mix.hist(sample, model, ylim = c(0, .5), xlim = c(0, 18))

}



plotComponents <- function(sample, model, idmrVec = NULL, snpsVec = NULL, chip = "", channel = "", type = "", legend = FALSE){

	naVec <- NULL
	if(legend) sub <- type else sub <- NULL
	#if(is.matrix(sample) || !is(model, 'Skew.normal')) stop('Function only accepts a single sample and a single model')

	x <- seq(0, 18, length = 1000000)
	pii <- model$pii
	mu <- model$mu
	sigma2 <- model$sigma2
	shape <- model$shape
	modes <- model$modes
	means <- model$means
	dens.list <- model$dens.list
	order <- order(mu)
	mu <- mu[order]
	pii <- pii[order]
	sigma2 <- sigma2[order]
	shape <- shape[order]
	modes <- modes[order]
	means <- means[order]
	dens.list <- dens.list[order]
	g <- length(pii)
	myColors <- brewer.pal(9, "Set1")[1:g]
	tickColors <- brewer.pal(9, "Set1")[4:6]
	name <- colnames(sample)

	dens.mix <- d.mixedSN(x, pii, mu, sigma2, shape) #xx to x
	y.max <- max(dens.mix)+0.2   #########################
	x.y.max <- x[which(dens.mix == y.max)]



	hist(sample, breaks = 60,probability=TRUE, border = 'grey', ylim = c(0, y.max),
		xlim = c(4,16),xlab = paste(channel, "Log2 Intensity"), main = NULL, sub = sub)
      	lines(x, dens.mix, lwd = 1) #xx to x


## finding the mean and dropping the line from the curve to the x-axis
	for(i in 1:g){
		lines(x, pii[i] * dens.list[[i]]$y, type = 'l', col = myColors[i], lwd = 1)
#		means[i] <- SNmean(mu[i], sigma2[i], shape[i])
		dens.max <- which.max(dens.list[[i]]$y)
#		modes[i] <- dens.list[[i]]$x[dens.max[i]]
		lines(segments(means[i], 0, means[i], pii[i] * dSN(means[i], mu[i], sigma2[i], shape[i]), lty = 2, lwd = 1, col = myColors[i]))
		lines(segments(modes[i], 0, modes[i], pii[i] * dens.list[[i]]$y[dens.max], col = myColors[i], lwd = 1))
	}
	lines(density(sample, na.rm = TRUE), lty = 3, lwd = 1)
	means <- round(means, 3)
	modes <- round(modes, 3)


## plot segments for points in pointsVec
	if(!is.null(idmrVec) && !is.null(snpsVec) && !is.null(naVec)){
		for(i in 1:length(naVec)){
			lines(segments(naVec[i], y.max * 0.10, naVec[i], y.max * 0.15, col = tickColors[3]))
		}
		for(i in 1:length(idmrVec)){
			lines(segments(idmrVec[i], y.max * .050, idmrVec[i], y.max * .099, col = tickColors[2]))
		}
		for(i in 1:length(snpsVec)){
			lines(segments(snpsVec[i], 0, snpsVec[i], y.max * .049, col = tickColors[1]))
		}
	}
	else if(!is.null(idmrVec) && !is.null(snpsVec)){
		for(i in 1:length(idmrVec)){
			lines(segments(idmrVec[i], y.max * .050, idmrVec[i], y.max * .099, col = tickColors[2]))
		}
		for(i in 1:length(snpsVec)){
			lines(segments(snpsVec[i], 0, snpsVec[i], y.max * .049, col = tickColors[1]))
		}
	}
	else if(!is.null(idmrVec) && !is.null(naVec)){
		for(i in 1:length(naVec)){
			lines(segments(naVec[i], y.max * .050, naVec[i], y.max * .099, col = tickColors[3]))
		}
		for(i in 1:length(idmrVec)){
			lines(segments(idmrVec[i], 0, idmrVec[i], y.max * .049, col = tickColors[2]))
		}
	}
	else if(!is.null(snpsVec) && !is.null(naVec)){
		for(i in 1:length(naVec)){
			lines(segments(naVec[i], y.max * .050, naVec[i], y.max * .099, col = tickColors[3]))
		}
		for(i in 1:length(snpsVec)){
			lines(segments(snpsVec[i], 0, snpsVec[i], y.max * .049, col = tickColors[1]))
		}
	}
	else if(!is.null(idmrVec)){
		for(i in 1:length(idmrVec)){
			lines(segments(idmrVec[i], 0, idmrVec[i], y.max * .05, col = tickColors[2]))
		}
	}
	else if(!is.null(snpsVec)){
		for(i in 1:length(snpsVec)){
			lines(segments(snpsVec[i], 0, snpsVec[i], y.max * .05, col = tickColors[1]))
		}
	}
	else if(!is.null(naVec)){
		for(i in 1:length(naVec)){
			lines(segments(naVec[i], 0, naVec[i], y.max * .05, col = tickColors[3]))
		}
	}
## create legend
## to alternate means and modes:
##
##	as.vector(rbind(means, modes))
	#if(x.y.max < 10) x.coord <- 4 #13
	#else x.coord <- 4 #5
  x.coord <- 4
	y.coord <- y.max #+ .05
	if(legend){
		legend(x.coord, y.coord, c(paste(rep(c("mean =", "mode ="), each = g), c(means, modes)), "dens sum","non-par dens") , col = c(rep(myColors,2), 'black', 'black'),
        		lty = c(rep(c(2, 1), each = g), 1, 3), lwd = c(rep(2, 2*g+1), 1), xpd = NA, cex = 0.9) #TRUE)
	}
}


panelPlots <- function(MethyLumiSet, typeIRedModels, typeIGreenModels, typeIIModels, 
                       plot = c('panel', 'frames'), samp.num = NULL, frame.nums = 1:9, norm = '',
					   idmr = TRUE, snps = TRUE){
	if(!any(grepl('^rs', rownames(MethyLumiSet)))) snps <- FALSE
	plot <- match.arg(plot)

	if (plot == 'frames'){
		if(is.null(samp.num)){
			if(ncol(MethyLumiSet) > 1) stop('To see single plots, a sample number must be specified')
			start <- 1
			end <- 1
			legend <- TRUE
		}
		else{
			if(length(samp.num) > 1) stop('Only one sample number can be specified at this time')
			start <- samp.num
			end <- samp.num
			legend <- TRUE
		}
		
	}
	else{
		if(!is.null(samp.num)){
			if(length(samp.num) > 1) stop('Only one sample number can be specified at this time')
			start <- samp.num
			end <- samp.num
			legend <- FALSE
		}
		else{
			start <- 1
			end <- ncol(MethyLumiSet)
			legend <- FALSE
		}
	}
	frame.nums <- 1:9 %in% frame.nums

	tickColors <- brewer.pal(9, "Set1")[4:6]
	all.betas <- betas(MethyLumiSet)
	typeI.red.meth <- prepSamp(subsetProbes(MethyLumiSet, 'M', 'I-red', idmr = idmr, snps = snps))
	typeI.red.unmeth <- prepSamp(subsetProbes(MethyLumiSet, 'U', 'I-red', idmr = idmr, snps = snps))
	typeI.red.betas <- cbind(all.betas[rownames(typeI.red.meth),], NA)
	typeI.green.meth <- prepSamp(subsetProbes(MethyLumiSet, 'M', 'I-green', idmr = idmr, snps = snps))
	typeI.green.unmeth <- prepSamp(subsetProbes(MethyLumiSet, 'U', 'I-green', idmr = idmr, snps = snps))
	typeI.green.betas <- cbind(all.betas[rownames(typeI.green.meth),], NA)
	typeII.meth <- prepSamp(subsetProbes(MethyLumiSet, 'M', 'II', idmr = idmr, snps = snps))
	typeII.unmeth <- prepSamp(subsetProbes(MethyLumiSet, 'U', 'II', idmr = idmr, snps = snps))
	typeII.betas <- cbind(all.betas[rownames(typeII.meth),], NA)

	if(snps){
		snp.names <- grep('rs', rownames(MethyLumiSet), value = TRUE)
		snpI.red.meth <- cbind(typeI.red.meth[rownames(typeI.red.meth) %in% snp.names,], NA)
		snpI.red.unmeth <- cbind(typeI.red.unmeth[rownames(typeI.red.unmeth) %in% snp.names,], NA)
		snpI.red.betas <- cbind(all.betas[rownames(snpI.red.meth),], NA)
		snpI.green.meth <- cbind(typeI.green.meth[rownames(typeI.green.meth) %in% snp.names,], NA)
		snpI.green.unmeth <- cbind(typeI.green.unmeth[rownames(typeI.green.unmeth) %in% snp.names,], NA)
		snpI.green.betas <- cbind(all.betas[rownames(snpI.green.meth),], NA)
		snpII.meth <- cbind(typeII.meth[rownames(typeII.meth) %in% snp.names,], NA)
		snpII.unmeth <- cbind(typeII.meth[rownames(typeII.unmeth) %in% snp.names,], NA)
		snpII.betas <- cbind(all.betas[rownames(snpII.meth),], NA)
	}
	if(idmr){
		data(iDMR, envir = environment())
		iDMR.typeI.red.meth <- cbind(typeI.red.meth[rownames(typeI.red.meth) %in% iDMR,], NA)
		iDMR.typeI.red.unmeth <- cbind(typeI.red.unmeth[rownames(typeI.red.unmeth) %in% iDMR,], NA)
		iDMR.typeI.red.betas <- cbind(all.betas[rownames(iDMR.typeI.red.meth),], NA)
		iDMR.typeI.green.meth <- cbind(typeI.green.meth[rownames(typeI.green.meth) %in% iDMR,], NA)
		iDMR.typeI.green.unmeth <- cbind(typeI.green.unmeth[rownames(typeI.green.unmeth) %in% iDMR,], NA)
		iDMR.typeI.green.betas <- cbind(all.betas[rownames(iDMR.typeI.green.meth),], NA)
		iDMR.typeII.meth <- cbind(typeII.meth[rownames(typeII.meth) %in% iDMR,], NA)
		iDMR.typeII.unmeth <- cbind(typeII.unmeth[rownames(typeII.unmeth) %in% iDMR,], NA)
		iDMR.typeII.betas <- cbind(all.betas[rownames(iDMR.typeII.meth),], NA)
	}


	finishPlot <- function(){
	  mtext(paste(colnames(MethyLumiSet)[i], norm), 3, adj = 0.5, line = 2, cex = 1.2, outer = TRUE)
	  oldPar <- par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
	  plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
	  legend("bottom", c("SNP", "iDMR"), xpd = TRUE, horiz = TRUE, inset = c(0, 0),
	         bty = "n", pch = 15, col = tickColors[1:2], cex = 1.2)
    par(oldPar)
	}

	for(i in start:end){
		if (plot == 'panel') par(mfcol = c(3,3), mgp = c(2, 1, 0), mar=c(4,4,0.5,0.5), oma = c(4,0,4,0))
		else par(oma = c(4,0,4,0), mar=c(4,4,2,2), mgp = c(2,1,0))

		if(frame.nums[1]){
			plotComponents(typeI.red.meth[,i], typeIRedModels[[1]][[i]], if(idmr) iDMR.typeI.red.meth[,i] else NULL, if(snps) snpI.red.meth[,i] else NULL, NULL, "Methylated", "I-red", legend)
			if (plot == 'frames') finishPlot()
		}
		if(frame.nums[2]){ 
			plotComponents(typeI.red.unmeth[,i], typeIRedModels[[2]][[i]], if(idmr) iDMR.typeI.red.unmeth[,i] else NULL, if(snps) snpI.red.unmeth[,i] else NULL, NULL, "Unmethylated", "I-red", legend)
			if (plot == 'frames') finishPlot()
		}
		if(frame.nums[3]){
			if (plot == 'frames') sub <- "Type I-Red" else sub <- NULL
			plot(density(typeI.red.betas[,i], na.rm = TRUE), xaxt = 'n', bty = 'n', main = NA, xlab = "Beta values", sub = sub)
			axis(1, at = c(0.0, 0.2, 0.4, 0.6, 0.8, 1.0), labels = c('0.0', '0.2', '0.4', '0.6', '0.8', '1.0'))
			axis(1, at = 0.5, labels = FALSE, col.ticks = 'red', tck = -.08)
			h <- max(density(typeI.red.betas[, i], na.rm = TRUE)$y)
			if(idmr && snps){
				idmrBetas <- iDMR.typeI.red.betas[,i]
				idmrBetas <- na.omit(idmrBetas)
				snpBetas <- snpI.red.betas[,i]
				snpBetas <- na.omit(snpBetas)
				idmrHgt <- h * .1
				idmrBtm <- h * .051
				for(j in 1:length(idmrBetas)){
					lines(segments(idmrBetas[j], idmrBtm , idmrBetas[j], idmrHgt, col = tickColors[2]))
				}
				snpHgt <- h * .049
				for(j in 1:length(snpBetas)){
					lines(segments(snpBetas[j], 0, snpBetas[j], snpHgt, col = tickColors[1]))
				}
			}
			else if(idmr && !snps){
				idmrBetas <- iDMR.typeI.red.betas[,i]
				idmrBetas <- na.omit(idmrBetas)
				height <- h * .05
				for(j in 1:length(idmrBetas)){
					lines(segments(idmrBetas[j], 0, idmrBetas[j], height, col = tickColors[2]))
				}
			}
			else if(!idmr && snps){
				snpBetas <- snpI.red.betas[,i]
				snpBetas <- na.omit(snpBetas)
				height <- h * .05
				for(j in 1:length(snpBetas)){
					lines(segments(snpBetas[j], 0, snpBetas[j], height, col = tickColors[1]))
				}
			}
			if (plot == 'frames') finishPlot()
		}
		if(frame.nums[4]){
			plotComponents(typeI.green.meth[,i], typeIGreenModels[[1]][[i]], if(idmr) iDMR.typeI.green.meth[,i] else NULL, if(snps) snpI.green.meth[,i] else NULL, NULL, "Methylated", "I-green", legend)
			if (plot == 'frames') finishPlot()
		}
		if(frame.nums[5]){
			plotComponents(typeI.green.unmeth[,i], typeIGreenModels[[2]][[i]], if(idmr) iDMR.typeI.green.unmeth[,i] else NULL, if(snps) snpI.green.unmeth[,i] else NULL, NULL, "Unmethylated", "I-green", legend)
			if (plot == 'frames') finishPlot()
		}
		if(frame.nums[6]){
			if (plot == 'frames') sub <- "Type I-Green" else sub <- NULL
			plot(density(typeI.green.betas[,i], na.rm = TRUE), xaxt = 'n', bty = 'n', main = NA, xlab = "Beta values", sub = sub)
			axis(1, at = c(0.0, 0.2, 0.4, 0.6, 0.8, 1.0), labels = c('0.0', '0.2', '0.4', '0.6', '0.8', '1.0'))
			axis(1, at = 0.5, labels = FALSE, col.ticks = 'red', tck = -.08)
			h <- max(density(typeI.green.betas[,i], na.rm = TRUE)$y)
			if(idmr && snps){
				idmrBetas <- iDMR.typeI.green.betas[,i]
				idmrBetas <- na.omit(idmrBetas)
				snpBetas <- snpI.green.betas[,i]
				snpBetas <- na.omit(snpBetas)
				idmrHgt <- h * .1
				idmrBtm <- h * .051
				for(j in 1:length(idmrBetas)){
					lines(segments(idmrBetas[j], idmrBtm , idmrBetas[j], idmrHgt, col = tickColors[2]))
				}
				snpHgt <- h * .049
				for(j in 1:length(snpBetas)){
					lines(segments(snpBetas[j], 0, snpBetas[j], snpHgt, col = tickColors[1]))
				}
			}
			else if(idmr && !snps){
				idmrBetas <- iDMR.typeI.green.betas[,i]
				idmrBetas <- na.omit(idmrBetas)
				height <- h * .05
				for(j in 1:length(idmrBetas)){
					lines(segments(idmrBetas[j], 0, idmrBetas[j], height, col = tickColors[2]))
				}
			}
			else if(!idmr && snps){
				snpBetas <- snpI.green.betas[,i]
				snpBetas <- na.omit(snpBetas)
				height <- h * .05
				for(j in 1:length(snpBetas)){
					lines(segments(snpBetas[j], 0, snpBetas[j], height, col = tickColors[1]))
				}
			}
			if (plot == 'frames') finishPlot()
		}
		if(frame.nums[7]){
			plotComponents(typeII.meth[,i], typeIIModels[[1]][[i]], if(idmr) iDMR.typeII.meth[,i] else NULL, if(snps) snpII.meth[,i] else NULL, NULL, "Methylated (green)", "II", legend)
			if (plot == 'frames') finishPlot()
		}
		if(frame.nums[8]){
			plotComponents(typeII.unmeth[,i], typeIIModels[[2]][[i]], if(idmr) iDMR.typeII.unmeth[,i] else NULL, if(snps) snpII.unmeth[,i] else NULL, NULL, "Unmethylated (red)", "II", legend)
			if (plot == 'frames') finishPlot()
		}
		if(frame.nums[9]){
			if (plot == 'frames') sub <- "Type II" else sub <- NULL
			plot(density(typeII.betas[,i], na.rm = TRUE), xaxt = 'n', bty = 'n', main = NA, xlab = "Beta values", sub = sub)
			axis(1, at = c(0.0, 0.2, 0.4, 0.6, 0.8, 1.0), labels = c('0.0', '0.2', '0.4', '0.6', '0.8', '1.0'))
			axis(1, at = 0.5, labels = FALSE, col.ticks = 'red', tck = -.08)
			h <- max(density(typeII.betas[,i], na.rm = TRUE)$y)
			if(idmr && snps){
				idmrBetas <- iDMR.typeII.betas[,i]
				idmrBetas <- na.omit(idmrBetas)
				snpBetas <- snpII.betas[,i]
				snpBetas <- na.omit(snpBetas)
				idmrHgt <- h * .1
				idmrBtm <- h * .051
				for(j in 1:length(idmrBetas)){
					lines(segments(idmrBetas[j], idmrBtm , idmrBetas[j], idmrHgt, col = tickColors[2]))
				}
				snpHgt <- h * .049
				for(j in 1:length(snpBetas)){
					lines(segments(snpBetas[j], 0, snpBetas[j], snpHgt, col = tickColors[1]))
				}
			}
			else if(idmr && !snps){
				idmrBetas <- iDMR.typeII.betas[,i]
				idmrBetas <- na.omit(idmrBetas)
				height <- h * .05
				for(j in 1:length(idmrBetas)){
					lines(segments(idmrBetas[j], 0, idmrBetas[j], height, col = tickColors[2]))
				}
			}
			else if(!idmr && snps){
				snpBetas <- snpII.betas[,i]
				snpBetas <- na.omit(snpBetas)
				height <- h * .05
				for(j in 1:length(snpBetas)){
					lines(segments(snpBetas[j], 0, snpBetas[j], height, col = tickColors[1]))
				}
			}
			if (plot == 'frames') finishPlot()
		}
	
		if(plot == 'panel'){
			mtext("Type I-Red", 3, adj = 0.15, outer = TRUE)
			mtext("Type I-Green", 3, adj = 0.5, outer = TRUE)
			mtext("Type II", 3, adj = 0.85, outer = TRUE)
			finishPlot()
		}

	}

}
