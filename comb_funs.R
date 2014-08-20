######################
#functions to combine estimates across different sites
######################


####################
##Combined Estimates##
####################

#Takes each region and calculates combined 
#estimate for region for a component

#Only for cities with monitors 
healthest.allcities <- function(cities, lag,
	outcome, reg.formula, estimates) {
	

	lcit <- length(cities)
	
	
	#Find health risk estimates for each city for the component
	#Five columns for overall est and 4 season-specific estimates
	regcoefs <- matrix(ncol = 5, nrow = lcit)
	regses <- matrix(ncol = 5, nrow = lcit)
	regcoefs.seas <- regcoefs
	regses.seas <- regses

	#variance covariance of season estimates for each city
    varall <- array(dim = c(4, 4, lcit))
    #sample size (number of days) for each city
 	sampsize <- matrix(ncol = 3, nrow = lcit)
 	
	for(j in 1:length(cities)){
		# print(cities[j])  
        
		cityrisk <- try(healthrisk.city(city = cities[j], 
				lag = lag,
				outcome = outcome, 
				reg.formula = reg.formula,
				estimates = estimates), silent = T)
				
		if(class(cityrisk) != "try-error") {	

			regcoefs[j, ] <- cityrisk[[1]][, 1]
			regses[j, ] <- cityrisk[[1]][, 2]
					
			regcoefs.seas[j, ] <- cityrisk[[2]][, 1]
			regses.seas[j, ] <- cityrisk[[2]][, 2]
	
	        varall[,, j] <- cityrisk[[3]]
	        sampsize[j, ] <- cityrisk[[4]]
			}

		}#end loop over cities


		rownames(regcoefs) <- cities
		rownames(regses) <- cities
		rownames(regcoefs.seas) <- cities
		rownames(regses.seas) <- cities


		#results
		list(list(regcoefs, regses),
			list(regcoefs.seas, regses.seas),
			varall, sampsize)

	}




comb.overall <- function(res.healthest, iqr, type,
	regions = NULL, print1 = F, print2 = F) {
	
	regcoefs <- res.healthest[[1]][[1]]
	regses <- res.healthest[[1]][[2]]
	varcovar <- res.healthest[[3]]
	
	if(type == "regional") {
		out <- comb.reg(regcoefs, regses, regions, iqr, print1)
		
	}else if(type == "natsea") {
		out <- comb.sea(regcoefs, regses, varcovar, iqr, print1, print2)
	}
	
	
	out
}




comb.sea <- function(regcoefs, regses, varcovar, 
	iqr, print1 = F, print2 = F){

          
    #find complete for annual      
	ann <- complete.cases(regcoefs[, 1])
	coef.ann <- regcoefs[ann, 1]
	ses.ann <- regses[ann, 1]

	#find complete for season
	season <- complete.cases(regcoefs[, -1])
	coef.seas <- regcoefs[season, -1]
	var.seas <- varcovar[,, season]

	#get rid of NA var terms
	keeps <- 0
	for(i in 1 : dim(var.seas)[3]) {
		# print(i)
		lna <- length(which(is.na(var.seas[, , i])))
		ln1 <- length(which(abs(var.seas[,, i]) > 10000))
		
		if(lna != 0 | ln1 != 0) {
			# browser()
			keeps <- c(keeps, i)
		}
		
	}
	keeps <- keeps[-1]
	if(length(keeps) > 0) {
		var.seas <- var.seas[,, -keeps]
		coef.seas <- coef.seas[-keeps, ]
		print(paste("remove large var/NA var:", cities[keeps]))
	}
	
	if(dim(var.seas)[3] < 10) {
		stop("error: less than 10 cities")
	}
	
	
	whl <- which(abs(coef.ann) > 1000)
	if(length(whl) > 0) {
		coef.ann <- coef.ann[-whl]
		ses.ann <- ses.ann[-whl]
		print(paste("remove: too large", cities[whl]))
		
	}
          
	tln.ann <- tlnise(Y = coef.ann, V = ses.ann^2, 
		brief = 0, maxiter = 5000, Tol = 10^(-7), prnt = print1)
	tln.sea <- tlnise(Y = coef.seas, V = var.seas,
		brief = 2, maxiter = 10000, Tol = 10^(-7), prnt = print2)
	var.post <- tln.sea$Dgamma
	

	combine.est <- matrix(nrow = 5, ncol = 2)
	combine.est[1, ] <- tln.ann$gamma[1 : 2]
	combine.est[-1, ] <- tln.sea$gamma[, 1 : 2]
	
	#get season-specific estimates
	combine.sea <- seas.ests(combine.est, var.post)
	
	rownames(combine.est)<-c("ann", "wint", "cspri", 
		"csumm", "cfall")
	rownames(combine.sea) <- c("ann", "wint", "spri",
		"summ", "fall")
	
	
	iqr.est <- 100 * (exp(combine.sea * iqr) - 1)
	est.sea <- list(combine.sea, iqr.est)
  	names(est.sea)<-c("Estimate(sea)","IQR increase")

	list(combine.est, est.sea, list(tln.ann$A, tln.sea$A),
		var.post)

}




comb.reg <- function(ests, ses, regions, iqr, print1 = F) {
		
	lreg <- length(regions)
	reg.ests <- matrix(nrow = lreg, ncol = 2)
	
	#Find combined estimate for each region
	for(i in 1 : lreg) {
		
		#find rows of coefficients for region
		cits <- names(areawh)[which(areawh == regions[i])]
        whr <- which(rownames(ests) %in% cits)
		regcoefs.reg <- ests[whr, 1]
		ses.reg <- ses[whr, 1]
		
		
		#check NAs
		whc <- which(is.na(regcoefs.reg) == TRUE)
		if(length(whc) != 0){
			print("regionNA")

		}
		
		tln.reg <- tlnise(Y = regcoefs.reg, V = ses.reg^2,
			brief = 0, maxiter = 5000, Tol = 10^(-7), prnt = print1)
		reg.ests[i, ] <- tln.reg$gamma[1 : 2]

	}

	reg.ests
	
	iqr.est <- 100 * (exp(reg.ests * iqr) - 1)
	est.reg <- list(reg.ests, iqr.est)
  	names(est.reg)<-c("Estimate(reg)", "IQR increase")
  	
  	list(est.reg)

}





