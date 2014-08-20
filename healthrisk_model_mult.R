#####################################
# Functions to find city-specific (including seasonal) 
# estimates of health risks 
# for multiple pollutant model, 1 lag
#
# 10/11/13
#
#
######################################




#This function takes each city and component and 
#calculates the health estimate using the regression


healthest.city.outm <- function(lags, outcome, reg.formula, 
	cons.names, cities, ests, iqrs = NULL, 
	seed = NULL, print1 = F) {

	set.seed(seed)
	
	healthest.city.lags <- list()
	healthest.overall <- list()
	
	
	#loop over lags
	for(j in 1 : length(lags)) {
		print(paste0("lag", lags[j]))
		
		cityest <- list()
		overallest <- list()
		
		
		#get city specific health effect estimates
		temp.est <- healthest.allcitiesm(cities = cities, 
			lag = lags[j], outcome = outcome, reg.formula = reg.formula,
			estimates = ests, cons = cons.names)

		
		#get national or regional estimates
		overallest <- comb.overallm(temp.est, iqrs[cons.names], cons.names,
				print1)
		
		healthest.city.lags[[j]] <- temp.est
		healthest.overall[[j]] <- overallest
            
            

		}#end loop over lags	
		
        list(healthest.city.lags, healthest.overall)
}





comb.overallm <- function(res.healthest, iqr, cons,
	print1 = F) {
	
	regcoefs <- res.healthest[[1]][[1]]
	regses <- res.healthest[[1]][[2]]

	combine.est <- matrix(nrow = length(cons), ncol = 2)
	iqr.est <- matrix(nrow = length(cons), ncol = 2)
	
	 #find complete for annual  
	for(i in 1 : length(cons)) {    
		ann <- complete.cases(regcoefs[, i])
		coef.ann <- regcoefs[ann, i]
		ses.ann <- regses[ann, i]
      
		tln.ann <- tlnise(Y = coef.ann, V = ses.ann^2, 
			brief = 0, maxiter = 5000, Tol = 10^(-7), prnt = print1)
			
		
		combine.est[i, ] <- tln.ann$gamma[1 : 2]
		iqr.est[i, ] <- 100 * (exp(combine.est[i, ] * iqr[i]) - 1)

		
		
	}

	rownames(combine.est)<- cons
	rownames(iqr.est)<- paste0("iqr", cons)
	

	list(combine.est, iqr.est, tln.ann$A)

}





#Only for cities with monitors 
healthest.allcitiesm <- function(cities, lag,
	outcome, reg.formula, estimates, cons) {
	

	lcit <- length(cities)
	
	
	#Find health risk estimates for each city for the component
	#Five columns for overall est and 4 season-specific estimates
	regcoefs <- matrix(ncol = length(cons), nrow = lcit)
	regses <- matrix(ncol = length(cons), nrow = lcit)

    #sample size (number of days) for each city
 	sampsize <- matrix(ncol = 2, nrow = lcit)
 	
	for(j in 1:length(cities)){
		# print(cities[j])  
        
		cityrisk <- try(healthrisk.citym(city = cities[j], 
				lag = lag,
				outcome = outcome, 
				reg.formula = reg.formula,
				estimates = estimates), silent = T)

		if(class(cityrisk) != "try-error") {
			regcoefs[j, ] <- cityrisk[[1]][, 1]
			regses[j, ] <- cityrisk[[1]][, 2]
	        sampsize[j, ] <- cityrisk[[2]]
        }


		}#end loop over cities


		rownames(regcoefs) <- cities
		rownames(regses) <- cities
		
		colnames(regcoefs) <- cons
		colnames(regses) <- cons

		#results
		list(list(regcoefs, regses),
					sampsize)

	}








#
#namecomponents<-c("SULFATE", "NITRATE")
#city is four letter community [e.g. city <- "pitt"]
#lag is single day lag [e.g. lag <- 0]
#outcome is from mortality data [e.g. outcome <- "death"]
#reg.formula corresponds to regression formula, including * for number of years
#[e.g. reg.formula <- "factor(agecat) + factor(dow) + ns(tmpd, df = 3) + ns(Lag(tmpd, 1, group = agecat),df = 3) + ns(date, df = 8 * ")  ]
#estimates are list of community specific estimates for constituents of choice
healthrisk.citym <- function(city, lag,
	outcome, reg.formula, estimates) {

	#find community average constituent concentrations
	comm.avg <- estimates[[1]][[city]]
	colnames(comm.avg) <- c("date", paste0("cons", 1))
	for(i in 2 : length(estimates)) {
		comm.avg2 <- estimates[[i]][[city]]
		colnames(comm.avg2) <- c("date", paste0("cons", i))
		
		comm.avg <- merge(comm.avg, comm.avg2, by = "date")
	}

	#Find city mortality data
	cityvar1 <- cityvar[[city]]
	#order by date
	cityvar1 <- cityvar1[ order(cityvar1[, "date"]), ]
	#after year 2000
	years <- as.numeric(substr(cityvar1[, 5], 1, 4))
	cityvar1 <- cityvar1[which(years >= 2000), ]
	
	#merge PM data with health data
	data <- merge(cityvar1, comm.avg,
		by = "date", all.x = TRUE)
		        
   	#if merge contains at least 1 day
	if(length(data[,1]) >= 1 ) {

		#create lag
		###HERE, fix lag
		lags <- matrix(nrow = nrow(data), ncol = length(estimates))
		for(i in 1 : length(estimates)) {
			lags[, i] <- Lag(data[, paste0("cons", i)], 
				k = lag, group = data[, "agecat"])
		}
		colnames(lags) <- paste0("lcons", seq(1, length(estimates)))
		data <-data.frame(data, lags)	


		#determine which days fall in which season
		seasons <- get.season(data)
		data <- data.frame(data, seasons)
		#rename last column
		colnames(data) <- c(colnames(data)[-ncol(data)], "season")

		#GET RESULTS
		

		out <- reg.healthm(data = data, outcome = outcome,
			reg.formula = reg.formula, city = city, nums = length(estimates))
		out
		
	}else{
		print("no overlap dates")
	}


}




reg.healthm <- function(data, outcome, reg.formula, city, nums) {
	
	nyears <- length(unique(substr(data[,"date"], 1, 4)))
	
	#set up formulas
	reg.form <- paste0(outcome,"~ ",paste(paste0("lcons", seq(1, nums)), collapse = "+"),
		" + ", reg.formula,
		nyears, ")")

	#run regression, no interaction, find coefficients		
	glm.ann <- try(glm(formula = eval(reg.form),
		data = data,
		family = quasipoisson,
		control = glm.control(maxit=1000)))


	#set up results
	#rows for annual, and 4 season estimates
	regcoef<-matrix(nrow = nums, ncol = 2)
				
	#save annual estimates to regcoef		
	if(class(glm.ann[1]) != "character") {
		glmann <- summary(glm.ann)[["coefficients"]]
		regcoef <- glmann[which(substr(rownames(glmann), 
			1, 5) == "lcons"), c(1, 2)]

	#no estimates for annual model
	}else{
		print("error,annual")
	}	

	###SAMPLE SIZE
	#find sample size as degrees of freedom
    sampsize.ann <- (glm.ann$df.null + 1) / 3	
	#find sample size as length of cons data
	sampsize.dat <- length(unique(data[complete.cases(data), 1]))
	samp.size <- c(sampsize.ann, sampsize.dat)
	names(samp.size) <- c("ann nulldf", "noNA dat")
	
	
	
	out <- list(regcoef, samp.size)
	names(out) <- c("Coefficient",  "sampsize")
	out
}

