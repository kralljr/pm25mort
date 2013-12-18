##########
# Functions to run all
# This is the outer function to run the health effects analysis
# 8/23/13
##########



##############
# Function to get health effect estimates for all constituents and all lags of interest
#lags is vector of single-day lags [e.g. lags = c(0, 1, 2)]
#outcome corresponds to column heading [e.g. outcome = "death"]
#reg.formula corresponds to regression formula, including * for number of years
#[e.g. reg.formula <- "factor(agecat) + factor(dow) + ns(tmpd, df = 3) + ns(Lag(tmpd, 1, group = agecat),df = 3) + ns(date, df = 8 * ")  ]
#cons.names is vector of constituents
#cities is list of cities
#ests is list of community average estimates for each city 
##where each element has two columns: date and est
#type indicates if want to do national/seasonal analysis or regional analysis
#regions is vector of regions if want to do regional analysis
########################
healthest.city.out <- function(lags, outcome, reg.formula, 
	cons.names, cities, ests, iqrs = NULL, 
	type = NULL, regions = NULL, seed = NULL, print1 = F, print2 = F) {

	set.seed(seed)
	
	healthest.city.lags <- list()
	healthest.overall <- list()
	
	
	#loop over lags
	for(j in 1 : length(lags)) {
		print(paste0("lag", lags[j]))
		
		cityest <- list()
		overallest <- list()
		
		
		#loop over constituetns
		for(i in 1 : length(cons.names)) { 
			print(cons.names[i])
			
			#get community average constituent concentration
			estimates <- ests[[cons.names[i]]]
			
			#get city specific health effect estimates
			temp.est <- healthest.allcities(cities = cities, 
				lag = lags[j], outcome = outcome, reg.formula = reg.formula,
				estimates = estimates)
			cityest[[i]] <- temp.est	
			
			#get national or regional estimates
			if(!is.null(type)) {
				overallest[[i]] <- comb.overall(temp.est, iqrs[cons.names[i]], 
					type, regions, print1, print2)
			}

			
			}#end loop over cons
			
			healthest.city.lags[[j]] <- list(cityest)
			healthest.overall[[j]] <- list(overallest)
            
            

		}#end loop over lags	
		
        list(healthest.city.lags, healthest.overall)
}



createci <- function(ests) {
	ci <- ests[1] + c(-1.96, 1.96) * ests[2]
	c(round(ests[1], 2), round(ci, 2))
}


tableest <- function(ests, cons) {
	n <- length(ests[[2]][[1]][[1]])
	out <- matrix(nrow = n, ncol = 4)
	rownames(out) <- cons
	for(i in 1 : n) {
		est2 <- ests[[2]][[1]][[1]][[i]][[2]][[2]][1, ]
		out[i, ] <- c(createci(est2), postprob(est2[1], est2[2]))
	}


	out
}



postprob <- function(pe, se){
	pn <- vector()
	for(i in 1 : length(pe)) {
		pn[i] <- pnorm(0, mean = pe[i], sd = se[i],
			lower.tail = F)
	}
	round(pn, 2)
}
