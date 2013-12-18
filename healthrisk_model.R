#####################################
# Functions to find city-specific (including seasonal) 
# estimates of health risks
# for single pollutant model, 1 lag
#
# 8/23/13
#
#
######################################




#This function takes each city and component and 
#calculates the health estimate using the regression

#
#namecomponent<-"SULFATE"
#city is four letter community [e.g. city <- "pitt"]
#lag is single day lag [e.g. lag <- 0]
#outcome is from mortality data [e.g. outcome <- "death"]
#reg.formula corresponds to regression formula, including * for number of years
#[e.g. reg.formula <- "factor(agecat) + factor(dow) + ns(tmpd, df = 3) + ns(Lag(tmpd, 1, group = agecat),df = 3) + ns(date, df = 8 * ")  ]
#estimates are community specific estimates for constituent of choice
healthrisk.city <- function(city, lag,
	outcome, reg.formula, estimates) {

	#find community average constituent concentrations
	comm.avg <- estimates[[city]]
	colnames(comm.avg) <- c("date", "cons")


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
		data <-data.frame(data, 
			Lag(data[, "cons"], k = lag, group = data[, "agecat"]))	
		#rename lagged column	
		colnames(data)<-c(colnames(data)[-ncol(data)], "lcons")

		#determine which days fall in which season
		seasons <- get.season(data)
		data <- data.frame(data, seasons)
		#rename last column
		colnames(data) <- c(colnames(data)[-ncol(data)], "season")

		#GET RESULTS
		

		out <- reg.health(data = data, outcome = outcome,
			reg.formula = reg.formula, city = city)
		out
		
	}else{
		print("no overlap dates")
	}


}




reg.health <- function(data, outcome, reg.formula, city) {
	
	nyears <- length(unique(substr(data[,"date"], 1, 4)))
	
	#set up formulas
	reg.form.int <- paste0(outcome, "~lcons * season + ",
		reg.formula, nyears, ")")
	reg.form <- paste0(outcome,"~lcons + ", reg.formula,
		nyears, ")")

	#run regression interacting cons with season, find coefficients
	glm.int <- try(glm(formula = eval(reg.form.int),
		data = data,
		family = quasipoisson,
		control = glm.control(maxit = 1000)))
	glmcoef.int <- summary(glm.int)[["coefficients"]]
	#find interaction terms
	wh.int <- which(substring(rownames(glmcoef.int), 1, 5) == "lcons")
	seas.int.names <- paste0(rep("lcons", 4), c("", rep(":season", 3)), 
		c("", seq(2, 4)))
	match <- match(seas.int.names, rownames(glmcoef.int))
	
	if(length(wh.int) != 4) {
		print(c("Error: cannot fit seasons for", city))
		print(rownames(glmcoef.int)[wh.int])
	}
	
	#run regression, no interaction, find coefficients		
	glm.ann <- try(glm(formula = eval(reg.form),
		data = data,
		family = quasipoisson,
		control = glm.control(maxit=1000)))


	#set up results
	#rows for annual, and 4 season estimates
	regcoef<-matrix(nrow = 5, ncol = 2)
	rownames(regcoef) <- c("annual", paste0("wint", c("", 
			":chspri", ":chsumm", ":chfall")))
			
	vcovar.seas <- matrix(nrow = 4, ncol = 4)	
		
		
	#if no error in interaction model	
	if(class(glm.int[1]) != "character") {
		
		#get estimate and standard error for season coef
		regcoef[-1, ] <- glmcoef.int[match, c(1, 2)]
		vcovar.seas <- vcov(glm.int)[match, match]
			
		#save annual estimates to regcoef		
		if(class(glm.ann[1]) != "character") {
			glmann <- summary(glm.ann)[["coefficients"]]
			regcoef[1,] <- glmann[which(substr(rownames(glmann), 
				1, 5) == "lcons"), c(1, 2)]

		#no estimates for annual model
		}else{
			print("error,annual")
		}	

		#get estimates for seasons
		regcoef.seas <- seas.ests(regcoef, vcovar.seas)
		rownames(regcoef.seas) <- c("annual", "wint", "spri", "summ", "fall")
		
		}#no estimates for seasonal model

	###SAMPLE SIZE
	#find sample size as degrees of freedom
	sampsize.int <- (glm.int$df.null + 1) / 3
	#find sample size as degrees of freedom
    sampsize.ann <- (glm.ann$df.null + 1) / 3	
	#find sample size as length of cons data
	sampsize.dat <- length(unique(data[which(is.na(data$cons) == FALSE), 1]))
	samp.size <- c(sampsize.int, sampsize.ann, sampsize.dat)
	names(samp.size) <- c("int nulldf", "ann nulldf", "noNA dat")
	
	
	
	out <- list(regcoef, regcoef.seas, vcovar.seas, samp.size)
	names(out) <- c("Coefficient", "Season est", 
			"vcovar.seas", "sampsize")
	out
}


###function to calculate seasonal estimates 
# coefficients are in change from winter
seas.ests <- function(regcoef, vcovar.seas) {
	regcoef.seas <- regcoef
	for(k in 3 : 5) {
		#get regression coef for season
		regcoef.seas[k, 1] <- regcoef.seas[2, 1] + regcoef.seas[k, 1]
		#get standard error
		regcoef.seas[k, 2] <- sqrt(vcovar.seas[1, 1] + 
			vcovar.seas[(k - 1), (k - 1)] + 
			2 * vcovar.seas[1, (k - 1)])
		}
	regcoef.seas	
}
