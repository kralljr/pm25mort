


get.season <- function(data) {

	seasall <- rep(seq(1, 12), 2)
	seas.mat <- matrix(nrow = 4, ncol = 4)
	start <- 12
	for(i in 1 : 4) {
		seas.mat[i, ] <- seasall[start : (start + 3)]
		start <- start + 3
	}
	
	#find months and days
	months <- as.numeric(substr(data[, "date"], 6, 7))
	days <- as.numeric(substr(data[, "date"], 9, 10))
	
	seas <- vector(, length = nrow(data))
	for(i in 1:4){
	
		#if in middle two months
		whMIDDLE <- which(months %in% seas.mat[i,2:3])
		#if in first month, after the 20th
		whFIRST <- which((months %in% seas.mat[i,1]) & (days >= 21))
		#if in last month, and before the 20th
		whLAST <- which((months %in% seas.mat[i,4]) & (days < 21))
		
		seas[c(whMIDDLE, whFIRST, whLAST)] <- i
	}
	factor(seas)
	
}

	
