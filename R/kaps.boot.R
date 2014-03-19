##########################################################################
####
#### Multiway-splits adaptive partitioning with distribution free approach 
####
####		Soo-Heang EO and HyungJun CHO 
####
####		Version 1.0.0
####
####		30 Oct, 2013
####
##########################################################################

############################################################
# Fit KAPS with Boostrapping
############################################################
kaps.boot <- function(formula, data, type.boot = c("Boot", "random"), K = 2:5, B = 200, minors = kaps.control()){

	#options(warn = -1)
	#cat("Plz, wait some minutes.... \n")
	if(any(K == 1)) stop("the number of subgroups (K) is greater than 1.")
	type.boot <- match.arg(type.boot)
	
	n = nrow(data) # total number of observations concerning with time 
	rownames(data) = N <- 1:n

    if(n == 0L) stop("0 (non-NA) cases.")
	if(length(K) < 1) stop("the minimum number of subgroups (K) is greater than 1.")
	
	#call = match.call()

	Boot.over.stat = Boot.pair.stat = matrix(0, nrow = B, ncol = length(K))
	result = matrix(0.000, nrow = length(K), ncol = 4)
	colnames(result) <- c("over_mean", "over_se", "pair_mean", "pair_se")
	rownames(result) <- paste("K=",K,sep="")

	### generate Bootstrap samples
	#index.Boot = replicate(B, sample(1:n, n, replace = TRUE))
	#dat.index <- matrix(NA, nrow = n, ncol = B)
	#dat.index <- apply(dat.index, 2, function(x,n) sample(n, replace = TRUE), n = n)
	#index.OOB = apply(index,Boot, 2, function(x, n) n[-unique(a)], n=n)
	#dat.oob.index <- dat.index
	#cat("Now, KAPS Bootstrap is working. Please, wait a minute.^^\n")
	if(minors@ncl >= 2){
		#cl = makeCluster(minors@ncl)
		registerDoMC(minors@ncl)
		cat("Number of cores:", getDoParWorkers(), "\n")
		boot.res = foreach(i = 1:B, .combine = rbind) %dopar% {
			if(type.boot == "Boot"){
				index.boot = sample(n, n, replace = TRUE)
				learning = data[index.boot, , drop = FALSE]
				test = data[N[-unique(index.boot)], , drop = FALSE]
				mindat.boot = floor(nrow(learning) * 0.05)
			} else if(type.boot == "random"){
				index.random = sample(n, floor(n * 0.7), replace = FALSE)
				learning = data[index.random, , drop = FALSE]
				test = data[-index.random, , drop = FALSE]
				mindat.boot = floor(nrow(learning) * 0.05)
			}

			fit = try(lapply(K, kaps.fit, formula = formula, data = learning, mindat = mindat.boot, minors = minors), silent = TRUE)
			if(class(fit) != "try-error"){
				test.stat = sapply(fit, kaps.perm, newdata = test, permute = FALSE)	
				return(test.stat)		
			} 
		}
		#stopCluster(cl)
	} else{
		boot.res = foreach(i = 1:B, .combine = rbind) %do% {
			if(type.boot == "Boot"){
				index.boot = sample(n, n, replace = TRUE)
				learning = data[index.boot, , drop = FALSE]
				test = data[N[-unique(index.boot)], , drop = FALSE]
				mindat.boot = floor(nrow(learning) * 0.05)
			} else if(type.boot == "random"){
				index.random = sample(n, floor(n * 0.7), replace = FALSE)
				learning = data[index.random, , drop = FALSE]
				test = data[-index.random, , drop = FALSE]
				mindat.boot = floor(nrow(learning) * 0.05)
			}
		
			fit = try(lapply(K, kaps.fit, formula = formula, data = learning, mindat = mindat.boot, minors = minors), silent = TRUE)
			if(class(fit) != "try-error"){
				test.stat = sapply(fit, kaps.perm, newdata = test, permute = FALSE)
				return(test.stat)
			} 
		}
	}


	### output
	ind.row = 1:nrow(boot.res) %% 2
	Boot.over.stat[1:B,] <- boot.res[ind.row == 1,]
	Boot.pair.stat[1:B,] <- boot.res[ind.row == 0,]
	if(minors@sel == "mean"){
		# mean
		result[,1] <- apply(Boot.over.stat, 2, mean, na.rm = TRUE)
		result[,3] <- apply(Boot.pair.stat, 2, mean, na.rm = TRUE)
		#result[,2] <- apply(Boot.over.stat, 2, sd, na.rm = TRUE) / sqrt(B)
		#result[,4] <- apply(Boot.pair.stat, 2, sd, na.rm = TRUE) / sqrt(B)
	} else if(minors@sel == "median"){
		# median
		result[,1] <- apply(Boot.over.stat, 2, median, na.rm = TRUE)
		result[,3] <- apply(Boot.pair.stat, 2, median, na.rm = TRUE)
		#result[,2] <- apply(Boot.over.stat, 2, sd, na.rm = TRUE) / sqrt(B)
		#result[,4] <- apply(Boot.pair.stat, 2, sd, na.rm = TRUE) / sqrt(B)
	} else if(minors@sel == "trim"){
		# censoring-related trimmed mean
		CR = table(data[,all.vars(formula)[2]])[1]
		result[,1] <- apply(Boot.over.stat, 2, mean, na.rm = TRUE, trim = CR)
		result[,3] <- apply(Boot.pair.stat, 2, mean, na.rm = TRUE, trim = CR)
		#result[,2] <- apply(Boot.over.stat, 2, sd, na.rm = TRUE) / sqrt(B)
		#result[,4] <- apply(Boot.pair.stat, 2, sd, na.rm = TRUE) / sqrt(B)
	} else if(minors@sel == "test"){
		# test of proportion
		# this algorithm is used for testing the null that the proportion are the same
		# H0: p = 0.05
		# H1: p < 0.05
		over.sum = pair.sum = c()
		for( i in 1:ncol(Boot.over.stat) ){
			over.sum[i] <- sum(Boot.over.stat[,i] > minors@alpha)
			test.tmp <- prop.test(over.sum[i], n = B, alternative = "less")
			result[i,1] <- test.tmp$p.value
			result[i,2] <- test.tmp$statistic

			pair.sum[i] <- sum(Boot.pair.stat[,i] > minors@alpha)
			test.tmp <- prop.test(pair.sum[i], n = B, alternative = "less")
			result[i,3] <- test.tmp$p.value
			result[i,4] <- test.tmp$statistic
		}
	} 

	res = list()
	res$stat =result
	res$boot.over.stat = Boot.over.stat
	res$boot.pair.stat = Boot.pair.stat
	return(res)
}
# END @ 09 Nov 2013
