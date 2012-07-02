kaps <- function(formula, data, K = 2:5, V = 5, mindat  , ...){
###################################################################
####
####		Multiway-splits by adaptive partitioning
####
###################################################################
#####           pre-processing step
	options(warn = -1)
	if(missing(mindat)) mindat = floor(nrow(data) * 0.05)
	minors = kaps.control(...)
	if(any(K == 1)) stop("the number of subgroups (K) must be greater than .")
	n <- nrow(data) # total number of observations concerning with time 
	rownames(data) <- 1:n
    if(n == 0L) stop("0 (non-NA) cases.")
	if(length(K) < 1) stop("the minimum number of subgroups (K) is greater than 1.")
	if(length(K) > 10) stop("the maximum number of subgroups (K) is too large.")
	call <- match.call()
	# Obtain candidate log-rank statistics at each group.
	## CHECKME: modify by parallel computing
	### parallel computing in order to find optimal k subgroups
	if(minors@ncl ==1){
		cat("Now, selecting a set of cut-off points...\n")
		aps <- lapply(K, apss, formula = formula, data = data, mindat= mindat, minors = minors)
		if(minors@fold & length(K) >= 2){
			cat("Now, finding an optimal K...\n")
			elbow.stat <- K.apss(K, V, formula, data, mindat, minors )
		}
	}
	else{
		ncl <- makeCluster(minors@ncl)
		Export.func <- c("as.Formula", "model.part", "survdiff", "Surv")
		clusterExport(ncl, Export.func)
		cat("Now, running kaps algorithm for your data\n")
		aps <- parLapply(ncl,K, apss, formula = formula, data = data, mindat= mindat, minors = minors)
		stopCluster(ncl)
		if(minors@fold & length(K) >= 2) {
			cat("Now, finding optimal K. \nPlease, wait...\n")
			elbow.stat <- K.apss(K, V, formula, data, mindat, minors )
		}
	}
	# Declare Group k as the number of terminal nodes if k is maximum among WH statistic.
	test.stat <- sapply(aps, adj.test)
	## output of test.stat
	# test.stat[1,] = overall statistic
	# test.stat[2,] = optimal pair test staitstic
	# test.stat[3,] = cube-root transformation
	# test.stat[4,] = WH approximation statistic
	if(minors@fold & length(K) >= 2){
		if(length(K) > 1) {
			test.tmp <- elbow.stat[,elbow.stat[2,] >= 3.8, drop = FALSE] # optimal pair statistic
			if( ncol(test.tmp) >= 1) {
				Xvalue <- round(1 - pchisq(q = test.tmp[2,],df=1),4) # optimal pair p-value
				index <- which(Xvalue <= 0.05)
				index <- which.max(test.tmp[4,index]) # WH statistic
				if(length(index) >= 2) index <- index[length(index)]
			}
			else index <- 1
		}
		else index <- 1
	}
	else{
		if(length(K) > 1) {
			test.tmp <- test.stat[,test.stat[2,] >= 3.8, drop = FALSE]
			if( ncol(test.tmp) >= 1){
				Xvalue <- round(1 - pchisq(q = test.tmp[2,],df=1),4)
				index <- which(Xvalue <= 0.05)
				index <- which.max(test.tmp[4,index])
				if(length(index) >= 2) index <- index[length(index)]
			}
			else index <- 1
		}
		else index <- 1
	}

	result <- aps[[index]]
	result@index <- as.integer(index)
	result@pvalue <- 1 - pchisq(q = test.stat[2,], df = 1) # optimal pair p-value
	result@groups <- K
	attr(result@groups,"names") <- paste("K=",K,sep="")
	result@Z <- test.stat[1,] #overall statistic for selection candidate
	result@X <- test.stat[2,] #optimal pair p-value
	result@WH <- test.stat[3,] #cube-root transformation
	result@t <- test.stat[4,] #WH approximation statistic
	result@results <- aps
	if(minors@fold & length(K) >= 2) result@elbow <- elbow.stat
	result@call <- call
	return(result)
}
