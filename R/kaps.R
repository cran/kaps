kaps <- function(formula, data, K = 2:5, V = 5, mindat  , ...){
###################################################################
####
####		Multiway-splits by adaptive partitioning
####
###################################################################
#####           pre-processing step
	options(warn = -1)
	if(missing(mindat)) mindat = floor(nrow(data) * 0.05)
	minors = apss.control(...)
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
		if(minors@fold){
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
		if(minors@fold) {
			cat("Now, finding optimal K. \nPlease, wait...\n")
			elbow.stat <- K.apss(K, V, formula, data, mindat, minors )
		}
	}
	## Do permutation test
	# Declare Group k as the number of terminal nodes if k is maximum among WH statistic.
	test.stat <- sapply(aps, adj.test)
	if(minors@fold ){
		if(length(K) > 1) {
			test.tmp <- elbow.stat[,elbow.stat[2,] >= 3.8, drop = FALSE]
			if( ncol(test.tmp) >= 0) {
				pvalue <- round(1 - pnorm(q = test.tmp[4,]),4)
				zvalue <- round(1 - pchisq(q = test.tmp[2,],df=1),4)
				value.std <- (pvalue + zvalue) / 2
				index <- which(value.std == min(value.std))
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
				pvalue <- round(1 - pnorm(q = test.tmp[4,]),4)
				zvalue <- round(1 - pchisq(q = test.tmp[2,],df=1),4)
				value.std <- (pvalue + zvalue) / 2
				index <- which(value.std == min(value.std))
				if(length(index) >= 2) index <- index[length(index)]
			}
			else index <- 1
		}
		else index <- 1
	}

	result <- aps[[index]]
	result@index <- as.integer(index)
	#result@pvalue <- 1 - pt(q = test.stat[4,], df = 1 )
	result@pvalue <- 1 - pnorm(q = test.stat[4,])
	result@groups <- K
	attr(result@groups,"names") <- paste("K=",K,sep="")
	if(minors@fold) result@elbow <- elbow.stat
 	result@Chisq <- test.stat[1,]
	result@Z <- test.stat[2,]
	result@WH <- test.stat[3,]
	result@t <- test.stat[4,]
	#result@g.ID <- sapply(aps,function(x) x@where)
	result@candid <- aps
	result@call <- call
	return(result)
}
