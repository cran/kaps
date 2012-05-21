apss <- function(formula, data, K = 3, mindat = floor(nrow(data) * 0.01), minors = apss.control()) {
###################################################################
####
####		Adaptive Partitioning and Substaging for survival data
####
###################################################################

#####           pre-processing step
	options(warn = -1)

	n <- nrow(data) # total number of observations concerning with time 
	rownames(data) <- 1:n
    if(n == 0L) stop("0 (non-NA) cases.")

	# original formula
    mf <- match.call() # save a copy of the call
	formula <- as.Formula(formula)
	#mf$formula <- f
	#mf <- eval(mf, parent.frame())
	vars <- all.vars(formula)
	
	if(!is.data.frame(data)) data <- as.data.frame(data)
	## Set the object of result
    result <- new("apss")
	result@formula <- formula
	
	##### Model Fitting
	X <- model.part(formula, data = data, rhs = 1, drop = FALSE)
	f <- update(formula, . ~ K)
	## treat pre-determined split points
	if(is.null(attr(minors@pre.pt,"names"))) pt.set <- lapply(X, function(x) sort(unique(x)) )
	else{
		pre.name <- colnames(X) != attr(minors@pre.pt, "names")
		pt.set <- lapply(X[,pre.name, drop = FALSE], function(x) sort(unique(x)) )
		pt.set <- c(pt.set, minors@pre.pt)
		X <- X[,attr(pt.set, "names"), drop = FALSE]
	}
	## treat pre-determined ranges
	scope <- minors@scope
	if(!is.null(attr(scope, "names"))) {
		scope.vars <- which(names(pt.set) == attr(scope, "names"))
		pt.set.name <- names(pt.set)
	}
	
	result@Chisq <- 0
	pt.set <- lapply(pt.set, function(x,upper, lower) x[x <= upper & x >= lower], upper = minors@upper.limit, lower = minors@lower.limit)
	
	v <- K-1
	for(i in 1: length(pt.set)){
		if(!is.null(attr(scope, "names"))){
			if(i %in% scope.vars) {
				rngs <- scope[[pt.set.name[i]]]  
				pt.set[[i]] <- pt.set[[i]][ pt.set[[i]] >= rngs[1] & pt.set[[i]] <= rngs[2] ]
			}
		}
		ord <- order(X[,i])
		data <- data[ord,]
		X <- X[ord, , drop = FALSE]
		test.where <- group.sel(X[,i], pt.set[[i]], K, mindat,data, f, minors)
		x.test <- test.where$test
		
		##### Result
		if(x.test[1,test.where$index] >= result@Chisq) {
			index <- test.where$index
			result@Chisq <- x.test[1,index]
			result@WH <- (x.test[1,index] / v)^(1/3)
			result@t <- (result@WH - (1- (2 / (9*v)))) / sqrt(2 / (9 * v))
			result@pvalue <-  1 - pt(q = result@t, df = 1 )
			result@pair <- x.test[3,index]
			result@where <- test.where$where
			result@split.pt <- sapply(unique(result@where), function(x,y,where) max(y[where == x]), y = X[,i], where = result@where) 
			result@split.pt <- result@split.pt[-length(result@split.pt)]
			result@split.var <- colnames(X[,i,drop =FALSE])
			result@data <- data
		}
	}
	result@groups <- K
	attr(result@groups,"names") <- paste("G=",K,sep="")
	result@mindat <- mindat
	result@Options <- minors
	return(result)
}

K.apss <- function(K, V, formula, data, mindat, minors){
	## Finding optimal K by cross-validation
	mf <- match.call() 
	formula <- as.Formula(formula)
	## apss by V-fold Cross-Validation
	if(V < 2) stop("V must be an integer which is greater than and equal to 2.")
	fold.no <- rep(1:V, nrow(data) %/% V) # generate new var.: fold number
	fold.no <- c(fold.no, (1:V)[nrow(data) %% V])
	fold.chi <- matrix(NA, nrow = V, ncol = length(K))
	fold.adj.chi <- matrix(NA, nrow = V, ncol = length(K))
	fold.WH <- matrix(NA, nrow = V, ncol = length(K))
	fold.t <- matrix(NA, nrow = V, ncol = length(K))
	result <- matrix(NA, nrow = 4, ncol = length(K))
	## CHECKME: modify by parallel computing
	if(minors@ncl == 1){
		for(i in 1:V){	
			train <- data[fold.no != i,, drop = FALSE]
			test <- data[fold.no == i,, drop = FALSE]
			sapss <- lapply(K, apss, formula = formula, data = train, mindat= mindat, minors = minors)
			tmp <- sapply(sapss, predict, newdata = test)
			fold.chi[i,] <- tmp[1,]
			fold.adj.chi[i,] <- tmp[2,]
			fold.WH[i,] <- tmp[3,]
			fold.t[i,] <- tmp[4,]
			rm(sapss, train, test, tmp)
		}
	}
	else{
		ncl <- makeCluster(minors@ncl)
		Export.func <- c("as.Formula", "model.part", "survdiff", "Surv")
		clusterExport(ncl, Export.func)
		for(i in 1:V){	
			train <- data[fold.no != i,, drop = FALSE]
			test <- data[fold.no == i,, drop = FALSE]
			sapss <- parLapply(ncl, K, apss, formula = formula, data = train, mindat= mindat, minors = minors)
			tmp <- parSapply(ncl, sapss, predict, newdata = test)
			fold.chi[i,] <- tmp[1,]
			fold.adj.chi[i,] <- tmp[2,]
			fold.WH[i,] <- tmp[3,]
			fold.t[i,] <- tmp[4,]
			rm(sapss, train, test, tmp)
		}
		stopCluster(ncl)
	}
	### output
	result[1,] <- apply(fold.chi,2, mean, na.rm = TRUE)
	result[2,] <- apply(fold.adj.chi,2, mean, na.rm = TRUE)
	result[3,] <- apply(fold.WH,2, mean, na.rm = TRUE)
	result[4,] <- apply(fold.t,2, mean, na.rm = TRUE)
	return(result)
}


apss.boot <- function(fit, B = 200){
### Adaptive Partitioning by Bootstrap
	if(!inherits(fit, "apss")) stop("This function requires the object of apss class as a main object.")
	data <- fit@data
	formula <- fit@formula
	K <- length(unique(fit@where))
	minors <- fit@Options
	mindat <- fit@mindat
	
	### generate Bootstrap samples
	n <- nrow(data)
	data.index <- matrix(NA, nrow = n, ncol = B)
	data.index <- apply(data.index, 2, function(x,n) sample(n, replace = TRUE), n = n)
	apss.match <- function(x, data, ...){
		data.tmp <- data[x,]
		res <- apss(formula = formula, data = data.tmp, K = K, mindat = mindat, minors = minors)
		res <- res@split.pt
	}
	cat("Now, Boostrap APSS working. Please, wait a minute.^^\n")
	split.pt <- apply(data.index, 2, apss.match, data = data, formula = formula, 
		K = K,  mindat = mindat, minors = minors)  
	#split.pt <- apply(split.pt, 1, function(x) names(table(x))[table(x) == max(table(x))])
	split.pt <- apply(split.pt, 1, mean, na.rm = TRUE, trim = 0.25) 
	return(as.numeric(split.pt))
}
