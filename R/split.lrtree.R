split.lrtree <-
function(dataSet, controls, duplicated.data){
#################################################################
## Binary Splits for Censored Data
## Select the growing method: logrank (Segal, 1988) 
#################################################################
## Add bootstrapping to find ideal split point.
	result <- new("SplitRes")
	if(controls@TGCtrl@tree.size == "Boot"){
		n <- nrow(dataSet@X)
		tree.boot <-function(n, dataSet, controls){
			ind<-sample(n,replace=TRUE)
			dataSet@Y <- dataSet@Y[ind,]
			dataSet@X <- dataSet@X[ind,,drop=FALSE]
			rownames(dataSet@Y) <- 1:n
			rownames(dataSet@X) <- 1:n
			tree.imp <- sapply(dataSet@X[duplicated.data,, drop = FALSE], logrank.tree,
				y = dataSet@Y[duplicated.data,],
				controls = controls
			)
			imp.which <- which(tree.imp[[2]] == max(tree.imp[[2]], na.rm = TRUE))
			if(length(imp.which) >= 2) imp.which <- imp.which[sample(length(imp.which), 1)]
			#res <- matrix(NA, nrow = 1, ncol= 3)
			res <- list()
			res[[1]] <- tree.imp[,imp.which][[2]] 
			res[[2]] <- tree.imp[,imp.which][[1]]
			res[[3]] <- colnames(tree.imp)[imp.which]
			return(res)
		}
		imp.pt <- replicate(controls@TGCtrl@B, tree.boot(n, dataSet, controls))
		#imp.pt <- replicate(20, tree.boot(n, dataSet, controls))
		# Majority voting to select split variable.
		tmp.var <- unlist(imp.pt[3,])
		result@Var <- names(which.max(table(tmp.var)))
		tmp.pt <- unlist(imp.pt[2,])
		#result@Point <- mean(tmp.pt[tmp.var == result@Var], na.rm = TRUE, trim = 0.1)
		tmp.pt <- names(which.max(table(tmp.pt)) ) # mode of candidate points
		result@Point <- as.numeric(tmp.pt)
		gr <- ifelse(dataSet@X <= result@Point, 1, 2)
		if((length(gr[gr == 1]) >= controls@SplitCtrl@minbucket) & 
			(length(gr[gr == 2]) >= controls@SplitCtrl@minbucket) ) {
			fit <- survdiff(dataSet@Y ~ gr, rho = controls@VarCtrl@rho)
			result@Impurity <- fit$chisq
		}
		else result@Impurity <- 0
		result@type <- ifelse(is.factor(dataSet@X[,result@Var]),"Factor","Numeric") 
		result@Obs <- which(dataSet@X[, result@Var] <= result@Point)
		return(result)
	}
	else{
		tree.imp <- sapply(dataSet@X[duplicated.data,, drop = FALSE], logrank.tree,
			y = dataSet@Y[duplicated.data,],
			controls = controls
		)
		imp.which <- which(tree.imp[[2]] == max(tree.imp[[2]], na.rm = TRUE))
		if(length(imp.which) >= 2) imp.which <- imp.which[sample(length(imp.which), 1)]
		# numeric: continuous predictor, factor: categorical predictor
		result@Impurity <- tree.imp[,imp.which][[2]] 
		result@Var <- colnames(tree.imp)[imp.which]
		result@type <- ifelse(is.factor(dataSet@X[,result@Var]),"Factor","Numeric") 
		# two-sample statistic or R(t)
		if(result@type == "Numeric") {
			result@Point <- as.numeric(tree.imp[,imp.which][[1]])
			result@Obs <- which(dataSet@X[, result@Var] <= result@Point)
		}
		else if(result@type == "Factor"){
			result@Set <- tree.imp[,imp.which][[2]]
			result@Obs <- which(dataSet@X[, result@Var] %in% result@Set)
		}
		return(result)
	}
}
