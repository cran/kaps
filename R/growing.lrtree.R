growing.lrtree <- function(Formula, dataSet, frame = NULL, nodenum = 1, controls, where){
#####################################################################################
### Main function for recursive partitioning
### Growing trees until at least one stopping rule is met
### Stopping rule
### 1. the number of sample in a node t (minsplit & minbucket)
### 2. statistic criterion (eps)
### 1&2 are requsted in a function split.lrtree()
### 3. the maximum depth (maxdepth)
###
#####################################################################################

	if(controls@VarCtrl@random){
		# For frailty models in the near future^^
		id.var <- all.vars(Formula)[length(all.vars(Formula))]
		id.unique <- unique(dataSet@Z[,id.var])
		duplicated.data <- !duplicated(dataSet@Z[,id.var])
	} 
	else{
		duplicated.data <- rep(TRUE, nrow(dataSet@Y))
	}
	tmp  <- coxph(dataSet@Y ~ 1) # fit a constant cox ph model.
	tmp1 <- survfit(tmp)
	
	if( ( length(duplicated.data[duplicated.data]) <= controls@SplitCtrl@minsplit ) | 
		(nodenum >= 2^(controls@TGCtrl@maxdepth-1) & nodenum < 2^controls@TGCtrl@maxdepth )){
		## Stop growing due to the sample size or the tree depth
		result <- list()
		result$frame <- rbind(frame, c(
			Node = nodenum , 
			var = "<Leaf>", 
			point = 999, 
			set = "", 
			n = tmp$n, 
			nLogLik = round(tmp$loglik * (-1),2), 
			yval = survmed(tmp1$surv, tmp1$time), 
			terminal = 1)
		)

		where[match(rownames(dataSet@X),names(where))] <- nodenum
		result$ID <- where

		return(result)
	}# No split condition END
	
	## Binary Splits for Censored Data
	spt.info <- split.lrtree(dataSet, controls, duplicated.data )
	#names(spt.info@Impurity) <- NULL		

	if( spt.info@Impurity <= controls@VarCtrl@eps | spt.info@Var == "" ){
		# stop growing tree structure
		result <- list()
		result$frame <- rbind(frame, c(
			Node = nodenum , 
			var = "<Leaf>", 
			point = 999, 
			set = "", 
			n = tmp$n, 
			nLogLik = round(tmp$loglik * (-1),2), 
			yval = survmed(tmp1$surv, tmp1$time), 
			terminal = 1)
		)

		where[match(rownames(dataSet@X),names(where))] <- nodenum
		result$ID <- where
		return(result)
	}
	#else if( spt.info@WH.stat <= controls@SplitCtrl@WH.strd | spt.info@Impurity <= controls@VarCtrl@eps ){
		# use Multi-Level stopping rule
	#}
	else {
		# still grow tree size
		result <- list()
		result$frame <- rbind(frame, c(
			Node = nodenum , 
			var = spt.info@Var, 
			point = round(spt.info@Point,3), 
			set = spt.info@Set, 
			n = tmp$n, 
			nLogLik = round(tmp$loglik * (-1),2), 
			yval = survmed(tmp1$surv, tmp1$time), 
			terminal = 0)
		)

		where[match(rownames(dataSet@X),names(where))] <- nodenum
		result$ID <- where
	 
		leftnode <- new("dataset")
		leftnode@Y <- dataSet@Y[spt.info@Obs,]
		leftnode@X <- dataSet@X[spt.info@Obs,, drop = FALSE]
		leftnode@Z <- dataSet@Z[spt.info@Obs,, drop = FALSE]
		
		result <- growing.lrtree(Formula, leftnode, frame = result$frame, 
			nodenum = (2*nodenum), controls = controls, where = result$ID)

		rightnode <- new("dataset")
		rightnode@Y <- dataSet@Y[-spt.info@Obs,]
		rightnode@X <- dataSet@X[-spt.info@Obs,,drop = FALSE]
		rightnode@Z <- dataSet@Z[-spt.info@Obs,,drop = FALSE]
		result <- growing.lrtree(Formula, rightnode, frame = result$frame, 
			nodenum = (2*nodenum + 1), controls = controls, where = result$ID)
		return(result)
	}
}

