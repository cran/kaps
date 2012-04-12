### Set S4 class
## Top-level functions for lrtree algorithms
lrtree <- function(Formula, data, subset = NULL, weights = NULL, ...) {
#########################################################
#####           pre-processing step
#########################################################
	options(warn = -1)
	controls = lrtree.control(...)
	#controls = lrtree.control()
	n <- nrow(data) # total number of observations concerning with time 
	rownames(data) <- 1:n
    if(n == 0L) stop("0 (non-NA) cases.")
	# original formula
    mf <- match.call() # save a copy of the call
	Formula <- as.Formula(Formula)
	#mf[[1]] <- as.name("model.frame")
	#mf$formula <- f
	#mf <- eval(mf, parent.frame())
	vars <- all.vars(Formula)
	
	if(!is.data.frame(data)) data <- as.data.frame(data)
	## Set the object of result
    result <- new("Tree")
	result@call <- mf
	result@formula <- Formula
	result@data <- data

######################################
#####           Model Fitting
######################################
## select fitting model
## Missing handling    
	dataSet <- new("dataset")
	dataSet@Y <- Surv(data[,vars[1]], data[,vars[2]])
	dataSet@X <- model.part(Formula, data = data, rhs = 1)
	dataSet@Z <- model.part(Formula, data = data, rhs = 2)

	where <- rep(1, nrow(dataSet@Y)) 
	attr(where,"names") <- rownames(dataSet@X)
	cat("Now, recursive partitioning for your data. Plz, wait a minute^^. \n")
	frameID <- growing.lrtree(Formula, dataSet, frame = NULL, nodenum = 1, controls, where)
	result@frame <- as.data.frame(frameID$frame)
	result@where <- frameID$ID
    result@controls <- controls
    result@weights <- weights

		
######################################
#####   for all of the pruning 
######################################

    candidate.tree <- candidate.prune(result)
    result@alpha <- candidate.tree$alpha
 	result@alpha.prime <- sqrt(cumprod(result@alpha))
	result@alpha.prime[1] <- 0
    result@frame$complexity <- as.numeric(0)
    result@frame$complexity[pmatch(candidate.tree$node.change, result@frame$Node)]<- round(result@alpha,3)
    return(result)
}
