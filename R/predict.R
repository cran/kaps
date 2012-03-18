### predict functions for pass class
setGeneric("predict")
setMethod("predict","apss", function(object, newdata){
	## predict the ID for terminal subgroups 
	if(missing(newdata)) where <- object@where
	else where <- pred.apss(object@split.pt, object@formula, newdata)
	
	newdata$subgroups <- where
	f <- update(object@formula, . ~  subgroups)
	test.chi <- survdiff(f, data = newdata)$chisq
	test.adj.chi <- pairwise.test(where, data = newdata, formula = object@formula, rho = object@Options@rho, adj = object@Options@p.adjust.methods, splits = object@Options@splits)[1]
	v <- length(object@split.pt)
	x.mean <- 1 - 2 / (9 * v)
	x.std <- sqrt(2 / (9 * v))
	WH <- (test.chi / v)^(1/3)
	t <-(WH - x.mean)/(x.std)
	if(t <=0) t <- 0
	return(pred.stat = c(test.chi, test.adj.chi, WH, t))
	}
)

pred.apss <- function(split.pt, f, newdata){
#### find ID number for new data set
##Input
# split points by train data
# new covariates data with the type of data.frame
## output
# ID number (the observations that assins each terminal group)
## FIXME by more efficient way in the near future
	X <- model.part(f, data = newdata, rhs = 1, drop = FALSE)
	nc <- length(split.pt)
	gClass <- matrix(NA, ncol = nc, nrow = nrow(newdata))
	gClass <- sapply(split.pt, function(x,y) y > x, y = X)
	#for(i in 1:nc) gClass[,i] <- X > split.pt[i]
	where <- apply(gClass, 1, sum)
	where <- where + 1
	return(where)
}

