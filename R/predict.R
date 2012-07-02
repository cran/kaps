### predict functions for apss class
setGeneric("predict")
setMethod("predict","apss", function(object, newdata, type = c("predict", "apss")){
	## predict the ID for terminal subgroups 
	type <- match.arg(type)
	if(missing(newdata)) {
		where <- object@groupID
		newdata <- object@data
	}
	else where <- pred.apss(object@split.pt, object@formula, newdata)
	
	if( type == "predict"){
		#tmps <- range(object@data[,object@split.var])
		#match(, )
		#pts <- c(tmps[1], round(object@split.pt,2), tmps[2])
		result <- data.frame(Group = where, newdata)
		rownames(result) <- "The estimated group="
		#result$Range <- paste("(", "<X<=",")", sep = "")
		#colnames(result) <- c("Newdata", "Group", "Range")
		colnames(result) <- c("Group", colnames(newdata))
		return(result)
	}
	else if(type == "apss"){
		newdata$subgroups <- where
		f <- update(object@formula, . ~  subgroups)
		test.chi <- survdiff(f, data = newdata)$chisq
		test.adj.chi <- pairwise.test(where, data = newdata, formula = object@formula, rho = object@Options@rho, adj = object@Options@p.adjust.methods, splits = object@Options@splits, shortcut = object@Options@shortcut)[1]
		v <- length(object@split.pt)
		x.mean <- 1 - 2 / (9 * v)
		x.std <- sqrt(2 / (9 * v))
		WH <- (test.chi / v)^(1/3)
		t <-abs((WH - x.mean)/(x.std))
		return(pred.stat = c(test.chi, test.adj.chi, WH, t))
	}
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
	if(is.vector(gClass)) gClass <- t(gClass)
	where <- apply(gClass, 1, sum)
	where <- where + 1
	return(where)
}

### predict functions for Tree class
setMethod("predict","Tree", function(object, newdata){
	## predict the ID for terminal nodes
	if(missing(newdata)) where <- object@where
	else {
		if(!is.data.frame(newdata)) newdata <- as.data.frame(newdata)
		where <- pred.Tree(object@frame, object@formula, labels(object), newdata)
	}
	frame <- object@frame
	node <- as.numeric(format(frame$Node))
	yval <- as.numeric(format(frame$yval))
	pred <- yval[match(where, node)]
	return(pred)
})

pred.Tree <- function(frame, f, splits, newdata){
### find ID number for new data set 
## Input
# tree frame by train data
# formula with the S4 class "Formula"
# new covariates data with the type of data.frame
# where: ID
	if(nrow(frame) == 1L) return(where = rep(1, nrow(newdata)))
	else{
		nr <- nrow(newdata)
		rownames(newdata)<- 1:nr
		n <- format(frame$Node)
		node <- as.numeric(n)
		terminal <- node[frame$terminal == 1L]
		path.node <- path.Tree(frame, nodes = terminal, splits, print.it = FALSE)
		where <- rep(NA, nrow(newdata))
		attr(where, "names") <- 1:nr
		for(i in 1:length(terminal)){
			x <- path.node[[i]][-1]
			path.sel <- paste(x, collapse = " & ")
			tmp <- subset(newdata, eval(parse(text = path.sel)))
			rownames(tmp)
			where[match(rownames(tmp),names(where))] <- terminal[i]
		}
		return(where)
	}
}

path.Tree <- function(frame, nodes, splits, print.it = TRUE){
### find the path to go to a terminal node
### This function is based on the path.rpart on the rpart package
### We appreciate the authors who wrote the rpart and party packages
### Modifited by Soo-Heang Eo, 1st April 2012	
	frame <- frame
	n <- format(frame$Node)
	node <- as.numeric(n)
	which <- desc.mat(node)
	path <- list()
	
	if(missing(nodes)){
		for(i in 1:ncol(which)){
			path[[n[i]]] <- path.i <- splits[which[,i]]
			if(print.it){
				cat("\n","node number:",n[i],"\n")
				cat(paste("   ",path.i), sep ="\n")
			}
		}
	}
	else {
		if(length(nodes <- node.match(nodes, node)) == 0L)
			return(invisible())
		for(i in nodes){
			path[[n[i]]] <- path.i <- splits[which[,i]]
			if(print.it){
				cat("\n","node number:",n[i],"\n")
				cat(paste("   ",path.i), sep ="\n")
			}
		}
	}
	invisible(path)
}

desc.mat <- function (nodes, include = TRUE) {
	## Borrowed from descendants() in the rpart package!
    n <- length(nodes)
    if (n == 1) 
        return(matrix(TRUE, 1, 1))
    ind <- 1:n
    desc <- matrix(FALSE, n, n)
    if (include) 
        diag(desc) <- TRUE
    parents <- match((nodes%/%2), nodes)
    lev <- floor(log(nodes, base = 2))
    desc[1, 2:n] <- TRUE
    for (i in max(lev):2) {
        desc[cbind(ind[parents[lev == i]], ind[lev == i])] <- TRUE
        parents[lev == i] <- parents[parents[lev == i]]
        lev[lev == i] <- i - 1
    }
    return(desc)
}

node.match <- function(nodes, nodelist, leaves, print.it = TRUE){
	## Borrowed from node.match() in the rpart package!
    node.index <- match(nodes, nodelist, nomatch = 0)
    bad <- nodes[node.index == 0]
    if (length(bad) > 0 & print.it) 
        warning("supplied nodes ", paste(bad, collapse = ","), 
            " are not in this tree")
    good <- nodes[node.index > 0]
    if (!missing(leaves) && any(leaves <- leaves[node.index])) {
        warning("supplied nodes ", paste(good[leaves], collapse = ","), 
            " are leaves")
        node.index[node.index > 0][!leaves]
    }
    else node.index[node.index > 0]
}
