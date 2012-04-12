subtree <- function(frame, k, terminal.node,res = NULL)
{
    if(any(terminal.node == k)){
        res <- rbind(res, frame[frame$Node == k,c("Node","nLogLik","terminal")])
        return(res)
    }
    else{
        res <- rbind(res, frame[frame$Node == k,c("Node","nLogLik","terminal")])
        res <- subtree(frame, 2 * k, terminal.node,res)
        res <- subtree(frame, (2*k)+1, terminal.node,res)
        return(res)
    }
}

####################################################################################################################
cp.prune <- function(frame, res.view = FALSE){
    #Step 1 : obtaining complexity parameter
    frame$complexity <- 0
    frame[,1] <- as.numeric(format(frame[,1]))
    frame[,6] <- as.numeric(format(frame[,6]))
    max.node <- max(frame$Node)
    terminal.node <- frame[frame$terminal == 1L,1]
    int.node <- frame[frame$terminal == 0L,1]
    T.num = R.t = R.T = alpha <- NA
    
    for(i in 1:length(int.node)){
        subtree.tmp <- subtree(frame, k= int.node[i], terminal.node)
        R.t[i] <- subtree.tmp[1,2]
        T.num[i] <- length(subtree.tmp[subtree.tmp$terminal == 1L,3])            
        R.T[i] <- sum(subtree.tmp[subtree.tmp$terminal == 1L,2])  
        alpha[i] <- (R.t[i] - R.T[i]) / (T.num -1)
        frame$complexity[frame$Node == int.node[i]] <- alpha[i]
        rm(subtree.tmp)
    }
    
    if(res.view) {
        result <- data.frame( node = int.node,
                              R.t = R.t,
                              R.T = R.T,
                              ter.num = T.num,
                              cp = alpha)
    }
    else{
        result <- list()
        result$alpha <- min(alpha)
        result$change <- int.node[alpha ==  result$alpha]
        pruned.node <- int.node[alpha == result$alpha]
        subtree.node <- subtree(frame, k = pruned.node, terminal.node)[-1,1]
        result$frame <- frame[is.na(match(frame$Node, subtree.node)),]
        result$frame[result$frame[,1] == pruned.node,8] <- 1
        result$frame[result$frame[,1] == pruned.node,2] <- "<Leaf>"
        rownames(result$frame) <- 1:length(result$frame[,1])
    }
    return(result)
}
####################################################################################################################
candidate.prune <- function(treeres, ...){
    if(nrow(treeres@frame) == 1L){
        result <- list(alpha = NA, node.change = NA, sub.tree = NA)
        return(result)
    }
    else{
        frame <- treeres@frame   
        candidate.tree <- list()
        candidate.tree[[1]] <- treeres@frame
        
        frame.tmp <- frame
        i <- 0
        alpha = node.change <- NA
        repeat{
            i <- i+1
            tmp <- cp.prune(frame.tmp)
            candi.tmp <- tmp$frame
            alpha[i] <- tmp$alpha
            node.change[i] <- tmp$change
            candidate.tree[[i+1]] <- candi.tmp
            frame.tmp <- candi.tmp
            if(all(frame.tmp$terminal ==1)) break
        }
        result <- list(alpha = alpha, node.change = node.change, sub.tree = candidate.tree)
        return(result)
    }
}
####################################################################################################################
alpha.prune <- function(treeres, alpha){
	if(nrow(treeres@frame) == 1L) return(treeres@frame)
	else{
		frame <- treeres@frame
		frame[,1] <- as.numeric(format(frame[,1]))
		frame[,6] <- as.numeric(format(frame[,6]))
		int.node <- frame[frame$terminal == 0L,1]
		ter.node <- frame[frame$terminal == 1L,1]
		alpha.int <- treeres@frame$complexity[frame$Node %in% int.node] 
		pruned.node <- int.node[alpha.int <= alpha]
		if( length(pruned.node) == 0L) return(frame)
		else{
			subtree.tmp <- lapply(pruned.node, subtree, frame = frame, terminal.node = ter.node)
			subtree.node <- lapply(subtree.tmp, function(x) x[-1,1])
			subtree.node <- unique(unlist(subtree.node))
			frame[frame[,1] %in% pruned.node,8] <- 1
			frame[frame[,1] %in% pruned.node,2] <- "<Leaf>"
			result <- frame[is.na(match(frame$Node, subtree.node)),]
			rownames(result) <- 1:length(result[,1])
			return(result)
		}
	}
}
####################################################################################################################
## create a new generic function, with a default method
prune <- function(treeres, ...) attributes(treeres, ...)
setGeneric("prune")
setMethod("prune","Tree", function(treeres, K = 10,
	cost = function(y, yhat, status) sum(abs(y - yhat) * status)
	){
    # This function is modified function of both cv.glm in boot package and xpred.rpart in rpart package
    # by Soo-heang Eo in March, 29, 2012. 
    # Thanks for the authors who write cv.glm and xpred.rpart function.
		
##################################################
####     initial part for general data set
##################################################
		
   		if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
    	    runif(1)
	    }
	    seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
	    var.list <- all.vars(treeres@formula)
	    resp.var <- var.list[1]
	    status <- var.list[2]
	    alpha <- treeres@alpha.prime
		nc <- length(alpha)
	    frame <- format(treeres@frame)
	    for(i in c(1,3,5,6,7,9)) frame[,i] <- as.numeric(frame[,i])
	    pruned.table <- matrix(NA, nrow = nc, ncol = 5)
	    colnames(pruned.table) <- c("Tree", "node", "T.tilda", "alpha", "Rcv(T)")
		n <- frame$n[1]
		
		K.o <- K
		K <- trunc(K)
		out <- NULL
		
		if(K <2) stop("K should be greater than or equal to 2.")
		if(K >n) stop("K should be less than or equal to the number of observations.")
		
		kvals <- unique(round(n / (1:floor(n/2))))
		temp <- abs(kvals - K)
    	if (!any(temp == 0)) K <- kvals[temp == min(temp)][1]
	    if (K != K.o) warning("K has been set to ", K)
    
	    ff <- ceiling(n/K)
    	s <- sample(rep(1:K, ff), n)
	    n.s <- table(s)
    	ms <- max(s)
	    #names(s) <- unique(treeres@where)
		
	    cost.0 <- cost(y = treeres@data[,resp.var],
	    	 yhat= frame$yval[match(treeres@where, frame$Node)],
	    	 status = treeres@data[,status])
		R.CV <- matrix(NA, nrow = K, ncol = nc)
   		 
   		for(i in 1:ms){
   			j.out <- c(1:n)[s == i]
   			j.in <- c(1:n)[s != i]
   			
   			# define a new learning
   			cat("K =", i,"\n")
			new.learning <- treeres@data[j.in,, drop = FALSE]
   			new.test <- treeres@data[j.out, , drop = FALSE]
   			
   			# obtain the maximal tree using a new training data set
			Call <- treeres@call
			Call$data <- as.name("new.learning")
			train.tree <- eval(Call) 
		    
		    for(j in 1:nc){
		    	train.tree@frame <- alpha.prune(train.tree, alpha[j])
		    	fitted.time <- predict(train.tree, newdata = new.test)
				# must treat NA value.
				na.pos <- which(!is.na(fitted.time))
				T.tilda <- sum(train.tree@frame$terminal == 1)
		    	cost.tmp <- cost(y = new.test[na.pos,resp.var], yhat = fitted.time[na.pos],
		    		status = new.test[na.pos,status])
				R.CV[i,j] <- cost.tmp + alpha[j] * T.tilda
			}
   		}
		sel.cv <- apply(R.CV, 2, sum)
		sel.cv <- which(sel.cv == min(sel.cv))
		## pruned tree by CV
		result <- new("Tree")
		result@call <- treeres@call
		result@formula <- treeres@formula
		result@data <- treeres@data
		
		result@frame <- alpha.prune(treeres, alpha[sel.cv])
		result@where <- pred.Tree(result@frame, result@formula, labels(result), result@data)
		
		result@controls <- treeres@controls
		result@weights <- treeres@weights
		result@alpha.prime <- alpha[sel.cv]
		return(result)
	}
)
        
