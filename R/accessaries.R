kapsNews <- function(){
	file.locate <- file.path(system.file(package = "kaps"), "NEWS")
	file.show(file.locate)
}

## Count minimum sample size
count.mindat <- function(formula, data, part = 10){
## output: the minimum number of samle size.

	formula <- as.Formula(formula)
	X <- model.part(formula, data, rhs = 1)
	subgr <- apply(X, 2, function(x) {
		res <- table(x) / length(x)
		res <- res[order(res, decreasing = FALSE)]
		return(res)
		}
	    )
	res <- apply(subgr, 2, function(x, data, part){
		dat.prop <- x[length(x) - part]
		res <- floor(nrow(data) * dat.prop)
		return(res)}, 
		data =  data, part = part
		)
	names(res) <- NULL
	return(res)
}

## Calculate all possible subset
# from gRbase package
combnat <- function(x,m){
  if (length(x)==1 && is.numeric(x))
    x <- seq(x)
  if (length(x) < m)
    stop("Error in combnat: n < m\n")
  NCAND = as.integer(length(x))
  NSEL  = as.integer(m)
  NSET <- as.integer(choose(NCAND,NSEL))
  ANS  <- as.integer(rep.int(0, NSET*NSEL))
  res <- .C("combnC", NSEL, NCAND, NSET, ANS, DUP=FALSE
            ,PACKAGE="kaps" )[[4]]
  res <- x[res]
  dim(res) <- c(NSEL, NSET)  
  return(res)
}

## Calculate Chi-square and Adj. chi-square test statistics among subgroups.
adj.test <- function(x){
	data <- x@data
	data$subgroups <- x@where
	#data$group <- gr.sel[,index]
	f <- update(x@formula, . ~ subgroups)		
	test.stat <- survdiff(f, data = data)$chisq
	v <- length(x@split.pt)
	x.mean <- 1 - (2 / (9 * v))
	x.std <- sqrt(2 / (9 * v))
	WH <- (test.stat / v)^(1/3)
	t <-(WH-x.mean)/(x.std)
	#WH <- ( (7/9) + sqrt(nu)*((test.stat/ nu )^(1/3) - 1 + (2 / (9 * nu))))^3
	#WH <- max(0, WH)
	# return(c(Chisq = test.stat, adj.Chisq = WH))
	return(c(Chisq = test.stat, adj.Chisq = x@Chisq, WH = WH, t = t))
}

## Median Survival Time, refered to the print.survfit() in the survival package
survmed <- function(surv,time, tol= 1.0e-5) {
	keep <- (!is.na(surv) & surv <(.5 + tol))
	if (!any(keep)) NA
	else {
		time <- time[keep]
		surv <- surv[keep]
		if (abs(surv[1]-.5) <tol  && any(surv< surv[1])) 
			(time[1] + time[min(which(surv<surv[1]))])/2
		else time[1]
	}
}

## summary functions for adaptive staging algorithms
surv.yrs <- function(pt, surv, time){
	if(any(time == pt)) return(surv[time == pt])
	else {
		mod <- min(abs(time - pt))
		if(mod > 7) {
			return(0)
		}
		if(any(time == (pt + mod))) return(surv[time ==(pt + mod)])
		else return(surv[time ==(pt - mod)])
	}
}

## summary functions for apss package
setGeneric("summary")
setMethod("summary","apss", function(object,K){
	if(!missing(K)) object <- object@candid[[which(object@groups == K)]]

	data <- object@data
	f <- update(object@formula, . ~ 1)
	surv.root <- survfit(f, data = data)
	rootS <- summary(surv.root)
	
	data$Group <- object@where 
	f <- update(object@formula, . ~ Group)
	surv.all <- survfit(f, data = data)
	objS <- summary(surv.all)
	 
	subgr.surv <- list()
	subgr.time <- list()
	level.objS <- levels(objS$strata)
	for(i in 1:length(level.objS)) {
		subgr.surv[[i]] <- objS$surv[objS$strata == level.objS[i]] 
		subgr.time[[i]] <- objS$time[objS$strata == level.objS[i]] 
	}
	gr.med <- objS$table[,"median"]
	root.med <-  rootS$table["median"]
	root.1 <- surv.yrs(pt = 12, surv = rootS$surv, time = rootS$time)
	root.3 <- surv.yrs(pt = 36, surv = rootS$surv, time = rootS$time)
	root.5 <- surv.yrs(pt = 60, surv = rootS$surv, time = rootS$time)
	root <- round(c(nrow(data), root.med, root.1, root.3, root.5), 3)
	names(root) <- NULL
	gr.1yrs <- mapply( surv.yrs, surv = subgr.surv, time = subgr.time, MoreArgs = list(pt = 12))
	gr.3yrs <- mapply( surv.yrs, surv = subgr.surv, time = subgr.time, MoreArgs = list(pt = 36))
	gr.5yrs <- mapply( surv.yrs, surv = subgr.surv, time = subgr.time, MoreArgs = list(pt = 60))
	obs <- table(object@where)
	obs <- as.vector(obs)
	names(obs) <- names(gr.med)
	res <- data.frame( 
		N = obs,
		Med = gr.med,
		yrs.1 = round(gr.1yrs,3),
		yrs.3 = round(gr.3yrs,3),
		yrs.5 = round(gr.5yrs,3)
	)
	rownames(res) <- names(gr.med)
	res <- res[order(res$Med, decreasing = TRUE),]
	res <- rbind(root = root, res)
	return(res)
	}
)

setMethod("summary","Tree", function(object){
	data <- object@data
	f <- update(object@formula, . ~ 1)
	surv.root <- survfit(f, data = data)
	rootS <- summary(surv.root)
	
	data$Node <- object@where 
	f <- update(object@formula, . ~ Node)
	surv.all <- survfit(f, data = data)
	objS <- summary(surv.all)
	 
	 
	subgr.surv <- list()
	subgr.time <- list()
	level.objS <- levels(objS$strata)
	for(i in 1:length(level.objS)) {
		subgr.surv[[i]] <- objS$surv[objS$strata == level.objS[i]] 
		subgr.time[[i]] <- objS$time[objS$strata == level.objS[i]] 
	}
	gr.med <- objS$table[,"median"]
	root.med <-  rootS$table["median"]
	root.1 <- surv.yrs(pt = 12, surv = rootS$surv, time = rootS$time)
	root.3 <- surv.yrs(pt = 36, surv = rootS$surv, time = rootS$time)
	root.5 <- surv.yrs(pt = 60, surv = rootS$surv, time = rootS$time)
	root <- round(c(nrow(data), root.med, root.1, root.3, root.5), 3)
	names(root) <- NULL
	gr.1yrs <- mapply( surv.yrs, surv = subgr.surv, time = subgr.time, MoreArgs = list(pt = 12))
	gr.3yrs <- mapply( surv.yrs, surv = subgr.surv, time = subgr.time, MoreArgs = list(pt = 36))
	gr.5yrs <- mapply( surv.yrs, surv = subgr.surv, time = subgr.time, MoreArgs = list(pt = 60))
	obs <- table(object@where)
	obs <- as.vector(obs)
	names(obs) <- names(gr.med)
	res <- data.frame( 
		N = obs,
		Med = gr.med,
		yrs.1 = round(gr.1yrs,3),
		yrs.3 = round(gr.3yrs,3),
		yrs.5 = round(gr.5yrs,3)
	)
	rownames(res) <- names(gr.med)
	res <- res[order(res$Med, decreasing = TRUE),]
	res <- rbind(root = root, res)
	return(res)
	}
)

####################################################################
## Labels for tree diagram and plot
setGeneric("labels")
setMethod("labels","Tree", function(object, type = "diagram") {
	
	ff <- object@frame
	ff <- format(ff)
	n <- nrow(ff)
	
	res <- "Root"
	names(res) <- 1
    if (n==1) return(res)  #special case of no splits
	
    is.leaf <- (ff$var == "<Leaf>")
    whichrow <- !is.leaf

    index <- cumsum(c(1, 1*(!is.leaf)))
    irow  <- index[c(whichrow, FALSE)]     #we only care about the primary split

    if(type == "diagram"){
	lsplit <- rsplit <- vector(mode='character', length= length(irow))
	lsplit <- ifelse(ff[whichrow, 4] == "", 
			paste(ff[whichrow, 2],  "<=", ff[whichrow, 3],sep = ""),
			paste(ff[whichrow, 2],  "==", ff[whichrow, 4],sep = "")
		)	
	names(lsplit) <- as.numeric(ff[whichrow, 1]) * 2
	rsplit <- ifelse(ff[whichrow, 4] == "", 
		paste(ff[whichrow, 2],  "> ", ff[whichrow, 3],sep = ""),
		paste(ff[whichrow, 2],  "!=", ff[whichrow, 4],sep = "")
	)	
	names(rsplit) <- as.numeric(ff[whichrow, 1])  * 2 + 1
	}
	else if(type == "plot"){
		lsplit <- rsplit <- vector(mode='character', length= length(irow))
		lsplit <- ifelse(ff[whichrow, 4] == "", 
				paste(  "<=", ff[whichrow, 3],sep = ""),
				paste(  "==", ff[whichrow, 4],sep = "")
			)	
		names(lsplit) <- as.numeric(ff[whichrow, 1]) * 2
		rsplit <- ifelse(ff[whichrow, 4] == "", 
			paste( "> ", ff[whichrow, 3],sep = ""),
			paste("!=", ff[whichrow, 4],sep = "")
		)	
		names(rsplit) <- as.numeric(ff[whichrow, 1])  * 2 + 1
	}
	else if(type == "var"){
		lsplit <- rsplit <- vector(mode='character', length= length(irow))
		lsplit <- ff[whichrow, 2]
		names(lsplit) <- as.numeric(ff[whichrow, 1]) * 2
		rsplit <- ff[whichrow, 2]
		names(rsplit) <- as.numeric(ff[whichrow, 1])  * 2 + 1
	}
	res <- c(res, lsplit, rsplit)
	res <- res[pmatch(as.numeric(ff$Node),as.numeric(names(res)))]
	return(res)
	}
)

####################################################################
## print furnctions for an apss object used in apss()
setGeneric("print")
setMethod("print","apss", function(x,K){
		if(!missing(K)) x <- x@candid[[which(x@groups == K)]]
		x@Options@fold <- FALSE
		show(x)
	}
)

setGeneric("show")
setMethod("show", "apss", function(object){
	cat("Call:\n")
	print(object@call)
	cat("\n")
	cat("	      K-Adaptive Partitioning for Survival Data\n\n")
	cat(" Sample: ", nrow(object@data), "\n")
	
	pts <- round(object@split.pt,2)
	cutpoint <- matrix(object@split.pt, nrow = 1, ncol = length(pts))
	rownames(cutpoint) <- object@split.var
	colnames(cutpoint) <- c("1st pt", "2nd pt", "3rd pt", paste(4:20,"th pt", sep = ""))[1:length(object@split.pt)]
	Signif <- symnum(object@pvalue, corr = FALSE, na = FALSE, 
                  cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1), 
                  symbols = c("***", "**", "*", ".", " "))

	cat("\n")
	if(length(object@groups) >= 2){
		test.stat <- data.frame(
			Z = object@Z,
			X2 = object@Chisq,
			df = as.integer(object@groups) - 1,
			WH = object@WH,
			t = object@t)
		test.stat<-	zapsmall(test.stat,digits = 3)
		test.stat <- cbind(Z = test.stat$Z, 
			Zp = round(1 - pchisq(q = object@Z,df=1),3),
			test.stat[,-1],
			pvalue = round(object@pvalue,3),
			sig = format(Signif))
		dimnames(test.stat) <- list(attr(object@groups,"names"), c("Z","Pr(>|Z|)","x2(s)","df", "WH","t value","Pr(>|t|)", ""))
		print(test.stat)
		cat("---\nSignif. codes: ", attr(Signif, "legend"), "\n")
		if(object@Options@fold){
			cat("\n By Cross-Validation\n")
			cv.pvalue <- round(1 - pnorm(q = object@elbow[4,]),3)
			Signif <- symnum(cv.pvalue, corr = FALSE, na = FALSE, 
			  cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1), 
			  symbols = c("***", "**", "*", ".", " "))
			cv.stat <- data.frame(
				Z = object@elbow[2,],
				X2 = object@elbow[1,],
				df = as.integer(object@groups) - 1,
				WH = object@elbow[3,],
				t = object@elbow[4,])
			cv.stat<- zapsmall(cv.stat,digits = 3)
			cv.stat <- cbind(cv.stat[,1], 
				Zp = round(1 - pchisq(q = object@elbow[2,],df=1),3),
				cv.stat[,-1],
				pvalue = cv.pvalue,sig = format(Signif))
			dimnames(cv.stat) <- list(attr(object@groups,"names"), c("Z","Pr(>|Z|)","x2(s)","df", "WH","t value","Pr(>|t|)", ""))
			print(cv.stat)
			cat("---\nSignif. codes: ", attr(Signif, "legend"), "\n")

		}
		cat("\n Estimated Cut-off Points for K =",object@groups[object@index],":\n")
	}
	else cat("\n Estimated Cut-off Points for K =",object@groups,":\n")
	print(round(cutpoint,2))

	cat("\n Pairwise comparisons using log-rank tests \n")
	data <- object@data
	data$group <- object@where
	pair <- combn(unique(object@where),2)
	f <- update(object@formula, . ~ group)
	res <- c()
	for(i in 1:ncol(pair)){
		data.tmp <- data[data$group %in% pair[,i],]
		tmp <- survdiff(f, data = data.tmp, rho = object@Options@rho)
		res[i] <- 1 - pchisq(tmp$chisq, 1)
	}
	# Adjsut P-values for Multiple Comparisons 
	res <- p.adjust(res, method = object@Options@p.adjust.methods)
	pair.mat <- matrix(NA, ncol = length(object@split.pt), nrow = length(pts))
	pair.mat[!upper.tri(pair.mat)] <- res
	pair.mat <- round(pair.mat, 5)
	## adjust split points
	tmps <- range(data[,object@split.var])
	pts <- c(tmps[1], pts, tmps[2])
	n.pts <- length(pts)
	colnames(pair.mat) <- paste(pts[-(n.pts:(n.pts-1))],"<x<=", sep = "")
	colnames(pair.mat) <- paste(colnames(pair.mat),pts[-c(1,length(pts))], sep = "")
	rownames(pair.mat) <- paste(pts[-c(1,length(pts))],"<x<=", sep = "")
	rownames(pair.mat) <- paste(rownames(pair.mat), pts[-(1:2)],sep="")
	pair.mat <- ifelse(pair.mat < 1.0e-6, "<.0000", pair.mat)
	pair.mat <- ifelse(!is.na(pair.mat), pair.mat, "-")
	pair.mat <- as.data.frame(pair.mat)
	print(pair.mat)
	cat("\n P-value adjustment method:" ,object@Options@p.adjust.methods, "\n")
	invisible(object)
	}
)

setMethod("show", "Tree", function(object){
	options(warn = -1)
	spaces = 3
	frame <- object@frame
	var <- as.character(frame$var)
    n <- frame$n
    node <- as.numeric(format(frame$Node))
    depth <- tree.depth(node)
    indent <- paste(rep(" ", spaces * 32), collapse = "")
    if (length(node) > 1) {
        indent <- substring(indent, 1, spaces * seq(depth))
        indent <- paste(c("", indent[depth]), format(node), ")", sep = "")
    }
    else indent <- paste(format(node), ")", sep = "")
	
    term <- rep(" ", length(depth))
    term[frame$terminal == 1L] <- "*"
    z <- labels(object)
    z <- paste(indent, z, n, frame$nLogLik, frame$yval, term)

	cat("Call:\n")
	print(object@call)
	cat("\n")
	cat("Tree-structured Models for Censored Data\n\n")
	cat("Sample: ", format(n[1]),"          Methods: ", object@controls@SplitCtrl@fitted.model, "\n\n", sep = "")
        
    #This is stolen, unabashedly, from print.tree
    #if (x$method=="class")
    #     cat("node), split, n, loss, yval, (yprob)\n")
    #else cat("node), split, n, deviance, yval\n")
    cat("node), split, n, Loss, yval\n")
    cat("                * denotes terminal node\n\n")
    cat(z, sep = "\n")
    return(invisible(object))
	}
)

####################################################################
### panel function for Kaplan-Meier curves in nodes

## A standard tree plot with KM curves at the terminal nodes
## Soo-Heang Eo
setGeneric("plot")
setMethod("plot","apss",function(x, K, ...){
		if(!missing(K) | !x@Options@fold) {
			if (!missing(K)) x <- x@candid[[which(x@groups == K)]]
			else {
				par(mfrow = c(1,2))

				fvars <- all.vars(x@formula)
				plot(x@data[,fvars[3]], x@data[,fvars[1]], pch = c(1,3)[x@data[,fvars[2]]+1] , 
				axes=FALSE, col= ifelse(x@data[,fvars[2]] == 1, "red2","blue"), xlab="", ylab="")
				mtext(fvars[3], side=1, line=3, cex=1.2)
				mtext("Survival months", side=2, line=2, cex=1.2)
				axis(1)
				axis(2, at= c(12,24,48,72,96,120, 144, 168, 192, 216, 240), labels=c(12,24,48,72,96,120, 144, 168, 192, 216, 240))
				legend("topright",c("Event","Censored"), col=c("blue","red2"), cex=1, pch=c(1,3), bty = "n")
				r <- locfit.censor(x = x@data[,fvars[3]], y = x@data[,fvars[1]], cens= (1- x@data[,fvars[2]])) 
				lines(r, lwd=2, lty=1)
			}
			km.curve(x)
			legend("topright", legend = paste("G", 1:(length(x@split.pt)+1),sep = ""), 
				bty = "n", col = unique(x@where), lty = unique(x@where))
			return(invisible(x))
		}
		else K <- x@groups[x@index]
		
		if( (length(x@groups) == 1) ){
			km.curve(x, main = paste("# of subgroups = ",K))
			legend("topright", legend = paste("G", 1:(length(x@split.pt)+1),sep = ""), 
				bty = "n", col = unique(x@where), lty = unique(x@where))
			return(invisible(x))
		}
	
		fvars <- all.vars(x@formula)
		par(mfrow = c(2,2))
		
		plot(x@data[,fvars[3]], x@data[,fvars[1]], pch = c(1,3)[x@data[,fvars[2]]+1] , 
			axes=FALSE, col= ifelse(x@data[,fvars[2]] == 1, "red2","blue"), xlab="", ylab="")
		mtext(fvars[3], side=1, line=3, cex=1.2)
		mtext("Survival months", side=2, line=2, cex=1.2)
		axis(1)
		axis(2, at= c(12,24,48,72,96,120, 144, 168, 192, 216, 240), labels=c(12,24,48,72,96,120, 144, 168, 192, 216, 240))
		legend("topright",c("Event","Censored"), col=c("blue","red2"), cex=1, pch=c(1,3), bty = "n")
		r <- locfit.censor(x = x@data[,fvars[3]], y = x@data[,fvars[1]], cens= (1-x@data[,fvars[2]])) 
		lines(r, lwd=2, lty=1)
		
		km.curve(x)
		legend("topright", legend = paste("G", 1:(length(x@split.pt)+1),sep = ""), 
			bty = "n", col = unique(x@where), lty = unique(x@where))

		plot(x@elbow[2,], type = "b", col = "blue", xlab = "", ylab = "",axes = FALSE, ylim = c(0,max(x@elbow[2,])),...)
		mtext("# of subgroups (K)", side=1, line=3, cex=1.2)
		mtext(expression(bar(Z)^cv), side=2, line=2, cex=1.2)
		abline(h = 3.84, col = "gray", lty = 2)
		axis(side = 1, at = 1:ncol(x@elbow), labels =  x@groups)
		axis(side = 2)

		plot(x@elbow[4,], type = "b", col = colors()[630], xlab = "", ylab = "", ylim = c(0,ceiling(max(x@elbow[4,]))), axes = FALSE,...)
		mtext("# of subgroups (K)", side=1, line=3, cex=1.2)
		mtext(expression(bar(t)^cv), side=2, line=2, cex=1.2)
		abline(h = 1.64, col = "gray", lty = 2)
		axis(side = 1, at = 1:ncol(x@elbow), labels =  x@groups)
		axis(side = 2)

	}
)

## plot Kaplan-Meire survival curves for termninal nodes
km.curve <- function(object, 
	x.lab = c(0,24,48,72,96,120, 144, 168, 192, 216, 240), ...){
	object@data$where <- object@where 
	f <- update(object@formula, . ~ where)
	surv.all <- survfit(f, data = object@data)
	id.n <- length(unique(object@where))
	plot(surv.all, col= 1:id.n, lty = 1:id.n, axes=FALSE, cex=1, lwd=1.5, ...)
	mtext("Survival months", side=1, line=3, cex=1.2)
	mtext("Survival probability", side=2, line=3, cex=1.2)
	axis(2)
	axis(1, at= x.lab, labels=x.lab)
}

setMethod("plot", "Tree", 
	function(x, uniform=TRUE, branch= .2, margin= .1, newpage = TRUE, type = "naive", 
		compress = TRUE, nspace, minbranch = .5, ...) {
		if(type == "grid"){
			require(grid)
			### total number of terminal nodes
			nx <- sum(as.numeric(format(x@frame$terminal)))
			### maximal depth of the tree
			ny <- tree.depth(node)
			ny <- max(ny) +1

			## setup newpage
			if (newpage) grid.newpage()

			## setup root viewport
			root_vp <- viewport(layout = grid.layout(3, 3, 
						heights = unit(c(0, 1, 1), c("lines", "null", "lines")),
						widths = unit(c(1, 1, 1), c("lines", "null", "lines"))), 
						name = "root")       
			pushViewport(root_vp)

			## setup viewport for tree
			tree_vp <- viewport(layout.pos.col = 2, layout.pos.row = 2, 
						xscale = c(0, nx), yscale = c(0, ny + 1), 
								name = "tree")
			pushViewport(tree_vp)

			if((nx <= 1 & ny <= 1)) {
			  pushViewport(plotViewport(margins = rep(1.5, 4), name = paste("Node", ptr$nodeID, sep = "")))
			  terminal_panel(ptr)    
			} else {
				## call the workhorse
			upViewport()
			if (pop) popViewport() else upViewport()
			}
		}
		else if(type %in% c("km", "naive")){
			oldpar <- par(mar=c(2,2,4,2)+.1)
			on.exit(invisible(par(oldpar)))

			if (nrow(x@frame) <= 1L) stop("fit is not a tree, just a root")

			if (compress & missing(nspace)) nspace <- branch
			if (!compress) nspace <- -1L     #means no compression
			if (dev.cur() == 1L) dev.set()
			assign(paste(".studi.parms", dev.cur(), sep = "."),
			list(uniform=uniform, branch=branch, nspace=nspace,
			 minbranch=minbranch), envir=.GlobalEnv)

			#define the plot region
			temp <- gtlco(x)
			xx <- temp$x
			yy <- temp$y
			temp1 <- range(xx) + diff(range(xx))*c(-margin, margin)
			temp2 <- range(yy) + diff(range(yy))*c(-margin, margin)
			plot(temp1, temp2, type='n', axes=FALSE, xlab='', ylab='', ...)

			# Draw a series of horseshoes or V's, left son, up, down to right son
			# NA's in the vector cause lines() to "lift the pen"
			node <- as.numeric(format(x@frame$Node))
			temp <- studi.branch(xx, yy, node, branch)

			if (branch > 0) text(xx[1L], yy[1L], '|', col = "gray")
			lines(c(temp$x), c(temp$y), col = "gray")
			text(x, type = type, x.range = temp1[2], y.range = temp2[2])
			invisible(list(x=xx, y=yy))
		}
	}
)
####################################################################
setGeneric("text")
setMethod("text", "Tree",
	function(x, splits = TRUE, all=FALSE,
		digits = getOption("digits") - 3, fwidth=.7, fheight =.95, type = "naive",
		x.range, y.range, ...){
		# SCCS @(#)text.rpart.s 1.12 06/06/01
		# This is a modified version of text.tree function.
		# Again!!
		# This is a modified version of text.rpart function.
		# Modified by Soo-heang, Eo
		# 2011-12-22

		if (nrow(x@frame) == 1)	stop("fit is not a tree, just a root")

		frame <- x@frame
		frame <- format(frame)
		
		cxy <- par("cxy")                   #character width and height
		if(!is.null(srt <- list(...)$srt) && srt == 90)  cxy <- rev(cxy)

		xy <- gtlco(x)

		node <- as.numeric(frame$Node)
		is.left <- (node%%2 == 0)            #left hand sons
		node.left <- node[is.left]
		parent <- match(node.left/2, node)

		left.child <- match(2 * node, node)
		right.child <- match(node * 2 + 1, node)
		vars <- labels(x, type = "var")
		rows <- labels(x,type = "plot")
		#text(xy$x, xy$y + .5 * cxy[2L], vars[left.child])
		#text(xy$x, xy$y - .5 * cxy[2L], rows[left.child])

		leaves <- if(all) rep(TRUE, nrow(frame)) else frame$terminal == 1L
		nodename <- paste("Node ", frame$Node, sep = "")
		stat <- paste("n= ", frame$n, sep = "")
		stat.n <- paste("med= ", frame$yval,sep = "")
		##find maximum length of stat
		maxlen <- max(string.bounding.box(stat.n)$columns) +1L
		maxht <- max(string.bounding.box(stat.n)$rows) + 1L

		if(fwidth<1)  a.length <- fwidth* cxy[1L] * maxlen
		else a.length <- fwidth*cxy[1L]

		if(fheight<1) b.length <- fheight* cxy[2L] * maxht
		else b.length <- fheight*cxy[2L]

		child <- match(node[frame$terminal == 1L],node)
		## parameters for type = km
		if(type == "km"){
			f <- update(x@formula, .~1)
			old.dir <- getwd()	
			setwd(system.file("data", package = "rap"))
			for(i in parent) oval(xy$x[i],xy$y[i], a= a.length/2, b=b.length)
			for(i in child) {
				## Save KM curve for a specific terminal node in temporary folder.
				trmnl<- x@data[x@where == node[i],]
				trmnl.fit <- survfit(f, data = trmnl)
				new.coords <- dostep(trmnl.fit$time, trmnl.fit$surv)
				#postscript(paste("node",i,".ps",sep=""))
					#plot(trmnl.fit, col = 4, lwd = 2,  conf.int= FALSE)
				#dev.off()
				## Step 1. transforming from a PostScript image format to a specialized RGML format
				#PostScriptTrace(paste("node",i,".ps",sep=""))
				## Step 2. reading RGML formats into R
				###FIXME: change the angle to 90 degree
				#petal <- readPicture(paste("node",i,".ps.xml",sep=""))
				## Step 3. R to grid
				###FIXME: have to be modified xy coordinate system
				pushViewport(viewport(x = unit( (xy$x[i] - 2 * cxy[1]) / (x.range + 7*cxy[1]), "npc"), y = unit( (xy$y[i] )/ (y.range + 1.5*cxy[2]), "npc"),
					width = unit(a.length / (x.range), "npc"), 
					height = unit(a.length / ( x.range), "npc"),
					just = c("right","top")
					))
					grid.rect(gp = gpar(col = "red", fill = colors()[606]))
					grid.lines(new.coords$x/max(new.coords$x), new.coords$y, gp = gpar(col ="blue", lwd = 2))
					grid.text(paste("Node ",i,sep = ""), y = unit(4.5, "lines"),just = c("center","top"))
					grid.text(stat.n[i], y = unit(3.5, "lines"),just = c("center","top"))
					#grid.symbols(petal)
					#grid.lines(trmnl.fit$time, trmnl.fit$surv, gp = gpar(lwd = 2, col = "red"))
				upViewport()
			}
			#setwd(old.dir)
		}
		else if(type == "naive"){
			for(i in parent) oval(xy$x[i],xy$y[i], a= a.length/2, b=b.length)
			for(i in child) rectangle(xy$x[i],xy$y[i], a=a.length/2,b=b.length)
			text(xy$x[leaves], xy$y[leaves] + .8 * cxy[2], nodename[leaves], adj=.5, ...)
			text(xy$x[leaves], xy$y[leaves] , stat[leaves], adj=.5, ...)
			text(xy$x[leaves], xy$y[leaves] - .8 * cxy[2], stat.n[leaves], adj=.5, ...)	
		}
		
		text(xy$x, xy$y + .5 * cxy[2L], vars[left.child], adj = .5, font = 2, cex= 1.2, ...)
		text(xy$x, xy$y - .5 * cxy[2L], rows[left.child], adj = .5, ...)
		text(xy$x[!leaves] + 2.8*cxy[1L], xy$y[!leaves] + .8 * cxy[2L], nodename[!leaves], pos = 4, cex=.9, ...)
		text(xy$x[!leaves] + 2.8*cxy[1L], xy$y[!leaves] , stat[!leaves],  pos = 4, cex=.9,...)
		text(xy$x[!leaves] + 2.8*cxy[1L], xy$y[!leaves] - .8 * cxy[2L], stat.n[!leaves], pos = 4, cex=.9,...)
	
		## stick values on nodes
		invisible()
	}
)

####################################################################
dostep <- function(x, y) {

    ### create a step function based on x, y coordinates
    ### borrowed from `party:dostep'
    if (is.na(x[1] + y[1])) {
        x <- x[-1]
        y <- y[-1]
    }
    n <- length(x)
    if (n > 2) {  
        # replace verbose horizonal sequences like
        # (1, .2), (1.4, .2), (1.8, .2), (2.3, .2), (2.9, .2), (3, .1)
        # with (1, .2), (3, .1).  They are slow, and can smear the looks
        # of the line type.
        dupy <- c(TRUE, diff(y[-n]) !=0, TRUE)
        n2 <- sum(dupy)

        #create a step function
        xrep <- rep(x[dupy], c(1, rep(2, n2-1)))
        yrep <- rep(y[dupy], c(rep(2, n2-1), 1))
        RET <- list(x = xrep, y = yrep)
    } else {
        if (n == 1) {
            RET <- list(x = x, y = y)
        } else {
            RET <- list(x = x[c(1,2,2)], y = y[c(1,1,2)])
        }
    }
    return(RET)
}

tree.depth <- function (nodes){
    depth <- floor(log(nodes, base = 2) + 1e-7)
    as.vector(depth - min(depth))
}

studi.branch <- function(x, y, node, branch){
    if (missing(branch)) {
		if (exists(parms <-paste(".studi.parms", dev.cur(), sep="." ),
					   envir=.GlobalEnv)) {
				parms <- get(parms, envir=.GlobalEnv)
				branch <- parms$branch
			}
		else branch <- 0
    }

    # Draw a series of horseshoes, left son, up, over, down to right son
    #   NA's in the vector cause lines() to "lift the pen"
    is.left <- (node%%2 ==0)        #left hand sons
    node.left <- node[is.left]
    parent <- match(node.left/2, node)
    sibling <- match(node.left+1, node)
    temp <- (x[sibling] - x[is.left])*(1-branch)/2
    xx <- rbind(x[is.left], x[is.left] + temp,
                x[sibling]- temp, x[sibling], NA)
    yy <- rbind(y[is.left], y[parent], y[parent], y[sibling], NA)
    list(x=xx, y=yy)
}

# SCCS  @(#)formatg.s	1.3 06/06/01
# format a set of numbers using C's "g" format
#  It is applied on an element by element basis, which is more
#  appropriate for rpart output than the standard Splus format()
#  command.
# For instance if x=(123, 1.23, .00123)
#	  format(x) = "123.00000", "1.23000", "0.00123"
#  but formatg does not add all of those zeros to the first two numbers
#
formatg <- function(x, digits= unlist(options('digits')),
		         format= paste("%.", digits, "g", sep='')) {
    if (!is.numeric(x)) stop("x must be a numeric vector")

    n <- length(x)
    temp <- sprintf(format, x)
    if (is.matrix(x)) matrix(temp, nrow=nrow(x))
    else temp
}

#SCCS @(#)rpartco.s 1.7 02/07/00
# Compute the x-y coordinates for a tree
# This is a modification of rpartco.
# Modified by Soo-heang, Eo

gtlco <- function(tree, parms =  paste(".studi.parms", dev.cur(), sep = ".")){
    frame <- tree@frame
	frame <- format(frame)
    node <- as.numeric(frame$Node)
    depth <- tree.depth(node)
    is.leaf <- (frame$terminal == 1L)
    if (exists(parms, envir=.GlobalEnv)){
        parms <- get(parms, envir=.GlobalEnv)
        uniform <- parms$uniform
        nspace <-parms$nspace
        minbranch <- parms$minbranch
    }
    else {
        uniform <- FALSE
        nspace <- -1
        minbranch <- .4
    }

    if(uniform) y <- (1 + max(depth) - depth) / max(depth,4)
    else {                    
        y <- dev <- as.numeric(frame$nLogLik) * 2
        temp <- split(seq(node), depth)     
        parent <- match(floor(node/2), node)
        sibling <- match(ifelse(node %% 2, node - 1, node + 1), node)

        # assign the depths
        #for(i in temp[-1]) {
        #    temp2 <- dev[parent[i]] - (dev[i] + dev[sibling[i]])
        #    y[i] <- y[parent[i]] - temp2
        #}

        # For some problems, classification & loss matrices in particular
        #   the gain from a split may be 0.  This is ugly on the plot.
        # Hence the "fudge" factor of  .3* the average step
        
        fudge <-  minbranch * diff(range(y)) / max(depth)
        for(i in temp[-1]) {
            temp2 <- dev[parent[i]] - (dev[i] + dev[sibling[i]])
            haskids <- !(is.leaf[i] & is.leaf[sibling[i]])
            y[i] <- y[parent[i]] - ifelse(temp2<=fudge & haskids, fudge, temp2)
        }
        y <- y / (max(y))
    }

    # Now compute the x coordinates, by spacing out the leaves and then
    #   filling in
	x   <-  double(length(node))         #allocate, then fill it in below
    x[is.leaf] <- seq(sum(is.leaf))      # leaves at 1, 2, 3, ....
    left.child <- match(node * 2, node)
    right.child <- match(node * 2 + 1, node)

    # temp is a list of non-is.leaf, by depth
    temp <- split(seq(node)[!is.leaf], depth[!is.leaf])
    for(i in rev(temp))  x[i] <- 0.5 * (x[left.child[i]] + x[right.child[i]])

    if (nspace < 0) return(list(x=x, y=y))

    #
    # Now we get fancy, and try to do overlapping
    #
    #  The basic algorithm is, at each node:
    #      1: get the left & right edges, by depth, for the left and
    #           right sons, of the x-coordinate spacing.
    #      2: find the minimal free spacing.  If this is >0, slide the
    #           right hand son over to the left
    #      3: report the left & right extents of the new tree up to the
    #           parent
    #   A way to visualize steps 1 and 2 is to imagine, for a given node,
    #      that the left son, with all its descendants, is drawn on a
    #      slab of wood.  The left & right edges, per level, give the
    #      width of this board.  (The board is not a rectangle, it has
    #      'stair step' edges). Do the same for the right son.  Now
    #      insert some spacers, one per level, and slide right hand
    #      board over until they touch.  Glue the boards and spacer
    #      together at that point.
    #
    #  If a node has children, its 'space' is considered to extend left
    #    and right by the amount "nspace", which accounts for space
    #    used by the arcs from this node to its children.  For
    #    horseshoe connections nspace usually is 1.
    #
    #  To make it global for a recursive function, the x coordinate list
    #    is written into frame 0.
    #
    compress <- function(me, depth) 
    {
        lson <- me +1
        x <- x
        if (is.leaf[lson]) left <- list(left=x[lson], right=x[lson],depth=depth+1, sons=lson)
        else               left <- compress(me+1, depth+1)

        rson <- me + 1 + length(left$sons)        #index of right son
        if (is.leaf[rson]) right<- list(left=x[rson], right=x[rson],depth=depth+1, sons=rson)
        else               right<- compress(rson, depth+1)

        maxd <- max(left$depth, right$depth) - depth
        mind <- min(left$depth, right$depth) - depth

        # Find the smallest distance between the two subtrees
        #   But only over depths that they have in common
        # 1 is a minimum distance allowed
        slide <- min(right$left[1:mind] - left$right[1:mind]) -1
        if (slide >0) 
        { # slide the right hand node to the left
            x[right$sons] <- x[right$sons] - slide;
            x[me] <- (x[right$sons[1]] + x[left$sons[1]])/2
#       assign("x", x)
            x <<- x
        }
        else slide <- 0

        # report back
        if (left$depth > right$depth) 
        {
            templ <- left$left
            tempr <- left$right
            tempr[1:mind] <- pmax(tempr[1:mind], right$right -slide)
        }
        else 
        {
            templ <- right$left  - slide
            tempr <- right$right - slide
            templ[1:mind] <- pmin(templ[1:mind], left$left)
        }

        list(left = c(x[me]- nspace*(x[me] -x[lson]), templ),
             right= c(x[me]- nspace*(x[me] -x[rson]), tempr),
             depth= maxd+ depth, sons=c(me, left$sons, right$sons))
    }
#    assign('compress', compress)
#    assign('x', x)
#    assign('is.leaf', is.leaf)
#    assign('nspace', nspace)

#    temp <- compress(1, 1)
#    x <- get('x')
#    remove(c('compress', 'x', 'is.leaf', 'nspace'))
    list(x = x, y = y)
}


string.bounding.box <- function(s) {
	s2 <- strsplit(s, "\n")
	rows <- sapply(s2, length)
	columns <- sapply(s2, function(x) max(nchar(x)))
	list(columns=columns, rows=rows)
}

oval <- function(middlex,middley,a,b) {
	theta <- seq(0, 2*pi, pi/30)
	newx <- middlex + a * cos(theta)
	newy <- middley + b * sin(theta)
	polygon(newx,newy,border=TRUE,col=colors()[361])
}

rectangle <- function(middlex, middley,a,b) {
	newx <- middlex + c(a,a,-a,-a)
	newy <- middley + c(b,-b,-b,b)
	polygon(newx,newy,col=colors()[606], border = "red")
}

