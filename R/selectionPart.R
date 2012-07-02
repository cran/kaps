logrank.tree <- function(x, y, controls, verbose = FALSE){
################################################
### splitting algorithms for logrank and median survival tree
### Soo-Heang Eo
### x: design vector
### y: survival object
### controls: hyper parameters for studi object
################################################
	options(warn = -1)
	cut.pt <- 0
	test.stat <- 0
	test.tmp <- 0

	gr <- rep(NA, length(y))
	xx <- c()
	yy <- c()
	if(is.factor(x)){
		
	}
	else{
		if(is.na(any(controls@SplitCtrl@pre.pt))) pt <- sort(unique(x))
		else pt <- controls@SplitCtrl@pre.pt
		for(i in 1:(length(pt) - 1)){
			gr <- ifelse(x <= pt[i], 1, 2)
			if(length(unique(gr)) >= 2) {
				if((length(gr[gr == 1]) >= controls@SplitCtrl@minbucket) & 
				(length(gr[gr == 2]) >= controls@SplitCtrl@minbucket) ) {
					if(controls@SplitCtrl@fitted.model == "logrank") {
						fit <- survdiff(y ~ gr, rho = controls@VarCtrl@rho)
						test.tmp <- fit$chisq
					}
					else if(controls@SplitCtrl@fitted.model == "cart") {
						fit <- survfit(y ~ gr)
						test.tmp <- survmed(fit$surv, fit$time)
					}
					xx <- c(xx, pt[i])
					yy <- c(yy, test.tmp)
					if(test.tmp > test.stat){
						test.stat <- test.tmp
						cut.pt <- pt[i]
					}
				}
			}
		}
	}
	
	if(length(xx) >= 2 & controls@VarCtrl@plot) {
		#readline("Plot Log-rank satatistics at a node:")
		par(ask = TRUE)
		plot(xx,yy, type = "b", ylim = c(0,max(yy, na.rm = TRUE)), xlab="Cutoff point", ylab="Log-rank statistic",axes=FALSE)
		axis(2)
		axis(1, at=pt, labels=pt)
		box()
		abline(h=3.84, lty=2, col = 4)
	}	
	#cat("cut.pt:", cut.pt, "test-statiatic: ", Chisq, "\n")
	return(list(cut.pt = cut.pt, test.stat = test.stat))
}
