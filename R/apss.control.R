apss.control <- function(pre.pt = list(), scope = list(), rho = 0, ncl = 1, 
    fold = TRUE, lower.limit = 0, upper.limit = 12, 
	p.adjust.methods = c("holm", "hochberg", "hommel", "bonferroni",  "BH", "BY", "fdr", "none")){
	p.adjust.methods <- match.arg(p.adjust.methods)
	#method <- match.arg(method)
	res <- new("apssOptions")
	#res@method <- method
	res@lower.limit <- lower.limit
	res@upper.limit <- upper.limit
	res@fold <- fold
	res@pre.pt = pre.pt
	res@scope = scope
	res@rho = rho
	res@ncl = as.integer(ncl)
	res@p.adjust.methods = p.adjust.methods
	return(res)
}

