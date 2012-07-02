kaps.control <- function(pre.pt = list(), scope = list(), rho = 0, ncl = 1, 
    fold = TRUE, lower.limit = 0, upper.limit = 12,
	shortcut = TRUE, splits = c("logrank", "maxstat"),
	p.adjust.methods = c("none", "holm", "hochberg", "hommel", "bonferroni",  "BH", "BY", "fdr")){
	p.adjust.methods <- match.arg(p.adjust.methods)
	splits <- match.arg(splits)
	res <- new("kapsOptions")
	#res@method <- method
	res@lower.limit <- lower.limit
	res@upper.limit <- upper.limit
	res@fold <- fold
	res@pre.pt = pre.pt
	res@scope = scope
	res@rho = rho
	res@shortcut = shortcut
	res@splits = splits
	res@ncl = as.integer(ncl)
	res@p.adjust.methods = p.adjust.methods
	return(res)
}

