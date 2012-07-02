## Control the hyer parameters
lrtree.control <-
function(minsplit = 20, minbucket = round(minsplit/2), maxdepth = 5, eps = 3.84, WH.strd = 1,
		fitted.model = c("logrank","cart"),
		tree.size = c("direct", "Boot"), B = 200,
		random = FALSE, plot = FALSE, 
		rho = 0, pre.pt = NA,
        m.zero = 3, Jt = 2, classwt = FALSE, L.split = FALSE, time.varying = FALSE, ncl = 2
		) {
	## This function consists of 3 parts, Variable Control, Split Control, and Tree Growing Control.
	## Each of parts has S4 class.

    if (missing(minsplit) && !missing(minbucket))
        minsplit <- ceiling(minbucket * 2)
	fitted.model <- match.arg(fitted.model)
	#con.spt.sel <- match.arg(con.spt.sel)
	#cat.spt.sel <- match.arg(cat.spt.sel)
	#resi.type <- match.arg(resi.type)
	#imp.type <- match.arg(imp.type)
	tree.size <- match.arg(tree.size)
	
    res <- new("HyperParaControl")
	res@VarCtrl@eps <- eps
	res@VarCtrl@m.zero <- m.zero
	res@VarCtrl@Jt <- as.integer(Jt)
	res@VarCtrl@classwt <-classwt
	res@VarCtrl@random <-random
	res@VarCtrl@plot <-plot
	res@VarCtrl@rho <-rho
	
	res@SplitCtrl@fitted.model <- fitted.model
	res@SplitCtrl@minsplit <- as.integer(minsplit)
	res@SplitCtrl@minbucket <- as.integer(minbucket)
	res@SplitCtrl@WH.strd <- WH.strd
	#res@SplitCtrl@con.spt.sel <- con.spt.sel
	#res@SplitCtrl@cat.spt.sel <- cat.spt.sel
	res@SplitCtrl@L.split <- L.split
	res@SplitCtrl@time.varying <- time.varying
	#res@SplitCtrl@MultiLevel <- MultiLevel
	res@SplitCtrl@pre.pt <- pre.pt
	
	res@TGCtrl@maxdepth <- as.integer(maxdepth)
	res@TGCtrl@tree.size <- tree.size
	res@TGCtrl@B <- B
	res@TGCtrl@ncl <- as.integer(ncl)
	
    if (!validObject(res))
        stop("res is not a valid object of class", class(res))
    res
}
