## S4 classes for apss package
## Soo-Heang Eo, 2011-01-31
setOldClass("Formula")
setOldClass("apssOptions")
setOldClass("Surv")
##########################################
### A class for adaptive partitioning
setClass(Class = "apssOptions",
	representation = representation(
		pre.pt 		=	"list",
		scope		=	"list",
		lower.limit	=	"numeric",
		upper.limit	=	"numeric",
		rho			=	"numeric",
		fold		=	"logical",
		ncl			=	"integer",
		p.adjust.methods	= 	"character"
	)
)

setClass(Class = "apss", 
	representation = representation(
		call				= "language",
		formula				= "Formula",
		data				= "data.frame",
		where  				= "vector",  
		index				= "integer",
		Chisq 		 		= "numeric",
		pvalue				= "numeric",
		WH					= "numeric",
		t					= "numeric",
		Z					= "numeric",
		pair				= "numeric",
		split.var		 	= "character",
		split.pt		 	= "numeric",
		#g.ID				= "matrix",
		mindat				= "numeric",
		elbow				= "matrix",
		groups				= "vector",
		candid				= "list",
		Options				= "apssOptions"
	)
)


###########################################################
##############  Simple Tree Classes  ######################
###########################################################

#####################################
	## 1) Variable Control arguments
	### a. eps
	### c. m.zero
	### d. Jt
	### e. random: use random effects
	### f. model.type
	### g. resi.type
	### h. resi.type2
	### i. imp.type
	### j. na.type
	### k. opt
	### l. classwt

setClass(Class = "VarCtrl",
    representation = representation(
        eps = "numeric",
		m.zero = "numeric",
        Jt      = "integer",
        random   = "logical",
        na.type = "function",
		opt = "character",
		classwt = "logical",
		plot = "logical",
		rho = "numeric"
    ),
    prototype = list(
        eps = 1e-4,
		m.zero = 2,
        Jt      = as.integer(2),
        random   = FALSE,
        na.type = na.omit,
		opt = "optim",
		classwt =  FALSE,
		plot = FALSE,
		rho = 0
    )
)

	#####################################
	## 2) Split Control arguments
	### a. minsplit
	### b-1. WH.strd
	### b-2. minbucket
	### c. con.spt.sel
	### d. cat.spt.sel
	### f. L.split
	### g. time.varying
	
setClass(Class = "SplitCtrl",
	representation = representation(
		minsplit = "integer",
		minbucket = "integer",
        WH.strd   = "numeric",
		fitted.model = "character",
 		con.spt.sel = "character",
		cat.spt.sel = "character",
		L.split = "logical",
		time.varying = "logical",
		MultiLevel = "logical",
		pre.pt = "vector"
	),
	prototype = list(
		minsplit = as.integer(20),
		minbucket = as.integer(7),
        WH.strd   = 1,
		fitted.model = "studi",
		con.spt.sel = "Greedy",
		cat.spt.sel = "ESA",
		L.split = TRUE,
		time.varying = FALSE,
		MultiLevel = FALSE,
		pre.pt = NA
		),
	validity = function(object) {
        if (any(c(object@minsplit, object@minbucket, object@WH.strd) < 0)) {
            warning("no negative values allowed in objects of class ", 
                    sQuote("SplitCtrl"))
            return(FALSE)
        }
        return(TRUE)
    }
)

	#####################################
	## 3) Tree Growing Control arguments
	### a. maxdepth
	### b. pruning
	### c. ncl

setClass(Class = "TGCtrl",
	representation = representation(
		maxdepth = "integer",
		tree.size = "character",
		B = "numeric",
		ncl = "integer"
	),
	prototype = list(
		maxdepth = as.integer(5),
		pruning = "growing",
		B = 200,
		ncl = as.integer(2)
	)
)
	
	#####################################
	### Top-level hyper parameter arguments
	
setClass(Class = "HyperParaControl",
    representation = representation(
        VarCtrl   = "VarCtrl",
        SplitCtrl = "SplitCtrl",
        TGCtrl    = "TGCtrl"
    ),
    prototype = list(VarCtrl = new("VarCtrl"),
                     SplitCtrl = new("SplitCtrl"),
                     TGCtrl = new("TGCtrl")
    ),
    validity = function(object) {
        (validObject(object@VarCtrl) && 
        validObject(object@SplitCtrl)) &&
        validObject(object@TGCtrl)
    }
)


	#####################################
	### A class for data by recursive binary splits
setClass(Class = "dataset",
	representation = representation(
		Y = "Surv",
		X = "data.frame",
		Z = "data.frame",
		resid = "vector",
		resid.sign = "vector"
	)
)

	#####################################
	### S4 class for node result by recursive partitioning
setClass(Class = "SplitRes",
	representation = representation(
		Var = "character",
		WH.stat = "numeric",
		Point = "numeric",
		Set = "character",
		Impurity = "numeric",
		type = "character",
		Obs = "vector"
	), prototype = list(
		Var = "",
		WH.stat = 0,
		Point = 0,
		Set = "",
		Impurity = 0,
		type = "",
		Obs = vector(mode = "integer")
		)
)

	#####################################
	### S4 class using frame work for showing
#setClass(Class = "frame",
#	representation = representation(
#		Node = "numeric",
#		var = "character",
#		point = "numeric",
#		set = "character",
#		n = "numeric",
#		N = "numeric",
#		nLogLik = "numeric",
#		yval = "numeric",
#		slope = "numeric",
#		terminal = "logical"
#	)
#)


	#####################################
	### A class for binary trees   
setClass(Class = "Tree", 
	representation = representation(
		call				= "language",
		formula				= "Formula",
		data				= "data.frame",
		frame 				= "data.frame",
		controls			= "HyperParaControl",
		where  				= "vector",     
		subset 				= "NULL",
		weights 			= "NULL",
		alpha           	= "numeric" ,
		alpha.prime 		= "numeric"
	)
)


