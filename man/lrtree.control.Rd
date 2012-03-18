\name{lrtree.control}
\alias{lrtree.control}
\title{Control for lrtree Models}
\description{
    various parameters that control aspects of the \code{lrtree} fit. 
}
\usage{
	lrtree.control(minsplit = 20, minbucket = round(minsplit/2), maxdepth = 5, eps = 3.84, WH.strd = 1,
		fitted.model = c("logrank","cart"),
		tree.size = c("Boot", "ML"), B = 200,
		random = FALSE, plot = FALSE, 
		rho = 0, pre.pt = NA,
        m.zero = 3, Jt = 2, classwt = FALSE, L.split = FALSE, time.varying = FALSE, ncl = 2
		)
}
\arguments{
    \item{minsplit}{
        the minimum number of subjects that must exist in a node, in order for 
        a split to be attempted. 
    }
    \item{minbucket}{
      the minimum number of observations in any terminal node.
      If only one of \code{minbucket} or \code{minsplit} is specified,
      the code either sets \code{minsplit} to \code{minbucket*3}
      or \code{minbucket} to \code{minsplit/3}, as appropriate.
    }
    \item{maxdepth}{
	  the maximum of tree depth. Default is 5.
    }
  \item{eps}{
    a numerical criteria in split variable statistic in order to stop growing tree. 
  }
  \item{WH.strd}{
    a numerical criteria in WH approximation in order to stop growing tree .
  }
  \item{fitted.model}{
     a parameter for a splitting and growing altorithm: \cr
	 \code{"logrank"} - fit Segal's log rank survival tree\cr
	 \code{"cart"} - fit survival cart algrotihm. \cr
  }
  \item{tree.size}{
    a parameter to determine the size of the tree:\cr
    \code{"Boot"} - use a bootstrapping split rule at each split.\cr
    \code{"ML"} - use Multi-Level direct stopping rule (Eo and Cho, 2012+). \cr
	\code{"BU"} - use the Bottom-Up approach (Segal, 1988). \cr
  }
  \item{B}{
   a scalar parameter to choose the number of bootrapping at each split.\cr
  }
  \item{random}{
    a logical parameter to treat multiple responses censored data.
  }
  \item{plot}{
    a logical parameter to draw test statistic at split node t.
  }
  \item{rho}{
    a numeric parameter to choose the type of methos in survfit function.
  }
  \item{pre.pt}{predetermined split point for split point selection. 
 }
  \item{m.zero}{
    a parameter to "RG" methods. Only splits yielding at least m.zero cases in each subnode.
  }
  \item{Jt}{
    numerical values to set criterion for constructing variable selection table.
  }
  \item{classwt}{
    Priors of the classes. Need not add up to one. Ignored for regression.
  }
  \item{L.split}{ 
    Logical parameter to use Linear splits. If TRUE, run linear split test.
  }
  \item{time.varying}{ 
    Logical parameter to treat time varying predictors. 
  }
  \item{ncl}{ 
    the number of cores to use for cost-complexity pruning.
  }
}
\seealso{
    \code{\link{lrtree}}
}
\keyword{tree}
