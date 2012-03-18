\name{plot}
\docType{methods}
\alias{plot}
\alias{plot,Tree-method}
\alias{plot,apss-method}
\title{Methods for Function \code{plot} of an S4 class object "Tree" and "apss".}
\description{
  Plot a tree diagram for \code{"Tree"} object.
  Plot an elbow polot for subgroups of the \code{"apss"} object.
}
\usage{
\S4method{plot}{apss}(x, K, ...)

\S4method{plot}{Tree}(x, uniform=TRUE, branch=0.2, margin = .1, newpage = TRUE, 
	type = "naive", compress=TRUE, nspace, minbranch=.5, ...)
}
\arguments{
 \item{x}{an object of the S4 class \code{"apss"} and \code{"Tree"} }
 \item{K}{a scalar object that plots the kaplan meire survival curves for the proposed number of K. If missing, it works with selected K in the fitting model.}
 \item{uniform}{if TRUE, uniform vertical spacing of the nodes is used; this may be less cluttered 
	when fitting a large plot onto a page. The default is to use a non-uniform spacing 
	proportional to the error in the fit.
 }
 \item{branch}{controls the shape of the branches from parent to child node. Any number 
	from 0 to 1 is allowed. A value of 1 gives square shouldered branches, a value of 0 
	give V shaped branches, with other values being intermediate.
}
 \item{margin}{an extra percentage of white space to leave around the borders of the tree.
	(Long labels sometimes get cut off by the default computation).
 }
 \item{newpage}{a logical parameter to select newpage.
 }
 \item{type}{a logical parameter to select the type of tree diagram.
 }
 \item{compress}{if FALSE, the leaf nodes will be at the horizontal plot coordinates of 1:nleaves. 
	If TRUE, the routine attempts a more compact arrangement of the tree. The compaction algorithm 
	assumes uniform=TRUE; surprisingly, the result is usually an improvement even when that is not the case.
}
 \item{nspace}{the amount of extra space between a node with children and a leaf, as compared to 
	the minimal space between leaves. Applies to compressed trees only. The default is the value of branch.
}
 \item{minbranch}{set the minimum length for a branch to minbranch times the average branch length. This 
	parameter is ignored if uniform=TRUE. Sometimes a split will give very little improvement, or even (in
	the classification case) no improvement at all. A tree with branch lengths strictly proportional to 
	improvement leaves no room to squeeze in node labels.
}
 \item{\dots}{other arguments to the \code{\link[=graphics]{plot.default}} function can be passed here.}
}
\details{ 
	This function provides elbow plot and Kaplan-Meire survival curves for S4 formulation \code{"apss"} and tree diagram for S4 formulation \code{"Tree"}, respectively. 
}
\seealso{
	\code{\link{lrtree}} and \code{\linkS4class{Tree}} class \cr
    \code{\link{kaps}} and \code{\linkS4class{apss}} class
}
\keyword{methods}
