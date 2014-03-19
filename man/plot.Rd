\name{plot}
\docType{methods}
\alias{plot}
\alias{plot,kaps-method}
\alias{plot,kaps,missing-method}
\title{Visualize an object "kaps"}
\usage{
	\S4method{plot}{kaps}(x, y, ...) 
}
\description{
	\code{plot} method for "kaps" with extended facilities. It provides four panels consisting of the scatter plot, the Kaplan-Meier survival curves, overall p-values, and the worst-pair p-values against K.
}
\arguments{
	\item{x}{an object from \code{kaps}}
	\item{y}{missing object}
	\item{\dots}{other arguments to the \code{\link[=graphics]{plot.default}} function can be passed here. Among them, the argument \code{K} is a scalar object that plots the Kaplan-Meire survival curves for the K. If missing, it works with selected K in the model fitting.}
}
\details{ 
	This function generates four plots. The top left panel is the scatterplot of survival times against the selected prognostic factor with the line fitted by local censored regression using \code{locfit}. The top right panel is the Kaplan-Meier survival curves for the subgroups selected with the optimal \emph{K}. At the bottom are displayed the plots of the overall and worst-pair p-values against K. The dotted lines indicate thresholds for significance 0.05.The outputs for a specific K can also be printed out with the argument K.

	For the sake of the Kaplan-Meire curve with estimated subgroups, the function \code{km.curve} is provided. 
}
\seealso{
    \code{\link{kaps}} \cr
    \code{\link{km.curve}}
}
\keyword{methods}
