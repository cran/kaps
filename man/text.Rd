\name{text}
\docType{methods}
\alias{text}
\alias{text,Tree-method}
\title{Methods for Function \code{text}}
\description{ Labels the current graphic device of the tree dendrogram with text.}
\usage{
\S4method{text}{Tree}(x, splits = TRUE, all=FALSE,
		 digits = getOption("digits") - 3, fwidth=.7, fheight =.95,
		 type = "naive", x.range, y.range,...)
}
\arguments{
  \item{x}{object of class \code{Tree}.}
 \item{splits}{logical flag. If TRUE (default), then the splits in the tree are labeled with the criterion for the split.}
 \item{all}{Logical. If TRUE, all nodes are labeled, otherwise just terminal nodes.}
 \item{digits}{number of significant digits to include in numerical labels.}
 \item{fwidth}{Relates to option fancy and the width of the ellipses and rectangles. If fwidth < 1 then it is a scaling factor (default = 0.8). If fwidth > 1 then it represents the number of character widths (for current graphical device) to use.}
 \item{fheight}{Relates to option fancy and the height of the ellipses and rectangles. If fheight <1 then it is a scaling factor (default = 0.8). If fheight > 1 then it represents the number of character heights (for current graphical device) to use.}
 \item{type}{ an argument to define the structure of terminal node. If type is "naive", print node number, number of sample, and median survival time. If type if "km", plot Kaplan-Meire curves for each terminal node. }
 \item{x.range}{Relates to viewport object}
 \item{y.range}{Relates to viewport object}
 \item{\ldots}{other arguments to the \code{text} function can be passed here.}
}
\keyword{methods}
