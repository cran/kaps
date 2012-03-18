\name{summary}
\docType{methods}
\alias{summary}
\alias{summary,apss-method}
\alias{summary,Tree-method}
\title{Methods for Function \code{summary} for the S4 classes "apss" and "Tree"}
\usage{
\S4method{summary}{apss}(object, K)
\S4method{summary}{Tree}(object)
}
\description{ This function provides tabloid information such as survival median, 1-, 3-, and 5 years actual suvival time for each terminal group as a table format.}
\arguments{
  \item{object}{an object of S4 class \code{apss} or \code{Tree}.}
  \item{K}{a scalar object to determine the number of subgroups K. If missing, the estimated subgroup K is selected.}
}
\seealso{
  \code{\link{kaps}} and \code{\linkS4class{apss}} class \cr
  \code{\link{lrtree}} and \code{\linkS4class{apss}} class
}
\keyword{methods}
