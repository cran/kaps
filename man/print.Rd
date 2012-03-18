\name{print}
\docType{methods}
\alias{print}
\alias{print,apss-method}
\alias{print,Tree-method}
\title{Methods for Function \code{print} for the "apss" class. }
\usage{
\S4method{print}{apss}(x, K)
}
\description{Print simple summary for \code{apss} calss with the specified number of subgroups K}
\arguments{
  \item{x}{an object of S4 class \code{apss}.}
  \item{K}{a scalar object to determine the number of subgroups K. If missing, the estimated subgroup K is selected.}
}
\seealso{
  \code{\link{kaps}} and \code{\linkS4class{apss}} class
}
\keyword{methods}
