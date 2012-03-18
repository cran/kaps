\name{predict}
\docType{methods}
\alias{predict}
\alias{predict,apss-method}
\title{Predict new data witgh the fitted model \code{kaps} .}
\usage{
\S4method{predict}{apss}(object, newdata)
}
\description{
This function produces predicted values, obtained by evaluationg the kaps in the frame newdata predict. 
}
\arguments{
  \item{object}{object of S4 class \code{apss}.}
  \item{newdata}{An optinal data frame in which to look for variables with which to predict. If omitted, the fitted values are used.}
}
\seealso{
  \code{\link{kaps}}
}
\keyword{methods}