\name{predict}
\docType{methods}
\alias{predict}
\alias{predict,apss-method}
\alias{predict,Tree-method}
\title{Predict new data witgh the fitted model \code{kaps} .}
\usage{
\S4method{predict}{apss}(object, newdata, type = c("predict","apss"))
\S4method{predict}{Tree}(object, newdata)
}
\description{
This function produces predicted values, obtained by evaluationg the kaps in the frame newdata predict. 
}
\arguments{
  \item{object}{object of S4 class \code{apss}.}
  \item{newdata}{An optinal data frame in which to look for variables with which to predict. If omitted, the fitted values are used.}
  \item{type}{A logical argument. If "predict", predict an apss class. The type of "apss" is used for the internal useness.}
}
\seealso{
  \code{\link{kaps}}, \code{\link{lrtree}}
}
\keyword{methods}
