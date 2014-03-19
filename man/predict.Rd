\name{predict}
\docType{methods}
\alias{predict}
\alias{predict,kaps-method}
\title{Predict new values using the fitted object "kaps".}
\description{
	This function provides the predicted subgroup or test statistics.
}
\usage{
\S4method{predict}{kaps}(object, newdata, type = c("predict","kaps"))
}
\arguments{
  \item{object}{object from \code{kaps}.}
  \item{newdata}{An optinal argument in which the name of predicted object is located. If omitted, the dataset used in the model fitting is utilized.}
  \item{type}{a type of preddiction. If "predict", predict subgroups based on the fitted model. If a type is "kaps", it returns the overall and worst-pair test staiststics for estimated subgroups.}
}
\seealso{
  \code{\link{kaps}}
}
\keyword{methods}
