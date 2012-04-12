\name{prune}
\docType{methods}
\alias{prune}
\alias{prune,Tree-method}
\title{Conduct cost-complexity pruining to obtain ideal tree size}
\usage{
\S4method{prune}{Tree}(treeres, K = 10, 
	cost = function(y, yhat, status) sum(abs(y - yhat) * status))
}
\description{
This function conduct a cost-complexity pruning in order to obtain ideal tree size. 
}
\arguments{
  \item{treeres}{object of S4 class \code{Tree}.}
  \item{K}{The number of partitions}
  \item{cost}{cost function}
}
\seealso{
  \code{\link{lrtree}}
}
\keyword{methods}
