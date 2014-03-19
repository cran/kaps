\name{summary}
\docType{methods}
\alias{summary}
\alias{summary,kaps-method}
\title{Summarize an object "kaps"}
\description{ This function provides tabloid information such as survival median, 1-, 3-, and 5 years actual suvival time for each partition.}
\usage{
	\S4method{summary}{kaps}(object, K)
}
\arguments{
  \item{object}{an object with the class \code{kaps}.}
  \item{K}{a scalar object to determine the number of subgroups K. If missing, the estimated subgroup K is selected.}
}
\seealso{
  \code{\link{kaps}}
}
\keyword{methods}
