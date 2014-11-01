\name{show}
\docType{methods}
\alias{show}
\alias{show,kaps-method}
\title{Show an object "kaps"}
\usage{
\S4method{show}{kaps}(object)
}
\description{
	It returns the outputs of the object "kaps" consisting of three parts. The first part displays the model formula with a dataset and the selected number for K. Next, the information regarding the selection of an optimal set of cut-off points is provided. Lastly, the p-values of pairwise two-sample test comparisons among all the pairs of subgroups are provided.
}
\arguments{
  \item{object}{object from \code{kaps}.}
}
\seealso{
  \code{\link{kaps}}
}
\keyword{methods}
